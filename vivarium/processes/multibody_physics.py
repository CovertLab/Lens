"""
==========================
Multibody physics process
==========================
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import argparse

import random
import math

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# vivarium imports
from vivarium.library.pymunk_multibody import PymunkMultibody
from vivarium.library.units import units, remove_units
from vivarium.core.emitter import timeseries_from_data
from vivarium.core.process import Process
from vivarium.core.composition import (
    process_in_experiment,
    simulate_experiment,
    PROCESS_OUT_DIR,
)
from vivarium.processes.coarse_motor import run, tumble
from vivarium.processes.derive_globals import volume_from_length

# plotting
from vivarium.plots.multibody_physics import (
    check_plt_backend,
    plot_agent,
    plot_agents,
    plot_snapshots,
    plot_temporal_trajectory,
    plot_motility,
)


NAME = 'multibody'

DEFAULT_BOUNDS = [10, 10]

# constants
PI = math.pi



def random_body_position(body):
    # pick a random point along the boundary
    width, length = body.dimensions
    if random.randint(0, 1) == 0:
        # force along ends
        if random.randint(0, 1) == 0:
            # force on the left end
            location = (random.uniform(0, width), 0)
        else:
            # force on the right end
            location = (random.uniform(0, width), length)
    else:
        # force along length
        if random.randint(0, 1) == 0:
            # force on the bottom end
            location = (0, random.uniform(0, length))
        else:
            # force on the top end
            location = (width, random.uniform(0, length))
    return location


def daughter_locations(parent_location, parent_values):
    parent_length = parent_values['length']
    parent_angle = parent_values['angle']
    pos_ratios = [-0.25, 0.25]
    daughter_locations = []
    for daughter in range(2):
        dx = parent_length * pos_ratios[daughter] * math.cos(parent_angle)
        dy = parent_length * pos_ratios[daughter] * math.sin(parent_angle)
        location = [parent_location[0] + dx, parent_location[1] + dy]
        daughter_locations.append(location)
    return daughter_locations



class Multibody(Process):
    """Simulates collisions and forces between agent bodies with a multi-body physics engine.

    :term:`Ports`:
    * ``agents``: The store containing all agent sub-compartments. Each agent in
      this store has values for location, angle, length, width, mass, thrust, and torque.

    Arguments:
        initial_parameters(dict): Accepts the following configuration keys:

        * **jitter_force**: force applied to random positions along agent
          bodies to mimic thermal fluctuations. Produces Brownian motion.
        * **agent_shape** (:py:class:`str`): agents can take the shapes
          ``rectangle``, ``segment``, or ``circle``.
        * **bounds** (:py:class:`list`): size of the environment in
          micrometers, with ``[x, y]``.
        * **mother_machine** (:py:class:`bool`): if set to ``True``, mother
          machine barriers are introduced.
        * ***animate*** (:py:class:`bool`): interactive matplotlib option to
          animate multibody. To run with animation turned on set True, and use
          the TKAgg matplotlib backend:

          .. code-block:: console

              $ MPLBACKEND=TKAgg python vivarium/processes/multibody_physics.py

    Notes:
        * rotational diffusion in liquid medium with viscosity = 1 mPa.s: :math:`Dr = 3.5 \pm0.3 rad^{2}/s`
          (Saragosti, et al. 2012. Modeling E. coli tumbles by rotational diffusion.)
        * translational diffusion in liquid medium with viscosity = 1 mPa.s: :math:`Dt = 100 um^{2}/s`
          (Saragosti, et al. 2012. Modeling E. coli tumbles by rotational diffusion.)
    """

    name = NAME
    defaults = {
        'agents': {},
        'jitter_force': 1e-3,  # pN
        'agent_shape': 'segment',
        'bounds': DEFAULT_BOUNDS,
        'mother_machine': False,
        'animate': False,
        'time_step': 2,
    }

    def __init__(self, initial_parameters=None):
        if initial_parameters is None:
            initial_parameters = {}

        # multibody parameters
        jitter_force = self.or_default(
            initial_parameters, 'jitter_force')
        self.agent_shape = self.or_default(
            initial_parameters, 'agent_shape')
        self.bounds = self.or_default(
            initial_parameters, 'bounds')
        self.mother_machine = self.or_default(
            initial_parameters, 'mother_machine')

        # make the multibody object
        self.time_step = self.or_default(
            initial_parameters, 'time_step')
        multibody_config = {
            'agent_shape': self.agent_shape,
            'jitter_force': jitter_force,
            'bounds': self.bounds,
            'barriers': self.mother_machine,
            'physics_dt': self.time_step / 10,
        }
        self.physics = PymunkMultibody(multibody_config)

        # interactive plot for visualization
        self.animate = initial_parameters.get('animate', self.defaults['animate'])
        if self.animate:
            plt.ion()
            self.ax = plt.gca()
            self.ax.set_aspect('equal')

        parameters = {'time_step': self.defaults['time_step']}
        parameters.update(initial_parameters)

        super(Multibody, self).__init__(parameters)

    def ports_schema(self):
        glob_schema = {
            '*': {
                'boundary': {
                    'location': {
                        '_emit': True,
                        '_default': [0.5 * bound for bound in self.bounds],
                        '_updater': 'set',
                        '_divider': {
                            'divider': daughter_locations,
                            'topology': {
                                'length': ('length',),
                                'angle': ('angle',)}}},
                    'length': {
                        '_emit': True,
                        '_default': 2.0,
                        '_divider': 'split',  # TODO -- might want this to be set by agent
                        '_updater': 'set'},
                    'width': {
                        '_emit': True,
                        '_default': 1.0,
                        '_updater': 'set'},
                    'angle': {
                        '_emit': True,
                        '_default': 0.0,
                        '_updater': 'set'},
                    'mass': {
                        '_emit': True,
                        '_default': 1339 * units.fg,
                        '_updater': 'set'},
                    'thrust': {
                        '_default': 0.0,
                        '_updater': 'set'},
                    'torque': {
                        '_default': 0.0,
                        '_updater': 'set'},
                }
            }
        }
        schema = {'agents': glob_schema}
        return schema

    def next_update(self, timestep, states):
        agents = states['agents']

        # animate before update
        if self.animate:
            self.animate_frame(agents)

        # update multibody with new agents
        self.physics.update_bodies(remove_units(agents))

        # run simulation
        self.physics.run(timestep)

        # get new agent positions
        agent_positions = self.physics.get_body_positions()
        update = {'agents': agent_positions}

        # for mother machine configurations, remove cells above the channel height
        if self.mother_machine:
            channel_height = self.mother_machine['channel_height']
            delete_agents = []
            for agent_id, position in agent_positions.items():
                location = position['boundary']['location']
                y_loc = location[1]
                if y_loc > channel_height:
                    # cell has moved past the channels
                    delete_agents.append(agent_id)
            if delete_agents:
                update['agents'] = {
                    agent_id: position
                    for agent_id, position in agent_positions.items()
                    if agent_id not in delete_agents}

                update['agents']['_delete'] = [
                    (agent_id,)
                    for agent_id in delete_agents]

        return update

    ## matplotlib interactive plot
    def animate_frame(self, agents):
        plt.cla()
        for agent_id, data in agents.items():
            # location, orientation, length
            data = data['boundary']
            x_center = data['location'][0]
            y_center = data['location'][1]
            angle = data['angle'] / PI * 180 + 90  # rotate 90 degrees to match field
            length = data['length']
            width = data['width']

            # get bottom left position
            x_offset = (width / 2)
            y_offset = (length / 2)
            theta_rad = math.radians(angle)
            dx = x_offset * math.cos(theta_rad) - y_offset * math.sin(theta_rad)
            dy = x_offset * math.sin(theta_rad) + y_offset * math.cos(theta_rad)

            x = x_center - dx
            y = y_center - dy

            if self.agent_shape == 'rectangle' or self.agent_shape == 'segment':
                # Create a rectangle
                rect = patches.Rectangle((x, y), width, length, angle=angle, linewidth=1, edgecolor='b')
                self.ax.add_patch(rect)

            elif self.agent_shape == 'circle':
                # Create a circle
                circle = patches.Circle((x, y), width, linewidth=1, edgecolor='b')
                self.ax.add_patch(circle)

        plt.xlim([0, self.bounds[0]])
        plt.ylim([0, self.bounds[1]])
        plt.draw()
        plt.pause(0.01)


# configs
def make_random_position(bounds):
    return [
        np.random.uniform(0, bounds[0]),
        np.random.uniform(0, bounds[1])]

def single_agent_config(config):
    # cell dimensions
    width = 1
    length = 2
    volume = volume_from_length(length, width)
    bounds = config.get('bounds', DEFAULT_BOUNDS)
    location = config.get('location')
    if location:
        location = [loc * bounds[n] for n, loc in enumerate(location)]
    else:
        location = make_random_position(bounds)

    return {'boundary': {
        'location': location,
        'angle': np.random.uniform(0, 2 * PI),
        'volume': volume,
        'length': length,
        'width': width,
        'mass': 1339 * units.fg,
        'thrust': 0,
        'torque': 0}}

def agent_body_config(config):
    agent_ids = config['agent_ids']
    agent_config = {
        agent_id: single_agent_config(config)
        for agent_id in agent_ids}
    return {'agents': agent_config}

def get_baseline_config(config={}):
    animate = config.get('animate', False)
    bounds = config.get('bounds', [500, 500])
    jitter_force = config.get('jitter_force', 0)
    n_agents = config.get('n_agents', 1)
    initial_location = config.get('initial_location')

    # agent settings
    agent_ids = [str(agent_id) for agent_id in range(n_agents)]
    motility_config = {
        'animate': animate,
        'jitter_force': jitter_force,
        'bounds': bounds}
    body_config = {
        'bounds': bounds,
        'agent_ids': agent_ids,
        'location': initial_location}
    motility_config.update(agent_body_config(body_config))
    return motility_config

# tests and simulations
def test_multibody(config={'n_agents':1}, time=10):
    n_agents = config.get('n_agents',1)
    agent_ids = [str(agent_id) for agent_id in range(n_agents)]

    body_config = agent_body_config({'agent_ids': agent_ids})
    multibody = Multibody(body_config)

    # initialize agent's boundary state
    initial_agents_state = body_config['agents']
    experiment = process_in_experiment(multibody)
    experiment.state.set_value({'agents': initial_agents_state})

    # run experiment
    settings = {
        'timestep': 1,
        'total_time': time,
        'return_raw_data': True}
    return simulate_experiment(experiment, settings)

def simulate_growth_division(config, settings):
    initial_agents_state = config['agents']

    # make the process
    multibody = Multibody(config)
    experiment = process_in_experiment(multibody)
    experiment.state.update_subschema(
        ('agents',), {
            'boundary': {
                'mass': {
                    '_divider': 'split'},
                'length': {
                    '_divider': 'split'}}})
    experiment.state.apply_subschemas()

    # get initial agent state
    experiment.state.set_value({'agents': initial_agents_state})
    agents_store = experiment.state.get_path(['agents'])

    ## run simulation
    # get simulation settings
    growth_rate = settings.get('growth_rate', 0.0006)
    growth_rate_noise = settings.get('growth_rate_noise', 0.0)
    division_volume = settings.get('division_volume', 0.4)
    total_time = settings.get('total_time', 10)
    timestep = 1

    time = 0
    while time < total_time:
        experiment.update(timestep)
        time += timestep
        agents_state = agents_store.get_value()

        agent_updates = {}
        remove_agents = []
        add_agents = {}
        for agent_id, state in agents_state.items():
            state = state['boundary']
            location = state['location']
            angle = state['angle']
            length = state['length']
            width = state['width']
            mass = state['mass'].magnitude

            # update
            growth_rate2 = (growth_rate + np.random.normal(0.0, growth_rate_noise)) * timestep
            new_mass = mass + mass * growth_rate2
            new_length = length + length * growth_rate2
            new_volume = volume_from_length(new_length, width)

            if new_volume > division_volume:
                daughter_ids = [str(agent_id) + '0', str(agent_id) + '1']

                daughter_updates = []
                for daughter_id in daughter_ids:
                    daughter_updates.append({
                        'daughter': daughter_id,
                        'path': (daughter_id,),
                        'processes': {},
                        'topology': {},
                        'initial_state': {}})

                # initial state will be provided by division in the tree
                update = {
                    '_divide': {
                        'mother': agent_id,
                        'daughters': daughter_updates}}
                invoked_update = InvokeUpdate({'agents': update})
                experiment.send_updates([invoked_update])
            else:
                agent_updates[agent_id] = {
                    'boundary': {
                        'volume': new_volume,
                        'length': new_length,
                        'mass': new_mass * units.fg}}

        # update experiment
        invoked_update = InvokeUpdate({'agents': agent_updates})
        experiment.send_updates([invoked_update])

    return experiment.emitter.get_data()


class InvokeUpdate(object):
    def __init__(self, update):
        self.update = update
    def get(self, timeout=0):
        return self.update

def simulate_motility(config, settings):
    # time of motor behavior without chemotaxis
    run_time = 0.42  # s (Berg)
    tumble_time = 0.14  # s (Berg)

    total_time = settings['total_time']
    timestep = settings['timestep']
    initial_agents_state = config['agents']

    # make the process
    config['time_step'] = timestep
    multibody = Multibody(config)
    experiment = process_in_experiment(multibody)
    experiment.state.update_subschema(
        ('agents',), {
            'boundary': {
                'thrust': {
                    '_emit': True,
                    '_updater': 'set'},
                'torque': {
                    '_emit': True,
                    '_updater': 'set'}},
            'cell': {
                'motor_state': {
                    '_value': 0,
                    '_updater': 'set',
                    '_emit': True,
                }}})
    experiment.state.apply_subschemas()

    # get initial agent state
    experiment.state.set_value({'agents': initial_agents_state})
    agents_store = experiment.state.get_path(['agents'])

    # initialize hidden agent motile states, and update agent motile forces in agent store
    agent_motile_states = {}
    motile_forces = {}
    for agent_id, specs in agents_store.get_value().items():
        [thrust, torque] = run()
        agent_motile_states[agent_id] = {
            'motor_state': 1,  # 0 for run, 1 for tumble
            'time_in_motor_state': 0}
        motile_forces[agent_id] = {
            'boundary': {
                'thrust': thrust,
                'torque': torque},
            'cell': {
                'motor_state': 1}}

    invoked_update = InvokeUpdate({'agents': motile_forces})
    experiment.send_updates([invoked_update])

    ## run simulation
    # test run/tumble
    time = 0
    while time < total_time:
        experiment.update(timestep)
        time += timestep

        # update motile force and apply to state
        motile_forces = {}
        for agent_id, motile_state in agent_motile_states.items():
            motor_state = motile_state['motor_state']
            time_in_motor_state = motile_state['time_in_motor_state']
            thrust = None

            if motor_state == 1:  # tumble
                if time_in_motor_state < tumble_time:
                    # [thrust, torque] = tumble()
                    time_in_motor_state += timestep
                else:
                    # switch
                    [thrust, torque] = run()
                    motor_state = 0
                    time_in_motor_state = 0

            elif motor_state == 0:  # run
                if time_in_motor_state < run_time:
                    # [thrust, torque] = run()
                    time_in_motor_state += timestep
                else:
                    # switch
                    [thrust, torque] = tumble()
                    motor_state = 1
                    time_in_motor_state = 0

            agent_motile_states[agent_id] = {
                'motor_state': motor_state,  # 0 for run, 1 for tumble
                'time_in_motor_state': time_in_motor_state}
            motile_forces[agent_id] = {
                'cell': {
                    'motor_state': motor_state}}

            if thrust:
                motile_forces[agent_id]['boundary'] = {
                    'thrust': thrust,
                    'torque': torque}

        invoked_update = InvokeUpdate({'agents': motile_forces})
        experiment.send_updates([invoked_update])

    return experiment.emitter.get_data()

def run_jitter(config={}, out_dir='out', filename='jitter'):
    total_time = config.get('total_time', 30)
    timestep = config.get('timestep', 0.05)
    motility_config = get_baseline_config({
        'animate': False,
        'jitter_force': 1e0,
        'bounds': [50, 50],
        'n_agents': 8,
    })

    # make the process
    multibody = Multibody(motility_config)
    experiment = process_in_experiment(multibody)
    experiment.state.update_subschema(
        ('agents',), {
            'cell': {
                'motor_state': {
                    '_value': 0,
                    '_updater': 'set',
                    '_emit': True,
                }}})
    experiment.state.apply_subschemas()

    time = 0
    while time < total_time:
        experiment.update(timestep)
        time += timestep
    data = experiment.emitter.get_data()

    # make trajectory plot
    timeseries = timeseries_from_data(data)
    plot_temporal_trajectory(timeseries, motility_config, out_dir, filename + '_trajectory')

def run_motility(config={}, out_dir='out', filename='motility'):
    total_time = config.get('total_time', 30)
    timestep = config.get('timestep', 0.05)
    config['initial_location'] = [0.5, 0.5]
    motility_config = get_baseline_config(config)

    # simulation settings
    motility_sim_settings = {
        'timestep': timestep,
        'total_time': total_time}

    # run motility sim
    motility_data = simulate_motility(motility_config, motility_sim_settings)
    motility_timeseries = timeseries_from_data(motility_data)

    # make motility plot
    plot_motility(motility_timeseries, out_dir, filename + '_analysis')
    plot_temporal_trajectory(motility_timeseries, motility_config, out_dir, filename + '_trajectory')

def run_growth_division(config={}):
    n_agents = 1
    agent_ids = [str(agent_id) for agent_id in range(n_agents)]

    bounds = [20, 20]
    settings = {
        'growth_rate': 0.02,
        'growth_rate_noise': 0.02,
        'division_volume': 2.6,
        'total_time': 140}

    gd_config = {
        'agent_shape': config.get('agent_shape', 'segment'),
        'animate': True,
        'jitter_force': 1e-3,
        'bounds': bounds}
    body_config = {
        'bounds': bounds,
        'agent_ids': agent_ids}
    gd_config.update(agent_body_config(body_config))
    gd_data = simulate_growth_division(gd_config, settings)

    # snapshots plot
    agents = {
        time: time_data['agents']
        for time, time_data in gd_data.items()
        if bool(time_data)}
    data = {
        'agents': agents,
        'config': gd_config}
    plot_config = {
        'out_dir': out_dir,
        'filename': 'growth_division_snapshots'}
    plot_snapshots(data, plot_config)

if __name__ == '__main__':
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    parser = argparse.ArgumentParser(description='multibody')
    parser.add_argument('--circles', '-c', action='store_true', default=False)
    parser.add_argument('--motility', '-m', action='store_true', default=False)
    parser.add_argument('--growth', '-g', action='store_true', default=False)
    parser.add_argument('--scales', '-s', action='store_true', default=False)
    parser.add_argument('--jitter', '-j', action='store_true', default=False)
    args = parser.parse_args()
    no_args = (len(sys.argv) == 1)

    if args.motility or no_args:
        run_motility({'animate': True}, out_dir)
    if args.growth or no_args:
        run_growth_division()
    if args.circles:
        run_growth_division({'agent_shape': 'circle'})
    if args.jitter:
        run_jitter({}, out_dir, 'jitter')
    if args.scales:
        bounds = [1000, 1000]
        jitter_force = 0

        ts_0p1 = {
            'timestep': 0.1,
            'bounds': bounds,
            'jitter_force': jitter_force}
        run_motility(ts_0p1, out_dir, 'ts_0p1')

        ts_0p01 = {
            'timestep': 0.01,
            'bounds': bounds,
            'jitter_force': jitter_force}
        run_motility(ts_0p01, out_dir, 'ts_0p01')

        # bounds_500 = {'bounds': [500, 500]}
        # run_motility(bounds_500, out_dir, 'bounds_500')
        #
        # bounds_5000 = {'bounds': [5000, 5000]}
        # run_motility(bounds_5000, out_dir, 'bounds_5000')
