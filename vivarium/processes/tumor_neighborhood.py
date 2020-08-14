"""
================================================
Multibody physics process with neighbor tracking
================================================
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
from vivarium.processes.multibody_physics import (
    random_body_position,
    daughter_locations,
    make_random_position,
)


# plotting
from vivarium.plots.multibody_physics import (
    plot_agent,
    plot_agents,
    plot_snapshots,
)


NAME = 'multibody_neighbors'

DEFAULT_BOUNDS = [10, 10]

# constants
PI = math.pi


class MultibodyNeighbors(Process):
    """Simulates collisions and forces between agent bodies with a multi-body physics engine.

    :term:`Ports`:
    * ``agents``: The store containing all agent sub-compartments. Each agent in
      this store has values for location, angle, length, width, mass, thrust, and torque.

    Arguments:
        parameters(dict): Accepts the following configuration keys:

        * **jitter_force**: force applied to random positions along agent
          bodies to mimic thermal fluctuations. Produces Brownian motion.
        * **agent_shape** (:py:class:`str`): agents can take the shapes
          ``rectangle``, ``segment``, or ``circle``.
        * **bounds** (:py:class:`list`): size of the environment in
          micrometers, with ``[x, y]``.
        * ***animate*** (:py:class:`bool`): interactive matplotlib option to
          animate multibody. To run with animation turned on set True, and use
          the TKAgg matplotlib backend:

          .. code-block:: console

              $ MPLBACKEND=TKAgg python vivarium/processes/tumor_neighborhood.py
    """

    name = NAME
    defaults = {
        'agents': {},
        'jitter_force': 1e-3,  # pN
        'agent_shape': 'circle',
        'bounds': DEFAULT_BOUNDS,
        'animate': False,
        'time_step': 2,
    }

    def __init__(self, parameters=None):
        super(MultibodyNeighbors, self).__init__(parameters)

        # multibody parameters
        jitter_force = self.parameters['jitter_force']
        self.agent_shape = self.parameters['agent_shape']
        self.bounds = self.parameters['bounds']

        # make the multibody object
        time_step = self.parameters['time_step']
        multibody_config = {
            'agent_shape': self.agent_shape,
            'jitter_force': jitter_force,
            'bounds': self.bounds,
            'physics_dt': time_step / 10,
        }
        self.physics = PymunkMultibody(multibody_config)

        # interactive plot for visualization
        self.animate = self.parameters['animate']
        if self.animate:
            plt.ion()
            self.ax = plt.gca()
            self.ax.set_aspect('equal')

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
                    'diameter': {
                        '_emit': True,
                        '_default': 2.0,
                        '_divider': 'split',  # TODO -- might want this to be set by agent
                        '_updater': 'set'},
                    'mass': {
                        '_emit': True,
                        '_default': 1 * units.fg,
                        '_updater': 'set'},
                    'neighbors': {
                        # TODO -- tumors required cytotoxic_packets and PD1
                        # TODO -- t cells require PDL1 and MHCI
                    }
                }
            }
        }
        schema = {'cells': glob_schema}
        return schema

    def next_update(self, timestep, states):
        cells = states['cells']

        # animate before update
        if self.animate:
            self.animate_frame(cells)

        # update multibody with new agents
        self.physics.update_bodies(remove_units(cells))

        # run simulation
        self.physics.run(timestep)

        # get new agent positions
        cell_positions = self.physics.get_body_positions()


        import ipdb; ipdb.set_trace()
        # TODO -- get neighbors


        return {
            'cells': cell_positions}

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
class InvokeUpdate(object):
    def __init__(self, update):
        self.update = update
    def get(self, timeout=0):
        return self.update

def simulate_multibody_neighbors(config, settings):
    initial_agents_state = config['agents']

    # make the process
    multibody = MultibodyNeighbors(config)
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

                import ipdb;
                ipdb.set_trace()

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

def multibody_neighbors_workflow(config={}, out_dir='out', filename='neighbors'):
    total_time = config.get('total_time', 30)
    timestep = config.get('timestep', 0.05)
    config['initial_location'] = [0.5, 0.5]
    motility_config = get_baseline_config(config)

    # simulation settings
    motility_sim_settings = {
        'timestep': timestep,
        'total_time': total_time}

    # run motility sim
    motility_data = simulate_multibody_neighbors(motility_config, motility_sim_settings)
    motility_timeseries = timeseries_from_data(motility_data)


if __name__ == '__main__':
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    multibody_neighbors_workflow({'animate': True}, out_dir)
