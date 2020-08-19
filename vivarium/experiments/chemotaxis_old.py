'''
====================
Chemotaxis Experiments
====================

Chemotaxis provides several pre-configured :py:class:`Experiments`
with different chemotactic agents and environments.
'''

from __future__ import absolute_import, division, print_function

import os
import copy
import sys
import argparse

import numpy as np
from arpeggio import (
    RegExMatch,
    ParserPython,
    OneOrMore,
)

from vivarium.library.dict_utils import deep_merge
from vivarium.core.emitter import (
    time_indexed_timeseries_from_data,
    timeseries_from_data
)
from vivarium.core.composition import (
    agent_environment_experiment,
    make_agents,
    simulate_experiment,
    plot_agents_multigen,
    process_in_compartment,
    EXPERIMENT_OUT_DIR,
)

# compartments
from vivarium.compartments.static_lattice import StaticLattice
from vivarium.compartments.chemotaxis_minimal import ChemotaxisMinimal
from vivarium.compartments.chemotaxis_master import ChemotaxisMaster
from vivarium.compartments.chemotaxis_flagella import (
    ChemotaxisVariableFlagella,
    ChemotaxisExpressionFlagella,
    ChemotaxisODEExpressionFlagella,
)

# processes
from vivarium.processes.coarse_motor import MotorActivity
from vivarium.processes.multibody_physics import agent_body_config
from vivarium.processes.static_field import make_field

# plots
from vivarium.plots.multibody_physics import (
    plot_temporal_trajectory,
    plot_agent_trajectory,
    plot_motility,
)
from vivarium.plots.coarse_motor import plot_variable_receptor

# make an agent from a lone MotorActivity process
# MotorActivityAgent = MotorActivity
MotorActivityAgent = process_in_compartment(
    MotorActivity,
    topology={
        'external': ('boundary',),
        'internal': ('cell',)})

# environment defaults
DEFAULT_ENVIRONMENT_TYPE = StaticLattice
DEFAULT_LIGAND_ID = 'glc__D_e'  # BiGG id for external glucose
DEFAULT_BOUNDS = [1000, 3000]
DEFAULT_AGENT_LOCATION = [0.5, 0.1]

# exponential field parameters
FIELD_SCALE = 4.0
EXPONENTIAL_BASE = 1.3e2
FIELD_CENTER = [0.5, 0.0]

# get initial ligand concentration based on location
LOC_DX = (DEFAULT_AGENT_LOCATION[0] - FIELD_CENTER[0]) * DEFAULT_BOUNDS[0]
LOC_DY = (DEFAULT_AGENT_LOCATION[1] - FIELD_CENTER[1]) * DEFAULT_BOUNDS[1]
DIST = np.sqrt(LOC_DX ** 2 + LOC_DY ** 2)
INITIAL_LIGAND = FIELD_SCALE * EXPONENTIAL_BASE ** (DIST / 1000)


def get_exponential_env_config():
    # multibody process config
    multibody_config = {
        'bounds': DEFAULT_BOUNDS}

    # static field config
    field_config = {
        'molecules': [DEFAULT_LIGAND_ID],
        'gradient': {
            'type': 'exponential',
            'molecules': {
                DEFAULT_LIGAND_ID: {
                    'center': FIELD_CENTER,
                    'scale': FIELD_SCALE,
                    'base': EXPONENTIAL_BASE}}},
        'bounds': DEFAULT_BOUNDS}

    return {
        'multibody': multibody_config,
        'field': field_config}

def get_linear_env_config():
    # field parameters
    slope = 1.0
    base = 1e-1
    field_center = [0.5, 0.0]

    # multibody process config
    multibody_config = {
        'animate': False,
        'bounds': DEFAULT_BOUNDS}

    # static field config
    field_config = {
        'molecules': [DEFAULT_LIGAND_ID],
        'gradient': {
            'type': 'linear',
            'molecules': {
                DEFAULT_LIGAND_ID: {
                    'base': base,
                    'center': field_center,
                    'slope': slope,
                }
            }
        },
        'bounds': DEFAULT_BOUNDS}

    return {
        'multibody': multibody_config,
        'field': field_config}

DEFAULT_ENVIRONMENT_CONFIG = {
    'type': DEFAULT_ENVIRONMENT_TYPE,
    'config': get_exponential_env_config()}

DEFAULT_AGENT_CONFIG = {
    'ligand_id': DEFAULT_LIGAND_ID,
    'initial_ligand': INITIAL_LIGAND,
    'external_path': ('global',),
    'agents_path': ('..', '..', 'agents'),
    'daughter_path': tuple()}

def set_environment_config(config={}):
    # override default environment config
    return deep_merge(dict(DEFAULT_ENVIRONMENT_CONFIG), config)

def set_agent_config(config={}):
    # override default agent config
    return deep_merge(dict(DEFAULT_AGENT_CONFIG), config)


# configs with faster timescales, to support close agent/environment coupling
FAST_TIMESCALE = 0.001
tumble_jitter = 4000
# fast timescale minimal agents
FAST_MOTOR_CONFIG = set_agent_config({
        'tumble_jitter': tumble_jitter,
        'time_step': FAST_TIMESCALE})
FAST_MINIMAL_CHEMOTAXIS_CONFIG = set_agent_config({
    'receptor': {
        'time_step': FAST_TIMESCALE},
    'motor': {
        'tumble_jitter': tumble_jitter,
        'time_step': FAST_TIMESCALE}})

# fast timescale environment
FAST_TIMESCALE_ENVIRONMENT_CONFIG = set_environment_config({
    'config': {
        'multibody': {'time_step': FAST_TIMESCALE},
        'field': {'time_step': FAST_TIMESCALE} }})

# agent types
agents_library = {
    'motor': {
        'name': 'motor',
        'type': MotorActivityAgent,
        'config': FAST_MOTOR_CONFIG,
    },
    'minimal': {
        'name': 'minimal',
        'type': ChemotaxisMinimal,
        'config': FAST_MINIMAL_CHEMOTAXIS_CONFIG,
    },
    'variable': {
        'name': 'variable',
        'type': ChemotaxisVariableFlagella,
        'config': DEFAULT_AGENT_CONFIG
    },
    'expression': {
        'name': 'expression',
        'type': ChemotaxisExpressionFlagella,
        'config': DEFAULT_AGENT_CONFIG
    },
    'ode': {
        'name': 'ode_expression',
        'type': ChemotaxisODEExpressionFlagella,
        'config': DEFAULT_AGENT_CONFIG
    }
}

# preset experimental configurations
preset_experiments = {
    'minimal': {
        'agents_config': [
            {
                'number': 6,
                'name': 'minimal',
                'type': ChemotaxisMinimal,
                'config': DEFAULT_AGENT_CONFIG,
            }
        ],
        'environment_config': DEFAULT_ENVIRONMENT_CONFIG,
        'simulation_settings': {
            'total_time': 30,
            'emit_step': 0.1,
        },
    },
    'ode': {
        'agents_config': [
            {
                'number': 1,
                'name': 'ode_expression',
                'type': ChemotaxisODEExpressionFlagella,
                # 'config': DEFAULT_AGENT_CONFIG,
                'config': deep_merge(
                    dict(DEFAULT_AGENT_CONFIG),
                    {'growth_rate': 0.0005})  # fast growth
            }
        ],
        'environment_config': DEFAULT_ENVIRONMENT_CONFIG,
        'simulation_settings': {
            'total_time': 50,
            'emit_step': 1.0,
        },
    },
    'master': {
        'agents_config': [
            {
                'number': 1,
                'name': 'master',
                'type': ChemotaxisMaster,
                'config': DEFAULT_AGENT_CONFIG
            }
        ],
        'environment_config': DEFAULT_ENVIRONMENT_CONFIG,
        'simulation_settings': {
            'total_time': 30,
            'emit_step': 0.1,
        },
    },
    'variable': {
        'agents_config': [
            {
                'number': 1,
                'name': '{}_flagella'.format(n_flagella),
                'type': ChemotaxisVariableFlagella,
                'config': set_agent_config({'n_flagella': n_flagella})
            }
            for n_flagella in [0, 3, 6, 9, 12]
        ],
        'environment_config': DEFAULT_ENVIRONMENT_CONFIG,
        'simulation_settings': {
            'total_time': 720,
            'emit_step': 0.1,
        },
    },
    'motor': {
        'agents_config': [
            {
                'type': MotorActivityAgent,
                'name': 'motor',
                'number': 1,
                'config': FAST_MOTOR_CONFIG,
            }
        ],
        'environment_config': FAST_TIMESCALE_ENVIRONMENT_CONFIG,
        'simulation_settings': {
            'total_time': 120,
            'emit_step': FAST_TIMESCALE,
        },
    },
    'fast_minimal': {
        'agents_config': [
            {
                'number': 1,
                'name': 'motor + receptor',
                'type': ChemotaxisMinimal,
                'config': FAST_MINIMAL_CHEMOTAXIS_CONFIG,
            }
        ],
        'environment_config': FAST_TIMESCALE_ENVIRONMENT_CONFIG,
        'simulation_settings': {
            'total_time': 60,
            'emit_step': FAST_TIMESCALE,
        },
    },
    'mixed': {
        'agents_config': [
            {
                'type': ChemotaxisMinimal,
                'name': 'motor + receptor',
                'number': 1,
                'config': FAST_MINIMAL_CHEMOTAXIS_CONFIG,
            },
            {
                'type': MotorActivityAgent,
                'name': 'motor',
                'number': 1,
                'config': FAST_MOTOR_CONFIG,
            }
        ],
        'environment_config': FAST_TIMESCALE_ENVIRONMENT_CONFIG,
        'simulation_settings': {
            'total_time': 300,
            'emit_step': FAST_TIMESCALE,
        },
    },
    'many_mixed': {
        'agents_config': [
            {
                'type': MotorActivityAgent,
                'name': 'motor',
                'number': 5,
                'config': FAST_MOTOR_CONFIG,
            },
            {
                'type': ChemotaxisMinimal,
                'name': 'motor + receptor',
                'number': 5,
                'config': FAST_MINIMAL_CHEMOTAXIS_CONFIG,
            },
        ],
        'environment_config': FAST_TIMESCALE_ENVIRONMENT_CONFIG,
        'simulation_settings': {
            'total_time': 300,
            'emit_step': FAST_TIMESCALE*10,
        },
    },
}


def run_chemotaxis_experiment(
    agents_config=None,
    environment_config=None,
    initial_state=None,
    simulation_settings=None,
    experiment_settings=None):

    if not initial_state:
        initial_state = {}
    if not experiment_settings:
        experiment_settings = {}

    total_time = simulation_settings['total_time']
    emit_step = simulation_settings['emit_step']

    # agents ids
    agent_ids = []
    for config in agents_config:
        number = config['number']
        if 'name' in config:
            name = config['name']
            if number > 1:
                new_agent_ids = [name + '_' + str(num) for num in range(number)]
            else:
                new_agent_ids = [name]
        else:
            new_agent_ids = [str(uuid.uuid1()) for num in range(number)]
        config['ids'] = new_agent_ids
        agent_ids.extend(new_agent_ids)

    initial_agent_body = agent_body_config({
        'bounds': DEFAULT_BOUNDS,
        'agent_ids': agent_ids,
        'location': DEFAULT_AGENT_LOCATION})
    initial_state.update(initial_agent_body)

    # make the experiment
    experiment = agent_environment_experiment(
        agents_config,
        environment_config,
        initial_state,
        experiment_settings)

    # simulate
    settings = {
        'total_time': total_time,
        'emit_step': emit_step,
        'return_raw_data': True}
    return simulate_experiment(experiment, settings)


# plotting
def plot_chemotaxis_experiment(
        data,
        field_config,
        out_dir,
        filename=''):

    # multigen agents plot
    plot_settings = {
        'agents_key': 'agents',
        'max_rows': 30,
        'skip_paths': [
            ('boundary', 'mass'),
            ('boundary', 'length'),
            ('boundary', 'width'),
            ('boundary', 'location'),
        ]}
    plot_agents_multigen(data, plot_settings, out_dir, 'agents')

    # trajectory and motility
    indexed_timeseries = time_indexed_timeseries_from_data(data)
    field = make_field(field_config)
    trajectory_config = {
        'bounds': field_config['bounds'],
        'field': field,
        'rotate_90': True}

    plot_temporal_trajectory(copy.deepcopy(indexed_timeseries), trajectory_config, out_dir, filename + 'temporal')
    plot_agent_trajectory(copy.deepcopy(indexed_timeseries), trajectory_config, out_dir, filename + 'trajectory')

    embdedded_timeseries = timeseries_from_data(data)
    plot_motility(embdedded_timeseries, out_dir, filename + 'motility_analysis')


# parsing expression grammar for agents
def agent_type(): return RegExMatch(r'[a-zA-Z0-9.\-\_]+')
def number(): return RegExMatch(r'[0-9]+')
def specification(): return agent_type, number
def rule(): return OneOrMore(specification)
agent_parser = ParserPython(rule)

def make_agent_config(agent_specs):
    agent_type = agent_specs[0].value
    number = int(agent_specs[1].value)
    config = agents_library[agent_type]
    config['number'] = number
    return config

def parse_agents_string(agents_string):
    all_agents = agent_parser.parse(agents_string)
    agents_config = []
    for idx, agent_specs in enumerate(all_agents):
        agents_config.append(make_agent_config(agent_specs))
    return agents_config

def make_dir(out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

def add_arguments():
    parser = argparse.ArgumentParser(description='chemotaxis control')
    parser.add_argument(
        '--agents', '-a',
        type=str,
        nargs='+',
        default='"minimal 1"',
        help='A list of agent types and numbers in the format "agent_type1 number1 agent_type2 number2"')
    parser.add_argument(
        '--environment', '-v',
        type=str,
        default='exponential',
        help='the environment type ("linear" or "exponential")')
    parser.add_argument(
        '--time', '-t',
        type=int,
        default=10,
        help='total simulation time, in seconds')
    parser.add_argument(
        '--emit', '-m',
        type=int,
        default=1,
        help='emit interval, in seconds')
    parser.add_argument(
        '--experiment', '-e',
        type=str,
        default=None,
        help='preconfigured experiments')
    return parser.parse_args()


def run_chemotaxis_simulation():
    """
    Execute a chemotaxis simulation with any number of chemotactic agent types
    """
    out_dir = os.path.join(EXPERIMENT_OUT_DIR, 'chemotaxis')
    make_dir(out_dir)

    args = add_arguments()

    if args.experiment:
        # get a preset experiment
        # make a directory for this experiment
        experiment_name = str(args.experiment)
        control_out_dir = os.path.join(out_dir, experiment_name)
        make_dir(control_out_dir)

        experiment_config = preset_experiments[experiment_name]
        agents_config = experiment_config['agents_config']
        environment_config = experiment_config['environment_config']
        simulation_settings = experiment_config['simulation_settings']

    else:
        # make a directory for this experiment
        experiment_name = '_'.join(args.agents)
        control_out_dir = os.path.join(out_dir, experiment_name)
        make_dir(control_out_dir)

        # configure the agents
        agents_config = []
        if args.agents:
            agents_string = ' '.join(args.agents)
            agents_config = parse_agents_string(agents_string)

        # configure the environment
        if args.environment == 'linear':
            env_config = get_linear_env_config()
        else:
            env_config = get_exponential_env_config()
        environment_config = {
            'type': DEFAULT_ENVIRONMENT_TYPE,
            'config': env_config,
        }

        # configure the simulation
        total_time = args.time
        emit_step = args.emit
        simulation_settings = {
            'total_time': total_time,
            'emit_step': emit_step,
        }

    # simulate
    data = run_chemotaxis_experiment(
        agents_config=agents_config,
        environment_config=environment_config,
        simulation_settings=simulation_settings,
    )

    # plot
    field_config = environment_config['config']['field']
    plot_chemotaxis_experiment(data, field_config, control_out_dir)


if __name__ == '__main__':
    run_chemotaxis_simulation()
