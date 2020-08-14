from __future__ import absolute_import, division, print_function

import os
import sys
import argparse
import uuid
import random
import pytest

import numpy as np

from vivarium.core.emitter import path_timeseries_from_data
from vivarium.library.timeseries import (
    process_path_timeseries_for_csv,
    save_flat_timeseries,
)
from vivarium.core.experiment import (
    Experiment
)
from vivarium.core.composition import (
    agent_environment_experiment,
    simulate_experiment,
    plot_agents_multigen,
    EXPERIMENT_OUT_DIR,
    REFERENCE_DATA_DIR,
    load_timeseries,
    assert_timeseries_close,
)
from vivarium.plots.multibody_physics import (
    plot_snapshots,
    plot_tags
)

# processes
from vivarium.processes.metabolism import (
    Metabolism,
    get_iAF1260b_config,
)

# compartments
from vivarium.compartments.lattice import Lattice
from vivarium.compartments.growth_division import GrowthDivision
from vivarium.compartments.growth_division_minimal import GrowthDivisionMinimal
from vivarium.compartments.transport_metabolism import TransportMetabolism
from vivarium.compartments.flagella_expression import (
    FlagellaExpressionMetabolism,
    get_flagella_metabolism_initial_state,
)



NAME = 'lattice'
OUT_DIR = os.path.join(EXPERIMENT_OUT_DIR, NAME)
DEFAULT_ENVIRONMENT_TYPE = Lattice
TIME_STEP = 1


# agents and their configurations
agents_library = {
    'growth_division': {
        'name': 'growth_division',
        'type': GrowthDivision,
        'config': {
            'agents_path': ('..', '..', 'agents'),
            'fields_path': ('..', '..', 'fields'),
            'dimensions_path': ('..', '..', 'dimensions'),
        }
    },
    'growth_division_minimal': {
        'name': 'growth_division_minimal',
        'type': GrowthDivisionMinimal,
        'config': {
            'agents_path': ('..', '..', 'agents'),
            'fields_path': ('..', '..', 'fields'),
            'dimensions_path': ('..', '..', 'dimensions'),
            'growth_rate': 0.001,
            'division_volume': 2.6
        }
    },
    'flagella_metabolism': {
        'name': 'flagella_metabolism',
        'type': FlagellaExpressionMetabolism,
        'config': {
            'agents_path': ('..', '..', 'agents'),
            'fields_path': ('..', '..', 'fields'),
            'dimensions_path': ('..', '..', 'dimensions'),
        }
    },
    'transport_metabolism': {
        'name': 'transport_metabolism',
        'type': TransportMetabolism,
        'config': {
            'agents_path': ('..', '..', 'agents'),
            'fields_path': ('..', '..', 'fields'),
            'dimensions_path': ('..', '..', 'dimensions'),
            'metabolism': {'time_step': 10},
            'transport': {'time_step': 10},
        }
    },
}


# environment config
def get_lattice_config(
    time_step=TIME_STEP,
    bounds=[20, 20],
    n_bins=[10, 10],
    jitter_force=1e-4,
    depth=3000.0,
    diffusion=1e-2,
    molecules=['glc__D_e', 'lcts_e'],
    gradient={},
    keep_fields_emit=[],
):

    environment_config = {
        'multibody': {
            'time_step': time_step,
            'bounds': bounds,
            'jitter_force': jitter_force,
            'agents': {}
        },
        'diffusion': {
            # 'time_step': time_step,
            'molecules': molecules,
            'n_bins': n_bins,
            'bounds': bounds,
            'depth': depth,
            'diffusion': diffusion,
            'gradient': gradient,
            '_schema': {
                'fields': {
                    field_id: {
                        '_emit': False}
                    for field_id in molecules
                    if field_id not in keep_fields_emit}}
        }
    }

    return environment_config

def get_iAF1260b_environment(
    time_step=TIME_STEP,
    bounds=[20,20],
    n_bins=[10, 10],
    jitter_force=1e-4,
    depth=3000.0,
    scale_concentration=1,  # scales minimal media
    diffusion=5e-3,
    override_initial={},
    keep_fields_emit=[],
):
    # get external state from iAF1260b metabolism
    config = get_iAF1260b_config()
    metabolism = Metabolism(config)
    molecules = {
        mol_id: conc * scale_concentration
        for mol_id, conc in metabolism.initial_state['external'].items()
    }
    for mol_id, conc in override_initial.items():
        molecules[mol_id] = conc
    gradient = {
        'type': 'uniform',
        'molecules': molecules}
    return get_lattice_config(
        time_step=time_step,
        bounds=bounds,
        molecules=list(molecules.keys()),
        n_bins=n_bins,
        jitter_force=jitter_force,
        depth=depth,
        diffusion=diffusion,
        gradient=gradient,
        keep_fields_emit=keep_fields_emit,
    )

environments_library = {
    'glc_lcts': {
        'type': DEFAULT_ENVIRONMENT_TYPE,
        'config': get_lattice_config(
            bounds=[30,30],
            jitter_force=1e-5,
        ),
    },
    'iAF1260b': {
        'type': DEFAULT_ENVIRONMENT_TYPE,
        'config': get_iAF1260b_environment(
            bounds=[17, 17],
        ),
    },
    'shallow_iAF1260b': {
        'type': DEFAULT_ENVIRONMENT_TYPE,
        'config': get_iAF1260b_environment(
            time_step=10,
            bounds=[30, 30],
            n_bins=[40, 40],
            jitter_force=2e-3,
            depth=1e1,
            scale_concentration=10000,
            diffusion=1e-1,
            override_initial={
                'glc__D_e': 0.2,
                'lcts_e': 8.0},
            keep_fields_emit=[
                'glc__D_e',
                'lcts_e'],
        ),
    }
}


# simulation settings
def get_experiment_settings(
        experiment_name='lattice',
        description='an experiment in the lattice environment',
        total_time=4000,
        emit_step=10,
        emitter='timeseries',
        agent_names=False,
        return_raw_data=True,
):
    return {
        'experiment_name': experiment_name,
        'description': description,
        'total_time': total_time,
        'emit_step': emit_step,
        'emitter': emitter,
        'agent_names': agent_names,
        'return_raw_data': return_raw_data
    }


# plot settings
def get_plot_settings(
    skip_paths=[],
    fields=[],
    tags=[],
    n_snapshots=6,
    background_color='black',
):
    settings = {
        'plot_types': {
            'agents': {
                'skip_paths': skip_paths,
                'remove_zeros': True,
            },
        }
    }
    if fields:
        settings['plot_types']['snapshots'] = {
            'fields': fields,
            'n_snapshots': n_snapshots,
        }
    if tags:
        settings['plot_types']['tags'] = {
            'tagged_molecules': tags,
            'n_snapshots': n_snapshots,
            'background_color': background_color,
        }
    return settings

def plot_experiment_output(
        data,
        plot_settings={},
        out_dir='out',
):
    environment_config = plot_settings['environment_config']
    agent_type = plot_settings.get('agent_type', 'agent')
    plot_types = plot_settings['plot_types']

    # extract data
    multibody_config = environment_config['config']['multibody']
    agents = {time: time_data['agents'] for time, time_data in data.items()}
    fields = {time: time_data['fields'] for time, time_data in data.items()}

    # pass to plots
    if 'agents' in plot_types:
        plot_settings = plot_types['agents']
        plot_settings['agents_key'] = 'agents'
        plot_agents_multigen(data, plot_settings, out_dir, agent_type)

    if 'snapshots' in plot_types:
        plot_config = plot_types['snapshots']
        field_ids = plot_types['snapshots']['fields']
        plot_fields = {
            time: {
                field_id: field_instance[field_id]
                for field_id in field_ids}
            for time, field_instance in fields.items()}
        data = {
            'agents': agents,
            'fields': plot_fields,
            'config': multibody_config}
        plot_config.update({
            'out_dir': out_dir,
            'filename': agent_type + '_snapshots',
        })
        plot_snapshots(data, plot_config)

    if 'tags' in plot_types:
        plot_config = plot_types['tags']
        data = {
            'agents': agents,
            'config': multibody_config}
        plot_config.update({
            'out_dir': out_dir,
            'filename': agent_type + '_tags',
        })
        plot_tags(data, plot_config)


# Experiment run function
def run_lattice_experiment(
        agents_config=None,
        environment_config=None,
        initial_state=None,
        initial_agent_state=None,
        experiment_settings=None,
):
    if experiment_settings is None:
        experiment_settings = {}
    if initial_state is None:
        initial_state = {}
    if initial_agent_state is None:
        initial_agent_state = {}

    # agents ids
    agent_ids = []
    for config in agents_config:
        number = config.get('number', 1)
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

    # make the experiment
    experiment = agent_environment_experiment(
        agents_config=agents_config,
        environment_config=environment_config,
        initial_state=initial_state,
        initial_agent_state=initial_agent_state,
        settings=experiment_settings,
    )

    # simulate
    return simulate_experiment(
        experiment,
        experiment_settings,
    )


def run_workflow(
        agent_type='growth_division_minimal',
        n_agents=1,
        environment_type='glc_lcts',
        initial_state=None,
        initial_agent_state=None,
        out_dir='out',
        experiment_settings=get_experiment_settings(),
        plot_settings=get_plot_settings()
):
    if initial_state is None:
        initial_state = {}
    if initial_agent_state is None:
        initial_agent_state = {}
    # agent configuration
    agent_config = agents_library[agent_type]
    agent_config['number'] = n_agents
    agents_config = [
        agent_config,
    ]

    # environment configuration
    environment_config = environments_library[environment_type]

    # simulate
    data = run_lattice_experiment(
        agents_config=agents_config,
        environment_config=environment_config,
        initial_state=initial_state,
        initial_agent_state=initial_agent_state,
        experiment_settings=experiment_settings,
    )

    plot_settings['environment_config'] = environment_config
    plot_settings['agent_type'] = agent_type
    plot_experiment_output(
        data,
        plot_settings,
        out_dir,
    )


def test_growth_division_experiment():
    '''test growth_division_minimal agent in lattice experiment'''
    growth_rate = 0.005  # fast!
    total_time = 150

    # get minimal agent config and set growth rate
    agent_config = agents_library['growth_division_minimal']
    agent_config['config']['growth_rate'] = growth_rate
    agent_config['number'] = 1
    agents_config = [agent_config]

    # get environment config
    environment_config = environments_library['glc_lcts']

    # simulate
    experiment_settings = get_experiment_settings(
        total_time=total_time,
        return_raw_data=True)

    data = run_lattice_experiment(
        agents_config=agents_config,
        environment_config=environment_config,
        experiment_settings=experiment_settings)

    # assert division
    time = list(data.keys())
    initial_agents = len(data[time[0]]['agents'])
    final_agents = len(data[time[-1]]['agents'])
    assert final_agents > initial_agents


@pytest.mark.slow
def test_transport_metabolism_experiment(seed=1):
    random.seed(seed)
    np.random.seed(seed)
    experiment_name = 'transport_metabolism_experiment'
    agent_config = agents_library['transport_metabolism']
    agents_config = [agent_config]
    environment_config = environments_library['shallow_iAF1260b']
    experiment_settings = get_experiment_settings(total_time=60)

    # simulate
    data = run_lattice_experiment(
        agents_config=agents_config,
        environment_config=environment_config,
        # initial_state=initial_state,
        # initial_agent_state=initial_agent_state,
        experiment_settings=experiment_settings,
    )
    path_ts = path_timeseries_from_data(data)
    del path_ts[('agents', 'transport_metabolism', 'boundary', 'divide')]  # these contain 'False'

    # save timeseries as csv
    processed_for_csv = process_path_timeseries_for_csv(path_ts)
    make_dir(OUT_DIR)
    save_flat_timeseries(
        processed_for_csv,
        OUT_DIR,
        experiment_name + '.csv')

    # compare timeseries to reference
    test_output = load_timeseries(os.path.join(OUT_DIR, experiment_name + '.csv'))
    expected = load_timeseries(os.path.join(REFERENCE_DATA_DIR, experiment_name + '.csv'))
    assert_timeseries_close(
        test_output, expected,
        default_tolerance=(1 - 1e-5),
    )


def make_dir(out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

def main():
    out_dir = OUT_DIR
    make_dir(out_dir)

    parser = argparse.ArgumentParser(description='lattice_experiment')
    parser.add_argument('--growth_division', '-g', action='store_true', default=False)
    parser.add_argument('--growth_division_minimal', '-m', action='store_true', default=False)
    parser.add_argument('--flagella_metabolism', '-f', action='store_true', default=False)
    parser.add_argument('--transport_metabolism', '-t', action='store_true', default=False)
    args = parser.parse_args()
    no_args = (len(sys.argv) == 1)

    if args.growth_division_minimal or no_args:
        minimal_out_dir = os.path.join(out_dir, 'minimal')
        make_dir(minimal_out_dir)
        run_workflow(
            agent_type='growth_division_minimal',
            experiment_settings=get_experiment_settings(
                total_time=6000
            ),
            out_dir=minimal_out_dir)

    elif args.growth_division:
        gd_out_dir = os.path.join(out_dir, 'growth_division')
        make_dir(gd_out_dir)
        run_workflow(
            agent_type='growth_division',
            environment_type='glc_lcts',
            experiment_settings=get_experiment_settings(
                total_time=18000
            ),
            plot_settings=get_plot_settings(
                skip_paths=[
                    ('boundary', 'location')
                ],
                fields=[
                    'glc__D_e',
                    'lcts_e',
                ],
                tags=[
                    ('internal', 'protein1'),
                    ('internal', 'protein2'),
                    ('internal', 'protein3'),
                ],
                n_snapshots=5,
            ),
            out_dir=gd_out_dir)

    elif args.flagella_metabolism:
        txp_mtb_out_dir = os.path.join(out_dir, 'flagella_metabolism')
        make_dir(txp_mtb_out_dir)
        run_workflow(
            agent_type='flagella_metabolism',
            environment_type='iAF1260b',
            initial_agent_state=get_flagella_metabolism_initial_state(),
            experiment_settings=get_experiment_settings(
                emit_step=120,
                agent_names=True,
                emitter='database',
                total_time=12000,
            ),
            plot_settings=get_plot_settings(
                skip_paths=[
                    ('boundary', 'external')
                ],
                fields=[
                    'glc__D_e',
                ],
                tags=[
                    ('proteins', 'flagella'),
                ],
                background_color='black',
                n_snapshots=5,
            ),
            out_dir=txp_mtb_out_dir)

    elif args.transport_metabolism:
        txp_mtb_out_dir = os.path.join(out_dir, 'transport_metabolism')
        make_dir(txp_mtb_out_dir)
        run_workflow(
            agent_type='transport_metabolism',
            n_agents=2,
            environment_type='shallow_iAF1260b',
            initial_agent_state={
                'boundary': {
                    'external': {
                        'glc__D_e': 1.0,
                        'lcts_e': 1.0}}},
            out_dir=txp_mtb_out_dir,
            experiment_settings=get_experiment_settings(
                experiment_name='glucose lactose diauxie',
                description='glucose-lactose diauxic shifters are placed in a shallow environment with glucose and '
                           'lactose. They start off with no internal LacY and uptake only glucose, but LacY is '
                           'expressed upon depletion of glucose they begin to uptake lactose. Cells have an iAF1260b '
                           'BiGG metabolism, kinetic transport of glucose and lactose, and ode-based gene expression '
                           'of LacY',
                total_time=6000,
                emit_step=200,
                emitter='database',
            ),
            plot_settings=get_plot_settings(
                fields=['glc__D_e', 'lcts_e'],
                tags=[
                    ('cytoplasm', 'LacY'),
                ],
            ),
        ),


if __name__ == '__main__':
    main()
