from __future__ import absolute_import, division, print_function

import copy
import os

from vivarium.core.process import Generator
from vivarium.core.composition import (
    compartment_in_experiment,
    COMPARTMENT_OUT_DIR,
)
from vivarium.library.dict_utils import deep_merge

# processes
from vivarium.processes.multibody_physics import (
    Multibody,
    agent_body_config,
)
from vivarium.plots.multibody_physics import plot_snapshots
from vivarium.processes.diffusion_field import (
    DiffusionField,
    get_gaussian_config,
)
from vivarium.processes.derive_colony_shape import ColonyShapeDeriver


NAME = 'lattice'


class Lattice(Generator):
    """
    Lattice:  A two-dimensional lattice environmental model with multibody physics and diffusing molecular fields.
    """

    defaults = {
        # To exclude a process, from the compartment, set its
        # configuration dictionary to None, e.g. colony_mass_deriver
        'multibody': {
            'bounds': [10, 10],
            'size': [10, 10],
            'agents': {}
        },
        'diffusion': {
            'molecules': ['glc'],
            'n_bins': [10, 10],
            'size': [10, 10],
            'depth': 3000.0,
            'diffusion': 1e-2,
        },
        'colony_shape_deriver': None,
        '_schema': {},
    }

    def __init__(self, config=None):
        super(Lattice, self).__init__(config)

    def generate_processes(self, config):
        processes = {
            'multibody': Multibody(config['multibody']),
            'diffusion': DiffusionField(config['diffusion'])
        }
        colony_shape_config = config['colony_shape_deriver']
        if colony_shape_config is not None:
            processes['colony_shape_deriver'] = ColonyShapeDeriver(
                colony_shape_config)
        return processes

    def generate_topology(self, config):
        topology = {
            'multibody': {
                'agents': ('agents',),
            },
            'diffusion': {
                'agents': ('agents',),
                'fields': ('fields',),
                'dimensions': ('dimensions',),
            },
            'colony_shape_deriver': {
                'colony_global': ('colony_global',),
                'agents': ('agents',),
            }
        }
        return {
            process: process_topology
            for process, process_topology in topology.items()
            if config[process] is not None
        }


def get_lattice_config(config=None):
    if config is None:
        config = {}
    bounds = config.get('bounds', [25, 25])
    molecules = config.get('molecules', ['glc'])
    n_bins = config.get('n_bins', tuple(bounds))
    center = config.get('center', [0.5, 0.5])
    deviation = config.get('deviation', 5)
    diffusion = config.get('diffusion', 1e0)
    n_agents = config.get('n_agents', 1)
    agent_ids = [str(agent_id) for agent_id in range(n_agents)]

    # multibody config
    mbp_config = {
        # 'animate': True,
        'jitter_force': 1e2,
        'bounds': bounds}
    body_config = {
        'bounds': bounds,
        'agent_ids': agent_ids}
    mbp_config.update(agent_body_config(body_config))

    # diffusion config
    dff_config = get_gaussian_config({
        'molecules': molecules,
        'n_bins': n_bins,
        'bounds': bounds,
        'diffusion': diffusion,
        'center': center,
        'deviation': deviation})

    return {
        'bounds': bounds,
        'multibody': mbp_config,
        'diffusion': dff_config,
    }


def test_lattice(config=None, end_time=10):
    if config is None:
        config = get_lattice_config()
    # configure the compartment
    compartment = Lattice(config)

    # configure experiment
    experiment_settings = {
        'compartment': config}
    experiment = compartment_in_experiment(
        compartment,
        experiment_settings)

    # run experiment
    timestep = 1
    time = 0
    while time < end_time:
        experiment.update(timestep)
        time += timestep
    data = experiment.emitter.get_data()
    return data


def main():
    out_dir = os.path.join(COMPARTMENT_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    config = get_lattice_config()
    data = test_lattice(config, 40)

    # make snapshot plot
    agents = {time: time_data['agents'] for time, time_data in data.items()}
    fields = {time: time_data['fields'] for time, time_data in data.items()}
    data = {
        'agents': agents,
        'fields': fields,
        'config': config}
    plot_config = {
        'out_dir': out_dir,
        'filename': 'snapshots'}
    plot_snapshots(data, plot_config)


if __name__ == '__main__':
    main()
