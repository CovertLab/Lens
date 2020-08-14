from __future__ import absolute_import, division, print_function

import os
import random

from vivarium.core.process import Generator
from vivarium.core.composition import (
    simulate_compartment_in_experiment,
    plot_simulation_output,
    COMPARTMENT_OUT_DIR
)

# processes
from vivarium.processes.chemoreceptor_cluster import (
    ReceptorCluster,
    get_exponential_random_timeline
)
from vivarium.processes.coarse_motor import MotorActivity



NAME = 'chemotaxis_minimal'

class ChemotaxisMinimal(Generator):

    defaults = {
        'ligand_id': 'MeAsp',
        'initial_ligand': 0.1,
        'boundary_path': ('boundary',),
        'receptor': {},
        'motor': {},
    }

    def __init__(self, config):
        super(ChemotaxisMinimal, self).__init__(config)

    def generate_processes(self, config):

        receptor_config = config['receptor']
        motor_config = config['motor']

        ligand_id = config['ligand_id']
        initial_ligand = config['initial_ligand']
        receptor_config.update({
            'ligand_id': ligand_id,
            'initial_ligand': initial_ligand})

        # declare the processes
        receptor = ReceptorCluster(receptor_config)
        motor = MotorActivity(motor_config)

        return {
            'receptor': receptor,
            'motor': motor}

    def generate_topology(self, config):
        boundary_path = config['boundary_path']
        external_path = boundary_path + ('external',)
        return {
            'receptor': {
                'external': external_path,
                'internal': ('cell',)},
            'motor': {
                'external': boundary_path,
                'internal': ('cell',)}}


def get_chemotaxis_config(config={}):
    ligand_id = config.get('ligand_id', 'MeAsp')
    initial_ligand = config.get('initial_ligand', 5.0)
    external_path = config.get('external_path', 'external')
    return {
        'external_path': (external_path,),
        'ligand_id': ligand_id,
        'initial_ligand': initial_ligand}


if __name__ == '__main__':
    out_dir = os.path.join(COMPARTMENT_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    environment_port = 'external'
    ligand_id = 'MeAsp'
    initial_conc = 0
    total_time = 60

    # configure timeline
    exponential_random_config = {
        'ligand': ligand_id,
        'environment_port': environment_port,
        'time': total_time,
        'timestep': 1,
        'initial_conc': initial_conc,
        'base': 1+4e-4,
        'speed': 14}

    # make the compartment
    config = {
        'ligand_id': ligand_id,
        'initial_ligand': initial_conc,
        'external_path': environment_port}
    compartment = ChemotaxisMinimal(get_chemotaxis_config(config))

    # run experiment
    experiment_settings = {
        'timeline': {
            'timeline': get_exponential_random_timeline(exponential_random_config),
            'ports': {'external': ('boundary', 'external')}},
        'timestep': 0.01,
        'total_time': 100}
    timeseries = simulate_compartment_in_experiment(compartment, experiment_settings)

    # plot settings for the simulations
    plot_settings = {
        'max_rows': 20,
        'remove_zeros': True,
        'overlay': {
            'reactions': 'flux'},
        'skip_ports': ['prior_state', 'null', 'global']}
    plot_simulation_output(
        timeseries,
        plot_settings,
        out_dir,
        'exponential_timeline')
