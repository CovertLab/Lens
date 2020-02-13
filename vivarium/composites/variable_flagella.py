from __future__ import absolute_import, division, print_function

import copy
import os
import random

from vivarium.actor.process import initialize_state
from vivarium.actor.composition import get_derivers

# processes
from vivarium.processes.Endres2006_chemoreceptor import ReceptorCluster
from vivarium.processes.Mears2014_flagella_activity import FlagellaActivity
from vivarium.processes.membrane_potential import MembranePotential
from vivarium.processes.division import Division, divide_condition, divide_state



# the composite function
def compose_variable_flagella(config):

    ## Declare the processes

    ## Chemotaxis
    # receptor
    receptor_parameters = copy.deepcopy(config)
    receptor_parameters.update({'ligand': 'glc__D_e'})
    receptor = ReceptorCluster(receptor_parameters)

    # flagella
    flagella_config = copy.deepcopy(config)
    flagella_range = [0, 1, 5]  #list(range(0, 4))
    flagella_config.update({'n_flagella': random.choice(flagella_range)})
    flagella = FlagellaActivity(flagella_config)

    # proton motive force
    PMF = MembranePotential(config)

    # Other processes
    division_config = copy.deepcopy(config)
    division = Division(division_config)

    # Place processes in layers
    processes = [
        {'PMF': PMF},
        {'receptor': receptor},
        {'flagella': flagella},
        {'division': division}]

    ## Make the topology
    # for each process, map process roles to compartment roles
    topology = {
        'receptor': {
            'internal': 'cell',
            'external': 'environment'},
        'flagella': {
            'internal': 'cell',
            'counts': 'counts',
            'external': 'environment',
            'membrane': 'membrane',
            'flagella': 'flagella'},
        'PMF': {
            'internal': 'cell',
            'external': 'environment',
            'membrane': 'membrane'},
        'division': {
            'global': 'global'}}

    # add derivers
    derivers = get_derivers(processes, topology)
    processes.extend(derivers['deriver_processes'])  # add deriver processes
    topology.update(derivers['deriver_topology'])  # add deriver topology

    # initialize the states
    states = initialize_state(processes, topology, config.get('initial_state', {}))

    options = {
        'environment_role': 'environment',
        # 'exchange_role': 'exchange',
        'topology': topology,
        'initial_time': config.get('initial_time', 0.0),
        'divide_condition': divide_condition,
        'divide_state': divide_state}

    return {
        'processes': processes,
        'states': states,
        'options': options}


if __name__ == '__main__':
    from vivarium.actor.process import load_compartment
    from vivarium.actor.composition import simulate_with_environment, convert_to_timeseries, plot_simulation_output

    out_dir = os.path.join('out', 'tests', 'variable_flagella_composite')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    compartment = load_compartment(compose_variable_flagella)

    # settings for simulation and plot
    options = compartment.configuration

    # define timeline
    timeline = [(5.0, {})]

    settings = {
        'environment_role': options['environment_role'],
        # 'exchange_role': options['exchange_role'],
        'environment_volume': 1e-12,  # L
        'timeline': timeline}

    plot_settings = {
        'max_rows': 20,
        'skip_roles': ['prior_state']}

    # saved_state = simulate_compartment(compartment, settings)
    saved_data = simulate_with_environment(compartment, settings)
    del saved_data[0]  # remove first state
    timeseries = convert_to_timeseries(saved_data)
    plot_simulation_output(timeseries, plot_settings, out_dir)
