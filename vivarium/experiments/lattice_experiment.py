from __future__ import absolute_import, division, print_function

import os
import uuid

from vivarium.compartment.tree import generate_state, Experiment
from vivarium.compartment.composition import (
    simulate_compartment,
    load_compartment,
    get_derivers
)

from vivarium.processes.multibody_physics import plot_snapshots
from vivarium.processes.diffusion_field import plot_field_output

# compartments
from vivarium.composites.lattice_environment import (
    make_lattice_environment,
)
from vivarium.composites.growth_division import growth_division




# # TODO -- this can be made into a general function
# def make_agents(settings):
#     n_agents = settings.get('n_agents', {})
#     compartment = settings.get('compartment', {})
#     config = settings.get('config', {})

#     processes = []
#     topologies = {}
#     agent_ids = []

#     for agent in range(n_agents):
#         agent_id = str(uuid.uuid1())
#         agent_ids.extend(agent_id)

#         # make the agent
#         agent_config = config.copy()
#         agent_config.update({'agent_id': agent_id})
#         agent = compartment(agent_config)  # TODO -- pass in compartment

#         # processes -- each process id is associated with its agent id with a tuple
#         a_processes = agent['processes']
#         a_processes = {
#             (agent_id, process_id): process
#             for process_id, process in a_processes.items()}

#         # topology
#         a_topology = {}
#         for process_id, topology in agent['topology'].items():
#             ports = {}
#             for port_id, store_id in topology.items():
#                 if store_id == BOUNDARY_STATE:
#                     ports[port_id] = store_id
#                 else:
#                     ports[port_id] = (agent_id, store_id)

#             a_topology[(agent_id, process_id)] = ports

#         # save processes and topology
#         processes.append(a_processes)
#         topologies[agent_id] = a_topology

#     return {
#         'agent_ids': agent_ids,
#         'processes': processes,
#         'topologies': topologies}


# TODO -- this can be made into a general function
def make_agents(count, compartment, config):
    processes = {}
    topology = {}

    for agent in range(count):
        # agent_id = str(uuid.uuid1())
        agent_id = str(agent)

        # make the agent
        agent_config = config.copy()
        agent_config.update({'agent_id': agent_id})
        agent = compartment(dict(
            agent_config,
            agent_id=agent_id))  # TODO -- pass in compartment

        # save processes and topology
        processes[agent_id] = {
            'cell': agent['processes']}
        topology[agent_id] = {
            'cell': agent['topology']}

    return {
        'processes': processes,
        'topology': topology}


# TODO -- this can move to a separate experiments directory
def lattice_experiment(config):
    # configure the experiment
    count = config.get('count')

    # get the environment
    environment = make_lattice_environment(config.get('environment', {}))
    environment_processes = environment['processes']
    environment_topology = environment['topology']
    inner_key = 'agents'  # TODO -- get this from config of each env process

    processes = {
        'environment': environment_processes}
    topology = {
        'environment': environment_topology}

    agents = make_agents(count, growth_division, {
        'cells_key': ['..', 'agents']})
    environment_processes['agents'] = agents['processes']
    environment_topology['agents'] = agents['topology']

    state = generate_state(
        environment_processes,
        environment_topology,
        config.get('initial_state', {}))

    experiment = Experiment({
        'processes': environment_processes,
        'topology': environment_topology,
        'initial_state': config.get('initial_state', {})})

    print('processes ------------------------')
    print(experiment.processes)

    print('topology ------------------------')
    print(experiment.topology)

    print('before ------------------------')
    print(experiment.state.get_value())

    experiment.update(10.0)

    print('after ----------------------------------')
    print(experiment.state.get_value())

    import ipdb; ipdb.set_trace()

    # ## add derivers
    # derivers = get_derivers(processes, topology)
    # deriver_processes = derivers['deriver_processes']
    # all_processes = processes + derivers['deriver_processes']
    # topology.update(derivers['deriver_topology'])  # add derivers to the topology

    # # initialize the states
    # # TODO -- pull out each agent_boundary, make a special initialize_state that can connect these up
    # stores = initialize_state(
    #     all_processes,
    #     topology,
    #     config.get('initial_state', {}))

    # print('state: '.format(stores[BOUNDARY_STATE].state))

    # options = {
    #     'name': config.get('name', 'lattice_environment'),
    #     'topology': topology,
    #     'initial_time': config.get('initial_time', 0.0)}

    # return {
    #     'processes': processes,
    #     'derivers': deriver_processes,
    #     'states': stores,
    #     'options': options}



# toy functions/ defaults
def get_lattice_config():

    environment_config = {
        'molecules': ['glc'],
        'bounds': [10, 10],
        'size': [10, 10],
    }

    agent_config = {}

    return {
        'count': 3,
        'environment': environment_config,
        'agents': agent_config
    }

def test_lattice_experiment(config=get_lattice_config(), time=10):
    lattice_environment = load_compartment(lattice_experiment, config)
    settings = {
        'return_raw_data': True,
        'total_time': time}
    return simulate_compartment(lattice_environment, settings)



if __name__ == '__main__':
    out_dir = os.path.join('out', 'tests', 'lattice_experiment')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    config = get_lattice_config()
    lattice_experiment(config)

    data = test_lattice_experiment(config, 10)

    # timeseries = get_timeseries(data)
    # plot_field_output(timeseries, config, out_dir, 'lattice_field')

    # make snapshot
    agents = {time: time_data['boundary'] for time, time_data in data.items()}
    fields = {}
    plot_snapshots(agents, fields, config, out_dir, 'lattice_bodies')
