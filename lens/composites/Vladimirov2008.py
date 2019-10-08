from __future__ import absolute_import, division, print_function

from lens.actor.process import State, merge_default_states, merge_default_updaters, dict_merge
from lens.actor.emitter import get_emitter, configure_emitter
from lens.environment.lattice_compartment import LatticeCompartment

# processes
from lens.processes.Endres2006_chemoreceptor import ReceptorCluster
from lens.processes.Vladimirov2008_motor import MotorActivity
from lens.processes.derive_volume import DeriveVolume


exchange_key = '__exchange'  # TODO -- this is declared in multiple locations

def initialize_vladimirov2008(config):
    config.update({
        'emitter': {
            'type': 'database',
            'url': 'localhost:27017',
            'database': 'simulations',
        }
    })

    # declare the processes
    receptor = ReceptorCluster(config)
    motor = MotorActivity(config)
    deriver = DeriveVolume(config)
    processes = {
        'receptor': receptor,
        'motor': motor,
        'deriver': deriver}

    # initialize the states
    default_states = merge_default_states(processes)
    default_updaters = merge_default_updaters(processes)

    # get environment ids, and make exchange_ids for external state
    environment_ids = []
    initial_exchanges = {}
    for process_id, process in processes.iteritems():
        roles = {role: {} for role in process.roles.keys()}
        initial_exchanges.update(roles)

    for role, state_ids in default_updaters.iteritems():
        for state_id, updater in state_ids.iteritems():
            if updater is 'accumulate':
                environment_ids.append(state_id)
                initial_exchanges[role].update({state_id + exchange_key: 0.0})

    # set states according to the compartment_roles mapping.
    # This will not generalize to composites with processes that have different roles
    compartment_roles = {
        'external': 'environment',
        'internal': 'cell'}

    states = {
        compartment_roles[role]: State(
            initial_state=dict_merge(default_states.get(role, {}), initial_exchanges.get(role, {})),
            updaters=default_updaters.get(role, {}))
        for role in default_states.keys()}
    # configure the states to the roles for each process
    topology = {
        'receptor': {
            'external': 'environment',
            'internal': 'cell'},
        'motor': {
            'external': 'environment',
            'internal': 'cell'},
        'deriver': {
            'internal': 'cell'},
        }

    # configure emitter
    emitter = configure_emitter(config, processes, topology)

    options = {
        'topology': topology,
        'emitter': emitter,
        'environment': 'environment',
        'compartment': 'cell',
        'exchange_key': exchange_key,
        'environment_ids': environment_ids,
    }

    # create the compartment
    return LatticeCompartment(processes, states, options)