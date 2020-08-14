'''
==================
Metabolism Process
==================

This module defines a :term:`process class` for modeling a cell's
metabolic processes with flux balance analysis (FBA). The cobrapy
FBA library is used for solving the problems. This supports metabolic
models from the `BiGG model database <http://bigg.ucsd.edu>`_,
and other configurations that can be passed to :py:class:`Metabolism`
to create models of metabolism.
'''

from __future__ import absolute_import, division, print_function

import os
import argparse
import logging as log

import numpy as np
import matplotlib.pyplot as plt

from vivarium.core.process import Process
from vivarium.core.composition import (
    simulate_process_in_experiment,
    save_timeseries,
    flatten_timeseries,
    load_timeseries,
    assert_timeseries_close,
    REFERENCE_DATA_DIR,
    PROCESS_OUT_DIR,
    plot_simulation_output
)
from vivarium.library.make_network import (
    get_reactions,
    make_network,
    save_network
)
from vivarium.data.synonyms import get_synonym
from vivarium.library.units import units
from vivarium.library.cobra_fba import CobraFBA
from vivarium.library.dict_utils import tuplify_port_dicts
from vivarium.library.regulation_logic import build_rule
from vivarium.plots.metabolism import plot_exchanges, BiGG_energy_carriers, energy_synthesis_plot
from vivarium.processes.derive_globals import AVOGADRO


NAME = 'metabolism'



def get_fg_from_counts(counts_dict, mw):
    composition_mass = sum([
        coeff / AVOGADRO * mw.get(mol_id, 0.0) * (units.g / units.mol)
        for mol_id, coeff in counts_dict.items()])  # g
    return composition_mass.to('fg')



class Metabolism(Process):
    """A general class that is configured to match specific models

    This metabolism process class models metabolism using flux balance
    analysis (FBA). The FBA problem is defined using the provided
    configuration parameters.

    To see how to configure the process manually, look at the source
    code for :py:func:`test_toy_metabolism`.

    :term:`Ports`:

    * **external**: Holds the state of molecules external to the FBA
      reactions. For a model of a cell's metabolism, this will likely
      hold metabolite concentrations in the extracellular space.
    * **internal**: Holds the state of molecules internal to the FBA.
      For a model of a cell's metabolism, this will probably be the
      cytosolic concentrations.
    * **reactions**: Holds the IDs of the modeled metabolic reactions.
      The linked :term:`store` does not need to be shared with any other
      processes.
    * **fields**: The environmental fields that will be updated with
      cell intake and uptake.
    * **dimensions**: Holds the dimensions of the environment.
    * **flux_bounds**: The bounds on the FBA, which are imposed by the
      availability of metabolites. For example, for the metabolism of a
      cell, the bounds represent the limits of transmembrane transport.
    * **global**: Should be linked to the ``global`` :term:`store`.

    Args:
        initial_parameters (dict): Configures the process with the
            following keys/values:

            * **initial_state** (:py:class:`dict`): the default state,
              with a dict for internal and external, like this:
              ``{'external': external_state, 'internal':
              internal_state}``
            * **stoichiometry** (:py:class:`dict`): a map from reaction
              ID to that reaction's stoichiometry dictionary, e.g.
              ``{reaction_id: stoichiometry_dict}``
            * **objective** (:py:class:`dict`): the stoichiometry dict
              to be optimized
            * **external_molecules** (:py:class:`list`): the external
              molecules
            * **reversible_reactions** (:py:class:`list`)

    """

    name = NAME
    defaults = {
        'constrained_reaction_ids': [],
        'model_path': 'models/iAF1260b.json',
        'default_upper_bound': 0.0,
        'regulation': {},
        'initial_state': {},
        'exchange_threshold': 1e-4,  # concentrations lower than exchange_threshold are considered depleted
        'initial_mass': 1339,  # fg
        'global_deriver_key': 'global_deriver',
        'mass_deriver_key': 'mass_deriver',
        'time_step': 1,
    }

    def __init__(self, initial_parameters=None):
        if initial_parameters == None:
            initial_parameters = {}
        self.nAvogadro = AVOGADRO

        time_step = self.or_default(
            initial_parameters, 'time_step')

        # initialize FBA
        if 'model_path' not in initial_parameters and 'stoichiometry' not in initial_parameters:
            initial_parameters['model_path'] = self.defaults['model_path']
        self.fba = CobraFBA(initial_parameters)
        self.reaction_ids = self.fba.reaction_ids()
        self.exchange_threshold = self.defaults['exchange_threshold']

        # additional FBA options
        self.constrained_reaction_ids = self.or_default(
            initial_parameters, 'constrained_reaction_ids')
        self.default_upper_bound = self.or_default(
            initial_parameters, 'default_upper_bound')

        # make the regulation functions
        regulation_logic = self.or_default(
            initial_parameters, 'regulation')
        self.regulation = {
            reaction: build_rule(logic)
            for reaction, logic in regulation_logic.items()}

        # get internal molecules from fba objective
        self.objective_composition = {}
        for reaction_id, coeff1 in self.fba.objective.items():
            for mol_id, coeff2 in self.fba.stoichiometry[reaction_id].items():
                if mol_id in self.objective_composition:
                    self.objective_composition[mol_id] += coeff1 * coeff2
                else:
                    self.objective_composition[mol_id] = coeff1 * coeff2

        ## Get initial internal state from initial_mass
        initial_metabolite_mass = self.or_default(
            initial_parameters, 'initial_mass')
        mw = self.fba.molecular_weights
        composition = {
            mol_id: (-coeff if coeff < 0 else 0)
            for mol_id, coeff in self.objective_composition.items()}
        composition_mass = get_fg_from_counts(composition, mw)
        scaling_factor = (initial_metabolite_mass / composition_mass).magnitude
        internal_state = {mol_id: int(coeff * scaling_factor)
            for mol_id, coeff in composition.items()}
        self.initial_mass = get_fg_from_counts(internal_state, mw)
        log.info('metabolism initial mass: {}'.format(self.initial_mass))

        ## Get external state from minimal_external fba solution
        external_state = {state_id: 0.0 for state_id in self.fba.external_molecules}
        external_state.update(self.fba.minimal_external)  # optimal minimal media from fba

        # save initial state
        self.initial_state = {
            'external': external_state,
            'internal': internal_state,
            'flux_bounds': {reaction_id: self.default_upper_bound
                            for reaction_id in self.constrained_reaction_ids},
        }

        # parameters
        parameters = {'time_step': time_step}
        parameters.update(initial_parameters)

        self.global_deriver_key = self.or_default(
            initial_parameters, 'global_deriver_key')
        self.mass_deriver_key = self.or_default(
            initial_parameters, 'mass_deriver_key')

        super(Metabolism, self).__init__(parameters)

    def ports_schema(self):
        ports = [
            'internal',
            'external',
            'fields',
            'reactions',
            'flux_bounds',
            'global',
            'dimensions',
        ]

        schema = {port: {} for port in ports}

        # internal
        for state in list(self.objective_composition.keys()):
            schema['internal'][state] = {
                '_value': self.initial_state['internal'].get(state, 0),
                '_divider': 'split',
                '_default': 0.0,
                '_emit': True,
                '_properties': {
                    'mw': self.fba.molecular_weights[state] * units.g / units.mol},
            }

        # external
        for state in self.fba.external_molecules:
            schema['external'][state] = {
                '_default': self.initial_state['external'].get(state, 0.0),
                '_emit': True,
            }

        # fields
        for state in self.fba.external_molecules:
            schema['fields'][state] = {
                '_default': np.zeros((1, 1)),
            }

        # reactions
        for state in self.reaction_ids:
            schema['reactions'][state] = {
                '_default': 0.0,
                '_emit': state in self.constrained_reaction_ids,
                '_updater': 'set',
            }

        # flux_bounds
        for state in self.constrained_reaction_ids:
            schema['flux_bounds'][state] = {
                '_default': self.initial_state['flux_bounds'].get(state, self.default_upper_bound),
                '_emit': True,
            }

        # globals
        schema['global']['mass'] = {
            '_default': self.initial_mass,
            '_emit': True}
        schema['global'] = {
            'mmol_to_counts': {
                '_default': 0.0 * units.L / units.mmol,
                '_emit': True,
            },
            'location': {
                '_default': [0.5, 0.5],
            }
        }

        # dimensions
        schema['dimensions'] = {
            'bounds': {
                '_default': [1, 1],
            },
            'n_bins': {
                '_default': [1, 1],
            },
            'depth': {
                '_default': 1,
            },
        }

        return schema

    def derivers(self):
        return {
            self.mass_deriver_key: {
                'deriver': 'mass_deriver',
                'port_mapping': {
                    'global': 'global',
                },
                'config': {
                    'from_path': ('..', '..'),
                },
            },
            self.global_deriver_key: {
                'deriver': 'globals_deriver',
                'port_mapping': {
                    'global': 'global',
                },
                'config': {
                    'initial_mass': self.initial_mass,
                },
            }
        }

    def next_update(self, timestep, states):
        ## get the state
        external_state = states['external']
        constrained_reaction_bounds = states['flux_bounds']  # (units.mmol / units.L / units.s)
        mmol_to_counts = states['global']['mmol_to_counts']

        ## get flux constraints
        # exchange_constraints based on external availability
        exchange_constraints = {mol_id: 0.0
            for mol_id, conc in external_state.items() if conc <= self.exchange_threshold}

        # get state of regulated reactions (True/False)
        flattened_states = tuplify_port_dicts(states)
        regulation_state = {}
        for reaction_id, reg_logic in self.regulation.items():
            regulation_state[reaction_id] = reg_logic(flattened_states)

        ## apply flux constraints
        # first, add exchange constraints
        self.fba.set_exchange_bounds(exchange_constraints)

        # next, add constraints coming from flux_bounds
        # to constrain exchange fluxes, add the suffix 'EX_' to the external molecule ID
        if constrained_reaction_bounds:
            self.fba.constrain_flux(constrained_reaction_bounds)

        # finally, turn reactions on/off based on regulation
        self.fba.regulate_flux(regulation_state)

        ## solve the fba problem
        objective_exchange = self.fba.optimize() * timestep  # (units.mmol / units.L / units.s)
        exchange_reactions = self.fba.read_exchange_reactions()
        exchange_fluxes = self.fba.read_exchange_fluxes()  # (units.mmol / units.L / units.s)
        internal_fluxes = self.fba.read_internal_fluxes()  # (units.mmol / units.L / units.s)

        # timestep dependence on fluxes
        exchange_fluxes.update((mol_id, flux * timestep) for mol_id, flux in exchange_fluxes.items())
        internal_fluxes.update((mol_id, flux * timestep) for mol_id, flux in internal_fluxes.items())

        # update internal counts from objective flux
        # calculate added mass from the objective molecules' molecular weights
        objective_count = (objective_exchange * mmol_to_counts).magnitude
        internal_state_update = {}
        for reaction_id, coeff1 in self.fba.objective.items():
            for mol_id, coeff2 in self.fba.stoichiometry[reaction_id].items():
                if coeff2 < 0:  # pull out molecule if it is USED to make biomass (negative coefficient)
                    added_count = int(-coeff1 * coeff2 * objective_count)
                    internal_state_update[mol_id] = added_count

        # convert exchange fluxes to counts
        field_updates = {
            reaction: {
                '_value': int((flux * mmol_to_counts).magnitude),
                '_updater': {
                    'updater': 'update_field_with_exchange',
                    'port_mapping': {
                        'global': 'global',
                        'dimensions': 'dimensions',
                    },
                },
            }
            for reaction, flux in exchange_fluxes.items()
        }

        all_fluxes = {}
        all_fluxes.update(internal_fluxes)
        all_fluxes.update(exchange_reactions)

        return {
            'fields': field_updates,
            'internal': internal_state_update,
            'reactions': all_fluxes,
        }



# configs
def get_e_coli_core_config():
    """Get an *E. coli* core metabolism model

    The model is the `e_coli_core model from BiGG
    <http://bigg.ucsd.edu/models/e_coli_core>`_.

    Returns:
        A configuration for the model that can be passed to the
        :py:class:`Metabolism` constructor.
    """
    metabolism_file = os.path.join('models', 'e_coli_core.json')
    return {'model_path': metabolism_file}

def get_iAF1260b_config():
    """Get the metabolism config for the iAF1260b BiGG model

    The metabolism model is the `iAF1260b model from BiGG
    <http://bigg.ucsd.edu/models/iAF1260b>`_.

    Returns:
        A configuration for the model that can be passed to the
        :py:class:`Metabolism` constructor.
    """
    metabolism_file = os.path.join('models', 'iAF1260b.json')
    return {'model_path': metabolism_file}

def get_toy_configuration():
    stoichiometry = {
        'R1': {'A': -1, 'ATP': -1, 'B': 1},
        'R2a': {'B': -1, 'ATP': 2, 'NADH': 2, 'C': 1},
        'R2b': {'C': -1, 'ATP': -2, 'NADH': -2, 'B': 1},
        'R3': {'B': -1, 'F': 1},
        'R4': {'C': -1, 'G': 1},
        'R5': {'G': -1, 'C': 0.8, 'NADH': 2},
        'R6': {'C': -1, 'ATP': 2, 'D': 3},
        'R7': {'C': -1, 'NADH': -4, 'E': 3},
        'R8a': {'G': -1, 'ATP': -1, 'NADH': -2, 'H': 1},
        'R8b': {'G': 1, 'ATP': 1, 'NADH': 2, 'H': -1},
        'Rres': {'NADH': -1, 'O2': -1, 'ATP': 1},
        'v_biomass': {'C': -1, 'F': -1, 'H': -1, 'ATP': -10}}

    external_molecules = ['A', 'F', 'D', 'E', 'H', 'O2']

    objective = {'v_biomass': 1.0}

    reversible = ['R6', 'R7', 'Rres']

    default_reaction_bounds = 1000.0

    exchange_bounds = {
        'A': -0.02,
        'D': 0.01,
        'E': 0.01,
        'F': -0.005,
        'H': -0.005,
        'O2': -0.1}

    initial_state = {
        'external': {
            'A': 21.0,
            'F': 5.0,
            'D': 12.0,
            'E': 12.0,
            'H': 5.0,
            'O2': 100.0}}

    # molecular weight units are (units.g / units.mol)
    molecular_weights = {
        'A': 500.0,
        'B': 500.0,
        'C': 500.0,
        'D': 500.0,
        'E': 500.0,
        'F': 50000.0,
        'H': 1.00794,
        'O2': 31.9988,
        'ATP': 507.181,
        'NADH': 664.425}

    config = {
        'stoichiometry': stoichiometry,
        'reversible': reversible,
        'external_molecules': external_molecules,
        'objective': objective,
        'initial_state': initial_state,
        'exchange_bounds': exchange_bounds,
        'default_upper_bound': default_reaction_bounds,
        'molecular_weights': molecular_weights}

    return config


# toy functions
def make_kinetic_rate(mol_id, vmax, km=0.0):
    def rate(state):
        flux = (vmax * state[mol_id]) / (km + state[mol_id])
        return flux
    return rate

def toy_transport():
    # process-like function for transport kinetics, used by simulate_metabolism
    transport_kinetics = {
        "EX_A": make_kinetic_rate("A", -1e-1, 5),  # A import
    }
    return transport_kinetics

# sim functions
def run_sim_save_network(config=get_toy_configuration(), out_dir='out/network'):
    metabolism = Metabolism(config)

    # initialize the process
    stoichiometry = metabolism.fba.stoichiometry
    reaction_ids = list(stoichiometry.keys())
    external_mol_ids = metabolism.fba.external_molecules
    objective = metabolism.fba.objective

    settings = {
        # 'environment_volume': 1e-6,  # L   # TODO -- bring back environment?
        'timestep': 1,
        'total_time': 10}

    timeseries = simulate_process_in_experiment(metabolism, settings)
    reactions = timeseries['reactions']

    # save fluxes as node size
    reaction_fluxes = {}
    for rxn_id in reaction_ids:
        if rxn_id in reactions:
            flux = abs(np.mean(reactions[rxn_id][1:]))
            reaction_fluxes[rxn_id] = np.log(1000 * flux + 1.1)
        else:
            reaction_fluxes[rxn_id] = 1

    # define node type
    node_types = {rxn_id: 'reaction' for rxn_id in reaction_ids}
    node_types.update({mol_id: 'external_mol' for mol_id in external_mol_ids})
    node_types.update({rxn_id: 'objective' for rxn_id in objective.keys()})
    info = {
        'node_types': node_types,
        'reaction_fluxes': reaction_fluxes}

    nodes, edges = make_network(stoichiometry, info)
    save_network(nodes, edges, out_dir)

def run_metabolism(metabolism, settings=None):
    if not settings:
        settings = {
            'total_time': 10}
    return simulate_process_in_experiment(metabolism, settings)

# tests
def test_toy_metabolism():
    regulation_logic = {
        'R4': 'if (external, O2) > 0.1 and not (external, F) < 0.1'}

    toy_config = get_toy_configuration()
    transport = toy_transport()

    toy_config['constrained_reaction_ids'] = list(transport.keys())
    toy_config['regulation'] = regulation_logic
    toy_metabolism = Metabolism(toy_config)

    # TODO -- add molecular weights!

    # simulate toy model
    timeline = [
        (5, {('external', 'A'): 1}),
        (10, {('external', 'F'): 0}),
        (15, {})]

    settings = {
        'environment': {
            'volume': 1e-8 * units.L,
        },
        'timestep': 1.0,
        'timeline': {
            'timeline': timeline}}
    return simulate_process_in_experiment(toy_metabolism, settings)

def test_BiGG_metabolism(config=get_iAF1260b_config(), settings={}):
    metabolism = Metabolism(config)
    run_metabolism(metabolism, settings)

reference_sim_settings = {
    'environment': {
        'volume': 1e-5 * units.L,
    },
    'timestep': 1,
    'total_time': 10}

def metabolism_similar_to_reference():
    config = get_iAF1260b_config()
    metabolism = Metabolism(config)
    timeseries = run_metabolism(metabolism, reference_sim_settings)

    reference = load_timeseries(
        os.path.join(REFERENCE_DATA_DIR, NAME + '.csv'))
    assert_timeseries_close(timeseries, reference)


if __name__ == '__main__':
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    parser = argparse.ArgumentParser(description='metabolism process')
    parser.add_argument('--bigg', '-b', action='store_true', default=False,)
    args = parser.parse_args()

    if args.bigg:
        # configure BiGG metabolism
        config = get_iAF1260b_config()
        metabolism = Metabolism(config)

        # simulation settings
        sim_settings = {
            'environment': {
                'volume': 1e-5 * units.L,
            },
            'total_time': 2520,  # 2520 sec (42 min) is the expected doubling time in minimal media
        }

        # run simulation
        timeseries = simulate_process_in_experiment(metabolism, sim_settings)

        save_timeseries(timeseries, out_dir)
        volume_ts = timeseries['global']['volume']
        mass_ts = timeseries['global']['mass']
        print('volume growth: {}'.format(volume_ts[-1] / volume_ts[0]))
        print('mass growth: {}'.format(mass_ts[-1] / mass_ts[0]))

        # plot settings
        plot_settings = {
            'max_rows': 30,
            'remove_zeros': True,
            'skip_ports': ['exchange', 'reactions']}

        # make plots from simulation output
        plot_simulation_output(timeseries, plot_settings, out_dir, 'BiGG_simulation')
        plot_exchanges(timeseries, sim_settings, out_dir)

        # # make plot of energy reactions
        # stoichiometry = metabolism.fba.stoichiometry
        # energy_carriers = [get_synonym(mol_id) for mol_id in BiGG_energy_carriers]
        # energy_reactions = get_reactions(stoichiometry, energy_carriers)
        # energy_plot_settings = {'reactions': energy_reactions}
        # energy_synthesis_plot(timeseries, energy_plot_settings, out_dir)

        # make a gephi network
        run_sim_save_network(get_iAF1260b_config(), out_dir)

    else:
        timeseries = test_toy_metabolism()
        plot_settings = {}
        plot_simulation_output(timeseries, plot_settings, out_dir, 'toy_metabolism')
