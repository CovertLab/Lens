from __future__ import absolute_import, division, print_function

import copy
import numpy as np
from arrow import StochasticSystem

from vivarium.data.molecular_weight import molecular_weight
from vivarium.utils.dict_utils import keys_list
from vivarium.core.tree import Process
from vivarium.data.chromosomes.flagella_chromosome import FlagellaChromosome

chromosome = FlagellaChromosome()

def build_complexation_stoichiometry(
        stoichiometry,
        rates,
        reaction_ids,
        monomer_ids,
        complex_ids):

    molecule_ids = monomer_ids + complex_ids
    matrix = np.zeros((len(stoichiometry), len(molecule_ids)), dtype=np.int64)
    rates_array = np.zeros(len(stoichiometry))

    reverse_index = {
        molecule_id: index
        for index, molecule_id in enumerate(molecule_ids)}

    for reaction_index, reaction_id in enumerate(reaction_ids):
        reaction = stoichiometry[reaction_id]
        rates_array[reaction_index] = rates[reaction_id]
        for molecule_id, level in reaction.items():
            matrix[reaction_index][reverse_index[molecule_id]] = level

    return matrix, rates_array


class Complexation(Process):
    defaults = {
        'monomer_ids': chromosome.complexation_monomer_ids,
        'complex_ids': chromosome.complexation_complex_ids,
        'stoichiometry': chromosome.complexation_stoichiometry,
        'rates': chromosome.complexation_rates,
        'mass_deriver_key': 'mass_deriver'}

    def __init__(self, initial_parameters={}):
        self.default_parameters = copy.deepcopy(self.defaults)
        self.derive_defaults(initial_parameters, 'stoichiometry', 'reaction_ids', keys_list)

        self.parameters = self.default_parameters
        self.parameters.update(initial_parameters)

        self.monomer_ids = self.parameters['monomer_ids']
        self.complex_ids = self.parameters['complex_ids']
        self.reaction_ids = self.parameters['reaction_ids']

        self.stoichiometry = self.parameters['stoichiometry']
        self.rates = self.parameters['rates']

        self.complexation_stoichiometry, self.complexation_rates = build_complexation_stoichiometry(
            self.stoichiometry,
            self.rates,
            self.reaction_ids,
            self.monomer_ids,
            self.complex_ids)

        self.complexation = StochasticSystem(self.complexation_stoichiometry)

        self.mass_deriver_key = self.or_default(initial_parameters, 'mass_deriver_key')

        ports = {
            'monomers': self.monomer_ids,
            'complexes': self.complex_ids,
            'global': []}

        parameters = {}
        parameters.update(initial_parameters)

        super(Complexation, self).__init__(ports, parameters)

    def ports_schema(self):
        return {
            'monomers': {
                monomer: {
                    '_default': 0,
                    '_emit': True,
                    '_properties': {
                        'mass': molecular_weight[
                            monomer]} if monomer in molecular_weight else {}}
                for monomer in self.monomer_ids},
            'complexes': {
                complex: {
                    '_default': 0,
                    '_emit': True,
                    '_properties': {
                        'mass': molecular_weight[
                            complex]} if complex in molecular_weight else {}}
                for complex in self.complex_ids},
            'global': {}}

    def derivers(self):
        return {
            self.mass_deriver_key: {
                'deriver': 'mass',
                'port_mapping': {
                    'global': 'global'}}}

    def next_update(self, timestep, states):
        monomers = states['monomers']
        complexes = states['complexes']

        substrate = np.zeros(len(self.monomer_ids) + len(self.complex_ids), dtype=np.int64)

        for index, monomer_id in enumerate(self.monomer_ids):
            substrate[index] = monomers[monomer_id]
        for index, complex_id in enumerate(self.complex_ids):
            substrate[index + len(self.monomer_ids)] = complexes[complex_id]

        result = self.complexation.evolve(
            timestep,
            substrate,
            self.complexation_rates)

        outcome = result['outcome'] - substrate

        monomers_update = {
            monomer_id: outcome[index]
            for index, monomer_id in enumerate(self.monomer_ids)}

        complexes_update = {
            complex_id: outcome[index + len(self.monomer_ids)]
            for index, complex_id in enumerate(self.complex_ids)}

        update = {
            'monomers': monomers_update,
            'complexes': complexes_update}

        return update

def test_complexation():
    complexation = Complexation()
    state = {
        'monomers': {
            monomer: 1000
            for monomer in complexation.monomer_ids},
        'complexes': {
            complex: 0
            for complex in complexation.complex_ids}}

    update = complexation.next_update(1.0, state)
    print('initial state: {}'.format(state))
    print('complexation update: {}'.format(update))

if __name__ == '__main__':
    test_complexation()
