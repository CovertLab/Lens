import copy
import numpy as np
from arrow import StochasticSystem

from vivarium.actor.process import Process
from vivarium.data.amino_acids import amino_acid_records
from vivarium.utils.datum import Datum
from vivarium.utils.polymerize import Elongation, Polymerase, Template, build_stoichiometry, build_rates, all_products

class Ribosome(Polymerase):
    pass

class Transcript(Template):
    pass

def generate_template(id, length, product):
    return Template({
        'id': id,
        'position': 0,
        'direction': 1,
        'sites': [],
        'terminators': [
            {'position': length,
             'strength': 1.0,
             'product': product}]})

def shuffle(l):
    l = [item for item in l]
    np.random.shuffle(l)
    return l

class Translation(Process):
    def __init__(self, initial_parameters={}):
        self.monomer_ids = [record['abbreviation'] for record in amino_acid_records]
        self.unbound_ribosomes_key = 'unbound_ribosomes'

        self.default_parameters = {
            'sequences': {
                'oA': shuffle(self.monomer_ids),
                'oAZ': shuffle(self.monomer_ids),
                'oB': shuffle(self.monomer_ids),
                'oBY': shuffle(self.monomer_ids)},
            'templates': {
                'oA': generate_template('oA', 20, ['eA']),
                'oAZ': generate_template('oAZ', 20, ['eA', 'eZ']),
                'oB': generate_template('oB', 20, ['eB']),
                'oBY': generate_template('oBY', 20, ['eB', 'eY'])},
            'transcript_affinities': {
                'oA': 1.0,
                'oAZ': 1.0,
                'oB': 1.0,
                'oBY': 1.0},
            'elongation_rate': 5.0,
            'advancement_rate': 1.0,
            'monomer_ids': self.monomer_ids}

        self.default_parameters['protein_ids'] = all_products(
            self.default_parameters['templates'])
        self.default_parameters['transcript_order'] = list(
            self.default_parameters['transcript_affinities'].keys())
        self.default_parameters['molecule_ids'] = self.monomer_ids + [
            self.unbound_ribosomes_key]

        parameters = copy.deepcopy(self.default_parameters)
        parameters.update(initial_parameters)

        self.transcript_affinities = parameters['transcript_affinities']
        self.transcript_order = parameters['transcript_order']
        self.transcript_count = len(self.transcript_order)

        self.monomer_ids = parameters['monomer_ids']
        self.molecule_ids = parameters['molecule_ids']
        self.protein_ids = parameters['protein_ids']
        self.elongation = 0
        self.elongation_rate = parameters['elongation_rate']
        self.advancement_rate = parameters['advancement_rate']

        self.sequences = parameters['sequences']
        self.templates = parameters['templates']

        self.affinity_vector = np.array([
            self.transcript_affinities[transcript_key]
            for transcript_key in self.transcript_order], dtype=np.float64)

        self.stoichiometry = build_stoichiometry(self.transcript_count)
        self.rates = build_rates(
            self.affinity_vector,
            self.advancement_rate)

        print('stoichiometry: {}'.format(self.stoichiometry))
        print('rates: {}'.format(self.rates))
        self.initiation = StochasticSystem(self.stoichiometry, self.rates)

        self.ribosome_id = 0

        self.roles = {
            'ribosomes': ['ribosomes'],
            'molecules': self.molecule_ids,
            'transcripts': self.transcript_order,
            'proteins': self.protein_ids}

        super(Translation, self).__init__(self.roles, parameters)

    def default_settings(self):
        default_state = {
            'ribosomes': {
                'ribosomes': []},
            'molecules': dict({
                self.unbound_ribosomes_key: 10}),
            'transcripts': {
                'oA': 1,
                'oAZ': 1,
                'oB': 1,
                'oBY': 1},
            'proteins': {
                protein_id: 0
                for protein_id in self.protein_ids}}

        default_state['molecules'].update({
            monomer_id: 50
            for monomer_id in self.monomer_ids})

        operons = list(default_state['transcripts'].keys())
        default_emitter_keys = {
            'ribosomes': ['ribosomes'],
            'molecules': self.monomer_ids + [self.unbound_ribosomes_key],
            'transcripts': operons,
            'proteins': self.protein_ids}

        default_updaters = {
            'ribosomes': {'ribosomes': 'set'},
            'molecules': {},
            'transcripts': {},
            'proteins': {}}

        return {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'updaters': default_updaters,
            'parameters': self.default_parameters}

    def next_update(self, timestep, states):
        ribosomes = list(map(Ribosome, states['ribosomes']['ribosomes']))
        molecules = states['molecules']
        transcripts = states['transcripts']
        transcript_counts = np.array([
            transcripts.get(transcript_key, 0)
            for transcript_key in self.transcript_order], dtype=np.int64)

        # Find out how many transcripts are currently blocked by a
        # newly initiated ribosome
        bound_transcripts = np.zeros(self.transcript_count, dtype=np.int64)
        ribosomes_by_transcript = {
            transcript_key: []
            for transcript_key in self.transcript_order}
        for ribosome in ribosomes:
            ribosomes_by_transcript[ribosome.template].append(ribosome)
        for index, transcript in enumerate(self.transcript_order):
            bound_transcripts[index] = len([
                ribosome
                for ribosome in ribosomes_by_transcript[transcript]
                if ribosome.is_bound()])

        # Make the state for a gillespie simulation out of total number of each
        # transcript not blocked by a bound ribosome, concatenated with the number
        # of each transcript that is bound by a ribosome.
        # These are the two states for each transcript the simulation
        # will operate on, essentially going back and forth between
        # bound and unbound states.

        original_unbound_ribosomes = states['molecules'][self.unbound_ribosomes_key]
        monomer_limits = {
            monomer: states['molecules'][monomer]
            for monomer in self.monomer_ids}
        unbound_ribosomes = original_unbound_ribosomes

        time = 0
        now = 0
        elongation = Elongation(
            self.sequences,
            self.templates,
            monomer_limits,
            self.elongation)

        while time < timestep:
            print('time: {} -----------======--------------------------'.format(time))

            # build the state vector for the gillespie simulation
            substrate = np.concatenate([
                transcript_counts - bound_transcripts,
                bound_transcripts,
                [unbound_ribosomes]])

            print('state: {}'.format(substrate))
            print('unbound ribosomes: {}'.format(unbound_ribosomes))
            print('bound transcripts: {}'.format(
                bound_transcripts))
            print('monomer limits: {}'.format(monomer_limits))

            # find number of monomers until next terminator
            # distance = chromosome.terminator_distance()
            distance = 1

            print('distance: {}'.format(distance))

            # find interval of time that elongates to the point of the next terminator
            interval = distance / self.elongation_rate

            print('interval: {}'.format(interval))
            print('substrates: {}'.format(substrate))

            # run simulation for interval of time to next terminator
            result = self.initiation.evolve(interval, substrate)

            # go through each event in the simulation and update the state
            ribosome_bindings = 0
            for now, event in zip(result['time'], result['events']):

                print('event {}: {}'.format(now, event))

                # perform the elongation until the next event
                terminations, monomer_limits, ribosomes = elongation.elongate(
                    time + now,
                    self.elongation_rate,
                    monomer_limits,
                    ribosomes)
                unbound_ribosomes += terminations

                # ribosome has bound the transcript
                if event < self.transcript_count:
                    transcript_key = self.transcript_order[event]
                    transcript = self.templates[transcript_key]
                    bound_transcripts[event] += 1
                    print('bound transcripts for {}: {}'.format(event, bound_transcripts[event]))

                    self.ribosome_id += 1
                    new_ribosome = Ribosome({
                        'id': self.ribosome_id,
                        'template': transcript_key,
                        'position': 0})
                    new_ribosome.bind()
                    ribosomes.append(new_ribosome)
                    ribosomes_by_transcript[transcript_key].append(new_ribosome)

                    print('{}: ribosome binding {}'.format(time, new_ribosome))

                    ribosome_bindings += 1
                    unbound_ribosomes -= 1
                # ribosome has begun polymerizing its protein
                else:
                    transcript_index = event - self.transcript_count
                    transcript_key = self.transcript_order[transcript_index]

                    bound_transcripts[transcript_index] -= 1
                    print('bound transcripts for {}: {}'.format(transcript_index, bound_transcripts[transcript_index]))

                    ribosome = ribosomes_by_transcript[transcript_key].pop()
                    ribosome.start_transcribing()

                    print('{}: ribosome commencing {}'.format(time, ribosome))

            print('bound transcripts: {}'.format(
                bound_transcripts))

            # now that all events have been accounted for, elongate
            # until the end of this interval.
            terminations, monomer_limits, ribosomes = elongation.elongate(
                time + interval,
                self.elongation_rate,
                monomer_limits,
                ribosomes)
            unbound_ribosomes += terminations

            print('bound ribosomes: {}'.format(ribosomes))
            print('complete transcripts: {}'.format(elongation.complete_polymers))
            print('monomer limits: {}'.format(monomer_limits))

            time += interval

        # track how far elongation proceeded to start from next iteration
        self.elongation = elongation.elongation - int(elongation.elongation)

        molecules = {
            self.unbound_ribosomes_key: unbound_ribosomes - original_unbound_ribosomes}

        molecules.update({
            key: count * -1
            for key, count in elongation.monomers.items()})

        update = {
            'ribosomes': {
                'ribosomes': [ribosome.to_dict() for ribosome in ribosomes]},
            'molecules': molecules,
            'proteins': elongation.complete_polymers}

        print('molecules update: {}'.format(molecules))

        return update


def test_translation():
    parameters = {
        'transcript_affinities': {
                'oA': 1.0,
                'oAZ': 1.0,
                'oB': 1.0,
                'oBY': 1.0},
        'elongation_rate': 10.0,
        'advancement_rate': 10.0}

    parameters = {}
    translation = Translation(parameters)

    states = {
        'ribosomes': {'ribosomes': []},
        'molecules': {translation.unbound_ribosomes_key: 10},
        'transcripts': {
            'oA': 10,
            'oAZ': 10,
            'oB': 10,
            'oBY': 10}}
    states['molecules'].update({
        molecule_id: 100
        for molecule_id in translation.monomer_ids})

    update = translation.next_update(10.0, states)
    
    print(update)
    print('complete!')



if __name__ == '__main__':
    test_translation()
        
