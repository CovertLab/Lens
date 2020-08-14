'''
===============================================
Ordinary Differential Equation Expression Model
===============================================
'''

from __future__ import absolute_import, division, print_function

import os
import argparse
import random
import math
import logging as log

from vivarium.core.process import Process
from vivarium.library.dict_utils import deep_merge, tuplify_port_dicts
from vivarium.core.composition import (
    simulate_process_in_experiment,
    plot_simulation_output,
    PROCESS_OUT_DIR,
)
from vivarium.library.regulation_logic import build_rule
from vivarium.library.units import units
from vivarium.processes.derive_globals import AVOGADRO


NAME = 'ode_expression'

class ODE_expression(Process):
    '''Models gene expression using ordinary differential equations

     This :term:`process class` models the rates of transcription,
     translation, and degradation using ordinary differential
     equations (ODEs). These ODEs are as follows:

     * Transcription: :math:`\\frac{dM}{dt} = k_M - d_M M`

         * :math:`M`: Concentration of mRNA
         * :math:`k_M`: Trancription rate
         * :math:`d_M`: Degradation rate

     * Translation: :math:`\\frac{dP}{dt} = k_P m_P - d_P P`

         * :math:`P`: Concentration of protein
         * :math:`m_P`: Concentration of transcript
         * :math:`k_P`: Translation rate
         * :math:`d_P`: Degradation rate

     The transcription rate abstracts over factors that include gene
     copy number, RNA polymerase abundance, promoter strengths, and
     nucleotide availability. The translation rate similarly
     abstracts over factors that include ribosome availability, the
     binding strength of the mRNA to the ribosome, the availability
     of tRNAs, and the availability of free amino acids.

     This process class also models genetic regulation with boolean
     logic statements. This restricts us to modeling binary
     regulation: reactions can be completely suppressed, but they
     cannot be only slowed. For example, we can model a statement
     like:

         If :math:`[glucose] > 0.1`, completely inhibit LacY
         expression.

     But we cannot model a statement like this:

         If :math:`[glucose] > 0.1`, reduce LacY expression by 50%.

     .. note::
         This model does **not** include:

             * Kinetic regulation
             * Cooperativity
             * Autoinhibition
             * Autoactivation

     :term:`Ports`:

     * **internal**: Expects a :term:`store` with all of the
       transcripts and proteins being modeled. This store should also
       include any internal regulators.
     * **external**: Expects a store with all external regulators.
       These variables are not updated; they are only read.

     Arguments:
         initial_parameters: A dictionary of configuration options.
             The following configuration options may be provided:

             * **transcription_rates** (:py:class:`dict`): Maps
               transcript names (the keys of the dict) to
               transcription rates (values of the dict).
             * **translation_rates** (:py:class:`dict`): Maps protein
               names (the keys of the dict) to translation rates (the
               values of the dict).
             * **degradation_rates** (:py:class:`dict`): Maps from
               names of substrates (transcripts or proteins) to
               degrade to their degradation rates.
             * **protein_map** (:py:class:`dict`): Maps from protein
               name to transcript name for each transcript-protein
               pair.
             * **initial_state** (:py:class:`dict`): Maps from
               :term:`port` names to the initial state of that port.
               Each initial state is specified as a dictionary with
               :term:`variable` names as keys and variable values as
               values.
             * **regulators** (:py:class:`list`): List of the
               molecules that regulate any of the modeled
               transcription or translation reactions. Each molecule
               is a tuple of the port and the molecule's variable
               name.
             * **regulation** (:py:class:`dict`): Maps from reaction
               product (transcript or protein) name to a boolean
               regulation statement. For example, to only transcribe
               ``lacy_RNA`` if external ``glc__D_e`` concentration is
               at most 0.1, we would pass:

               .. code-block:: python

                 {'lacy_RNA': 'if not (external, glc__D_e) > 0.1'}

     Example instantiation of a :term:`process` modeling LacY
     expression:

     >>> initial_state = {
     ...     'internal': {
     ...         'lacy_RNA': 0.0,
     ...         'LacY': 0.0,
     ...     },
     ...     'external': {
     ...         'glc__D_e': 2.0,
     ...     }
     ... }
     >>> config = {
     ...     'transcription_rates': {
     ...         'lacy_RNA': 3.0,
     ...     },
     ...     'translation_rates': {
     ...         'LacY': 4.0,
     ...     },
     ...     'degradation_rates': {
     ...         'lacy_RNA': 1.0,
     ...         'LacY': 0.0,
     ...     },
     ...     'protein_map': {
     ...         'LacY': 'lacy_RNA',
     ...     },
     ...     'initial_state': initial_state,
     ...     'regulators': [('external', 'glc__D_e')],
     ...     'regulation': {
     ...         'lacy_RNA': 'if (external, glc__D_e) > 1.0',
     ...     },
     ... }
     >>> expression_process = ODE_expression(config)
     >>> state = expression_process.default_state()
     >>> state == initial_state
     True
     >>> # When external glc__D_e present, no transcription
     >>> expression_process.next_update(1, state)
     {'internal': {'lacy_RNA': 0.0, 'LacY': 0.0}}
     >>> # But translation and degradation still occur
     >>> state['internal']['lacy_RNA'] = 1.0
     >>> expression_process.next_update(1, state)
     {'internal': {'lacy_RNA': -1.0, 'LacY': 4.0}}
     >>> # When external glc__D_e low, LacY transcribed
     >>> state['external']['glc__D_e'] = 0.5
     >>> expression_process.next_update(1, state)
     {'internal': {'lacy_RNA': 2.0, 'LacY': 4.0}}

     '''

    name = NAME
    defaults = {
        'time_step': 1.0,
        'transcription_rates': {},
        'translation_rates': {},
        'degradation_rates': {},
        'protein_map': {},
        'transcription_leak': {
            'rate': 0.0,  # probability of leak in 1 second
            'magnitude': 0.0,  # the amount of leak
        },
        'regulation': {},
        'regulators': [],
        'initial_state': {},
        'counts_deriver_key': 'expression_counts',
    }

    def __init__(self, parameters=None):
        super(ODE_expression, self).__init__(parameters)

        # ode gene expression
        self.transcription = self.parameters.get('transcription_rates')
        self.translation = self.parameters.get('translation_rates')
        self.degradation = self.parameters.get('degradation_rates')
        self.protein_map = self.parameters.get('protein_map')
        transcription_leak = self.parameters.get('transcription_leak')
        self.transcription_leak_rate = transcription_leak['rate']
        self.transcription_leak_magnitude = transcription_leak['magnitude']

        # boolean regulation
        regulation_logic = self.parameters.get('regulation')
        self.regulation = {
            gene_id: build_rule(logic) for gene_id, logic in regulation_logic.items()}
        regulators = self.parameters.get('regulators')
        self.internal_regulators = [state_id for port_id, state_id in regulators if port_id == 'internal']
        self.external_regulators = [state_id for port_id, state_id in regulators if port_id == 'external']

        # get initial state
        states = list(self.transcription.keys()) + list(self.translation.keys())
        null_states = {'internal': {
            state_id: 0 for state_id in states}}
        initialized_states = self.parameters.get('initial_state')
        self.initial_state = deep_merge(null_states, initialized_states)
        self.internal = list(self.initial_state.get('internal', {}).keys())
        self.external = list(self.initial_state.get('external', {}).keys())
        self.counts_deriver_key = self.parameters.get('counts_deriver_key')


    def ports_schema(self):
        ports = [
            'internal',
            'external',
            'counts',
            'global',
        ]
        schema = {port: {} for port in ports}

        # internal
        for state in self.internal + self.internal_regulators:
            schema['internal'][state] = {
                '_default': self.initial_state['internal'].get(state, 0.0),
                '_emit': True,
            }

        # external
        for state in self.external + self.external_regulators:
            schema['external'][state] = {
                '_default': self.initial_state['external'].get(state, 0.0),
                '_emit': True,
            }

        # counts
        for state in self.internal + self.internal_regulators:
            schema['counts'][state] = {
                '_divider': 'split',
                '_emit': True
            }

        # global
        schema['global'] = {}

        return schema

    def derivers(self):
        return {
            self.counts_deriver_key: {
                'deriver': 'counts_deriver',
                'port_mapping': {
                    'global': 'global',
                    'concentrations': 'internal',
                    'counts': 'counts'},
                'config': {
                    'concentration_keys': self.internal + self.internal_regulators}}}

    def next_update(self, timestep, states):
        internal_state = states['internal']

        # get state of regulated reactions (True/False)
        flattened_states = tuplify_port_dicts(states)
        regulation_state = {}
        for gene_id, reg_logic in self.regulation.items():
            regulation_state[gene_id] = reg_logic(flattened_states)

        internal_update = {}
        # transcription: dM/dt = k_M - d_M * M
        # M: conc of mRNA, k_M: transcription rate, d_M: degradation rate
        for transcript, rate in self.transcription.items():
            transcript_state = internal_state[transcript]

            # do not transcribe inhibited genes, except for transcription leaks
            if transcript in regulation_state and regulation_state.get(transcript):
                # leak probability for probability as function of the time step
                rate = -math.log(1 - self.transcription_leak_rate)
                leak_probability = 1 - math.exp(-rate * timestep)
                if random.uniform(0, 1) < leak_probability:
                    rate = self.transcription_leak_magnitude
                else:
                    rate = 0.0

            internal_update[transcript] = \
                (rate - self.degradation.get(transcript, 0) * transcript_state) * timestep

        # translation: dP/dt = k_P * m_P - d_P * P
        # P: conc of protein, m_P: conc of P's transcript, k_P: translation rate, d_P: degradation rate
        for protein, rate in self.translation.items():
            transcript = self.protein_map[protein]
            transcript_state = internal_state[transcript]
            protein_state = internal_state[protein]

            internal_update[protein] = \
                (rate * transcript_state - self.degradation.get(protein, 0) * protein_state) * timestep

        return {'internal': internal_update}


# test functions
# toy config
def get_lacy_config():
    '''LacY ODE expresion configuration

    This configuration can be passed to :py:class:`ODE_expression` like
    this:

    >>> config = get_lacy_config()
    >>> expression_process = ODE_expression(config)
    '''
    transcription_rates = {
        'lacy_RNA': 1e-7}

    translation_rates = {
        'LacY': 5e-3}

    protein_map = {
        'LacY': 'lacy_RNA'}

    degradation_rates = {
        'lacy_RNA': 3e-3,  # a single RNA lasts about 5 minutes
        'LacY': 3e-5}

    # define regulation
    regulators = [
        ('external', 'glc__D_e'),
        ('internal', 'lcts_p')]
    regulation = {
        'lacy_RNA': 'if (external, glc__D_e) > 0.1 and (internal, lcts_p) < 0.01'} # inhibited in this condition
    transcription_leak = {
        'rate': 1e-4,
        'magnitude': 1e-6,
        }

    # initial state
    initial_state = {
        'internal': {
            'lacy_RNA': 0.0,
            'LacY': 0.0},
        'external': {
            'glc__D_e': 8.0,
            'lcts_e': 8.0}}

    return {
        'transcription_rates': transcription_rates,
        'translation_rates': translation_rates,
        'degradation_rates': degradation_rates,
        'protein_map': protein_map,
        'regulators': regulators,
        'regulation': regulation,
        'transcription_leak': transcription_leak,
        'initial_state': initial_state}

def get_flagella_expression():
    '''Flagella ODE expresion configuration

    This configuration can be passed to :py:class:`ODE_expression` like
    this:

    >>> config = get_flagella_expression()
    >>> expression_process = ODE_expression(config)
    '''
    transcription = {
        'flag_RNA': 1e-6}

    translation = {
        'flagella': 8e-5}

    degradation = {
        'flag_RNA': 2e-2}  # 1e-23}

    protein_map = {
        'flagella': 'flag_RNA'}

    # get initial concentrations from counts
    volume = 1.2 * units.fL
    mmol_to_counts = (AVOGADRO * volume).to('L/mmol')
    counts = {
        'flagella': 0,
        'flag_RNA': 0}
    concentrations = {}
    for state_id, count in counts.items():
        concentrations[state_id] = (count / mmol_to_counts).magnitude

    initial_state = {
        'counts': counts,
        'internal': concentrations}

    return {
        'transcription_rates': transcription,
        'translation_rates': translation,
        'degradation_rates': degradation,
        'protein_map': protein_map,
        'initial_state': initial_state}


def test_expression(config=get_lacy_config(), timeline=[(100, {})]):
    expression = ODE_expression(config)
    settings = {'timeline': {'timeline': timeline}}
    return simulate_process_in_experiment(expression, settings)


if __name__ == '__main__':
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    parser = argparse.ArgumentParser(description='ODE expression')
    parser.add_argument('--lacY', '-l', action='store_true', default=False)
    parser.add_argument('--flagella', '-f', action='store_true', default=False)
    args = parser.parse_args()

    if args.flagella:
        timeline = [(2520, {})]
        timeseries = test_expression(get_flagella_expression(), timeline)
        plot_simulation_output(timeseries, {}, out_dir, 'flagella_expression')
    else:
        total_time = 5000
        shift_time1 = int(total_time / 5)
        shift_time2 = int(3 * total_time / 5)
        timeline = [
            (0, {('external', 'glc__D_e'): 10}),
            (shift_time1, {('external', 'glc__D_e'): 0}),
            (shift_time2, {('external', 'glc__D_e'): 10}),
            (total_time, {})]
        timeseries = test_expression(get_lacy_config(), timeline)
        plot_simulation_output(timeseries, {}, out_dir, 'lacY_expression')
