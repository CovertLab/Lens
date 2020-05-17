from __future__ import absolute_import, division, print_function

from vivarium.processes.derive_globals import AVOGADRO
from vivarium.compartment.process import Deriver
from vivarium.utils.units import units


def get_default_state():
    mass = 1339 * units.fg  # wet mass in fg
    density = 1100 * units.g / units.L
    volume = mass / density
    mmol_to_counts = (AVOGADRO * volume).to('L/mmol')

    return {
        'global': {
            'volume': volume.magnitude,
            'mmol_to_counts': mmol_to_counts.magnitude}}


class DeriveCounts(Deriver):
    """
    Process for deriving counts from concentrations
    """
    defaults = {
        'concentration_keys': []}

    def __init__(self, initial_parameters={}):

        self.initial_state = initial_parameters.get('initial_state', get_default_state())

        self.concentration_keys = self.or_default(
            initial_parameters, 'concentration_keys')

        ports = {
            'global': ['volume', 'mmol_to_counts'],
            'concentrations': self.concentration_keys,
            'counts': self.concentration_keys}

        parameters = {}
        parameters.update(initial_parameters)

        super(DeriveCounts, self).__init__(ports, parameters)

    def ports_schema(self):
        return {
            'global': {
                'volume': {
                    '_default': 0.0},
                'mmol_to_counts': {
                    '_default': 0.0}},
            'concentrations': {
                concentration: {
                    '_default': 0.0}
                for concentration in self.concentration_keys},
            'counts': {
                concentration: {
                    '_default': 0,
                    '_updater': 'set'}
                for concentration in self.concentration_keys}}

    def default_settings(self):

        # default emitter keys
        default_emitter_keys = {}

        # schema
        schema = {
            'counts': {
                state_id : {
                    'updater': 'set',
                    'divide': 'split'}
                for state_id in self.ports['counts']}}

        default_settings = {
            'state': self.initial_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema}

        return default_settings

    def next_update(self, timestep, states):
        mmol_to_counts = states['global']['mmol_to_counts']
        concentrations = states['concentrations']

        counts = {}
        for molecule, concentration in concentrations.items():
            counts[molecule] = int(concentration * mmol_to_counts)

        # concentrations = {port: state for port, state in states.items() if port not in ['counts', 'global']}

        # counts = {}
        # for port, states in concentrations.items():
        #     for state_id, conc in states.items():
        #         counts[state_id] = int(conc * mmol_to_counts)

        return {
            'counts': counts}
