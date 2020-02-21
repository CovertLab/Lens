from __future__ import absolute_import, division, print_function

from vivarium.processes.derive_globals import AVOGADRO
from vivarium.actor.process import Process
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


class DeriveCounts(Process):
    """
    Process for deriving counts from concentrations
    """
    def __init__(self, initial_parameters={}):

        self.initial_state = initial_parameters.get('initial_state', get_default_state())

        roles = initial_parameters.get('roles')
        roles.update({
            'global': ['volume', 'mmol_to_counts']})

        parameters = {}
        parameters.update(initial_parameters)

        super(DeriveCounts, self).__init__(roles, parameters)

    def default_settings(self):

        # default emitter keys
        default_emitter_keys = {}

        # schema
        schema = {
            'counts': {
                state_id : {
                    'updater': 'set',
                    'divide': 'split'}
                for state_id in self.roles['counts']}}

        default_settings = {
            'state': self.initial_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema}

        return default_settings

    def next_update(self, timestep, states):
        mmol_to_counts = states['global']['mmol_to_counts']
        concentrations = {role: state for role, state in states.items() if role not in ['counts', 'global']}

        counts = {}
        for role, states in concentrations.items():
            for state_id, conc in states.items():
                counts[state_id] = int(conc * mmol_to_counts)

        return {
            'counts': counts}