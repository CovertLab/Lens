from __future__ import absolute_import, division, print_function

from vivarium.processes.derive_globals import get_default_global_state
from vivarium.core.process import Deriver
from vivarium.library.units import units



class DeriveCounts(Deriver):
    """
    Process for deriving counts from concentrations
    """
    name = 'counts_deriver'
    defaults = {
        'concentration_keys': [],
        'initial_state': get_default_global_state(),
    }

    def __init__(self, parameters=None):
        super(DeriveCounts, self).__init__(parameters)

    def ports_schema(self):
        return {
            'global': {
                'volume': {
                    '_default': self.parameters[
                        'initial_state']['global']['volume'].to('fL')},
                'mmol_to_counts': {
                    '_default': self.parameters[
                        'initial_state']['global']['mmol_to_counts'].to('L/mmol')}},
            'counts': {
                molecule: {
                    '_default': 0,
                    '_divider': 'split',
                    '_updater': 'set'}
                for molecule in self.parameters['concentration_keys']},
            'concentrations': {
                molecule: {
                    '_default': 0.0}
                for molecule in self.parameters['concentration_keys']}}

    def next_update(self, timestep, states):
        mmol_to_counts = states['global']['mmol_to_counts'].to('L/mmol').magnitude
        concentrations = states['concentrations']

        counts = {}
        for molecule, concentration in concentrations.items():
            counts[molecule] = int(concentration * mmol_to_counts)

        return {
            'counts': counts}

# register process by invoking upon import
DeriveCounts()
