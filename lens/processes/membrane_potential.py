from __future__ import absolute_import, division, print_function

import os
import numpy as np
import scipy.constants as constants

from lens.actor.process import Process, deep_merge
from lens.utils.units import units

# (mmol) http://book.bionumbers.org/what-are-the-concentrations-of-different-ions-in-cells/
# Schultz, Stanley G., and A. K. Solomon. "Cation Transport in Escherichia coli" (1961)
# TODO -- add Mg2+, Ca2+
DEFAULT_STATE = {
    'internal': {
        'K': 300,  # (mmol) 30-300
        'Na': 10,  # (mmol) 10
        'Cl': 10},  # (mmol) 10-200 media-dependent
    'external': {
        'K': 5,
        'Na': 145,
        'Cl': 110}  # (mmol)
    }

# TODO -- get references on these
DEFAULT_PARAMETERS = {
    'p_K': 1,  # unitless, relative membrane permeability of K
    'p_Na': 0.05,  # unitless, relative membrane permeability of Na
    'p_Cl': 0.05,  # unitless, relative membrane permeability of Cl
    }

PERMEABILITY_MAP = {
    'K': 'p_K',
    'Na': 'p_Na',
    'Cl': 'p_Cl'
    }

CHARGE_MAP = {
    'K': '+',
    'Na': '+',
    'Cl': '-'
    }

class NoChargeError(Exception):
    pass

class MembranePotential(Process):
    '''
    Need to add a boot method for this process to lens/environment/boot.py for it to run on its own
    '''
    def __init__(self, config={}):

        # set states
        self.initial_states = config.get('states', DEFAULT_STATE)

        # set parameters
        parameters = {
            'R': constants.gas_constant,  # (J * K^-1 * mol^-1) gas constant
            'F': constants.physical_constants['Faraday constant'][0], # (C * mol^-1) Faraday constant
            'k': constants.Boltzmann, # (J * K^-1) Boltzmann constant
            }
        parameters.update(config.get('parameters', DEFAULT_PARAMETERS))

        self.permeability = config.get('permeability', PERMEABILITY_MAP)
        self.charge = config.get('charge', CHARGE_MAP)

        # set roles
        roles = {
            'internal': ['c_in'],
            'membrane': ['PMF', 'd_V', 'd_pH'],  # proton motive force (PMF), electrical difference (d_V), pH difference (d_pH)
            'external': ['c_out'],
        }

        super(MembranePotential, self).__init__(roles, parameters)

    def default_state(self):
        config = {'external': {'T': 310.15}}
        default_state = deep_merge((self.initial_states), config)
        return default_state

    def default_emitter_keys(self):
        keys = {
            'membrane': ['d_V', 'd_pH', 'PMF'],
        }
        return keys

    def default_updaters(self):
        keys = {'membrane': {
            'd_V': 'set',
            'd_pH': 'set',
            'PMF': 'set',
        }}
        return keys

    def next_update(self, timestep, states):
        internal_state = states['internal']
        external_state = states['external']

        # parameters
        R = self.parameters['R']
        F = self.parameters['F']
        k = self.parameters['k']

        # state
        T = external_state['T']  # temperature
        # e = 1 # proton charge # TODO -- get proton charge from state

        # Membrane potential.
        numerator = 0
        denominator = 0
        for ion_id, p_ion_id in self.permeability.iteritems():
            charge = self.charge[ion_id]
            p_ion = self.parameters[p_ion_id]

            # ions states
            internal = internal_state[ion_id]
            external = external_state[ion_id]

            if charge is '+':
                numerator += p_ion * external
                denominator += p_ion * internal
            elif charge is '-':
                numerator += p_ion * internal
                denominator += p_ion * external
            else:
                raise NoChargeError(
                    "No charge given for {}".format(ion_id))

        # Goldman equation for membrane potential
        # expected d_V = -120 mV
        d_V = (R * T) / (F) * np.log(numerator / denominator) * 1e3  # (mV). 1e3 factor converts from V

        # Nernst equation for pH difference
        # -2.3 * k * T / e  # -2.3 Boltzmann constant * temperature
        # expected d_pH = -50 mV
        d_pH = -50  # (mV) for cells grown at pH 7. (Berg, H. "E. coli in motion", pg 105)

        # proton motive force
        PMF = d_V + d_pH

        update = {
            'membrane': {
                'd_V': d_V,
                'd_pH': d_pH,
                'PMF': PMF,
            }}
        return update

def test_mem_potential():

    initial_parameters = {
        'states': DEFAULT_STATE,
        'parameters': DEFAULT_PARAMETERS,
        'permeability': PERMEABILITY_MAP,
        'charge': CHARGE_MAP,
    }

    # configure process
    mp = MembranePotential(initial_parameters)

    # get initial state and parameters
    state = mp.default_state()
    saved_state = {'internal': {}, 'external': {}, 'membrane': {}, 'time': []}

    ## Simulation
    timeline = [
        (0, {'external': {
            'Na': 1}
        }),
        (100, {'external': {
            'Na': 2}
        }),
        (500, {}),
    ]

    time = 0
    timestep = 1  # sec
    while time < timeline[-1][0]:
        time += timestep
        for (t, change_dict) in timeline:
            if time >= t:
                for key, change in change_dict.iteritems():
                    state[key].update(change)

        update = mp.next_update(timestep, state)
        saved_state['time'].append(time)

        # update external state
        for role in ['internal', 'external']:
            for state_id, value in state[role].iteritems():
                if state_id in saved_state[role].keys():
                    saved_state[role][state_id].append(value)
                else:
                    saved_state[role][state_id] = [value]

        # update membrane state from update
        for state_id, value in update['membrane'].iteritems():
            if state_id in saved_state['membrane'].keys():
                saved_state['membrane'][state_id].append(value)
            else:
                saved_state['membrane'][state_id] = [value]

    return saved_state

def plot_mem_potential(saved_state, out_dir='out'):
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt

    data_keys = [key for key in saved_state.keys() if key is not 'time']
    time_vec = [float(t) / 3600 for t in saved_state['time']]  # convert to hours

    # make figure, with grid for subplots
    n_data = [len(saved_state[key].keys()) for key in data_keys]
    n_rows = sum(n_data)
    fig = plt.figure(figsize=(8, n_rows * 2.5))
    grid = plt.GridSpec(n_rows + 1, 1, wspace=0.4, hspace=1.5)

    # plot data
    plot_idx = 0
    for key in data_keys:
        for mol_id, series in sorted(saved_state[key].iteritems()):
            ax = fig.add_subplot(grid[plot_idx, 0])  # grid is (row, column)

            ax.plot(time_vec, series)
            ax.title.set_text(str(key) + ': ' + mol_id)
            ax.set_xlabel('time (hrs)')

            # if key is 'internal':
            #     ax.set_yticks([0.0, 1.0])
            #     ax.set_yticklabels(["False", "True"])

            plot_idx += 1

    # save figure
    fig_path = os.path.join(out_dir, 'membrane_potential')
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.savefig(fig_path + '.pdf', bbox_inches='tight')


if __name__ == '__main__':
    saved_state = test_mem_potential()
    out_dir = os.path.join('out', 'tests', 'membrane_potential')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    plot_mem_potential(saved_state, out_dir)
