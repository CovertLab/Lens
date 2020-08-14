from __future__ import absolute_import, division, print_function

import os
import sys
import random
import math
import uuid
import argparse

import numpy as np
from numpy import linspace

from vivarium.library.dict_utils import deep_merge
from vivarium.core.process import Process
from vivarium.core.composition import (
    simulate_process_in_experiment,
    PROCESS_OUT_DIR
)
from vivarium.plots.flagella_activity import plot_activity

NAME = 'flagella_activity'



class FlagellaActivity(Process):
    ''' Model multi-flagellar motor activity and CheY-P fluctuations.

        :term:`Ports`:

        * **internal**: Expects a :term:`store` with 'chemoreceptor_activity',
            'motile_state', 'CheY', 'CheY_P', 'cw_bias'
        * **internal_counts**: Expects a :term:`store` with an int for 'flagella'
        * **flagella**: Expects a :term:`store` with a subschema for flagellar
          sub-compartments. Each one is assigned a uuid.
        * **membrane**: Expects a :term:`store` with 'PMF', 'protons_flux_accumulated'
        * **boundary**: Expects a :term:`store` with 'thrust', 'torque'

    References:
        * Mears, P. J., Koirala, S., Rao, C. V., Golding, I., & Chemla, Y. R. (2014).
        Escherichia coli swimming is robust against variations in flagellar number.

    TODO:
        - Mears et al. has flagella with 3 conformational states for flagella (normal (CCW), semi (CW), curly (CW)).
        - Flagella subcompartments should be nested dictionaries with multiple states (rotational state, # of motors engaged).
    '''

    name = NAME
    initial_pmf = 170 # PMF ~170mV at pH 7, ~140mV at pH 7.7 (Berg H, E. coli in motion, 2004, pg 113)
    defaults = {
        'n_flagella': 5,
        'initial_state': {
            'CheY': 2.59,
            'CheY_P': 2.59,  # (uM) mean concentration of CheY-P
            'cw_bias': 0.5,
            'motile_state': 0,  # 1 for tumble, -1 for run, 0 for none
            'PMF': initial_pmf,
        },
        # parameters from Mears, Koirala, Rao, Golding, Chemla (2014)
        'ccw_to_cw': 0.26,  # (1/s) Motor switching rate from CCW->CW
        'cw_to_ccw': 1.7,  # (1/s) Motor switching rate from CW->CCW
        'CB': 0.13,  # average CW bias of wild-type motors
        'lambda': 0.68,  # (1/s) transition rate from semi-coiled to curly-w state

        # CheY-P flucutations
        'YP_ss': 2.59,  # (uM) steady state concentration of CheY-P
        'sigma2_Y': 1.0,  # (uM^2) variance in CheY-P
        'tau': 0.2,  # (s) characteristic time-scale fluctuations in [CheY-P]

        # CW bias
        'K_d': 3.1,  # (uM) midpoint of CW bias vs CheY-P response curve
        'H': 10.3,  # Hill coefficient for CW bias vs CheY-P response curve

        # rotational state of individual flagella
        # parameters from Sneddon, Pontius, and Emonet (2012)
        'omega': 1.3,  # (1/s) characteristic motor switch time
        'g_0': 40,  # (k_B/T) free energy barrier for CCW-->CW
        'g_1': 40,  # (k_B/T) free energy barrier for CW-->CCW
        'K_D': 3.06,  # binding constant of Chey-P to base of the motor

        # motile force parameters
        'flagellum_thrust': 25,  # (pN) (Berg H, E. coli in motion, 2004, pg 113)
        'tumble_jitter': 120.0,
        'tumble_scaling': 1 / initial_pmf,
        'run_scaling': 1 / initial_pmf,
        'time_step': 0.01,  # 0.001
    }

    def __init__(self, parameters=None):
        super(FlagellaActivity, self).__init__(parameters)

    def ports_schema(self):
        ports = [
            'internal',
            'boundary',
            'internal_counts',
            'membrane',
            'flagella']
        schema = {port: {} for port in ports}

        # boundary
        for state in ['thrust', 'torque']:
            schema['boundary'][state] = {
                '_default': 0.0,
                '_emit': True,
                '_updater': 'set'}

        # membrane
        for state in ['PMF', 'protons_flux_accumulated']:
            schema['membrane'][state] = {'_default': self.parameters['initial_state'].get(state, 0.0)}

        # internal
        schema['internal']['chemoreceptor_activity'] = {}
        for state in ['motile_state', 'CheY', 'CheY_P', 'cw_bias']:
            schema['internal'][state] = {
                '_default': self.parameters['initial_state'].get(state, 0),
                '_emit': True,
                '_updater': 'set'}

        # internal_counts
        schema['internal_counts']['flagella'] = {
            '_value': self.parameters['n_flagella'],
            '_default': self.parameters['n_flagella'],
            '_emit': True}

        # flagella
        schema['flagella'] = {
            '_divider': 'split_dict',
            '*': {
                '_default': 1,
                '_updater': 'set',
                '_emit': True}}

        return schema

    def next_update(self, timestep, states):

        internal = states['internal']
        n_flagella = states['internal_counts']['flagella']
        flagella = states['flagella']
        PMF = states['membrane']['PMF']

        # states
        CheY = internal['CheY']
        CheY_P = internal['CheY_P']

        # parameters
        tau = self.parameters['tau']
        YP_ss = self.parameters['YP_ss']
        sigma = self.parameters['sigma2_Y']**0.5
        K_d = self.parameters['K_d']
        H = self.parameters['H']

        # update number of flagella
        flagella_update = {}
        new_flagella = int(n_flagella) - len(flagella)
        if new_flagella < 0:
            # remove flagella
            flagella_update['_delete'] = []
            remove = random.sample(list(flagella.keys()), abs(new_flagella))
            for flagella_id in remove:
                flagella_update['_delete'].append((flagella_id,))

        elif new_flagella > 0:
            # add flagella
            flagella_update['_add'] = []
            for index in range(new_flagella):
                flagella_id = str(uuid.uuid1())
                flagella_update['_add'].append({
                    'path': (flagella_id, ),
                    'state': random.choice([-1, 1])})

        ## update CheY-P
        # Mears et al. eqn S6
        dYP = -(1 / tau) * (CheY_P - YP_ss) * timestep + sigma * (2 * timestep / tau)**0.5 * random.normalvariate(0, 1)
        CheY_P = max(CheY_P + dYP, 0.0)  # keep value positive
        CheY = max(CheY - dYP, 0.0)  # keep value positive

        ## CW bias, Hill function
        # Cluzel, P., Surette, M., & Leibler, S. (2000).
        cw_bias = CheY_P**H / (K_d**H + CheY_P**H)

        ## update all flagella
        for flagella_id, motor_state in flagella.items():
            new_motor_state = self.update_flagellum(motor_state, cw_bias, CheY_P, timestep)
            flagella_update.update({flagella_id: new_motor_state})

        ## get cell motile state.
        # if any flagella is rotating CW, the cell tumbles.
        # flagella motor state: -1 for CCW, 1 for CW
        # motile state: -1 for run, 1 for tumble, 0 for no state
        if any(state == 1 for state in flagella_update.values()):
            motile_state = 1
            [thrust, torque] = self.tumble(n_flagella, PMF)
        elif len(flagella_update) > 0:
            motile_state = -1
            [thrust, torque] = self.run(n_flagella, PMF)
        else:
            motile_state = 0
            thrust = 0
            torque = 0

        return {
            'flagella': flagella_update,
            'internal': {
                'CheY': CheY,
                'CheY_P': CheY_P,
                'cw_bias': cw_bias,
                'motile_state': motile_state},
            'boundary': {
                'thrust': thrust,
                'torque': torque,
            }}

    def update_flagellum(self, motor_state, cw_bias, CheY_P, timestep):
        '''
        Rotational state of an individual flagellum from:
            Sneddon, M. W., Pontius, W., & Emonet, T. (2012).
            Stochastic coordination of multiple actuators reduces
            latency and improves chemotactic response in bacteria.

        # TODO -- normal, semi, curly states from Sneddon
        '''
        g_0 = self.parameters['g_0']  # (k_B/T) free energy barrier for CCW-->CW
        g_1 = self.parameters['g_1']  # (k_B/T) free energy barrier for CW-->CCW
        K_D = self.parameters['K_D']  # binding constant of Chey-P to base of the motor
        omega = self.parameters['omega']  # (1/s) characteristic motor switch time

        # free energy barrier
        delta_g = g_0 / 4 - g_1 / 2 * (CheY_P / (CheY_P + K_D))

        # switching frequency
        CW_to_CCW = omega * math.exp(delta_g)
        CCW_to_CW = omega * math.exp(-delta_g)
        # switch_freq = CCW_to_CW * (1 - cw_bias) + CW_to_CCW * cw_bias

        # flagella motor state: -1 for CCW, 1 for CW
        if motor_state == -1:
            prob_switch = CCW_to_CW * timestep
            if np.random.random(1)[0] <= prob_switch:
                new_motor_state = 1
            else:
                new_motor_state = -1

        elif motor_state == 1:
            prob_switch = CW_to_CCW * timestep
            if np.random.random(1)[0] <= prob_switch:
                new_motor_state = -1
            else:
                new_motor_state = 1

        return new_motor_state

    def tumble(self, n_flagella, PMF):
        thrust = self.parameters['tumble_scaling'] * PMF * self.parameters['flagellum_thrust'] * n_flagella
        torque = random.normalvariate(0, self.parameters['tumble_jitter'])
        return [thrust, torque]

    def run(self, n_flagella, PMF):
        thrust = self.parameters['run_scaling'] * PMF * self.parameters['flagellum_thrust'] * n_flagella
        torque = 0.0
        return [thrust, torque]


# testing functions
default_params = {'n_flagella': 5}
default_timeline = [(10, {})]
def test_activity(parameters=default_params, timeline=default_timeline):
    motor = FlagellaActivity(parameters)
    settings = {
        'timeline': {'timeline': timeline},
        'return_raw_data': True
    }
    return simulate_process_in_experiment(motor, settings)

def test_motor_PMF():

    # range of PMF value for test
    PMF_values = linspace(50.0, 200.0, 501).tolist()
    timestep = 1

    # initialize process and state
    motor = FlagellaActivity()
    state = motor.default_state()

    motor_state_vec = []
    thrust_vec = []
    torque_vec = []
    for PMF in PMF_values:
        state['membrane']['PMF'] = PMF
        update = motor.next_update(timestep, state)

        motile_state = update['internal']['motile_state']
        thrust = update['boundary']['thrust']
        torque = update['boundary']['torque']

        # save
        motor_state_vec.append(motile_state)
        thrust_vec.append(thrust)
        torque_vec.append(torque)

    return {
        'motile_state': motor_state_vec,
        'thrust': thrust_vec,
        'torque': torque_vec,
        'PMF': PMF_values,
    }


def run_variable_flagella(out_dir):
    # variable flagella
    init_params = {'n_flagella': 5}
    timeline = [
        (0, {}),
        (20, {('internal_counts', 'flagella'): 6}),
        (40, {('internal_counts', 'flagella'): 7}),
        (60, {})]
    output3 = test_activity(init_params, timeline)
    plot_activity(output3, out_dir, 'variable_flagella')


if __name__ == '__main__':
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    parser = argparse.ArgumentParser(description='flagella expression')
    parser.add_argument('--variable', '-v', action='store_true', default=False,)
    parser.add_argument('--zero', '-z', action='store_true', default=False, )
    args = parser.parse_args()
    no_args = (len(sys.argv) == 1)

    if args.variable:
        run_variable_flagella(out_dir)
    elif args.zero:
        zero_flagella = {'n_flagella': 0}
        timeline = [(10, {})]
        output1 = test_activity(zero_flagella, timeline)
        plot_activity(output1, out_dir, 'motor_control_zero_flagella')
    else:
        five_flagella = {'n_flagella': 5}
        timeline = [(60, {})]
        output2 = test_activity(five_flagella, timeline)
        plot_activity(output2, out_dir, 'motor_control')
