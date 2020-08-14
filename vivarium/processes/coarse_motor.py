'''
==========================
Vladimirov Motor Processes
==========================
'''

from __future__ import absolute_import, division, print_function

import os
import random
import copy
import math

import numpy as np
from numpy import linspace
import matplotlib.pyplot as plt

from vivarium.core.process import Process
from vivarium.core.composition import (
    simulate_process_in_experiment,
    PROCESS_OUT_DIR,
)
from vivarium.plots.coarse_motor import (
    plot_variable_receptor,
    plot_motor_control,
)

NAME = 'coarse_motor'

class MotorActivity(Process):
    ''' Model of motor activity

    Based on the model described in:
        Vladimirov, N., Lovdok, L., Lebiedz, D., & Sourjik, V. (2008).
        Dependence of bacterial chemotaxis on gradient shape and adaptation rate.

    CheY phosphorylation model from:
        Kollmann, M., Lovdok, L., Bartholome, K., Timmer, J., & Sourjik, V. (2005).
        Design principles of a bacterial signalling network. Nature.
    Motor switching model from:
        Scharf, B. E., Fahrner, K. A., Turner, L., and Berg, H. C. (1998).
        Control of direction of flagellar rotation in bacterial chemotaxis. PNAS.

    An increase of attractant inhibits CheA activity (chemoreceptor_activity),
    but subsequent methylation returns CheA activity to its original level.

    TODO -- add CheB phosphorylation
    '''

    name = NAME
    defaults = {
        'time_step': 0.1,
        # Vladimirov parameters
        # 'k_A': 5.0,  #
        'k_y': 100.0,  # 1/uM/s
        'k_z': 30.0,  # / CheZ,
        'gamma_Y': 0.1,  # rate constant
        'k_s': 0.45,  # scaling coefficient
        'adapt_precision': 3,  # scales CheY_P to cluster activity
        # motor
        'mb_0': 0.65,  # steady state motor bias (Cluzel et al 2000)
        'n_motors': 5,
        'cw_to_ccw': 0.83,  # 1/s (Block1983) motor bias, assumed to be constant
        'cw_to_ccw_leak': 0.25,  # rate of spontaneous transition to tumble
        # parameters for multibody physics
        'tumble_jitter': 120.0,
        # initial state
        'initial_state': {
            'internal': {
                # response regulator proteins
                'CheY_tot': 9.7,  # (uM) #0.0097,  # (mM) 9.7 uM = 0.0097 mM
                'CheY_P': 0.5,
                'CheZ': 0.01*100,  # (uM) #phosphatase 100 uM = 0.1 mM (0.01 scaling from RapidCell1.4.2)
                'CheA': 0.01*100,  # (uM) #100 uM = 0.1 mM (0.01 scaling from RapidCell1.4.2)
                # sensor activity
                'chemoreceptor_activity': 1/3,
                # motor activity
                'ccw_motor_bias': 0.5,
                'ccw_to_cw': 0.5,
                'motor_state': 1,  # motor_state 1 for tumble, 0 for run
            },
            'external': {
                'thrust': 0,
                'torque': 0,
            }}
    }

    def __init__(self, parameters=None):
        super(MotorActivity, self).__init__(parameters)

    def ports_schema(self):
        ports = ['internal', 'external']
        schema = {port: {} for port in ports}

        # external
        for state, default in self.parameters['initial_state']['external'].items():
            schema['external'][state] = {
                '_default': default,
                '_emit': True,
                '_updater': 'set'}

        # internal
        set_and_emit = [
                'ccw_motor_bias',
                'ccw_to_cw',
                'motor_state',
                'CheA',
                'CheY_P']
        for state, default in self.parameters['initial_state']['internal'].items():
            schema['internal'][state] = {'_default': default}
            if state in set_and_emit:
                schema['internal'][state].update({
                    '_emit': True,
                    '_updater': 'set'})

        return schema

    def next_update(self, timestep, states):
        internal = states['internal']
        P_on = internal['chemoreceptor_activity']
        motor_state_current = internal['motor_state']

        # parameters
        adapt_precision = self.parameters['adapt_precision']
        k_y = self.parameters['k_y']
        k_s = self.parameters['k_s']
        k_z = self.parameters['k_z']
        gamma_Y = self.parameters['gamma_Y']
        mb_0 = self.parameters['mb_0']
        cw_to_ccw = self.parameters['cw_to_ccw']

        ## Kinase activity
        # relative steady-state concentration of phosphorylated CheY.
        CheY_P = adapt_precision * k_y * k_s * P_on / (k_y * k_s * P_on + k_z + gamma_Y)  # CheZ cancels out of k_z

        ## Motor switching
        # CCW corresponds to run. CW corresponds to tumble
        ccw_motor_bias = mb_0 / (CheY_P * (1 - mb_0) + mb_0)  # (1/s)
        ccw_to_cw = cw_to_ccw * (1 / ccw_motor_bias - 1)  # (1/s)
        # don't let ccw_to_cw get under leak value
        if ccw_to_cw < self.parameters['cw_to_ccw_leak']:
            ccw_to_cw = self.parameters['cw_to_ccw_leak']

        if motor_state_current == 0:  # 0 for run
            # switch to tumble (cw)?
            rate = -math.log(1 - ccw_to_cw)  # rate for probability function of time
            prob_switch = 1 - math.exp(-rate * timestep)
            if np.random.random(1)[0] <= prob_switch:
                motor_state = 1
                thrust, torque = tumble(self.parameters['tumble_jitter'])
            else:
                motor_state = 0
                thrust, torque = run()

        elif motor_state_current == 1:  # 1 for tumble
            # switch to run (ccw)?
            rate = -math.log(1 - cw_to_ccw)  # rate for probability function of time
            prob_switch = 1 - math.exp(-rate * timestep)
            if np.random.random(1)[0] <= prob_switch:
                motor_state = 0
                [thrust, torque] = run()
            else:
                motor_state = 1
                [thrust, torque] = tumble()

        return {
            'internal': {
                'ccw_motor_bias': ccw_motor_bias,
                'ccw_to_cw': ccw_to_cw,
                'motor_state': motor_state,
                'CheY_P': CheY_P},
            'external': {
                'thrust': thrust,
                'torque': torque
            }}

def tumble(tumble_jitter=120.0):
    thrust = 100  # pN
    # average = 160
    # sigma = 10
    # torque = random.choice([-1, 1]) * random.normalvariate(average, sigma)
    torque = random.normalvariate(0, tumble_jitter)
    return [thrust, torque]

def run():
    # average thrust = 200 pN according to:
    # Berg, Howard C. E. coli in Motion. Under "Torque-Speed Dependence"
    thrust = 250  # pN
    torque = 0.0
    return [thrust, torque]


def test_motor_control(total_time=10):
    motor = MotorActivity({})
    experiment_settings = {
        'total_time': total_time,
        'timestep': 0.01}
    return simulate_process_in_experiment(motor, experiment_settings)

def test_variable_receptor():
    motor = MotorActivity()
    state = motor.default_state()
    timestep = 1
    chemoreceptor_activity = linspace(0.0, 1.0, 501).tolist()
    CheY_P_vec = []
    ccw_motor_bias_vec = []
    ccw_to_cw_vec = []
    motor_state_vec = []
    for activity in chemoreceptor_activity:
        state['internal']['chemoreceptor_activity'] = activity
        update = motor.next_update(timestep, state)
        CheY_P = update['internal']['CheY_P']
        ccw_motor_bias = update['internal']['ccw_motor_bias']
        ccw_to_cw = update['internal']['ccw_to_cw']
        motile_state = update['internal']['motor_state']

        CheY_P_vec.append(CheY_P)
        ccw_motor_bias_vec.append(ccw_motor_bias)
        ccw_to_cw_vec.append(ccw_to_cw)
        motor_state_vec.append(motile_state)

    # check ccw_to_cw bias is strictly increasing with increased receptor activity
    assert all(i <= j for i, j in zip(ccw_to_cw_vec, ccw_to_cw_vec[1:]))

    return {
        'chemoreceptor_activity': chemoreceptor_activity,
        'CheY_P': CheY_P_vec,
        'ccw_motor_bias': ccw_motor_bias_vec,
        'ccw_to_cw': ccw_to_cw_vec,
        'motor_state': motor_state_vec}


if __name__ == '__main__':
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    output1 = test_motor_control(200)
    plot_motor_control(output1, out_dir)

    output2 = test_variable_receptor()
    plot_variable_receptor(output2, out_dir)
