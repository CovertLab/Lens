from __future__ import absolute_import, division, print_function

import os
import numpy as np
import matplotlib.pyplot as plt


def plot_variable_receptor(output, out_dir='out', filename='motor_variable_receptor'):
    receptor_activities = output['chemoreceptor_activity']
    CheY_P_vec = output['CheY_P']
    ccw_motor_bias_vec = output['ccw_motor_bias']
    ccw_to_cw_vec = output['ccw_to_cw']

    # plot results
    cols = 1
    rows = 2
    plt.figure(figsize=(5 * cols, 2 * rows))

    ax1 = plt.subplot(rows, cols, 1)
    ax2 = plt.subplot(rows, cols, 2)

    ax1.scatter(receptor_activities, CheY_P_vec, c='b')
    ax2.scatter(receptor_activities, ccw_motor_bias_vec, c='b', label='ccw_motor_bias')
    ax2.scatter(receptor_activities, ccw_to_cw_vec, c='g', label='ccw_to_cw')

    ax1.set_xticklabels([])
    ax1.set_ylabel("CheY_P", fontsize=10)
    ax2.set_xlabel("receptor activity \n P(on) ", fontsize=10)
    ax2.set_ylabel("motor bias", fontsize=10)
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.7, hspace=0.1)
    plt.savefig(fig_path + '.png', bbox_inches='tight')


def plot_motor_control(timeseries, out_dir='out'):
    # TODO -- make this into an analysis figure
    expected_run = 0.42  # s (Berg) expected run length without chemotaxis
    expected_tumble = 0.14  # s (Berg)

    # receptor_activities = output['receptor_activities']
    CheY_P_vec = timeseries['internal']['CheY_P']
    ccw_motor_bias_vec = timeseries['internal']['ccw_motor_bias']
    ccw_to_cw_vec = timeseries['internal']['ccw_to_cw']
    motor_state_vec = timeseries['internal']['motor_state']
    time_vec = timeseries['time']

    # plot results
    cols = 1
    rows = 4
    plt.figure(figsize=(6 * cols, 1 * rows))

    ax1 = plt.subplot(rows, cols, 1)
    ax2 = plt.subplot(rows, cols, 2)
    ax3 = plt.subplot(rows, cols, 3)
    ax4 = plt.subplot(rows, cols, 4)

    ax1.plot(CheY_P_vec, 'b')
    ax2.plot(ccw_motor_bias_vec, 'b', label='ccw_motor_bias')
    ax2.plot(ccw_to_cw_vec, 'g', label='ccw_to_cw')

    # get length of runs, tumbles
    run_lengths = []
    tumble_lengths = []
    prior_state = 0
    state_start_time = 0
    for state, time in zip(motor_state_vec,time_vec):
        if state == 0:  # run
            if prior_state != 0:
                tumble_lengths.append(time - state_start_time)
                state_start_time = time
        elif state == 1:  # tumble
            if prior_state != 1:
                run_lengths.append(time - state_start_time)
                state_start_time = time
        prior_state = state

    avg_run_lengths = sum(run_lengths) / len(run_lengths)
    avg_tumble_lengths = sum(tumble_lengths) / len(tumble_lengths)

    # plot run distributions
    max_length = max(run_lengths + [1])
    bins = np.linspace(0, max_length, 30)
    ax3.hist([run_lengths], bins=bins, label=['run_lengths'], color=['b'])
    ax3.axvline(x=avg_run_lengths, color='k', linestyle='dashed', label='mean run')
    ax3.axvline(x=expected_run, color='r', linestyle='dashed', label='expected run')

    # plot tumble distributions
    max_length = max(tumble_lengths + [1])
    bins = np.linspace(0, max_length, 30)
    ax4.hist([tumble_lengths], bins=bins, label=['tumble_lengths'], color=['m'])
    ax4.axvline(x=avg_tumble_lengths, color='k', linestyle='dashed', label='mean tumble')
    ax4.axvline(x=expected_tumble, color='r', linestyle='dashed', label='expected tumble')

    # labels
    ax1.set_xticklabels([])
    ax1.set_ylabel("CheY_P", fontsize=10)

    ax2.set_xticklabels([])
    ax2.set_ylabel("motor bias", fontsize=10)
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax3.set_xlabel("motor state length (sec)")
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax4.set_xlabel("motor state length (sec)")
    ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # save the figure
    fig_path = os.path.join(out_dir, 'motor_control')
    plt.subplots_adjust(wspace=0.7, hspace=0.5)
    plt.savefig(fig_path + '.png', bbox_inches='tight')