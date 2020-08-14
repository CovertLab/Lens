from __future__ import absolute_import, division, print_function

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import Patch

from vivarium.core.emitter import timeseries_from_data



def plot_activity(data, out_dir='out', filename='motor_control'):
    timeseries = timeseries_from_data(data)

    CheY_vec = timeseries['internal']['CheY']
    CheY_P_vec = timeseries['internal']['CheY_P']
    cw_bias_vec = timeseries['internal']['cw_bias']
    motile_state_vec = timeseries['internal']['motile_state']
    thrust_vec = timeseries['boundary']['thrust']
    flagella_activity = timeseries.get('flagella', {})
    time_vec = timeseries['time']

    # make flagella activity grid
    flagella_ids = list(flagella_activity.keys())
    flagella_indexes = {}
    next_index = 0
    activity_grid = np.zeros((len(flagella_ids), len(time_vec)))
    total_CW = np.zeros((len(time_vec)))
    for time_index, (time, time_data) in enumerate(data.items()):
        time_flagella = time_data.get('flagella', {})

        for flagella_id, rotation_states in time_flagella.items():

            # get flagella_index by order of appearance
            if flagella_id not in flagella_indexes:
                flagella_indexes[flagella_id] = next_index
                next_index += 1
            flagella_index = flagella_indexes[flagella_id]

            modified_rotation_state = 0
            CW_rotation_state = 0
            if rotation_states == -1:
                modified_rotation_state = 1
            elif rotation_states == 1:
                modified_rotation_state = 2
                CW_rotation_state = 1

            activity_grid[flagella_index, time_index] = modified_rotation_state
            total_CW += np.array(CW_rotation_state)

    # grid for cell state
    motile_state_grid = np.zeros((1, len(time_vec)))
    motile_state_grid[0, :] = motile_state_vec

    # set up colormaps
    # cell motile state
    cmap1 = colors.ListedColormap(['steelblue', 'lightgray', 'darkorange'])
    bounds1 = [-1, -1/3, 1/3, 1]
    norm1 = colors.BoundaryNorm(bounds1, cmap1.N)
    motile_legend_elements = [
        Patch(facecolor='steelblue', edgecolor='k', label='Run'),
        Patch(facecolor='darkorange', edgecolor='k', label='Tumble'),
        Patch(facecolor='lightgray', edgecolor='k', label='N/A')]

    # rotational state
    cmap2 = colors.ListedColormap(['lightgray', 'steelblue', 'darkorange'])
    bounds2 = [0, 0.5, 1.5, 2]
    norm2 = colors.BoundaryNorm(bounds2, cmap2.N)
    rotational_legend_elements = [
        Patch(facecolor='steelblue', edgecolor='k', label='CCW'),
        Patch(facecolor='darkorange', edgecolor='k', label='CW'),
        Patch(facecolor='lightgray', edgecolor='k', label='N/A')]

    # plot results
    cols = 1
    rows = 5
    plt.figure(figsize=(3.0 * cols, 1.5 * rows))
    plt.rc('font', size=12)

    # define subplots
    ax1 = plt.subplot(rows, cols, 1)
    ax2 = plt.subplot(rows, cols, 2)
    ax3 = plt.subplot(rows, cols, 3)
    ax4 = plt.subplot(rows, cols, 4)
    ax5 = plt.subplot(rows, cols, 5)

    # plot Che-P state
    ax1.plot(time_vec, CheY_vec, label='CheY')
    ax1.plot(time_vec, CheY_P_vec, label='CheY_P')
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax1.set_xticks([])
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.set_xlim(time_vec[0], time_vec[-1])
    ax1.set_ylabel('concentration \n (uM)')

    # plot CW bias
    ax2.plot(time_vec, cw_bias_vec)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xticks([])
    ax2.set_xlim(time_vec[0], time_vec[-1])
    ax2.set_ylabel('CW bias')

    # plot flagella states in a grid
    if len(activity_grid) > 0:
        ax3.imshow(activity_grid,
                   interpolation='nearest',
                   aspect='auto',
                   cmap=cmap2,
                   norm=norm2,
                   # extent=[-1,1,-1,1]
                   extent=[time_vec[0], time_vec[-1], len(flagella_ids)+0.5, 0.5]
                   )
        plt.locator_params(axis='y', nbins=len(flagella_ids))
        ax3.set_yticks(list(range(1, len(flagella_ids) + 1)))
        ax3.set_xticks([])
        ax3.set_ylabel('flagella #')
        ax3.tick_params(axis='both', labelsize=10)

        # legend
        ax3.legend(
            title='activity',
            handles=rotational_legend_elements,
            loc='center left',
            bbox_to_anchor=(1, 0.5))
    else:
        # no flagella
        ax3.set_axis_off()

    # plot cell motile state
    ax4.imshow(motile_state_grid,
               interpolation='nearest',
               aspect='auto',
               cmap=cmap1,
               norm=norm1,
               extent=[time_vec[0], time_vec[-1], 0, 1])
    ax4.set_yticks([])
    ax4.set_xticks([])

    # legend
    ax4.legend(
        title='motile state',
        handles=motile_legend_elements,
        loc='center left',
        bbox_to_anchor=(1, 0.5))

    # plot motor thrust
    ax5.plot(time_vec, thrust_vec)
    ax5.set_xlim(time_vec[0], time_vec[-1])
    ax5.spines['right'].set_visible(False)
    ax5.spines['top'].set_visible(False)
    ax5.set_ylabel('thrust (pN)')
    ax5.set_xlabel('time (sec)')

    # save figure
    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.7, hspace=0.3)
    plt.savefig(fig_path + '.png', bbox_inches='tight')


def plot_motor_PMF(output, out_dir='out', figname='motor_PMF'):
    motile_state = output['motile_state']
    thrust = output['thrust']
    motile_torque = output['motile_torque']
    PMF = output['PMF']

    # plot results
    cols = 1
    rows = 1
    plt.figure(figsize=(10 * cols, 2 * rows))

    # define subplots
    ax1 = plt.subplot(rows, cols, 1)

    # plot motile_state
    ax1.plot(PMF, thrust)
    ax1.set_xlabel('PMF (mV)')
    ax1.set_ylabel('force')

    # save figure
    fig_path = os.path.join(out_dir, figname)
    plt.subplots_adjust(wspace=0.7, hspace=0.3)
    plt.savefig(fig_path + '.png', bbox_inches='tight')
