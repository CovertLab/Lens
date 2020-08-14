import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import Patch


def plot_signal_transduction(timeseries, out_dir='out', filename='signal_transduction'):
    ligand = timeseries['boundary']['external']
    chemoreceptor_activity = timeseries['internal']['chemoreceptor_activity']
    CheY_P = timeseries['internal']['CheY_P']
    cw_bias = timeseries['internal']['cw_bias']
    motile_state = timeseries['internal']['motile_state']
    time_vec = timeseries['time']

    # grid for cell state
    motile_state_grid = np.zeros((1, len(time_vec)))
    motile_state_grid[0, :] = motile_state

    # set up colormaps
    # cell motile state
    cmap1 = colors.ListedColormap(['steelblue', 'lightgray', 'darkorange'])
    bounds1 = [-1, -1/3, 1/3, 1]
    norm1 = colors.BoundaryNorm(bounds1, cmap1.N)
    motile_legend_elements = [
        Patch(facecolor='steelblue', edgecolor='k', label='Run'),
        Patch(facecolor='darkorange', edgecolor='k', label='Tumble'),
        Patch(facecolor='lightgray', edgecolor='k', label='N/A')]

    # plot results
    cols = 1
    rows = 5
    plt.figure(figsize=(3.0 * cols, 1.5 * rows))
    plt.rc('font', size=12)

    ax1 = plt.subplot(rows, cols, 1)
    ax2 = plt.subplot(rows, cols, 2)
    ax3 = plt.subplot(rows, cols, 3)
    ax4 = plt.subplot(rows, cols, 4)
    ax5 = plt.subplot(rows, cols, 5)

    for ligand_id, ligand_vec in ligand.items():
        ax1.plot(time_vec, ligand_vec, 'steelblue')
    ax2.plot(time_vec, chemoreceptor_activity, 'steelblue')
    ax3.plot(time_vec, CheY_P, 'steelblue')
    ax4.plot(time_vec, cw_bias, 'steelblue')
    # ax5.plot(time_vec, motile_state, 'steelblue')

    # plot cell motile state
    ax5.imshow(motile_state_grid,
               interpolation='nearest',
               aspect='auto',
               cmap=cmap1,
               norm=norm1,
               extent=[time_vec[0], time_vec[-1], 0, 1])

    ax1.set_xticklabels([])
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.tick_params(right=False, top=False)
    ax1.set_ylabel("external ligand \n (mM) ", fontsize=10)

    ax2.set_xticklabels([])
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.tick_params(right=False, top=False)
    ax2.set_ylabel("cluster activity \n P(on)", fontsize=10)

    ax3.set_xticklabels([])
    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.tick_params(right=False, top=False)
    ax3.set_ylabel("CheY-P", fontsize=10)

    ax4.set_xticklabels([])
    ax4.spines['right'].set_visible(False)
    ax4.spines['top'].set_visible(False)
    ax4.tick_params(right=False, top=False)
    ax4.set_ylabel("CW bias", fontsize=10)

    ax5.spines['right'].set_visible(False)
    ax5.spines['top'].set_visible(False)
    ax5.tick_params(right=False, top=False)
    ax5.set_xlabel("time (s)", fontsize=12)
    ax5.set_ylabel("motile state", fontsize=10)

    # legend
    ax5.legend(
        title='motile state',
        handles=motile_legend_elements,
        loc='center left',
        bbox_to_anchor=(1, 0.5))

    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.7, hspace=0.3)
    plt.savefig(fig_path + '.png', bbox_inches='tight')
