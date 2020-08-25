'''
=====================
Colony Plotting Tools
=====================

This module contains tools to help you plot colony data.
'''

from __future__ import division, absolute_import, print_function

from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

from vivarium.library.units import remove_units
from vivarium.processes.derive_colony_shape import Variables


INCH_PER_COL = 4
INCH_PER_ROW = 2
SUBPLOT_W_SPACE = 0.4
SUBPLOT_H_SPACE = 1.5


#: Key for circumference in path timeseries
CIRCUMFERENCE_PATH = (Variables.CIRCUMFERENCE,)
#: Key for surface area in path timeseries
AREA_PATH = (Variables.AREA,)
#: Key for major axis in path timeseries
MAJOR_AXIS_PATH = (Variables.MAJOR_AXIS,)
#: Key for minor axis in path timeseries
MINOR_AXIS_PATH = (Variables.MINOR_AXIS,)
#: Key to which the cirumference-to-area ratio will be written in the
#: path timeseries
CIRCUMFERENCE_AREA_RATIO_PATH = 'circumference / surface_area'
#: Key to which the number of colonies is written in the path timeseries
NUM_COLONIES_PATH = 'Number of Colonies'
#: Key to which the ratio of major axis to minor axis is written in the
#: path timeseries
AXIS_RATIO_PATH = '(Major Axis) / (Minor Axis)'


# METRIC DERIVERS

def _derive_circumference_area_ratio(path_ts):
    if not (CIRCUMFERENCE_PATH in path_ts and AREA_PATH in path_ts):
        return
    circumference = path_ts[CIRCUMFERENCE_PATH]
    area = path_ts[AREA_PATH]
    ratio = [
        [
            c / a
            for c, a in zip(circumference_list, area_list)
        ]
        for circumference_list, area_list in zip(circumference, area)
    ]
    path_ts[CIRCUMFERENCE_AREA_RATIO_PATH] = ratio

def _derive_axis_ratio(path_ts):
    if not (MAJOR_AXIS_PATH in path_ts and MINOR_AXIS_PATH in path_ts):
        return
    major = path_ts[MAJOR_AXIS_PATH]
    minor = path_ts[MINOR_AXIS_PATH]
    ratio = [
        [
            major_val / minor_val
            for major_val, minor_val in zip(major_list, minor_list)
        ]
        for major_list, minor_list in zip(major, minor)
    ]
    path_ts[AXIS_RATIO_PATH] = ratio


#: List of metric derivers that will be applied to the path timeseries
_METRIC_DERIVERS = [
    _derive_circumference_area_ratio,
    _derive_axis_ratio,
]


def plot_colony_metrics(
    path_ts, title_size=16, tick_label_size=12, max_cols=5
):
    '''Plot colony metrics over time.

    Metric mean is plotted with SEM error bands.

    Arguments:
        path_ts (dict): Path timeseries of the data to plot. Each item
            in the dictionary should have as its key the path and as its
            value a list of values for each timepoint. Each value should
            be a list of metric values, one entry per colony. The
            dictionary should have one additional key, ``time``, whose
            value is a list of times for each timepoint.
        title_size (int): Font size for the title of each plot
        tick_label_size (int): Font size for each plot's axis tick
            labels.
        max_cols (int): The maximum number of columns. We add columns
            until we hit this limit, and only then do we add rows.

    Returns:
        matplotlib.figure.Figure: The plot as a Figure object.
    '''
    path_ts = remove_units(path_ts)
    for deriver in _METRIC_DERIVERS:
        deriver(path_ts)
    times = path_ts['time']
    del path_ts['time']
    # path_ts has tuples for keys. Here we turn those into strings so
    # that numpy doesn't iterate through the path elements
    path_ts = {
        str(key): val for key, val in path_ts.items()
    }
    arbitrary_metric = list(path_ts.keys())[0]
    path_ts[NUM_COLONIES_PATH] = [
        len(timepoint) for timepoint in path_ts[arbitrary_metric]
    ]
    # Create Figure
    paths = sorted(path_ts.keys())
    n_cols = min(len(paths), max_cols)
    n_rows = int(np.ceil(len(paths) / n_cols))
    fig = plt.figure(
        figsize=(INCH_PER_COL * n_cols, INCH_PER_ROW * n_rows))
    grid = plt.GridSpec(
        ncols=n_cols, nrows=n_rows, wspace=SUBPLOT_W_SPACE,
        hspace=SUBPLOT_H_SPACE
    )

    # Assign paths to subplot coordinates
    padding = [None] * int(n_cols * n_rows - len(paths))
    paths += padding
    paths_grid = np.array(paths)
    paths_grid = paths_grid.reshape((n_rows, n_cols))

    # Create the subplots
    for i in range(n_rows):
        for j in range(n_cols):
            path = paths_grid[i, j]
            if path is None:
                continue
            ax = fig.add_subplot(grid[i, j])
            # Configure axes and titles
            for tick_type in ('major', 'minor'):
                ax.tick_params(
                    axis='both', which=tick_type,
                    labelsize=tick_label_size
                )
            ax.title.set_text(path)
            ax.title.set_fontsize(title_size)
            ax.set_xlim([times[0], times[-1]])
            ax.xaxis.get_offset_text().set_fontsize(tick_label_size)
            ax.yaxis.get_offset_text().set_fontsize(tick_label_size)
            ax.set_xlabel('time (s)', fontsize=title_size)
            # Plot data
            data = path_ts[path]
            if path == NUM_COLONIES_PATH:
                ax.plot(times, data)
            else:
                means = []
                sems = []
                plot_times = []
                for i_time, metrics_list in enumerate(data):
                    if not metrics_list:
                        continue
                    array = np.array(metrics_list)
                    means.append(np.mean(array))
                    sems.append(
                        stats.sem(array) if len(array) > 1 else 0)
                    plot_times.append(times[i_time])
                x = np.array(plot_times)
                y = np.array(means)
                yerr = np.array(sems)
                yerr[np.isnan(yerr)] = 0
                ax.plot(x, y)
                ax.fill_between(x, y - yerr, y + yerr, alpha=0.2)

    return fig


def plot_metric_across_experiments(
    path_ts_dict, path, title=None, xlabel='time (s)', ylabel=None,
    title_size=16, tick_label_size=12,
):
    '''Overlay plots of a single metric from different experiments.

    Parameters:
        path_ts_dict (dict): Map from the string to use as the label for
            the experiment in the legend to that experiment's path
            timeseries.
        path (tuple): Path to plot. Should be a key in each value of
            ``path_ts_dict``.
        title (str): Plot title. If None, no title is set.
        xlabel (str): X-axis label. If None, no label is set.
        ylabel (str): Y-axis label. If None, no label is set.
        title_size (float): Font size for plot and axis titles.
        tick_label_size (float): Font size for tick labels.

    Returns:
        The figure with the plot.
    '''
    fig, ax = plt.subplots()

    # Set labels and font sizes
    if title is not None:
        ax.set_title(title, fontsize=title_size)
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=title_size)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=title_size)
    ax.xaxis.get_offset_text().set_fontsize(tick_label_size)
    ax.yaxis.get_offset_text().set_fontsize(tick_label_size)

    # Plot data
    for label, path_ts in path_ts_dict.items():
        data = path_ts[path]
        times = path_ts['time']
        if path == NUM_COLONIES_PATH:
            ax.plot(times, data, label=label)
        else:
            plot_times = []
            means = []
            for i, metrics_list in enumerate(data):
                if not metrics_list:
                    continue
                means.append(np.mean(metrics_list))
                plot_times.append(times[i])
            ax.plot(plot_times, means, label=label)

    ax.legend()
    fig.tight_layout()
    return fig
