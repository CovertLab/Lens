'''
======================================
Expression Dotplot Colored by Survival
======================================
'''
from matplotlib import pyplot as plt
import numpy as np

from vivarium.core.experiment import get_in


PATH_TO_AGENTS = ('agents',)
PATH_TO_DEAD = ('boundary', 'dead')
LIVE_COLOR = 'red'
DEAD_COLOR = 'black'
ALPHA = 0.5


def plot_expression_survival(
    data, path_to_variable, xlabel, time_range=(0, 1)
):
    '''Create Expression Dotplot Colored by Survival

    Plot one dot for each cell along an axis to indicate that cell's
    average expression level for a specified protein. The dot color
    reflects whether the cell survived long enough to divide.

    Note that only the expression levels while the cell is alive are
    considered in the average.

    Parameters:
        data (dict): The raw data emitted from the simulation.
        path_to_variable (tuple): Path from the agent root to the
            variable that holds the protein's expression level. We do
            not adjust for cell volume, so this should be a
            concentration.
        xlabel (str): Label for x-axis.
        time_range (tuple): Tuple of two :py:class:`float`s that are
            fractions of the total simulated time period. These
            fractions indicate the start and end points (inclusive) of
            the time range to consider when calculating average
            expression level.

    Returns:
        plt.Figure: The finished figure.
    '''
    expression_levels = dict()
    die = set()
    end_time = max(data.keys())
    for time, time_data in data.items():
        if (time < time_range[0] * end_time
                or time > time_range[1] * end_time):
            continue
        agents_data = get_in(time_data, PATH_TO_AGENTS)
        for agent, agent_data in agents_data.items():
            lst = expression_levels.setdefault(agent, [])
            value = get_in(agent_data, path_to_variable)
            if value is not None:
                lst.append(value)
            if get_in(agent_data, PATH_TO_DEAD, False):
                die.add(agent)
            # Only count values when cell is alive
            elif value is not None:
                lst.append(value)

    live_averages = []
    dead_averages = []
    for agent, levels in expression_levels.items():
        if not levels:
            continue
        if agent in die:
            dead_averages.append(np.mean(levels))
        else:
            live_averages.append(np.mean(levels))

    fig, ax = plt.subplots(figsize=(6, 2))
    ax.scatter(
        live_averages, [0.1] * len(live_averages),
        label='Survive', color=LIVE_COLOR, alpha=ALPHA,
    )
    ax.scatter(
        dead_averages, [0.1] * len(dead_averages),
        label='Die', color=DEAD_COLOR, alpha=ALPHA,
    )
    ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylim([0, 1.25])
    ax.get_yaxis().set_visible(False)
    for spine_name in ('left', 'top', 'right'):
        ax.spines[spine_name].set_visible(False)
    ax.spines['bottom'].set_position('zero')
    fig.tight_layout()
    return fig
