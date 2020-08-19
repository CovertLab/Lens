import argparse
import os

from matplotlib import pyplot as plt
import numpy as np

from vivarium.core.experiment import get_in

from vivarium.analysis.analyze import Analyzer


PATH_TO_AGENTS = ('agents',)
PATH_TO_DEAD = ('boundary', 'dead')
LIVE_COLOR = 'red'
DEAD_COLOR = 'black'
ALPHA = 0.5
FILENAME = 'expression_survival.pdf'
OUT_DIR = 'out'


def plot_expression_survival(
    data, path_to_variable, xlabel, time_range=(0, 1)
):
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

    live_averages = []
    dead_averages = []
    for agent, levels in expression_levels.items():
        if agent in die:
            dead_averages.append(np.mean(levels))
        else:
            live_averages.append(np.mean(levels))

    fig, ax = plt.subplots(figsize=(6, 2))
    ax.scatter(
        live_averages, np.random.random(len(live_averages)),
        label='Survive', color=LIVE_COLOR, alpha=ALPHA,
    )
    ax.scatter(
        dead_averages, np.random.random(len(dead_averages)),
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


def main():
    parser = argparse.ArgumentParser()
    Analyzer.add_connection_args(parser)
    parser.add_argument(
        'experiment_id',
        type=str,
        help='ID of experiment to analyze.',
    )
    parser.add_argument(
        'path',
        type=str,
        help='Comma-separated path from agent root to plotted variable',
    )
    parser.add_argument(
        '--force', '-f',
        action='store_true',
        default=False,
        help=(
            'Write plots even if output directory already exists. This '
            'could overwrite your existing plots'
        ),
    )
    parser.add_argument(
        '--xlabel',
        default='',
        help='Label for x-axis.',
    )
    parser.add_argument(
        '--start',
        default=0,
        type=float,
        help='Start of times to consider as fraction of total time.',
    )
    parser.add_argument(
        '--end',
        default=1,
        type=float,
        help='End of times to consider as fraction of total time.',
    )

    args = parser.parse_args()
    path_to_variable = tuple(args.path.split(','))
    data, _ = Analyzer.get_data(
        args, args.experiment_id)

    out_dir = os.path.join(OUT_DIR, args.experiment_id)
    if os.path.exists(out_dir):
        if not args.force:
            raise IOError('Directory {} already exists'.format(out_dir))
    else:
        os.makedirs(out_dir)
    fig = plot_expression_survival(
        data, path_to_variable, args.xlabel, (args.start, args.end))
    fig.savefig(os.path.join(out_dir, FILENAME))


if __name__ == '__main__':
    main()
