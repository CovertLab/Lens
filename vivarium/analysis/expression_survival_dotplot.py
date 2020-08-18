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

    args = parser.parse_args()
    path_to_variable = tuple(args.path.split(','))
    data, _ = Analyzer.get_data(
        args, args.experiment_id)

    expression_levels = dict()
    die = set()
    for _, time_data in data.items():
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

    fig, ax = plt.subplots()
    ax.scatter(
        live_averages, np.random.random(len(live_averages)),
        color=LIVE_COLOR, alpha=ALPHA,
    )
    ax.scatter(
        dead_averages, np.random.random(len(dead_averages)),
        color=DEAD_COLOR, alpha=ALPHA,
    )
    ax.set_ylim([0, 1])
    ax.get_yaxis().set_visible(False)
    for spine_name in ('left', 'top', 'right'):
        ax.spines[spine_name].set_visible(False)
    ax.spines['bottom'].set_position('zero')

    out_dir = os.path.join(OUT_DIR, args.experiment_id)
    if os.path.exists(out_dir):
        if not args.force:
            raise IOError('Directory {} already exists'.format(out_dir))
    else:
        os.makedirs(out_dir)
    fig.savefig(os.path.join(out_dir, FILENAME))


if __name__ == '__main__':
    main()
