from __future__ import absolute_import, division, print_function

import argparse
import csv
import os

from vivarium.plots.multibody_physics import (
    plot_snapshots,
    plot_tags,
)
from vivarium.core.composition import plot_agents_multigen
from vivarium.core.emitter import (
    get_atlas_client,
    get_local_client,
    data_from_database,
    SECRETS_PATH,
    timeseries_from_data,
    path_timeseries_from_embedded_timeseries,
)
from vivarium.plots.colonies import plot_colony_metrics


OUT_DIR = 'out'


class Analyzer:

    def __init__(
        self,
        snapshots_config=None,
        tags_config=None,
        timeseries_config=None,
    ):
        self.parser = self._get_parser()
        if snapshots_config is None:
            snapshots_config = {}
        if tags_config is None:
            tags_config = {}
        if timeseries_config is None:
            timeseries_config = {}
        self.snapshots_config = snapshots_config
        self.tags_config = tags_config
        self.timeseries_config = timeseries_config

    def run(self):
        args = self.parser.parse_args()
        self.plot(args)

    @staticmethod
    def get_data(args, experiment_id):
        if args.atlas:
            client = get_atlas_client(SECRETS_PATH)
        else:
            client = get_local_client(
                args.host, args.port, args.database_name)
        data, environment_config = data_from_database(
            experiment_id, client)
        del data[0]
        return data, environment_config

    @staticmethod
    def plot_snapshots(data, environment_config, out_dir, settings):
        snapshots_data = Analyzer.format_data_for_snapshots(
            data, environment_config)
        plot_config = {
            'out_dir': out_dir,
        }
        plot_config.update(settings)
        plot_snapshots(snapshots_data, plot_config)

    @staticmethod
    def format_data_for_snapshots(data, environment_config):
        agents = {
            time: timepoint['agents']
            for time, timepoint in data.items()
        }
        fields = {
            time: timepoint['fields']
            for time, timepoint in data.items()
        }
        snapshots_data = {
            'agents': agents,
            'fields': fields,
            'config': environment_config,
        }
        return snapshots_data

    @staticmethod
    def plot_tags(
        data, environment_config, tagged_molecules, out_dir, settings
    ):
        tags_data = Analyzer.format_data_for_tags(
            data, environment_config)
        plot_config = {
            'out_dir': out_dir,
            'tagged_molecules': tagged_molecules,
        }
        plot_config.update(settings)
        plot_tags(tags_data, plot_config)

    @staticmethod
    def format_data_for_tags(data, environment_config):
        agents = {
            time: timepoint['agents']
            for time, timepoint in data.items()
        }
        tags_data = {
            'agents': agents,
            'config': environment_config,
        }
        return tags_data

    @staticmethod
    def plot_timeseries(data, out_dir, settings):
        plot_settings = {
            'agents_key': 'agents',
            'title_size': 10,
            'tick_label_size': 10,
        }
        plot_settings.update(settings)
        plot_agents_multigen(data, plot_settings, out_dir)

    @staticmethod
    def plot_colony_metrics(data, out_dir):
        path_ts = Analyzer.format_data_for_colony_metrics(data)
        fig = plot_colony_metrics(path_ts)
        fig.savefig(os.path.join(out_dir, 'colonies'))

    @staticmethod
    def format_data_for_colony_metrics(data):
        embedded_ts = timeseries_from_data(data)
        colony_metrics_ts = embedded_ts['colony_global']
        colony_metrics_ts['time'] = embedded_ts['time']
        path_ts = path_timeseries_from_embedded_timeseries(
            colony_metrics_ts)
        return path_ts

    def plot(self, args):
        data, environment_config = Analyzer.get_data(
            args, args.experiment_id)
        out_dir = os.path.join(OUT_DIR, args.experiment_id)
        if os.path.exists(out_dir):
            if not args.force:
                raise IOError('Directory {} already exists'.format(out_dir))
        else:
            os.makedirs(out_dir)

        if args.snapshots:
            Analyzer.plot_snapshots(
                data, environment_config, out_dir,
                self.snapshots_config
            )
        if args.tags is not None:
            with open(args.tags, 'r') as f:
                reader = csv.reader(f)
                tagged_molecules = [
                    tuple(path) for path in reader
                ]
            Analyzer.plot_tags(
                data, environment_config, tagged_molecules, out_dir,
                self.tags_config
            )
        if args.timeseries:
            Analyzer.plot_timeseries(
                data, out_dir, self.timeseries_config)
        if args.colony_metrics:
            Analyzer.plot_colony_metrics(data, out_dir)

    def _get_parser(self):
        parser = argparse.ArgumentParser()
        Analyzer.add_connection_args(parser)
        parser.add_argument(
            'experiment_id',
            help='Experiment ID as recorded in the database',
        )
        parser.add_argument(
            '--snapshots', '-s',
            action='store_true',
            default=False,
            help='Plot snapshots',
        )
        parser.add_argument(
            '--tags', '-g',
            default=None,
            help=(
                'A path to a CSV file that lists the tagged molecules to '
                'plot. The first column should contain the name of the store '
                'under each agent boundary where the molecule is reported, '
                'and the second column should contain the name of the '
                'molecule. Setting this parameter causes a plot of the tagged '
                'molecules to be produced.'
            ),
        )
        parser.add_argument(
            '--timeseries', '-t',
            action='store_true',
            default=False,
            help='Generate line plot for each variable over time',
        )
        parser.add_argument(
            '--colony_metrics', '-c',
            action='store_true',
            default=False,
            help='Plot colony metrics',
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
        return parser

    @staticmethod
    def add_connection_args(parser):
        parser.add_argument(
            '--atlas', '-a',
            action='store_true',
            default=False,
            help=(
                'Read data from an mongoDB Atlas instead of a local mongoDB. '
                'Credentials, cluster subdomain, and database name should be '
                'specified in {}.'.format(SECRETS_PATH)
            )
        )
        parser.add_argument(
            '--port', '-p',
            default=27017,
            type=int,
            help=(
                'Port at which to access local mongoDB instance. '
                'Defaults to "27017".'
            ),
        )
        parser.add_argument(
            '--host', '-o',
            default='localhost',
            type=str,
            help=(
                'Host at which to access local mongoDB instance. '
                'Defaults to "localhost".'
            ),
        )
        parser.add_argument(
            '--database_name', '-d',
            default='simulations',
            type=str,
            help=(
                'Name of database on local mongoDB instance to read from. '
                'Defaults to "simulations".'
            )
        )


if __name__ == '__main__':
    analyzer = Analyzer()
    analyzer.run()
