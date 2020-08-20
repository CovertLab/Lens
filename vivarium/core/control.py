from __future__ import absolute_import, division, print_function

import os
import sys
import argparse

from arpeggio import (
    RegExMatch,
    ParserPython,
    OneOrMore,
)

from vivarium.core.experiment import timestamp
from vivarium.core.composition import (
    embedded_compartment_experiment,
    agent_environment_experiment,
    simulate_experiment,
    plot_agents_multigen,
    ToyCompartment,
    ToyEnvironment,
    EXPERIMENT_OUT_DIR,
)
from vivarium.analysis.analyze import Analyzer



# parsing expression grammar for agents
def agent_type(): return RegExMatch(r'[a-zA-Z0-9.\-\_]+')
def number(): return RegExMatch(r'[0-9]+')
def specification(): return agent_type, number
def rule(): return OneOrMore(specification)
agent_parser = ParserPython(rule)
def parse_agents_string(agents_string):
    all_agents = agent_parser.parse(agents_string)
    agents_config = []
    for idx, agent_specs in enumerate(all_agents):
        agents_config.append(make_agent_config(agent_specs))
    return agents_config

def make_dir(out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

class Control():

    def __init__(
            self,
            name=None,
            compartment_library={},
            experiment_library={},
            simulation_settings={},
            plot_settings={},
    ):
        if name is None:
            name = timestamp()
        self.compartment_library = compartment_library
        self.experiment_library = experiment_library
        self.args = self.add_arguments()

        # TODO experiment settings
        # TODO plot settings

        self.out_dir = os.path.join(EXPERIMENT_OUT_DIR, name)
        make_dir(self.out_dir)

    def add_arguments(self):
        parser = argparse.ArgumentParser(
            description='command line control of experiments'
        )
        parser.add_argument(
            '--agents', '-a',
            type=str,
            nargs='+',
            default=argparse.SUPPRESS,
            help='A list of agent types and numbers in the format "agent_type1 number1 agent_type2 number2"'
        )
        parser.add_argument(
            '--environment', '-v',
            type=str,
            default=argparse.SUPPRESS,
            help='the environment type'
        )
        parser.add_argument(
            '--time', '-t',
            type=int,
            default=60,
            help='simulation time, in seconds'
        )
        parser.add_argument(
            '--emit', '-m',
            type=int,
            default=1,
            help='emit interval, in seconds'
        )
        parser.add_argument(
            '--experiment', '-e',
            type=str,
            default=argparse.SUPPRESS,
            help='preconfigured experiments'
        )

        return vars(parser.parse_args())

    def execute(self):
        if self.args['experiment']:
            experiment_name = self.args['experiment']
            experiment_out_dir = os.path.join(self.out_dir, experiment_name)
            make_dir(experiment_out_dir)
            experiment_config = self.experiment_library[experiment_name]
            hierarchy = experiment_config['hierarchy']
            simulation_settings = experiment_config['simulation_settings']

        # simulate
        data = self.run_experiment(
            hierarchy=hierarchy,
            # initial_state=initial_state,
            # initial_agent_state=initial_agent_state,
            simulation_settings=simulation_settings,
        )

        # TODO -- use Analyzer
        # plot_settings['environment_config'] = environment_config
        # plot_settings['agent_type'] = agent_type
        # plot_experiment_output(
        #     data,
        #     plot_settings,
        #     out_dir,
        # )


    def run_experiment(
            self,
            # agents_config=None,
            # environment_config=None,
            initial_state=None,
            # initial_agent_state=None,
            hierarchy=None,
            simulation_settings=None,
            experiment_settings=None
    ):
        if experiment_settings is None:
            experiment_settings = {}
        if initial_state is None:
            initial_state = {}
        # if initial_agent_state is None:
        #     initial_agent_state = {}


        # make the experiment
        experiment = embedded_compartment_experiment(hierarchy)






        # # agents ids
        # agent_ids = []
        # for config in agents_config:
        #     number = config['number']
        #     if 'name' in config:
        #         name = config['name']
        #         if number > 1:
        #             new_agent_ids = [name + '_' + str(num) for num in range(number)]
        #         else:
        #             new_agent_ids = [name]
        #     else:
        #         new_agent_ids = [str(uuid.uuid1()) for num in range(number)]
        #     config['ids'] = new_agent_ids
        #     agent_ids.extend(new_agent_ids)

        # experiment = agent_environment_experiment(
        #     agents_config=agents_config,
        #     environment_config=environment_config,
        #     initial_state=initial_state,
        #     initial_agent_state=initial_agent_state,
        #     settings=experiment_settings,
        # )

        # simulate
        settings = {
            'total_time': simulation_settings['total_time'],
            'emit_step': simulation_settings['emit_step'],
            'return_raw_data': simulation_settings['return_raw_data']}
        return simulate_experiment(
            experiment,
            settings,
        )

    def get_plot_settings(
            self,
            fields=[],
            tags=[]
    ):
        return {
            'plot_types': {
                'agents': {},
                'snapshots': {
                    'fields': fields
                },
                'tags': {
                    'tag_ids': tags
                }
            }
        }

def run_control_test():
    compartment_library = {
        'agent': {
            'name': 'agent',
            'type': ToyCompartment,
            'config': {
                'external_path': ('..', 'external')
            },
        },
        'environment': {
            'name': 'environment',
            'type': ToyEnvironment,
            'config': {}
        }
    }
    experiment_library = {
        'toy': {
            'hierarchy': {
                'environment': {
                    'inner': [
                        {
                            'number': 2,
                            'name': 'agent',
                        }
                    ],
                    'topology': {
                        'agents': {'agent'}
                    }}  # TODO -- is this the best way to do embedded of agent within universe?
            },
            'simulation_settings': {
                'total_time': 30,
                'emit_step': 0.1,
            },
        },
    }

    workflow = Control(
        compartment_library=compartment_library,
        experiment_library=experiment_library
        )

    # TODO add an experiment structure? with agent in environment
    experiment = {}
    workflow.execute()


if __name__ == '__main__':
    run_control_test()
