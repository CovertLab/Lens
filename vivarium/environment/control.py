from __future__ import absolute_import, division, print_function

import copy
import time
import uuid

from vivarium.actor.control import ActorControl, AgentCommand
from vivarium.environment.make_media import Media
from vivarium.utils.units import units


class ShepherdControl(ActorControl):
    """Send messages to agents in the system to control execution."""

    def __init__(self, actor_config):
        super(ShepherdControl, self).__init__(str(uuid.uuid1()), actor_config)

    def add_cell(self, agent_type, actor_config):
        # TODO(jerry): Bring back the --variant choice?
        self.add_agent(
            str(uuid.uuid1()),
            agent_type,
            actor_config)

    def init_experiment(self, args, exp_config):
        default_experiment_id = exp_config.get('default_experiment_id')
        lattice_config = exp_config.get('lattice_config')
        environment_type = exp_config.get('environment_type')
        actor_config = exp_config.get('actor_config')
        agent_type = exp_config.get('agent_type')
        num_cells = exp_config.get('num_cells')
        actor_config['boot_config'].update(lattice_config)

        # get from args
        experiment_id = args['experiment_id']
        number = args.get('number')
        if number == 0:
            number = num_cells
        if not experiment_id:
            experiment_id = self.get_experiment_id(default_experiment_id)

        print('Creating experiment id {}: {} {} agents in {} environment\n'.format(
            experiment_id, number, agent_type, environment_type))

        # boot environment
        self.add_agent(experiment_id, environment_type, actor_config)
        time.sleep(10) # wait for the environment to boot

        # boot agents
        for index in range(number):
            self.add_cell(agent_type or args['type'], dict(actor_config, **{
                'boot': 'vivarium.environment.boot',
                'outer_id': experiment_id,
                'working_dir': args['working_dir'],
                'seed': index}))

    def lattice_experiment(self, args, actor_config):
        # define experiment: environment type and agent type
        experiment_id = 'lattice_experiment'
        environment_type = 'lattice'
        agent_type = 'growth_division'

        # overwrite default environment config
        lattice_config = {
            'name': 'lattice_experiment',
            'description': 'minimal growth_division agents are placed in a lattice environment'}

        exp_config = {
            'default_experiment_id': experiment_id,
            'lattice_config': lattice_config,
            'environment_type': environment_type,
            'actor_config': actor_config,
            'agent_type': agent_type,
            'num_cells': 1}

        self.init_experiment(args, exp_config)

    def growth_division_experiment(self, args, actor_config):

        # define experimental environment and agents
        experiment_id = 'growth_division'
        environment_type = 'lattice'
        agent_type = 'growth_division'

        # overwrite default environment config
        lattice_config = {
            'name': 'growth_division_experiment',
            'description': 'minimal growth_division agents are placed in a lattice environment'}

        exp_config = {
            'default_experiment_id': experiment_id,
            'lattice_config': lattice_config,
            'environment_type': environment_type,
            'actor_config': actor_config,
            'agent_type': agent_type,
            'num_cells': 1}

        self.init_experiment(args, exp_config)

    def ecoli_core_experiment(self, args, actor_config):

        # define experiment: environment type and agent type
        experiment_id = 'glc-g6p'
        environment_type = 'ecoli_core_glc'
        agent_type = 'shifter'

        # overwrite default environment config
        lattice_config = {
            'name': 'ecoli_core_experiment',
            'description': 'ecoli_core_experiment'}

        exp_config = {
            'default_experiment_id': experiment_id,
            'lattice_config': lattice_config,
            'environment_type': environment_type,
            'actor_config': actor_config,
            'agent_type': agent_type,
            'num_cells': 1}

        self.init_experiment(args, exp_config)

    def chemotaxis_experiment(self, args, actor_config):
        # define experiment: environment type and agent type
        experiment_id = 'chemotaxis'
        environment_type = 'measp_long'
        agent_type = 'minimal_chemotaxis'

        # overwrite default environment config
        lattice_config = {
            'name': 'chemotaxis_experiment',
            'description': 'a long environment with a static gradient of glucose and a-methyl-DL-aspartic acid (MeAsp) '
               'for observing chemotactic cells in action. Optimal chemotaxis is observed in a narrow range '
               'of CheA activity, where concentration of CheY-P falls into the operating range of flagellar motors.',
        }

        exp_config = {
            'default_experiment_id': experiment_id,
            'lattice_config': lattice_config,
            'environment_type': environment_type,
            'actor_config': actor_config,
            'agent_type': agent_type,
            'num_cells': 1}

        self.init_experiment(args, exp_config)

    def swarm_experiment(self, args, actor_config):
        # define experiment: environment type and agent type
        experiment_id = 'swarm'
        environment_type = 'measp_large'
        agent_type = 'minimal_chemotaxis'

        # overwrite default environment config
        lattice_config = {
            'name': 'swarm_experiment',
            'description': 'a large experiment for running swarms of chemotactic cells',
        }

        exp_config = {
            'default_experiment_id': experiment_id,
            'lattice_config': lattice_config,
            'environment_type': environment_type,
            'actor_config': actor_config,
            'agent_type': agent_type,
            'num_cells': 1}

        self.init_experiment(args, exp_config)



class EnvironmentCommand(AgentCommand):
    """
    Extend `AgentCommand` with new commands related to the lattice and ecoli experiments
    """

    def __init__(self, choices=[], description=''):
        full_description = '''
    Run an agent for the environmental context simulation.
    
    The commands are:
    `add --id OUTER_ID [--type T] [--config C]` ask the Shepherd to add an agent of
        type T with JSON configuration C to the environment OUTER_ID,
    `remove --id AGENT_ID` ask all Shepherds to remove agent AGENT_ID,
    `remove --prefix ID` ask all Shepherds to remove agents "ID*",
    `run [--id OUTER_ID]` start or resume one or all simulations,
    `pause [--id OUTER_ID]` pause one or all simulations,
    `divide --id AGENT_ID` ask a cell agent to divide,
    `shutdown [--id OUTER_ID]` shut down one or all environment agents and their
         connected agents,
    `experiment [--number N] [--type T] [--working-dir D]` ask the Shepherd to run
        a lattice environment with N agents of type T,
    'glc-g6p-experiment [--number N] [--type T]` ask the Shepherd to run a
        chemotaxis environment with N agents of type T
    'chemotaxis-experiment [--number N] [--type T]` ask the Shepherd to run a
        chemotaxis environment with N agents of type T
    ''' + description

        full_choices = [
            'growth-division-experiment',
            'ecoli-core-experiment',
            'chemotaxis-experiment',
            'swarm-experiment',
            ] + choices

        super(EnvironmentCommand, self).__init__(
            full_choices,
            {}, # no additional args
            full_description)

    def experiment(self, args):
        self.require(args, 'number', 'working_dir')
        control = ShepherdControl({'kafka_config': self.get_kafka_config()})
        control.lattice_experiment(args, self.actor_config)
        control.shutdown()

    def growth_division_experiment(self, args):
        self.require(args, 'number')
        control = ShepherdControl({'kafka_config': self.get_kafka_config()})
        control.growth_division_experiment(args, self.actor_config)
        control.shutdown()

    def ecoli_core_experiment(self, args):
        self.require(args, 'number')
        control = ShepherdControl({'kafka_config': self.get_kafka_config()})
        control.ecoli_core_experiment(args, self.actor_config)
        control.shutdown()

    def chemotaxis_experiment(self, args):
        self.require(args, 'number')
        control = ShepherdControl({'kafka_config': self.get_kafka_config()})
        control.chemotaxis_experiment(args, self.actor_config)
        control.shutdown()

    def swarm_experiment(self, args):
        self.require(args, 'number')
        control = ShepherdControl({'kafka_config': self.get_kafka_config()})
        control.swarm_experiment(args, self.actor_config)
        control.shutdown()

    def add_arguments(self, parser):
        parser = super(EnvironmentCommand, self).add_arguments(parser)

        parser.add_argument(
            '-e', '--experiment_id',
            type=str,
            help='The experiment id')

        parser.add_argument(
            '-f', '--emit_field',
            type=str,
            default='minimal',
            help='emitted media field')

        parser.add_argument(
            '-m', '--media',
            type=str,
            default='minimal',
            help='The environment media')

        parser.add_argument(
            '-t', '--timeline',
            type=str,
            default=None,
            help='The timeline')

        return parser

def run():
    command = EnvironmentCommand()
    command.execute()

if __name__ == '__main__':
    run()
