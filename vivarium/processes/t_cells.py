from __future__ import absolute_import, division, print_function

import os
import sys
import copy
import math
import random
import argparse

from vivarium.library.units import units
from vivarium.library.dict_utils import deep_merge
from vivarium.core.process import Process
from vivarium.core.composition import (
    simulate_process_in_experiment,
    plot_simulation_output,
    PROCESS_OUT_DIR,
)
from vivarium.core.process import Generator
from vivarium.processes.meta_division import MetaDivision


NAME = 'T_cell'
TIMESTEP = 60 * 10  # seconds


class TCellProcess(Process):
    """T-cell process with 2 states

    States:
        - PD1p (PD1+)
        - PD1n (PD1-)

    TODOs
        - make this work!
    """

    name = NAME
    defaults = {
        'time_step': TIMESTEP,
        'diameter': 10 * units.um,
        'initial_PD1n': 0.8,
        'transition_PD1n_to_PD1p_10days': 0.95,  # nearly 95% by 10 days (Wherry 2007)

        # death rates (Petrovas 2007)
        'PDL1_critical_number': 1e4,  # threshold number of PDL1 molecules/cell to trigger death
        'death_PD1p_14hr': 0.7,  # 0.7 / 14 hrs
        'death_PD1n_14hr': 0.2,  # 0.2 / 14 hrs
        'death_PD1p_next_to_PDL1p_14hr': 0.95,  # 0.95 / 14 hrs

        # production rates
        'PD1n_IFNg_production': 1.6e4/3600,  # molecules/cell/second (Bouchnita 2017)
        'PD1p_IFNg_production': 0.0,  # molecules/cell/second
        'PD1p_PD1_equilibrium': 5e4,  # equilibrium value of PD1 for PD1p (TODO -- get reference)

        # division rate (Petrovas 2007)
        'PD1n_growth_8hr': 0.95,  # 95% division in 8 hours
        'PD1p_growth_8hr': 0.05,  # 5% division in 8 hours

        # migration
        'PD1n_migration': 10.0,  # um/minute (Boissonnas 2007)
        'PD1n_migration_MHCIp_tumor': 2.0,  # um/minute (Boissonnas 2007)
        'PD1n_migration_MHCIp_tumor_dwell_time': 25.0,  # minutes (Thibaut 2020)
        'PD1p_migration': 5.0,   # um/minute (Boissonnas 2007)
        'PD1p_migration_MHCIp_tumor': 1.0,   # um/minute (Boissonnas 2007)
        'PD1p_migration_MHCIp_tumor_dwell_time': 10.0,  # minutes (Thibaut 2020)

        # killing  ()
        # TODO - @Eran - This value needs to be multiplied by 100 to deal with timestep usually 0.4 packet/min
        # linear production over 4 hr up to a total of 102+-20 granules # (Betts, 2004), (Zhang, 2006)
        'PD1n_cytotoxic_packets': 40,  # number of packets/min to each contacted tumor cell

        # 1:10 fold reduction of PD1+ T cell cytotoxic production (Zelinskyy, 2005)
        'PD1p_cytotoxic_packets': 4,  # number of packets/min to each contacted tumor cell

        # settings
        'self_path': tuple(),
    }

    def __init__(self, parameters=None):
        super(TCellProcess, self).__init__(parameters)

        if random.uniform(0, 1) < self.parameters['initial_PD1n']:
            self.initial_state = 'PD1n'
        else:
            self.initial_state = 'PD1p'
        self.self_path = self.parameters['self_path']

    def ports_schema(self):
        return {
            'globals': {
                'death': {
                    '_default': False,
                    '_emit': True,
                    '_updater': 'set'},
                'divide': {
                    '_default': False,
                    '_emit': True,
                    '_updater': 'set'}
            },
            'internal': {
                'cell_state': {
                    '_default': self.initial_state,
                    '_emit': True,
                    '_updater': 'set'
                }
            },
            'boundary': {
                'diameter': {
                    '_default': self.parameters['diameter']
                },
                'IFNg': {
                    '_default': 0,
                    '_emit': True,
                    '_updater': 'accumulate',
                },
                'PD1': {
                    '_default': 0,
                    '_emit': True,
                    '_updater': 'set',
                },  # membrane protein, promotes T-cell death
                'cytotoxic_packets': {
                    '_default': 0,
                    '_emit': True,
                    '_updater': 'set',
                },  # release into the tumor cells
            },
            'neighbors': {
                'PDL1': {
                    '_default': 0,
                    '_emit': True,
                },
                'MHCI': {
                    '_default': 0,
                    '_emit': True,
                },
            },
        }

    def next_update(self, timestep, states):
        cell_state = states['internal']['cell_state']
        PDL1 = states['neighbors']['PDL1']
        MHCI = states['neighbors']['MHCI']

        # death
        if cell_state == 'PD1n':
            rate = -math.log(1 - self.parameters['death_PD1n_14hr'])  # rate for probability function of time
            timestep_fraction = timestep / 50400  # for what fraction of 14 hours? (14*60*60 seconds)
            prob_divide = 1 - math.exp(-rate * timestep_fraction)
            if random.uniform(0, 1) < prob_divide:
                print('DEATH PD1- cell!')
                return {
                    '_delete': {
                        'path': self.self_path
                    }
                }

        elif cell_state == 'PD1p':
            if PDL1 <= self.parameters['PDL1_critical_number']:
                rate = -math.log(1 - self.parameters['death_PD1p_next_to_PDL1p_14hr'])  # rate for probability function of time
                timestep_fraction = timestep / 50400  # for what fraction of 14 hours? (14*60*60 seconds)
                prob_death = 1 - math.exp(-rate * timestep_fraction)
                if random.uniform(0, 1) < prob_death:
                    print('DEATH PD1+ cell with PDL1!')
                    return {
                        '_delete': {
                            'path': self.self_path
                        }
                    }

            else:
                rate = -math.log(1 - self.parameters['death_PD1p_14hr'])  # rate for probability function of time
                timestep_fraction = timestep / 50400  # for what fraction of 14 hours? (14*60*60 seconds)
                prob_death = 1 - math.exp(-rate * timestep_fraction)
                if random.uniform(0, 1) < prob_death:
                    print('DEATH PD1+ cell without PDL1!')
                    return {
                        '_delete': {
                            'path': self.self_path
                        }
                    }

        # division
        if cell_state == 'PD1n':
            rate = -math.log(1 - self.parameters['PD1n_growth_8hr'])  # rate for probability function of time
            timestep_fraction = timestep / 28800  # for what fraction of 8 hours? (8*60*60 seconds)
            prob_divide = 1 - math.exp(-rate * timestep_fraction)
            if random.uniform(0, 1) < prob_divide:
                print('DIVIDE PD1- cell!')
                return {
                    'globals': {
                        'divide': True
                    }
                }

        elif cell_state == 'PD1p':
            rate = -math.log(1 - self.parameters['PD1p_growth_8hr'])  # rate for probability function of time
            timestep_fraction = timestep / 28800  # for what fraction of 8 hours? (8*60*60 seconds)
            prob_divide = 1 - math.exp(-rate * timestep_fraction)
            if random.uniform(0, 1) < prob_divide:
                print('DIVIDE PD1+ cell!')
                return {
                    'globals': {
                        'divide': True
                    }
                }

        # state transition
        new_cell_state = cell_state
        if cell_state == 'PD1n':
            rate = -math.log(1 - self.parameters['transition_PD1n_to_PD1p_10days'])  # rate for probability function of time
            timestep_fraction = timestep / 864000  # for what fraction of 10 days? (60*60*24*10 seconds)
            prob_transition = 1 - math.exp(-rate * timestep_fraction)
            if random.uniform(0, 1) < prob_transition:
                new_cell_state = 'PD1p'
        elif cell_state == 'PD1p':
            pass

        # behavior
        IFNg = 0
        PD1 = 0
        cytotoxic_packets = 0

        # TODO migration
        #  #also dependent on MHCI+/PDL1+
        #  #50% bound vs. 30% bound in flow cytometry experiment on low vs. high
        #  #Migration rates above in parameters

        # Killing by passing cytotoxic packets to tumor
        if new_cell_state == 'PD1n':

            # TODO - @Eran - how do we make T cells only interact with 1 of neighbors at a time?
            # and how do we reference this cell - see last elif
            # check conditional 4 fold reduction
            if MHCI >= 1e4: #molecules/neighbor cell

                # Max production happens for PD1- T cells in contact with MHCI+ tumor
                cytotoxic_packets = self.parameters['PD1n_cytotoxic_packets'] * timestep

                # produce IFNg  # rates are determined above based on 1 minute
                IFNg = self.parameters['PD1n_IFNg_production'] * timestep

            elif MHCI < 1e4:  # molecules/neighbor cell

                # 4 fold reduction in production in T cells in contact with MHCI- tumor
                # (Bohm, 1998), (Merritt, 2003)
                cytotoxic_packets = self.parameters['PD1n_cytotoxic_packets'] / 4 * timestep

                # produce IFNg  # rates are determined above based on 1 minute
                IFNg = self.parameters['PD1n_IFNg_production'] / 4 * timestep

            # target behavior 3 contacts required for cell death, 1-4 cells killed/day

            # TODO  - elif T cell is not in contact with tumor (no cytotoxic packets)
            #   continue

        elif new_cell_state == 'PD1p':
            PD1 = self.parameters['PD1p_PD1_equilibrium']

            if MHCI >= 1e4:  # molecules/neighbor cell

                # Max production happens for PD1- T cells in contact with MHCI+ tumor
                cytotoxic_packets = self.parameters['PD1p_cytotoxic_packets'] * timestep

                # produce IFNg  # rates are determined above based on 1 minute
                IFNg = self.parameters['PD1p_IFNg_production'] * timestep

            elif MHCI < 1e4:  # molecules/neighbor cell

                # 4 fold reduction in production in T cells in contact with MHCI- tumor
                # (Bohm, 1998), (Merritt, 2003)
                cytotoxic_packets = self.parameters['PD1p_cytotoxic_packets'] / 4 * timestep

                # produce IFNg  # rates are determined above based on 1 minute
                IFNg = self.parameters['PD1p_IFNg_production'] / 4 * timestep

            # target behavior 3 contacts required for cell death, 1-4 cells killed/day

            # TODO  - elif T cell is not in contact with tumor (no cytotoxic packets)
            #   continue

        return {
            'internal': {
                'cell_state': new_cell_state
            },
            'boundary': {
                'IFNg': IFNg,
                'PD1': PD1,
                'cytotoxic_packets': cytotoxic_packets,
            },
        }


class TCellCompartment(Generator):

    defaults = {
        'boundary_path': ('boundary',),
        'agents_path': ('..', '..', 'agents',),
        'daughter_path': tuple()}

    def __init__(self, config):
        self.config = config
        for key, value in self.defaults.items():
            if key not in self.config:
                self.config[key] = value

        # paths
        self.boundary_path = config.get('boundary_path', self.defaults['boundary_path'])
        self.agents_path = config.get('agents_path', self.defaults['agents_path'])

    def generate_processes(self, config):
        daughter_path = config['daughter_path']
        agent_id = config['agent_id']

        division_config = dict(
            config.get('division', {}),
            daughter_path=daughter_path,
            agent_id=agent_id,
            compartment=self)

        t_cell = TCellProcess(config.get('growth', {}))
        division = MetaDivision(division_config)

        return {
            't_cell': t_cell,
            'division': division}

    def generate_topology(self, config):
        return {
            't_cell': {
                'internal': ('internal',),
                'boundary': self.boundary_path,
                'global': self.boundary_path},
            'division': {
                'global': self.boundary_path,
                'cells': self.agents_path},
            }

def get_timeline():
    timeline = [
        (0, {
            ('neighbors', 'PDL1'): 0.0,
            ('neighbors', 'MHCI'): 0.0,
        }),
        (10 * TIMESTEP, {
            ('neighbors', 'PDL1'): 0.0,
            ('neighbors', 'MHCI'): 0.0,
        }),
        (20 * TIMESTEP, {
            ('neighbors', 'PDL1'): 0.0,
            ('neighbors', 'MHCI'): 0.0,
        }),
        (30 * TIMESTEP, {
            ('neighbors', 'PDL1'): 5e5,
            ('neighbors', 'MHCI'): 5e5,
        }),
        (40 * TIMESTEP, {
            ('neighbors', 'PDL1'): 5e5,
            ('neighbors', 'MHCI'): 5e5,
        }),
        (50 * TIMESTEP, {
            ('neighbors', 'PDL1'): 5e5,
            ('neighbors', 'MHCI'): 5e5,
        }),
        (60 * TIMESTEP, {
            ('neighbors', 'PDL1'): 5e5,
            ('neighbors', 'MHCI'): 5e5,
        }),
        (70 * TIMESTEP, {
            ('neighbors', 'PDL1'): 0.0,
            ('neighbors', 'MHCI'): 0.0,
        }),
        (80 * TIMESTEP, {
            ('neighbors', 'PDL1'): 0.0,
            ('neighbors', 'MHCI'): 0.0,
        }),
        (90 * TIMESTEP, {
            ('neighbors', 'PDL1'): 0.0,
            ('neighbors', 'MHCI'): 0.0,
        }),
        (100 * TIMESTEP, {
            ('neighbors', 'PDL1'): 5e5,
            ('neighbors', 'MHCI'): 5e5,
        }),
        (110 * TIMESTEP, {}),
    ]
    return timeline


#TODO - @Eran - how do I get more in the output like graphs/data?
def test_single_t_cell(
    total_time=1200,
    timeline=None,
    out_dir='out'):

    t_cell_process = TCellProcess({})
    if timeline is not None:
        settings = {
            'timeline': {
                'timeline': timeline}}
    else:
        settings = {'total_time': total_time}

    # run experiment
    timeseries = simulate_process_in_experiment(t_cell_process, settings)

    # plot
    plot_settings = {
        'remove_zeros': False
    }
    plot_simulation_output(timeseries, plot_settings, out_dir)



if __name__ == '__main__':
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    parser = argparse.ArgumentParser(description='ODE expression')
    parser.add_argument('--single', '-s', action='store_true', default=False)
    parser.add_argument('--timeline', '-t', action='store_true', default=False)
    args = parser.parse_args()
    no_args = (len(sys.argv) == 1)

    total_time = 12000
    if args.single:
        test_single_t_cell(
            total_time=total_time,
            out_dir=out_dir)

    if args.timeline or no_args:
        timeline = get_timeline()
        test_single_t_cell(
            timeline=timeline,
            out_dir=out_dir)
