from __future__ import absolute_import, division, print_function

import os
import sys
import copy
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


NAME = 'Tumor'


class TumorProcess(Process):
    """Tumor process with 2 states

    States:
        - PDL1p (PDL1+, MHCI+)
        - PDL1n (PDL1-, MHCI-)

    Required parameters:
        -

    Target behavior:

    TODOs

    """

    name = NAME
    defaults = {
        'time_step': 60,
        'diameter': 20 * units.um,
        'initial_PDL1n': 1.0, #all start out this way based on data

        # death rates
        'death_apoptosis': 1e-4,  # negligible compared to growth/killing 5 day mean life span (Gong, 2017)

        # division rate
        'PDL1n_growth': 7e-4,  # probability of division - 1/24 hr (Eden, 2011)
        #'PDL1p_growth': 0,  # Cells arrested - do not divide (data, Thibaut 2020, Hoekstra 2020)

        # migration
        'tumor_migration': 0.25,  # um/minute (Weigelin 2012)
        #TODO - @Eran - how to manage migration with square grids if migration is smaller than grid?

        #membrane equillibrium amounts #TODO
        'PDL1p_PDL1_equilibrium': 5e4, #TODO ref
        'PDL1p_MHCI_equilibrium': 5e4, #TODO ref

        # settings
        'self_path': tuple(),
    }

    def __init__(self, parameters=None):
        super(TumorProcess, self).__init__(parameters)

        if random.uniform(0, 1) < self.parameters['initial_PDL1n']:
            self.initial_state = 'PDL1n'
        else:
            self.initial_state = 'PDL1p'
        self.self_path = self.parameters['self_path']

    def ports_schema(self):
        return {
            'globals': {
                'death': {
                    '_default': False,
                    '_updater': 'set'
                },
                'divide': {
                    '_default': False,
                    '_updater': 'set'
                }
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
                'PDL1': {
                    '_default': 0,
                    '_emit': True,
                    '_updater': 'accumulate',
                },  # membrane protein, promotes T cell exhuastion and deactivation with PD1
                'MHCI': {
                    '_default': 0,
                    '_emit': True,
                    '_updater': 'set',
                },  # membrane protein, promotes Tumor death and T cell activation with TCR
                'IFNg': {
                    '_default': 0 * units.ng/units.mL,
                    '_emit': True,
                },  # cytokine changes tumor phenotype to MHCI+ and PDL1+
                'IFNg_timer': {
                    '_default': 0,
                    '_updater': 'accumulate',
                },  # cytokine changes tumor phenotype
            },
            'neighbors': {
                'PD1': {
                    '_default': 0,
                    '_emit': True,},
                'cytotoxic_packets': {
                    '_default': 0,
                    '_emit': True,
                },
            }
        }

    def next_update(self, timestep, states):
        cell_state = states['internal']['cell_state']
        cytotoxic_packets = states['neighbors']['cytotoxic_packets']
        IFNg = states['boundary']['IFNg']
        IFNg_timer = states['boundary']['IFNg_timer']

        # death by apoptosis
        if random.uniform(0, 1) < self.parameters['death_apoptosis'] * timestep:
            print('Apoptosis!')
            return {
                '_delete': {
                    'path': self.self_path},
                'globals': {
                    'death': 'apoptosis'
                }
            }

        # death by cytotoxic packets from T cells
        # should take about 120 min from start of T cell contact and about 2-3 contacts
        # need to multiply total number by 10 because multiplied T cell number by this amount
        #number needed for death refs: (Verret, 1987), (Betts, 2004), (Zhang, 2006)
        if cytotoxic_packets >= 128*10:
            print('Tcell_death!')
            return {
                '_delete': {
                    'path': self.self_path},
                'globals': {
                    'death': 'Tcell_death'}
            }

        # division
        if cell_state == 'PDL1n':
            if random.uniform(0, 1) < self.parameters['PDL1n_growth'] * timestep:
                print('PDL1n DIVIDE!')
                return {
                    'globals': {
                        'divide': True
                    }
                }
        elif cell_state == 'PDL1p':
            pass

        update = {}
        # state transition
        new_cell_state = cell_state
        if cell_state == 'PDL1n':
            if IFNg >= 1*units.ng/units.mL:
                if IFNg_timer > 6*60*60:
                    new_cell_state = 'PDL1p'
                else:
                    update.update({
                        'boundary': {
                            'IFNg_timer':timestep
                        }
                    })

            pass
        elif cell_state == 'PDL1p':
            pass

        # behavior
        MHCI = 0
        PDL1 = 0

        # TODO - @Eran - I do not think dynamics of these receptors matters as much, but I put in
        #  anyways so that we would have them

        if new_cell_state == 'PDL1p':

            PDL1 = self.parameters['PDL1p_PDL1_equilibrium']
            MHCI = self.parameters['PDL1p_MHCI_equilibrium']

        elif new_cell_state == 'PDL1n':
            pass

        # TODO migration - Do after the environment is set up

        return update



class TumorCompartment(Generator):

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

        Tumor = TumorProcess(config.get('growth', {}))
        division = MetaDivision(division_config)

        return {
            'Tumor': Tumor,
            'division': division}

    def generate_topology(self, config):
        return {
            'Tumor': {
                'internal': ('internal',),
                'boundary': self.boundary_path,
                'global': self.boundary_path},
            'division': {
                'global': self.boundary_path,
                'cells': self.agents_path},
            }


def get_PD1_timeline():
    timeline = [
        (0, {('neighbors', 'PD1'): 0.0}),
        (10, {('neighbors', 'PD1'): 5e4}),
        (20, {('neighbors', 'PD1'): 0.0}),
        (30, {('neighbors', 'PD1'): 5e4}),
        (40, {('neighbors', 'PD1'): 5e4}),
        (50, {('neighbors', 'PD1'): 5e4}),
        (60, {('neighbors', 'PD1'): 5e4}),
        (70, {('neighbors', 'PD1'): 5e4}),
        (80, {('neighbors', 'PD1'): 5e4}),
        (90, {}),
    ]
    return timeline

def get_IFNg_timeline():
    timeline = [
        (0, {('boundary', 'IFNg'): 0.0*units.ng/units.mL}),
        (100, {('boundary', 'IFNg'): 1.0*units.ng/units.mL}),
        (200, {('boundary', 'IFNg'): 2.0*units.ng/units.mL}),
        (300, {('boundary', 'IFNg'): 3.0*units.ng/units.mL}),
        (400, {('boundary', 'IFNg'): 4.0*units.ng/units.mL}),
        (500, {('boundary', 'IFNg'): 3.0*units.ng/units.mL}),
        (600, {('boundary', 'IFNg'): 2.0*units.ng/units.mL}),
        (700, {('boundary', 'IFNg'): 2.0*units.ng/units.mL}),
        (800, {('boundary', 'IFNg'): 2.0*units.ng/units.mL}),
        (900, {}),
    ]
    return timeline

def get_cytoxic_packets_timeline():
    timeline = [
        (0, {('neighbors', 'cytotoxic_packets'): 0.0}),
        (100, {('neighbors', 'cytotoxic_packets'): 10.0}),
        (200, {('neighbors', 'cytotoxic_packets'): 20.0}),
        (300, {('neighbors', 'cytotoxic_packets'): 30.0}),
        (400, {('neighbors', 'cytotoxic_packets'): 40.0}),
        (500, {('neighbors', 'cytotoxic_packets'): 70.0}),
        (600, {('neighbors', 'cytotoxic_packets'): 100.0}),
        (700, {('neighbors', 'cytotoxic_packets'): 150.0}),
        (800, {('neighbors', 'cytotoxic_packets'): 160.0}),
        (900, {}),
    ]
    return timeline

def get_combined_timeline():
    timeline = [
        (0, {
            ('neighbors', 'cytotoxic_packets'): 0.0,
            ('boundary', 'IFNg'): 0.0*units.ng/units.mL,
            ('neighbors', 'PD1'): 0.0,
        }),
        (100, {
            ('neighbors', 'cytotoxic_packets'): 10.0,
            ('boundary', 'IFNg'): 1.0*units.ng/units.mL,
            ('neighbors', 'PD1'): 5e4,
        }),
        (200, {
            ('neighbors', 'cytotoxic_packets'): 20.0,
            ('boundary', 'IFNg'): 2.0 * units.ng / units.mL,
            ('neighbors', 'PD1'): 0.0,
        }),
        (300, {
            ('neighbors', 'cytotoxic_packets'): 30.0,
            ('boundary', 'IFNg'): 3.0 * units.ng / units.mL,
            ('neighbors', 'PD1'): 5e4,
        }),
        (400, {
            ('neighbors', 'cytotoxic_packets'): 40.0,
            ('boundary', 'IFNg'): 4.0 * units.ng / units.mL,
            ('neighbors', 'PD1'): 5e4,
        }),
        (500, {
            ('neighbors', 'cytotoxic_packets'): 70.0,
            ('boundary', 'IFNg'): 3.0 * units.ng / units.mL,
            ('neighbors', 'PD1'): 5e4,
        }),
        (600, {
            ('neighbors', 'cytotoxic_packets'): 100.0,
            ('boundary', 'IFNg'): 2.0 * units.ng / units.mL,
            ('neighbors', 'PD1'): 5e4,
        }),
        (700, {
            ('neighbors', 'cytotoxic_packets'): 150.0,
            ('boundary', 'IFNg'): 2.0 * units.ng / units.mL,
            ('neighbors', 'PD1'): 5e4,
        }),
        (800, {
            ('neighbors', 'cytotoxic_packets'): 160.0,
            ('boundary', 'IFNg'): 2.0 * units.ng / units.mL,
            ('neighbors', 'PD1'): 5e4,
        }),
        (900, {}),
    ]
    return timeline


def test_single_Tumor(
        total_time=1200,
        timeline=None,
        out_dir='out'):

    Tumor_process = TumorProcess({})

    if timeline is not None:
        settings = {
            'timeline': {
                'timeline': timeline}}
    else:
        settings = {'total_time': total_time}

    # run experiment
    timeseries = simulate_process_in_experiment(Tumor_process, settings)

    # plot
    plot_settings = {'remove_zeros': False}
    plot_simulation_output(timeseries, plot_settings, out_dir)



if __name__ == '__main__':
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    parser = argparse.ArgumentParser(description='tumor cells')
    parser.add_argument('--single', '-s', action='store_true', default=False)
    parser.add_argument('--timeline', '-t', action='store_true', default=False)
    args = parser.parse_args()
    no_args = (len(sys.argv) == 1)

    #TODO - @Eran is time always in seconds?
    total_time = 12000
    if args.single:
        test_single_Tumor(
            total_time=total_time,
            out_dir=out_dir)

    if args.timeline or no_args:
        timeline = get_combined_timeline()
        test_single_Tumor(
            timeline=timeline,
            out_dir=out_dir)
