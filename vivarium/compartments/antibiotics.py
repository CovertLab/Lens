from __future__ import absolute_import, division, print_function

import copy
import math
import os

from vivarium.core.process import Generator
from vivarium.core.composition import (
    simulate_compartment_in_experiment,
    plot_simulation_output,
    flatten_timeseries,
    save_timeseries,
    load_timeseries,
    REFERENCE_DATA_DIR,
    COMPARTMENT_OUT_DIR,
    assert_timeseries_close,
)
from vivarium.library.dict_utils import deep_merge
from vivarium.processes.antibiotic_transport import AntibioticTransport
from vivarium.processes.death import DeathFreezeState
from vivarium.processes.diffusion_cell_environment import (
    CellEnvironmentDiffusion,
)
from vivarium.processes.division_volume import DivisionVolume
from vivarium.processes.growth import Growth
from vivarium.processes.ode_expression import ODE_expression


NUM_DIVISIONS = 3
DEFAULT_DIVISION_SECS = 2400  # seconds to divide
INITIAL_INTERNAL_ANTIBIOTIC = 0.5
INITIAL_EXTERNAL_ANTIBIOTIC = 1


class Antibiotics(Generator):

    defaults = {
        'fields_path': ('fields',),
        'dimensions_path': ('dimensions',),
        'ode_expression': {
            'transcription_rates': {
                'AcrAB-TolC_RNA': 1e-3,
            },
            'translation_rates': {
                'AcrAB-TolC': 1.0,
            },
            'degradation_rates': {
                # Set for on the order of 100 RNAs at equilibrium
                'AcrAB-TolC_RNA': 1.0,
                # Set so exporter concentration reaches equilibrium
                'AcrAB-TolC': 1e-3,
            },
            'protein_map': {
                'AcrAB-TolC': 'AcrAB-TolC_RNA'
            },
        },
        'antibiotic_transport': {
            'initial_pump': 0.0,
            'initial_internal_antibiotic': INITIAL_INTERNAL_ANTIBIOTIC,
            'initial_external_antibiotic': INITIAL_EXTERNAL_ANTIBIOTIC,
            'pump_kcat': 5e-4,
        },
        'death': {
            'detectors': {
                'antibiotic': {
                    'antibiotic_threshold': 0.95,
                },
            },
            'targets': [
                'antibiotic_transport', 'growth', 'expression', 'death',
                'division',
            ]
        },
        'growth': {
            # Growth rate calculated so that 2 = exp(DIVISION_TIME *
            # rate) because division process divides once cell doubles
            # in size
            'growth_rate': math.log(2) / DEFAULT_DIVISION_SECS
        },
        'division': {},
        'diffusion': {
            'default_state': {
                'membrane': {
                    'porin': 1e-20,
                },
                'external': {
                    'antibiotic': INITIAL_EXTERNAL_ANTIBIOTIC,
                },
                'internal': {
                    'antibiotic': INITIAL_INTERNAL_ANTIBIOTIC,
                },
            },
            'molecules_to_diffuse': ['antibiotic'],
            'permeability_per_porin': {
                'porin': 5e-3,
            },
        },
    }
    name = 'antibiotics_compartment'

    def __init__(self, config):
        merged_config = copy.deepcopy(Antibiotics.defaults)
        deep_merge(merged_config, config)
        super(Antibiotics, self).__init__(merged_config)

    def generate_processes(self, config):
        antibiotic_transport = AntibioticTransport(
            config['antibiotic_transport'])
        growth = Growth(config['growth'])
        expression = ODE_expression(config['ode_expression'])
        death = DeathFreezeState(config['death'])
        division = DivisionVolume(config['division'])
        diffusion = CellEnvironmentDiffusion(config['diffusion'])

        return {
            'antibiotic_transport': antibiotic_transport,
            'growth': growth,
            'expression': expression,
            'death': death,
            'division': division,
            'diffusion': diffusion,
        }

    def generate_topology(self, config):
        return {
            'antibiotic_transport': {
                'internal': ('cell',),
                'external': ('external',),
                'pump_port': ('cell',),
                'fields': config['fields_path'],
                'fluxes': ('fluxes',),
                'global': ('global',),
                'dimensions': config['dimensions_path'],
            },
            'growth': {
                'global': ('global',),
            },
            'expression': {
                'counts': ('cell_counts',),
                'internal': ('cell',),
                'external': ('external',),
                'global': ('global',),
            },
            'division': {
                'global': ('global',),
            },
            'death': {
                'global': ('global',),
                'internal': ('cell',),
            },
            'diffusion': {
                'membrane': ('cell',),
                'internal': ('cell',),
                'external': ('external',),
                'fields': config['fields_path'],
                'global': ('global',),
                'dimensions': config['dimensions_path'],
            },
        }


def run_antibiotics_composite():
    DIVISION_TIME = DEFAULT_DIVISION_SECS
    sim_settings = {
        'emit_step': 1,
        'total_time': DIVISION_TIME * NUM_DIVISIONS,
        'environment': {}
    }
    config = {}
    compartment = Antibiotics(config)
    return simulate_compartment_in_experiment(compartment, sim_settings)

def test_antibiotics_composite_similar_to_reference():
    timeseries = run_antibiotics_composite()
    flattened = flatten_timeseries(timeseries)
    reference = load_timeseries(
        os.path.join(REFERENCE_DATA_DIR, Antibiotics.name + '.csv'))
    assert_timeseries_close(
        flattened, reference,
        tolerances={
            'cell_counts_AcrAB-TolC': 99999,
            'cell_counts_antibiotic': 999,
            'cell_counts_AcrAB-TolC_RNA': 9,
            'cell_counts_porin': 9999,
        }
    )


def main():
    out_dir = os.path.join(COMPARTMENT_OUT_DIR, Antibiotics.name)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    plot_settings = {
        'max_rows': 25,
    }

    timeseries = run_antibiotics_composite()
    plot_simulation_output(timeseries, plot_settings, out_dir)
    save_timeseries(timeseries, out_dir)


if __name__ == '__main__':
    main()
