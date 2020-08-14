from __future__ import absolute_import, division, print_function

import os

from vivarium.core.process import Generator
from vivarium.core.composition import (
    simulate_compartment_in_experiment,
    plot_simulation_output,
    plot_compartment_topology,
    COMPARTMENT_OUT_DIR
)
from vivarium.compartments.gene_expression import plot_gene_expression_output
from vivarium.compartments.flagella_expression import get_flagella_expression_config

# processes
from vivarium.processes.metabolism import (
    Metabolism,
    get_iAF1260b_config
)
from vivarium.processes.convenience_kinetics import (
    ConvenienceKinetics,
    get_glc_lct_config
)
from vivarium.processes.transcription import Transcription
from vivarium.processes.translation import Translation
from vivarium.processes.degradation import RnaDegradation
from vivarium.processes.complexation import Complexation
from vivarium.processes.chemoreceptor_cluster import ReceptorCluster
from vivarium.processes.flagella_activity import FlagellaActivity
from vivarium.processes.membrane_potential import MembranePotential
from vivarium.processes.division_volume import DivisionVolume

# compartments
from vivarium.compartments.master import default_metabolism_config
from vivarium.compartments.flagella_expression import get_flagella_expression_config

NAME = 'chemotaxis_master'


def metabolism_timestep_config(time_step=1):
    config = default_metabolism_config()
    config.update({'time_step': time_step})
    return config

class ChemotaxisMaster(Generator):

    defaults = {
        'dimensions_path': ('dimensions',),
        'fields_path': ('fields',),
        'boundary_path': ('boundary',),
        'transport': get_glc_lct_config(),
        'metabolism': metabolism_timestep_config(10),
        'transcription': get_flagella_expression_config({})['transcription'],
        'translation': get_flagella_expression_config({})['translation'],
        'degradation': get_flagella_expression_config({})['degradation'],
        'complexation': get_flagella_expression_config({})['complexation'],
        'receptor': {'ligand': 'MeAsp'},
        'flagella': {'n_flagella': 5},
        'PMF': {},
        'division': {},
    }

    def __init__(self, config=None):
        super(ChemotaxisMaster, self).__init__(config)

    def generate_processes(self, config):
        # Transport
        transport = ConvenienceKinetics(config.get('transport'))

        # Metabolism
        # add target fluxes from transport
        target_fluxes = transport.kinetic_rate_laws.reaction_ids
        config['metabolism']['constrained_reaction_ids'] = target_fluxes
        metabolism = Metabolism(config['metabolism'])

        # flagella expression
        transcription = Transcription(config['transcription'])
        translation = Translation(config['translation'])
        degradation = RnaDegradation(config['degradation'])
        complexation = Complexation(config['complexation'])

        # chemotaxis -- flagella activity, receptor activity, and PMF
        receptor = ReceptorCluster(config['receptor'])
        flagella = FlagellaActivity(config['flagella'])
        PMF = MembranePotential(config['PMF'])

        # Division
        # get initial volume from metabolism
        if 'division' not in config:
            config['division'] = {}
        config['division']['initial_state'] = metabolism.initial_state
        division = DivisionVolume(config['division'])

        return {
            'metabolism': metabolism,
            'transport': transport,
            'transcription': transcription,
            'translation': translation,
            'degradation': degradation,
            'complexation': complexation,
            'receptor': receptor,
            'flagella': flagella,
            'PMF': PMF,
            'division': division,
        }

    def generate_topology(self, config):
        dimensions_path = config['dimensions_path']
        fields_path = config['fields_path']
        boundary_path = config['boundary_path']
        external_path = boundary_path + ('external',)
        return {
            'transport': {
                'internal': ('internal',),
                'external': external_path,
                'fields': ('null',),  # metabolism's exchange is used
                'fluxes': ('flux_bounds',),
                'global': boundary_path,
                'dimensions': dimensions_path,
            },
            'metabolism': {
                'internal': ('internal',),
                'external': external_path,
                'reactions': ('reactions',),
                'fields': fields_path,
                'flux_bounds': ('flux_bounds',),
                'global': boundary_path,
                'dimensions': dimensions_path
            },
            'transcription': {
                'chromosome': ('chromosome',),
                'molecules': ('internal',),
                'proteins': ('proteins',),
                'transcripts': ('transcripts',),
                'factors': ('concentrations',),
                'global': boundary_path,
            },
            'translation': {
                'ribosomes': ('ribosomes',),
                'molecules': ('internal',),
                'transcripts': ('transcripts',),
                'proteins': ('proteins',),
                'concentrations': ('concentrations',),
                'global': boundary_path,
            },
            'degradation': {
                'transcripts': ('transcripts',),
                'proteins': ('proteins',),
                'molecules': ('internal',),
                'global': boundary_path,
            },
            'complexation': {
                'monomers': ('proteins',),
                'complexes': ('proteins',),
                'global': boundary_path,
            },
            'receptor': {
                'external': external_path,
                'internal': ('internal',),
            },
            'flagella': {
                'internal': ('internal',),
                'membrane': ('membrane',),
                'internal_counts': ('proteins',),
                'flagella': ('flagella',),
                'boundary': boundary_path,
            },
            'PMF': {
                'external': external_path,
                'membrane': ('membrane',),
                'internal': ('internal',),
            },
            'division': {
                'global': boundary_path,
            }
        }

def run_chemotaxis_master(out_dir):
    total_time = 10

    # make the compartment
    compartment = ChemotaxisMaster({})

    # save the topology network
    settings = {'show_ports': True}
    plot_compartment_topology(
        compartment,
        settings,
        out_dir)

    # run an experinet
    settings = {
        'timestep': 1,
        'total_time': total_time}
    timeseries = simulate_compartment_in_experiment(compartment, settings)

    volume_ts = timeseries['boundary']['volume']
    print('growth: {}'.format(volume_ts[-1]/volume_ts[0]))

    # plots
    # simulation output
    plot_settings = {
        'max_rows': 40,
        'remove_zeros': True,
        'skip_ports': ['reactions', 'prior_state', 'null']}
    plot_simulation_output(timeseries, plot_settings, out_dir)

    # gene expression plot
    gene_exp_plot_config = {
        'name': 'flagella_expression',
        'ports': {
            'transcripts': 'transcripts',
            'proteins': 'proteins',
            'molecules': 'internal'}}
    plot_gene_expression_output(
        timeseries,
        gene_exp_plot_config,
        out_dir)

def test_chemotaxis_master(total_time=5):
    compartment = ChemotaxisMaster({})
    settings = {
        'timestep': 1,
        'total_time': total_time}
    return simulate_compartment_in_experiment(compartment, settings)


if __name__ == '__main__':
    out_dir = os.path.join(COMPARTMENT_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    run_chemotaxis_master(out_dir)
