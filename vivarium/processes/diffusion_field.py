'''
===============
Diffusion Field
===============
'''

from __future__ import absolute_import, division, print_function

import sys
import os
import argparse

import numpy as np
from scipy import constants
from scipy.ndimage import convolve
import matplotlib.pyplot as plt

from vivarium.core.process import Process
from vivarium.core.composition import (
    simulate_process,
    PROCESS_OUT_DIR
)
from vivarium.library.lattice_utils import (
    count_to_concentration,
    get_bin_site,
    get_bin_volume,
)
from vivarium.library.units import units

from vivarium.plots.multibody_physics import plot_snapshots

NAME = 'diffusion_field'

# laplacian kernel for diffusion
LAPLACIAN_2D = np.array([[0.0, 1.0, 0.0], [1.0, -4.0, 1.0], [0.0, 1.0, 0.0]])
AVOGADRO = constants.N_A


def gaussian(deviation, distance):
    return np.exp(-np.power(distance, 2.) / (2 * np.power(deviation, 2.)))

def make_gradient(gradient, n_bins, size):
    '''Create a gradient from a configuration

    **Uniform**

    A uniform gradient fills the field evenly with each molecule, at
    the concentrations specified.

    Example configuration:

    .. code-block:: python

        'gradient': {
            'type': 'uniform',
            'molecules': {
                'mol_id1': 1.0,
                'mol_id2': 2.0
            }},

    **Gaussian**

    A gaussian gradient multiplies the base concentration of the given
    molecule by a gaussian function of distance from center and
    deviation. Distance is scaled by 1/1000 from microns to millimeters.

    Example configuration:

    .. code-block:: python

        'gradient': {
            'type': 'gaussian',
            'molecules': {
                'mol_id1':{
                    'center': [0.25, 0.5],
                    'deviation': 30},
                'mol_id2': {
                    'center': [0.75, 0.5],
                    'deviation': 30}
            }},

    **Linear**

    A linear gradient sets a site's concentration (c) of the given
    molecule as a function of distance (d) from center and slope (b),
    and base concentration (a). Distance is scaled by 1/1000 from
    microns to millimeters.

    .. math::
        c = a + b * d

    Example configuration:

    .. code-block:: python

        'gradient': {
            'type': 'linear',
            'molecules': {
                'mol_id1':{
                    'center': [0.0, 0.0],
                    'base': 0.1,
                    'slope': -10},
                'mol_id2': {
                    'center': [1.0, 1.0],
                    'base': 0.1,
                    'slope': -5}
            }},

    **Exponential**

    An exponential gradient sets a site's concentration (c) of the given
    molecule as a function of distance (d) from center, with parameters
    base (b) and scale (a). Distance is scaled by 1/1000 from microns to
    millimeters. Note: base > 1 makes concentrations increase from the
    center.

    .. math::

        c=a*b^d.

    Example configuration:

    .. code-block:: python

        'gradient': {
            'type': 'exponential',
            'molecules': {
                'mol_id1':{
                    'center': [0.0, 0.0],
                    'base': 1+2e-4,
                    'scale': 1.0},
                'mol_id2': {
                    'center': [1.0, 1.0],
                    'base': 1+2e-4,
                    'scale' : 0.1}
            }},

    Parameters:
        gradient: Configuration dictionary that includes the ``type``
            key to specify the type of gradient to make.
        n_bins: A list of two elements that specify the number of bins
            to have along each axis.
        size: A list of two elements that specifies the size of the
            environment.
    '''
    bins_x = n_bins[0]
    bins_y = n_bins[1]
    length_x = size[0]
    length_y = size[1]
    fields = {}

    if gradient.get('type') == 'gaussian':
        for molecule_id, specs in gradient['molecules'].items():
            field = np.ones((bins_x, bins_y), dtype=np.float64)
            center = [specs['center'][0] * length_x,
                      specs['center'][1] * length_y]
            deviation = specs['deviation']

            for x_bin in range(bins_x):
                for y_bin in range(bins_y):
                    # distance from middle of bin to center coordinates
                    dx = (x_bin + 0.5) * length_x / bins_x - center[0]
                    dy = (y_bin + 0.5) * length_y / bins_y - center[1]
                    distance = np.sqrt(dx ** 2 + dy ** 2)
                    scale = gaussian(deviation, (distance/1000))
                    # multiply gradient by scale
                    field[x_bin][y_bin] *= scale
            fields[molecule_id] = field

    elif gradient.get('type') == 'linear':
        for molecule_id, specs in gradient['molecules'].items():
            field = np.zeros((bins_x, bins_y), dtype=np.float64)
            center = [specs['center'][0] * length_x,
                      specs['center'][1] * length_y]
            base = specs.get('base', 0.0)
            slope = specs['slope']

            for x_bin in range(bins_x):
                for y_bin in range(bins_y):
                    dx = (x_bin + 0.5) * length_x / bins_x - center[0]
                    dy = (y_bin + 0.5) * length_y / bins_y - center[1]
                    distance = np.sqrt(dx ** 2 + dy ** 2)
                    field[x_bin][y_bin] += base + slope * (distance/1000)
            fields[molecule_id] = field

    elif gradient.get('type') == 'exponential':
        for molecule_id, specs in gradient['molecules'].items():
            field = np.zeros((bins_x, bins_y), dtype=np.float64)
            center = [specs['center'][0] * length_x,
                      specs['center'][1] * length_y]
            base = specs['base']
            scale = specs.get('scale', 1)

            for x_bin in range(bins_x):
                for y_bin in range(bins_y):
                    dx = (x_bin + 0.5) * length_x / bins_x - center[0]
                    dy = (y_bin + 0.5) * length_y / bins_y - center[1]
                    distance = np.sqrt(dx ** 2 + dy ** 2)
                    field[x_bin][y_bin] = scale * base ** (distance/1000)
            fields[molecule_id] = field

    elif gradient.get('type') == 'uniform':
        for molecule_id, fill_value in gradient['molecules'].items():
            fields[molecule_id] = np.full((bins_x, bins_y), fill_value, dtype=np.float64)

    return fields


class DiffusionField(Process):
    '''
    Diffusion in 2-dimensional fields of molecules with agent exchange

    Agent uptake and secretion occurs at agent locations.

    Notes:

    * Diffusion constant of glucose in 0.5 and 1.5 percent agarose gel
      is around :math:`6 * 10^{-10} \\frac{m^2}{s}` (Weng et al. 2005.
      Transport of glucose and poly(ethylene glycol)s in agarose gels).
    * Conversion to micrometers:
      :math:`6 * 10^{-10} \\frac{m^2}{s}=600 \\frac{micrometers^2}{s}`.
    '''

    name = NAME
    defaults = {
        'time_step': 1,
        'molecules': ['glc'],
        'initial_state': {},
        'n_bins': [10, 10],
        'bounds': [10, 10],
        'depth': 3000.0,  # um
        'diffusion': 5e-1,
        'gradient': {},
        'agents': {},
    }

    def __init__(self, initial_parameters=None):
        if initial_parameters is None:
            initial_parameters = {}

        # initial state
        self.molecule_ids = initial_parameters.get('molecules', self.defaults['molecules'])
        self.initial_state = initial_parameters.get('initial_state', self.defaults['initial_state'])

        # parameters
        self.n_bins = initial_parameters.get('n_bins', self.defaults['n_bins'])
        self.bounds = initial_parameters.get('bounds', self.defaults['bounds'])
        depth = initial_parameters.get('depth', self.defaults['depth'])

        # diffusion
        diffusion = initial_parameters.get('diffusion', self.defaults['diffusion'])
        bins_x = self.n_bins[0]
        bins_y = self.n_bins[1]
        length_x = self.bounds[0]
        length_y = self.bounds[1]
        dx = length_x / bins_x
        dy = length_y / bins_y
        dx2 = dx * dy
        self.diffusion = diffusion / dx2
        self.diffusion_dt = 0.01
        # self.diffusion_dt = 0.5 * dx ** 2 * dy ** 2 / (2 * self.diffusion * (dx ** 2 + dy ** 2))

        # volume, to convert between counts and concentration
        self.bin_volume = get_bin_volume(self.n_bins, self.bounds, depth)

        # initialize gradient fields
        gradient = initial_parameters.get('gradient', self.defaults['gradient'])
        if gradient:
            gradient_fields = make_gradient(gradient, self.n_bins, self.bounds)
            self.initial_state.update(gradient_fields)

        # agents
        self.initial_agents = initial_parameters.get('agents', self.defaults['agents'])

        parameters = {}
        parameters.update(initial_parameters)

        super(DiffusionField, self).__init__(parameters)

    def ports_schema(self):
        local_concentration_schema = {
            molecule: {
                '_default': 0.0,
                '_updater': 'set'}
            for molecule in self.molecule_ids}

        schema = {'agents': {}}
        for agent_id, states in self.initial_agents.items():
            location = states['boundary'].get('location', [])
            schema['agents'][agent_id] = {
                'boundary': {
                    'location': {
                        '_value': location},
                }
            }
        glob_schema = {
            '*': {
                'boundary': {
                    'location': {
                        '_default': [0.5 * bound for bound in self.bounds],
                        '_updater': 'set'},
                    'external': local_concentration_schema}}}
        schema['agents'].update(glob_schema)

        # fields
        fields_schema = {
            'fields': {
                field: {
                    '_value': self.initial_state.get(field, self.ones_field()),
                    '_updater': 'accumulate',
                    '_emit': True,
                }
                for field in self.molecule_ids
            },
        }
        schema.update(fields_schema)

        # dimensions
        dimensions_schema = {
            'dimensions': {
                'bounds': {
                    '_value': self.parameters['bounds'],
                    '_updater': 'set',
                    '_emit': True,
                },
                'n_bins': {
                    '_value': self.parameters['n_bins'],
                    '_updater': 'set',
                    '_emit': True,
                },
                'depth': {
                    '_value': self.parameters['depth'],
                    '_updater': 'set',
                    '_emit': True,
                }
            },
        }
        schema.update(dimensions_schema)
        return schema

    def next_update(self, timestep, states):
        fields = states['fields']
        agents = states['agents']

        # diffuse field
        delta_fields = self.diffuse(fields, timestep)

        # get each agent's local environment
        local_environments = self.get_local_environments(agents, fields)

        update = {'fields': delta_fields}
        if local_environments:
            update.update({'agents': local_environments})

        return update

    def count_to_concentration(self, count):
        return count_to_concentration(
            count * units.count, self.bin_volume * units.L
        ).to(units.mmol / units.L).magnitude

    def get_bin_site(self, location):
        return get_bin_site(location, self.n_bins, self.bounds)

    def get_single_local_environments(self, specs, fields):
        bin_site = self.get_bin_site(specs['location'])
        local_environment = {}
        for mol_id, field in fields.items():
            local_environment[mol_id] = field[bin_site]
        return local_environment

    def get_local_environments(self, agents, fields):
        local_environments = {}
        if agents:
            for agent_id, specs in agents.items():
                local_environments[agent_id] = {'boundary': {}}
                local_environments[agent_id]['boundary']['external'] = \
                    self.get_single_local_environments(specs['boundary'], fields)
        return local_environments

    def ones_field(self):
        return np.ones((self.n_bins[0], self.n_bins[1]), dtype=np.float64)

    # diffusion functions
    def diffusion_delta(self, field, timestep):
        ''' calculate concentration changes cause by diffusion'''
        field_new = field.copy()
        t = 0.0
        dt = min(timestep, self.diffusion_dt)
        while t < timestep:
            field_new += self.diffusion * dt * convolve(field_new, LAPLACIAN_2D, mode='reflect')
            t += dt

        return field_new - field

    def diffuse(self, fields, timestep):
        delta_fields = {}
        for mol_id, field in fields.items():

            # run diffusion if molecule field is not uniform
            if len(set(field.flatten())) != 1:
                delta = self.diffusion_delta(field, timestep)
            else:
                delta = np.zeros_like(field)
            delta_fields[mol_id] = delta

        return delta_fields


# testing
def get_random_field_config(config={}):
    bounds = config.get('bounds', (20, 20))
    n_bins = config.get('n_bins', (10, 10))
    return {
        'molecules': ['glc'],
        'initial_state': {
            'glc': np.random.rand(n_bins[0], n_bins[1])},
        'n_bins': n_bins,
        'bounds': bounds}

def get_gaussian_config(config={}):
    molecules = config.get('molecules', ['glc'])
    bounds = config.get('bounds', (50, 50))
    n_bins = config.get('n_bins', (20, 20))
    center = config.get('center', [0.5, 0.5])
    deviation = config.get('deviation', 5)
    diffusion = config.get('diffusion', 5e-1)

    return {
        'molecules': molecules,
        'n_bins': n_bins,
        'bounds': bounds,
        'diffusion': diffusion,
        'gradient': {
            'type': 'gaussian',
            'molecules': {
                'glc': {
                    'center': center,
                    'deviation': deviation}}}}

def get_exponential_config(config={}):
    molecules = config.get('molecules', ['glc'])
    bounds = config.get('bounds', (40, 40))
    n_bins = config.get('n_bins', (20, 20))
    center = config.get('center', [1.0, 1.0])
    base = config.get('base', 1 + 2e-4)
    scale = config.get('scale', 0.1)
    diffusion = config.get('diffusion', 1e1)

    return {
        'molecules': molecules,
        'n_bins': n_bins,
        'bounds': bounds,
        'diffusion': diffusion,
        'gradient': {
            'type': 'exponential',
            'molecules': {
                'glc': {
                    'center': center,
                    'base': base,
                    'scale': scale}}}}

def get_secretion_agent_config(config={}):
    molecules = config.get('molecules', ['glc'])
    bounds = config.get('bounds', (20, 20))
    n_bins = config.get('n_bins', (10, 10))
    depth = config.get('depth', 30)
    n_agents = config.get('n_agents', 3)

    agents = {}
    for agent in range(n_agents):
        agent_id = str(agent)
        agents[agent_id] = {
            'boundary': {
                'location': [
                        np.random.uniform(0, bounds[0]),
                        np.random.uniform(0, bounds[1])],
                'external': {
                    mol_id: 0 for mol_id in molecules}}}
    return {
        'molecules': molecules,
        'n_bins': n_bins,
        'bounds': bounds,
        'depth': depth,
        'agents': agents,
        'initial_state': {
            mol_id: np.ones((n_bins[0], n_bins[1]))
            for mol_id in molecules}}

def test_diffusion_field(config=get_gaussian_config(), time=10):
    diffusion = DiffusionField(config)
    settings = {
        'return_raw_data': True,
        'total_time': time,
        'timestep': 1}
    return simulate_process(diffusion, settings)

def plot_fields(data, config, out_dir='out', filename='fields'):
    fields = {time: time_data['fields'] for time, time_data in data.items()}
    snapshots_data = {
        'fields': fields,
        'config': config}
    plot_config = {
        'out_dir': out_dir,
        'filename': filename}
    plot_snapshots(snapshots_data, plot_config)


if __name__ == '__main__':
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    parser = argparse.ArgumentParser(description='diffusion_field')
    parser.add_argument('--random', '-r', action='store_true', default=False)
    parser.add_argument('--gaussian', '-g', action='store_true', default=False)
    parser.add_argument('--exponential', '-e', action='store_true', default=False)
    parser.add_argument('--secretion', '-s', action='store_true', default=False)
    args = parser.parse_args()
    no_args = (len(sys.argv) == 1)

    if args.random or no_args:
        config = get_random_field_config()
        data = test_diffusion_field(config, 10)
        plot_fields(data, config, out_dir, 'random_field')

    if args.gaussian or no_args:
        config = get_gaussian_config()
        data = test_diffusion_field(config, 10)
        plot_fields(data, config, out_dir, 'gaussian_field')

    if args.exponential or no_args:
        config = get_exponential_config()
        data = test_diffusion_field(config, 20)
        plot_fields(data, config, out_dir, 'exponential_field')

    if args.secretion or no_args:
        config = get_secretion_agent_config()
        data = test_diffusion_field(config, 10)
        plot_fields(data, config, out_dir, 'secretion')
