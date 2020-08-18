from __future__ import absolute_import, division, print_function

import numpy as np
from scipy import constants

from vivarium.library.dict_utils import deep_merge
from vivarium.library.units import units, remove_units, Quantity
from vivarium.core.process import Process


AVOGADRO = constants.N_A * 1 / units.mol


class CellEnvironmentDiffusionFicks(Process):

    name = "cell_environment_diffusion_ficks"
    defaults = {
        'molecules_to_diffuse': [],
        'default_state': {
            'global': {
                'volume': 1.2 * units.fL,
            }
        },
        'default_default': 0,
        # This is an order-of-magnitude estimate from Nagano, K., &
        # Nikaido, H. (2009). Kinetic behavior of the major multidrug
        # efflux pump AcrB of Escherichia coli. Proceedings of the
        # National Academy of Sciences, 106(14), 5854–5858.
        # https://doi.org/10.1073/pnas.0901695106
        'permeability': 1e-6 * units.cm / units.sec,
        'surface_area_mass_ratio': 132 * units.cm**2 / units.mg,
        'volume_variable': 'volume',
    }

    def ports_schema(self):

        schema = {
            'internal': {
                # Molecule concentration in mmol/L
                molecule: {
                    '_default': self.parameters['default_default'],
                    '_divider': 'set',
                }
                for molecule in self.parameters['molecules_to_diffuse']
            },
            'external': {
                # Molecule concentration in mmol/L
                molecule: {
                    '_default': self.parameters['default_default'],
                    '_divider': 'set',
                }
                for molecule in self.parameters['molecules_to_diffuse']
            },
            'fields': {
                # Molecule count in mmol
                molecule: {
                    '_default': np.full(
                        (1, 1), self.parameters['default_default']),
                }
                for molecule in self.parameters['molecules_to_diffuse']
            },
            'global': {
                self.parameters['volume_variable']: {
                    '_default': self.parameters['default_default'],
                },
                'location': {
                    '_default': [0.5, 0.5],
                },
                'dry_mass': {
                    '_default': (
                        self.parameters['default_default'] * units.fg),
                },
            },
            'dimensions': {
                'bounds': {
                    '_default': [1, 1],
                },
                'n_bins': {
                    '_default': [1, 1],
                },
                'depth': {
                    '_default': 1,
                },
            },
        }

        for port, port_conf in self.parameters['default_state'].items():
            for variable, default in port_conf.items():
                if variable in schema[port]:
                    schema[port][variable]['_default'] = default

        return schema

    def next_update(self, timestep, states):
        permeability = self.parameters['permeability']
        area_mass = self.parameters['surface_area_mass_ratio']
        mass = states['global']['dry_mass']
        flux_mmol = {}
        for molecule in self.parameters['molecules_to_diffuse']:
            # Flux is positive when leaving the cell
            delta_concentration = (
                states['internal'][molecule]
                - states['external'][molecule]
            ) * units.mmol / units.L
            # Fick's first law of diffusion:
            rate = permeability * area_mass * delta_concentration
            flux = rate * mass * timestep * units.sec
            flux_mmol[molecule] = flux
        flux_counts = {
            molecule: flux * AVOGADRO
            for molecule, flux in flux_mmol.items()
        }
        volume_key = self.parameters['volume_variable']
        if not isinstance(states['global'][volume_key], Quantity):
            cell_volume = states['global'][volume_key] * units.fL
        else:
            cell_volume = states['global'][volume_key]
        update = {
            'fields': {
                molecule: {
                    '_value': mol_flux.to(units.count).magnitude,
                    '_updater': {
                        'updater': 'update_field_with_exchange',
                        'port_mapping': {
                            'global': 'global',
                            'dimensions': 'dimensions',
                        },
                    },
                }
                for molecule, mol_flux in flux_counts.items()
            },
            'internal': {
                molecule: - (
                    mol_flux / cell_volume
                ).to(units.mmol / units.L).magnitude
                for molecule, mol_flux in flux_mmol.items()
            },
        }
        return update
