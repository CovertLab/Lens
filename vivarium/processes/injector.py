from __future__ import absolute_import, division, print_function

import os

from vivarium.compartment.composition import (
    convert_to_timeseries,
    plot_simulation_output,
    simulate_process,
    TEST_OUT_DIR,
    save_timeseries,
)
from vivarium.compartment.process import Process

NAME = 'injector'


class Injector(Process):

    def __init__(self, initial_parameters=None):
        if initial_parameters is None:
            initial_parameters = {}

        self.substrate_rate_map = initial_parameters['substrate_rate_map']
        ports = {
            'internal': [
                substrate for substrate in self.substrate_rate_map
            ]
        }
        super(Injector, self).__init__(ports, initial_parameters)

    def default_settings(self):
        default_settings = {
            'state': {
                'internal': {
                    substrate: 0 for substrate in self.substrate_rate_map
                }
            },
            'emitter_keys': {'internal': self.substrate_rate_map.keys},
        }
        return default_settings

    def next_update(self, timestep, states):
        return {
            'internal': {
                substrate: timestep * rate
                for substrate, rate in self.substrate_rate_map.items()
            }
        }


def run_injector():
    parameters = {
        'substrate_rate_map': {'toy': 1.0},
        }
    injector = Injector(parameters)
    settings = {
        'total_time': 10,
    }
    data = simulate_process(injector, settings)
    timeseries = convert_to_timeseries(data)
    return timeseries


def test_injector():
    timeseries = run_injector()
    # Expect [0, 1, ..., 10] because 0 at start
    expected = [i for i in range(11)]
    assert expected == timeseries['internal']['toy']


def main():
    out_dir = os.path.join(TEST_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    timeseries = run_injector()
    plot_settings = {}
    plot_simulation_output(timeseries, plot_settings, out_dir)
    save_timeseries(timeseries, out_dir)


if __name__ == '__main__':
    main()
