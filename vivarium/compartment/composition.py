import copy
import os

import matplotlib.pyplot as plt

from vivarium.utils.dict_utils import deep_merge, deep_merge_check
from vivarium.compartment.process import (
    initialize_state,
    simulate_compartment,
    Compartment,
    COMPARTMENT_STATE,
    Process,
)
from vivarium.utils.units import units

# processes
from vivarium.processes.derive_globals import DeriveGlobals, AVOGADRO
from vivarium.processes.derive_counts import DeriveCounts
from vivarium.processes.derive_concentrations import DeriveConcs


deriver_library = {
    'globals': DeriveGlobals,
    'mmol_to_counts': DeriveCounts,
    'counts_to_mmol': DeriveConcs}



def get_derivers(process_list, topology):
    '''
    get the derivers for a list of processes

    requires:
        - process_list -- (list) with configured processes
        - topology -- (dict) with topology of the processes connected to compartment ports

    returns: (dict) with:
        {'deriver_processes': processes,
        'deriver_topology': topology}
    '''


    # get deriver configuration
    deriver_config = {}
    full_deriver_topology = {}
    deriver_topology = {}
    for level in process_list:
        for process_id, process in level.items():
            process_settings = process.default_settings()
            setting = process_settings.get('deriver_setting', [])
            try:
                port_map = topology[process_id]
            except:
                print('{} topology port mismatch'.format(process_id))
                raise

            for set in setting:
                type = set['type']
                keys = set['keys']
                source_port = set['source_port']
                target_port = set['derived_port']
                try:
                    source_compartment_port = port_map[source_port]
                    target_compartment_port = port_map[target_port]
                except:
                    print('{} source/target port mismatch'.format(process_id))
                    raise

                deriver_topology = {
                    type: {
                        source_port: source_compartment_port,
                        target_port: target_compartment_port,
                        'global': 'global'}}
                deep_merge(full_deriver_topology, deriver_topology)

                # TODO -- what if multiple different source/targets?
                # TODO -- merge overwrites them. need list extend
                # ports for configuration
                ports = {
                    source_port: keys,
                    target_port: keys}
                config = {type: {'ports': ports}}

                deep_merge(deriver_config, config)

    # configure the derivers
    deriver_processes = {}
    for type, config in deriver_config.items():
        deriver_processes[type] = deriver_library[type](config)

    # add global deriver
    # TODO -- configure global deriver to get mass
    global_deriver = {
        'global_deriver': DeriveGlobals({})}

    global_deriver_topology = {
        'global_deriver': {
            'global': 'global'}}

    deep_merge(deriver_topology, global_deriver_topology)

    # put the global deriver and additional derivers in two layers
    processes = [
        global_deriver,
        deriver_processes]

    return {
        'deriver_processes': processes,
        'deriver_topology': deriver_topology}

def get_schema(process_list, topology):
    schema = {}
    for level in process_list:
        for process_id, process in level.items():
            process_settings = process.default_settings()
            process_schema = process_settings.get('schema', {})
            try:
                port_map = topology[process_id]
            except:
                print('{} topology port mismatch'.format(process_id))
                raise

            # go through each port, and get the schema
            for process_port, settings in process_schema.items():
                compartment_port = port_map[process_port]
                compartment_schema = {
                    compartment_port: settings}

                ## TODO -- check for mismatch
                deep_merge_check(schema, compartment_schema)

    return schema

def process_in_compartment(process, settings={}):
    ''' put a process in a compartment, with all derivers added '''
    process_settings = process.default_settings()
    compartment_state_port = settings.get('compartment_state_port')

    processes = [{'process': process}]
    topology = {
        'process': {
            port: port for port in process.ports
            if (not compartment_state_port
                or port != compartment_state_port)
        }
    }

    if compartment_state_port:
        topology['process'][compartment_state_port] = COMPARTMENT_STATE

    derivers = get_derivers(processes, topology)
    deriver_processes = derivers['deriver_processes']
    deriver_topology = derivers['deriver_topology']

    # add deriver processes
    processes.extend(deriver_processes)

    # add deriver topology
    topology.update(deriver_topology)

    # make the state
    state_dict = process_settings['state']
    states = initialize_state(processes, topology, state_dict)

    options = {
        'topology': topology}

    return Compartment(processes, states, options)

def simulate_process_with_environment(process, settings={}):
    ''' simulate a process in a compartment with an environment '''
    compartment = process_in_compartment(process, settings)
    return simulate_with_environment(compartment, settings)

def simulate_process(process, settings={}):
    ''' simulate a process in a compartment with no environment '''
    compartment = process_in_compartment(process, settings)
    return simulate_compartment(compartment, settings)

def simulate_with_environment(compartment, settings={}):
    '''
    run a compartment simulation with an environment.
    requires processes made for LatticeCompartment, with environment_port and exchange_port
    '''

    # parameters
    nAvogadro = AVOGADRO

    # get environment configuration
    environment_port = settings['environment_port']
    env_volume = settings.get('environment_volume', 1e-12) * units.L
    exchange_port = settings.get('exchange_port')
    if exchange_port:
        exchange_ids = list(compartment.states[exchange_port].keys())
    else:
        print('no exchange port! simulate environment without exchange')
    environment = compartment.states.get(environment_port)
    exchange = compartment.states.get(exchange_port)

    # get timeline
    total_time = settings.get('total_time', 10)
    timeline = copy.deepcopy(settings.get('timeline', [(total_time, {})]))
    end_time = timeline[-1][0]
    timestep = compartment.time_step

    # initialize saved_state
    saved_state = {}

    ## run simulation
    time = 0
    saved_state[time] = compartment.current_state()
    while time < end_time:
        time += timestep
        for (t, change_dict) in timeline:
            if time >= t:
                for port_id, change in change_dict.items():
                    port = compartment.states.get(port_id)
                    port.assign_values(change)
                timeline.pop(0)

        # update compartment
        compartment.update(timestep)

        ## apply exchange to environment
        # get counts, convert to change in concentration
        if exchange:
            delta_counts = exchange.state_for(exchange_ids)
            mmol_to_counts = (nAvogadro.to('1/mmol') * env_volume).to('L/mmol').magnitude
            delta_concs = {mol_id: counts / mmol_to_counts for mol_id, counts in delta_counts.items()}
            environment.apply_update(delta_concs)

            # reset exchange
            reset_exchange = {key: 0 for key in exchange_ids}
            exchange.assign_values(reset_exchange)

        saved_state[time] = compartment.current_state()

    return saved_state

def convert_to_timeseries(sim_output):
    '''
    input:
        - saved_states (dict) with {timestep: state_dict}
    returns:
        - timeseries (dict) with timeseries in lists {'time': [], 'port1': {'state': []}}
    TODO --  currently assumes state is 1 dictionary deep. make a more general state embedding
    '''

    time_vec = list(sim_output.keys())
    initial_state = sim_output[time_vec[0]]
    timeseries = {port: {state: []
        for state, initial in states.items()}
        for port, states in initial_state.items()}
    timeseries['time'] = time_vec

    for time, all_states in sim_output.items():
        for port, states in all_states.items():
            for state_id, state in states.items():
                timeseries[port][state_id].append(state)

    return timeseries


def set_axes(ax, show_xaxis=False):
    ax.ticklabel_format(style='sci', axis='y', scilimits=(-5,5))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(right=False, top=False)
    if not show_xaxis:
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(bottom=False, labelbottom=False)


def plot_simulation_output(timeseries, settings={}, out_dir='out'):
    '''
    plot simulation output, with rows organized into separate columns.

    Requires:
        - timeseries (dict). This can be obtained from simulation output with convert_to_timeseries()
        - settings (dict) with:
            {
            'max_rows': (int) ports with more states than this number of states get wrapped into a new column
            'remove_zeros': (bool) if True, timeseries with all zeros get removed
            'remove_flat': (bool) if True, timeseries with all the same value get removed
            'skip_ports': (list) entire ports that won't be plotted
            'overlay': (dict) with
                {'bottom_port': 'top_port'}  ports plotted together by matching state_ids, with 'top_port' in red
            'show_state': (list) with [('port_id', 'state_id')]
                for all states that will be highlighted, even if they are otherwise to be removed
            }
    TODO -- some molecules have 'inf' concentrations for practical reasons. How should these be plotted?
    '''

    skip_keys = ['time']

    # get settings
    max_rows = settings.get('max_rows', 25)
    remove_zeros = settings.get('remove_zeros', False)
    remove_flat = settings.get('remove_flat', False)
    skip_ports = settings.get('skip_ports', [])
    overlay = settings.get('overlay', {})
    show_state = settings.get('show_state', [])
    top_ports = list(overlay.values())
    bottom_ports = list(overlay.keys())

    ports = [port for port in timeseries.keys() if port not in skip_keys + skip_ports]
    time_vec = timeseries['time']

    # remove selected states
    # TODO -- plot removed_states as text
    removed_states = []
    if remove_flat:
        # find series with all the same value
        for port in ports:
            for state_id, series in timeseries[port].items():
                if series.count(series[0]) == len(series):
                    removed_states.append((port, state_id))
    elif remove_zeros:
        # find series with all zeros
        for port in ports:
            for state_id, series in timeseries[port].items():
                if all(v == 0 for v in series):
                    removed_states.append((port, state_id))

    # if specified in show_state, keep in timeseries
    for port_state in show_state:
        if port_state in removed_states:
            removed_states.remove(port_state)

    # remove from timeseries
    for (port, state_id) in removed_states:
        del timeseries[port][state_id]

    # get the number of states in each port
    n_data = [len(timeseries[key]) for key in ports if key not in top_ports]
    if 0 in n_data:
        n_data.remove(0)

    # limit number of rows to max_rows by adding new columns
    columns = []
    for n_states in n_data:
        new_cols = n_states / max_rows
        if new_cols > 1:
            for col in range(int(new_cols)):
                columns.append(max_rows)

            mod_states = n_states % max_rows
            if mod_states > 0:
                columns.append(mod_states)
        else:
            columns.append(n_states)
    n_cols = len(columns)
    n_rows = max(columns)

    # make figure and plot
    fig = plt.figure(figsize=(n_cols * 6, n_rows * 1.5))
    grid = plt.GridSpec(n_rows, n_cols)

    row_idx = 0
    col_idx = 0
    for port in ports:
        top_timeseries = {}
        if port in bottom_ports:
            # get overlay
            top_port = overlay[port]
            top_timeseries = timeseries[top_port]
        elif port in top_ports + skip_ports:
            # don't give this row its own plot
            continue

        for state_id, series in sorted(timeseries[port].items()):
            ax = fig.add_subplot(grid[row_idx, col_idx])  # grid is (row, column)

            # check if series is a list of ints or floats
            # TODO -- plot non-numeric states as well (in particular dicts)
            if not all(isinstance(state, (int, float)) for state in series):
                break

            # plot line at zero if series crosses the zero line
            if any(x == 0.0 for x in series) or (any(x < 0.0 for x in series) and any(x > 0.0 for x in series)):
                zero_line = [0 for t in time_vec]
                ax.plot(time_vec, zero_line, 'k--')

            if (port, state_id) in show_state:
                ax.plot(time_vec, series, 'indigo', linewidth=2)
            else:
                ax.plot(time_vec, series)

            # overlay
            if state_id in top_timeseries.keys():
                ax.plot(time_vec, top_timeseries[state_id], 'm', label=top_port)
                ax.legend()

            ax.title.set_text(str(port) + ': ' + str(state_id))
            ax.title.set_fontsize(16)

            if row_idx == columns[col_idx]-1:
                # if last row of column
                set_axes(ax, True)
                ax.set_xlabel('time')
                row_idx = 0
                col_idx += 1
            else:
                set_axes(ax)
                row_idx += 1

    # save figure
    fig_path = os.path.join(out_dir, 'simulation')
    plt.subplots_adjust(wspace=0.3, hspace=0.5)
    plt.savefig(fig_path, bbox_inches='tight')


# TESTS


class ToyLinearGrowthDeathProcess(Process):

    GROWTH_RATE = 1.0
    THRESHOLD = 5.0

    def __init__(self, initial_parameters={}):
        ports = {
            'compartment': ['processes'],
            'global': ['mass'],
        }
        super(ToyLinearGrowthDeathProcess, self).__init__(
            ports, initial_parameters)

    def default_settings(self):
        default_settings = {
            'state': {
                'global': {
                    'mass': 0.0
                }
            },
        }
        return default_settings

    def next_update(self, timestep, states):
        mass = states['global']['mass']
        mass_grown = (
            ToyLinearGrowthDeathProcess.GROWTH_RATE * timestep)
        update = {
            'global': {'mass': mass_grown},
        }
        if mass > ToyLinearGrowthDeathProcess.THRESHOLD:
            update['compartment'] = {
                'processes': [],
            }
        return update


class TestSimulateProcess:

    def test_compartment_state_port(self):
        '''Check that compartment state ports are handled'''
        process = ToyLinearGrowthDeathProcess()
        settings = {
            'compartment_state_port': 'compartment',
        }
        saved_total_states = simulate_process(
            process, settings)
        timeseries = convert_to_timeseries(
            saved_total_states)
        expected_masses = [
            # Mass stops increasing the iteration after mass > 5 because
            # cell dies
            0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 7.0, 7.0, 7.0]
        masses = timeseries['global']['mass']
        assert masses == expected_masses