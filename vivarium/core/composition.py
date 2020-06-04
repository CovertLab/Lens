import copy
import csv
import os
import io

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

from vivarium.core.experiment import (
    Experiment,
    update_in,
    generate_derivers,
)
from vivarium.core.process import Process, Deriver
from vivarium.core.experiment import Compartment as TreeCompartment
from vivarium.core import emitter as emit
from vivarium.library.dict_utils import (
    deep_merge,
    deep_merge_check,
    flatten_timeseries,
)
from vivarium.library.units import units

# processes
from vivarium.processes.derive_globals import AVOGADRO
from vivarium.processes.timeline import (
    TimelineCompartment,
    TimelineProcess
)
from vivarium.processes.homogeneous_environment import HomogeneousEnvironment

REFERENCE_DATA_DIR = os.path.join('vivarium', 'reference_data')
TEST_OUT_DIR = os.path.join('out', 'tests')
PROCESS_OUT_DIR = os.path.join('out', 'processes')
COMPARTMENT_OUT_DIR = os.path.join('out', 'compartments')
EXPERIMENT_OUT_DIR = os.path.join('out', 'experiments')


# loading functions
def make_agents(agent_ids, compartment, config=None):

    if config is None:
        config = {}

    processes = {}
    topology = {}
    for agent_id in agent_ids:
        agent_config = copy.deepcopy(config)
        agent = compartment.generate(dict(
            agent_config,
            agent_id=agent_id))

        # save processes and topology
        processes[agent_id] = agent['processes']
        topology[agent_id] = agent['topology']

    return {
        'processes': processes,
        'topology': topology}

def process_in_experiment(process, settings={}):
    initial_state = settings.get('initial_state', {})
    emitter = settings.get('emitter', {'type': 'timeseries'})
    timeline = settings.get('timeline', [])
    environment = settings.get('environment', {})

    processes = {'process': process}
    topology = {
        'process': {
            port: (port,) for port in process.ports_schema().keys()}}

    if timeline:
        timeline_process = TimelineProcess({'timeline': timeline})
        processes.update({'timeline_process': timeline_process})
        topology.update({
            'timeline_process': {
                port: (port,) for port in timeline_process.ports}})

    if environment:
        environment_process = HomogeneousEnvironment(environment)
        processes.update({'environment_process': environment_process})
        topology.update({
            'environment_process': {
                port: (port,) for port in environment_process.ports}})

    # add derivers
    derivers = generate_derivers(processes, topology)
    processes = deep_merge(processes, derivers['processes'])
    topology = deep_merge(topology, derivers['topology'])

    return Experiment({
        'processes': processes,
        'topology': topology,
        'emitter': emitter,
        'initial_state': initial_state})

def compartment_in_experiment(compartment, settings={}):
    compartment_config = settings.get('compartment', {})
    emitter = settings.get('emitter', {'type': 'timeseries'})
    environment = settings.get('environment', {})
    outer_path = settings.get('outer_path', tuple())

    network = compartment.generate(compartment_config, outer_path)
    processes = network['processes']
    topology = network['topology']

    if environment:
        # TODO -- make a compartment_in_environment
        environment_process = HomogeneousEnvironment(environment)

        update_in(
            processes,
            outer_path,
            lambda existing: deep_merge(
                existing,
                {'environment_process': environment_process}))

        update_in(
            topology,
            outer_path,
            lambda existing: deep_merge(
                existing,
                {'environment_process': {
                    port: (port,) for port in environment_process.ports}}))

    return Experiment({
        'processes': processes,
        'topology': topology,
        'emitter': emitter,
        'initial_state': settings.get('initial_state', {})})

def add_timeline_to_compartment(compartment, settings={}):
    timeline = settings['timeline']
    path = settings['path']
    return TimelineCompartment({
        'timeline': timeline,
        'compartment': compartment,
        'path': path})


# simulation functions
def simulate_process(process, settings={}):
    experiment = process_in_experiment(process)
    return simulate_experiment(experiment, settings)

def simulate_process_in_experiment(process, settings={}):
    experiment = process_in_experiment(process, settings)
    return simulate_experiment(experiment, settings)

def simulate_compartment_in_experiment(compartment, settings={}):
    experiment = compartment_in_experiment(compartment, settings)
    return simulate_experiment(experiment, settings)

def simulate_experiment(experiment, settings={}):
    '''
    run an experiment simulation
        Requires:
        - a configured experiment

    Returns:
        - a timeseries of variables from all ports.
        - if 'return_raw_data' is True, it returns the raw data instead
    '''
    timestep = settings.get('timestep', 1)
    total_time = settings.get('total_time', 10)
    return_raw_data = settings.get('return_raw_data', False)

    if 'timeline' in settings:
        total_time = settings['timeline'][-1][0]

    # run simulation
    experiment.update_interval(total_time, timestep)

    if return_raw_data:
        return experiment.emitter.get_data()
    else:
        return experiment.emitter.get_timeseries()


# plotting functions
def plot_compartment_topology(compartment, settings, out_dir='out', filename='topology'):
    """
    Make a plot of the topology
     - compartment: a compartment
     - settings (dict): 'network_layout' can be 'bipartite' or 'process_layers'
    """
    store_rgb = [x/255 for x in [239,131,148]]
    process_rgb = [x / 255 for x in [249, 204, 86]]
    node_size = 2500
    node_distance = 1
    layer_distance = 10

    network = compartment.generate({})
    topology = network['topology']
    processes = network['processes']

    # get figure settings
    show_ports = settings.get('show_ports', True)

    # make graph from topology
    G = nx.Graph()
    process_nodes = []
    store_nodes = []
    edges = {}
    for process_id, connections in topology.items():
        process_nodes.append(process_id)
        G.add_node(process_id)

        for port, store_id in connections.items():
            if store_id not in store_nodes:
                store_nodes.append(store_id)
            if store_id not in list(G.nodes):
                G.add_node(store_id)

            edge = (process_id, store_id)
            edges[edge]= port

            G.add_edge(process_id, store_id)

    # are there overlapping names?
    overlap = [name for name in process_nodes if name in store_nodes]
    if overlap:
        print('{} shared by processes and stores'.format(overlap))


    # get positions
    pos = {}
    n_rows = max(len(process_nodes), len(store_nodes))
    plt.figure(3, figsize=(12, 1.2 * n_rows))

    for idx, node_id in enumerate(process_nodes, 1):
        pos[node_id] = np.array([-1, -idx*node_distance])
    for idx, node_id in enumerate(store_nodes, 1):
        pos[node_id] = np.array([1, -idx*node_distance])


    # plot
    nx.draw_networkx_nodes(G, pos,
                           nodelist=process_nodes,
                           with_labels=True,
                           node_color=process_rgb,
                           node_size=node_size,
                           node_shape='s')
    nx.draw_networkx_nodes(G, pos,
                           nodelist=store_nodes,
                           with_labels=True,
                           node_color=store_rgb,
                           node_size=node_size,
                           node_shape='o')

    # edges
    colors = list(range(1,len(edges)+1))
    nx.draw_networkx_edges(G, pos,
                           edge_color=colors,
                           width=1.5)

    # labels
    nx.draw_networkx_labels(G, pos,
                            font_size=8,
                            )
    if show_ports:
        nx.draw_networkx_edge_labels(G, pos,
                                 edge_labels=edges,
                                 font_size=6,
                                 label_pos=0.85)

    # save figure
    fig_path = os.path.join(out_dir, filename)
    plt.figure(3, figsize=(12, 12))
    plt.axis('off')
    plt.savefig(fig_path, bbox_inches='tight')

    plt.close()

def set_axes(ax, show_xaxis=False):
    ax.ticklabel_format(style='sci', axis='y', scilimits=(-5,5))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(right=False, top=False)
    if not show_xaxis:
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(bottom=False, labelbottom=False)

def get_plot_columns(ports, settings={}):
    top_ports = settings.get('top_ports', [])
    max_rows = settings.get('max_rows', 30)

    n_data = [len(ports[key]) for key in ports if key not in top_ports]
    if 0 in n_data:
        n_data.remove(0)

    # limit number of rows to max_rows by adding new columns
    columns = []
    for n_states in n_data:
        if n_states == 0:
            continue
        new_cols = n_states / max_rows
        if new_cols > 1:
            for col in range(int(new_cols)):
                columns.append(max_rows)
            mod_states = n_states % max_rows
            if mod_states > 0:
                columns.append(mod_states)
        else:
            columns.append(n_states)

    return columns

def plot_simulation_output(timeseries_raw, settings={}, out_dir='out', filename='simulation'):
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

    plot_fontsize = 8
    plt.rc('font', size=plot_fontsize)
    plt.rc('axes', titlesize=plot_fontsize)

    skip_keys = ['time']

    # get settings
    timeseries = copy.deepcopy(timeseries_raw)
    max_rows = settings.get('max_rows', 25)
    remove_zeros = settings.get('remove_zeros', True)
    remove_flat = settings.get('remove_flat', False)
    skip_ports = settings.get('skip_ports', [])
    overlay = settings.get('overlay', {})
    top_ports = list(overlay.values())
    bottom_ports = list(overlay.keys())

    time_vec = timeseries.pop('time')

    ports = {}
    for port_id, states in timeseries.items():
        if port_id in skip_keys + skip_ports:
            continue
        if port_id not in ports and len(states) != 0:
            ports[port_id] = []
        for state_id in list(states.keys()):
            if state_id not in ports[port_id]:
                ports[port_id].append(state_id)

    # remove selected states
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

    # remove from timeseries
    for (port, state_id) in removed_states:
        del timeseries[port][state_id]

    # limit number of rows to max_rows by adding new columns
    column_settings = {
        'top_ports': top_ports,
        'max_rows': max_rows}
    columns = get_plot_columns(timeseries, column_settings)

    # make figure and plot
    n_cols = len(columns)
    n_rows = max(columns)
    fig = plt.figure(figsize=(n_cols * 3, n_rows * 1))
    grid = plt.GridSpec(n_rows, n_cols)
    row_idx = 0
    col_idx = 0
    for port in ports:
        top_timeseries = {}

        # set up overlay
        if port in bottom_ports:
            top_port = overlay[port]
            top_timeseries = timeseries[top_port]
        elif port in top_ports + skip_ports:
            # don't give this row its own plot
            continue

        for state_id, series in sorted(timeseries[port].items()):
            ax = fig.add_subplot(grid[row_idx, col_idx])  # grid is (row, column)

            # check if series is a list of ints or floats
            if not all(isinstance(state, (int, float, np.int64, np.int32)) for state in series):
                ax.title.set_text(str(port) + ': ' + str(state_id) + ' (non numeric)')
            else:
                # plot line at zero if series crosses the zero line
                if any(x == 0.0 for x in series) or (any(x < 0.0 for x in series) and any(x > 0.0 for x in series)):
                    zero_line = [0 for t in time_vec]
                    ax.plot(time_vec, zero_line, 'k--')
                ax.plot(time_vec, series)

                # overlay
                if state_id in top_timeseries.keys():
                    ax.plot(time_vec, top_timeseries[state_id], 'm', label=top_port)
                    ax.legend()
                ax.title.set_text(str(port) + ': ' + str(state_id))

            if row_idx == columns[col_idx]-1:
                # if last row of column
                set_axes(ax, True)
                ax.set_xlabel('time (s)')
                row_idx = 0
                col_idx += 1
            else:
                set_axes(ax)
                row_idx += 1

    # save figure
    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.8, hspace=0.8)
    plt.savefig(fig_path, bbox_inches='tight')

def plot_agent_data(data, settings={}, out_dir='out', filename='agents'):
    '''
    Make a plot of all agent data
    TODO -- add agent color
    '''

    agents_key = settings.get('agents_key', 'cells')

    time_vec = list(data.keys())
    agents_timeseries = agent_timeseries_from_data(data, agents_key)

    # assume the initial agents have the same port schema as all subsequent agents
    initial_agents = data[time_vec[0]][agents_key]
    ports = {}
    for agent_id, state_data in initial_agents.items():
        for port_id, states in state_data.items():
            if port_id not in ports:
                ports[port_id] = []
            for state_id in states.keys():
                if state_id not in ports[port_id]:
                    ports[port_id].append(state_id)

    # get the column sizes
    columns = get_plot_columns(ports)

    # make figure
    n_rows = max(columns)
    n_cols = len(columns)
    fig = plt.figure(figsize=(4 * n_cols, 2 * n_rows))
    grid = plt.GridSpec(n_rows, n_cols, wspace=0.4, hspace=1.5)

    # set up the axes
    port_axes = {}
    for port_idx, (port_id, states) in enumerate(ports.items()):
        n_states = columns[port_idx]
        for state_idx, state_id in enumerate(states):
            ax = fig.add_subplot(grid[state_idx, port_idx])
            ax.title.set_text(str(port_id) + ': ' + str(state_id))
            ax.title.set_fontsize(16)
            if state_idx is not n_states-1:
                set_axes(ax)
            else:
                # if last state in this port, add time ticks
                set_axes(ax, True)
                ax.set_xlabel('time (s)')

            # save axis
            port_axes[(port_id, state_id)] = ax

    # plot the agents
    plotted_agents = []
    for time_idx, (time, time_data) in enumerate(data.items()):
        agents = time_data[agents_key]
        for agent_id, agent_data in agents.items():
            if agent_id not in plotted_agents:
                plotted_agents.append(agent_id)
                agent_ts = agents_timeseries[agent_id]
                for port_id, state_ts in agent_ts.items():
                    for state_id, state in state_ts.items():
                        if not isinstance(state[0], (float, int)):
                            continue
                        n_times = len(state)
                        plot_times = time_vec[time_idx:time_idx+n_times]
                        ax = port_axes[(port_id, state_id)]
                        ax.plot(plot_times, state)

    # save figure
    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.2, hspace=0.2)
    plt.savefig(fig_path, bbox_inches='tight')


# timeseries functions
def agent_timeseries_from_data(data, agents_key='cells'):
    timeseries = {}
    for time, all_states in data.items():
        agent_data = all_states[agents_key]
        for agent_id, ports in agent_data.items():
            if agent_id not in timeseries:
                timeseries[agent_id] = {}
            for port_id, states in ports.items():
                if port_id not in timeseries[agent_id]:
                    timeseries[agent_id][port_id] = {}
                for state_id, state in states.items():
                    if state_id not in timeseries[agent_id][port_id]:
                        timeseries[agent_id][port_id][state_id] = []
                    timeseries[agent_id][port_id][state_id].append(state)
    return timeseries

def save_timeseries(timeseries, out_dir='out'):
    '''Save a timeseries as a CSV in out_dir'''
    flattened = flatten_timeseries(timeseries)
    rows = np.transpose(list(flattened.values())).tolist()
    with open(os.path.join(out_dir, 'simulation_data.csv'), 'w') as f:
        writer = csv.writer(f)
        writer.writerow(flattened.keys())
        writer.writerows(rows)

def load_timeseries(path_to_csv):
    '''Load a timeseries saved as a CSV using save_timeseries.

    The timeseries is returned in flattened form.
    '''
    with io.open(path_to_csv, 'r', newline='') as f:
        reader = csv.DictReader(f)
        timeseries = {}
        for row in reader:
            for header, elem in row.items():
                timeseries.setdefault(header, []).append(float(elem))
    return timeseries

def timeseries_to_ndarray(timeseries, keys=None):
    if keys is None:
        keys = timeseries.keys()
    filtered = {key: timeseries[key] for key in keys}
    array = np.array(list(filtered.values()))
    return array

def _prepare_timeseries_for_comparison(
    timeseries1, timeseries2, keys=None,
    required_frac_checked=0.9,
):
    '''Prepare two timeseries for comparison

    Arguments:
        timeseries1: One timeseries. Must be flattened and include times
            under the 'time' key.
        timeseries2: The other timeseries. Same requirements as
            timeseries1.
        keys: Keys of the timeseries whose values will be checked for
            correlation. If not specified, all keys present in both
            timeseries are used.
        required_frac_checked: The required fraction of timepoints in a
            timeseries that must be checked. If this requirement is not
            satisfied, which might occur if the two timeseries share few
            timepoints, the test wll fail.

    Returns:
        A tuple of an ndarray for each of the two timeseries and a list of
        the keys for the rows of the arrays. Each ndarray has a row for
        each key, in the order of keys. The ndarrays have only the
        columns corresponding to the timepoints common to both
        timeseries.

    Raises:
        AssertionError: If a correlation is strictly below the
            threshold or if too few timepoints are common to both
            timeseries.
    '''
    if 'time' not in timeseries1 or 'time' not in timeseries2:
        raise AssertionError('Both timeseries must have key "time"')
    if keys is None:
        keys = set(timeseries1.keys()) & set(timeseries2.keys())
    else:
        if 'time' not in keys:
            keys.append('time')
    keys = list(keys)
    time_index = keys.index('time')
    shared_times = set(timeseries1['time']) & set(timeseries2['time'])
    frac_timepoints_checked = (
        len(shared_times)
        / min(len(timeseries1), len(timeseries2))
    )
    if frac_timepoints_checked < required_frac_checked:
        raise AssertionError(
            'The timeseries share too few timepoints: '
            '{} < {}'.format(
                frac_timepoints_checked, required_frac_checked)
        )
    array1 = timeseries_to_ndarray(timeseries1, keys)
    array2 = timeseries_to_ndarray(timeseries2, keys)
    shared_times_mask1 = np.isin(array1[time_index], list(shared_times))
    shared_times_mask2 = np.isin(array2[time_index], list(shared_times))
    return (
        array1[:, shared_times_mask1],
        array2[:, shared_times_mask2],
        keys,
    )

def assert_timeseries_correlated(
    timeseries1, timeseries2, keys=None,
    default_threshold=(1 - 1e-10), thresholds={},
    required_frac_checked=0.9,
):
    '''Check that two timeseries are correlated.

    Uses a Pearson correlation coefficient. Only the data from
    timepoints common to both timeseries are compared.

    Arguments:
        timeseries1: One timeseries. Must be flattened and include times
            under the 'time' key.
        timeseries2: The other timeseries. Same requirements as
            timeseries1.
        keys: Keys of the timeseries whose values will be checked for
            correlation. If not specified, all keys present in both
            timeseries are used.
        default_threshold: The threshold correlation coefficient to use
            when a threshold is not specified in thresholds.
        thresholds: Dictionary of key-value pairs where the key is a key
            in both timeseries and the value is the threshold
            correlation coefficient to use when checking that key
        required_frac_checked: The required fraction of timepoints in a
            timeseries that must be checked. If this requirement is not
            satisfied, which might occur if the two timeseries share few
            timepoints, the test wll fail.

    Raises:
        AssertionError: If a correlation is strictly below the
            threshold or if too few timepoints are common to both
            timeseries.
    '''
    array1, array2, keys = _prepare_timeseries_for_comparison(
        timeseries1, timeseries2, keys, required_frac_checked)
    for index, key in enumerate(keys):
        corrcoef = np.corrcoef(
            array1[index],
            array2[index],
        )[0][1]
        threshold = thresholds.get(key, default_threshold)
        if corrcoef < threshold:
            raise AssertionError(
                'The correlation coefficient for '
                '{} is too small: {} < {}'.format(
                    key, corrcoef, threshold)
            )

def assert_timeseries_close(
    timeseries1, timeseries2, keys=None,
    default_tolerance=(1 - 1e-10), tolerances={},
    required_frac_checked=0.9,
):
    '''Check that two timeseries are similar.

    Ensures that each pair of data points between the two timeseries are
    within a tolerance of each other, after filtering out timepoints not
    common to both timeseries.

    Arguments:
        timeseries1: One timeseries. Must be flattened and include times
            under the 'time' key.
        timeseries2: The other timeseries. Same requirements as
            timeseries1.
        keys: Keys of the timeseries whose values will be checked for
            correlation. If not specified, all keys present in both
            timeseries are used.
        default_tolerance: The tolerance to use when not specified in
            tolerances.
        tolerances: Dictionary of key-value pairs where the key is a key
            in both timeseries and the value is the tolerance to use
            when checking that key.
        required_frac_checked: The required fraction of timepoints in a
            timeseries that must be checked. If this requirement is not
            satisfied, which might occur if the two timeseries share few
            timepoints, the test wll fail.

    Raises:
        AssertionError: If a pair of data points have a difference
            strictly above the tolerance threshold or if too few
            timepoints are common to both timeseries.
    '''
    array1, array2, keys = _prepare_timeseries_for_comparison(
        timeseries1, timeseries2, keys, required_frac_checked)
    for index, key in enumerate(keys):
        tolerance = tolerances.get(key, default_tolerance)
        if not np.allclose(
            array1[index], array2[index], atol=tolerance
        ):
            raise AssertionError(
                'The data for {} differed by more than {}'.format(
                    key, tolerance)
            )


# TESTS
class ToyLinearGrowthDeathProcess(Process):

    GROWTH_RATE = 1.0
    THRESHOLD = 6.0

    def __init__(self, initial_parameters={}):
        self.targets = initial_parameters.get('targets')
        ports = {
            'global': ['mass'],
        }
        super(ToyLinearGrowthDeathProcess, self).__init__(
            ports, initial_parameters)

    def ports_schema(self):
        schema = {
            'global': {
                'mass': {
                    '_default': 0.0,
                    '_emit': True}}}

        schema['global'].update({
            target: {
                '_default': None}
            for target in self.targets})
        return schema


    def next_update(self, timestep, states):
        mass = states['global']['mass']
        mass_grown = (
            ToyLinearGrowthDeathProcess.GROWTH_RATE * timestep)
        update = {
            'global': {'mass': mass_grown},
        }
        if mass > ToyLinearGrowthDeathProcess.THRESHOLD:
            update['global'] = {
                '_delete': [(target,) for target in self.targets]}

        return update

class TestSimulateProcess:

    def test_process_deletion(self):
        '''Check that processes are successfully deleted'''
        process = ToyLinearGrowthDeathProcess({'targets': ['process']})
        settings = {}

        timeseries = simulate_process(process, settings)
        expected_masses = [
            # Mass stops increasing the iteration after mass > 5 because
            # cell dies
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 7.0, 7.0, 7.0]
        masses = timeseries['global']['mass']
        assert masses == expected_masses


# toy processes
class ToyMetabolism(Process):
    def __init__(self, initial_parameters={}):
        ports = {'pool': ['GLC', 'MASS']}
        parameters = {'mass_conversion_rate': 1}
        parameters.update(initial_parameters)

        super(ToyMetabolism, self).__init__(ports, parameters)

    def ports_schema(self):
        return {
            port_id: {
                key: {
                    '_default': 0.0,
                    '_emit': True}
                for key in keys}
            for port_id, keys in self.ports.items()}

    def next_update(self, timestep, states):
        update = {}
        glucose_required = timestep / self.parameters['mass_conversion_rate']
        if states['pool']['GLC'] >= glucose_required:
            update = {
                'pool': {
                    'GLC': -2,
                    'MASS': 1}}

        return update

class ToyTransport(Process):
    def __init__(self, initial_parameters={}):
        ports = {
            'external': ['GLC'],
            'internal': ['GLC']}
        parameters = {'intake_rate': 2}
        parameters.update(initial_parameters)

        super(ToyTransport, self).__init__(ports, parameters)

    def ports_schema(self):
        return {
            port_id: {
                key: {
                    '_default': 0.0,
                    '_emit': True}
                for key in keys}
            for port_id, keys in self.ports.items()}

    def next_update(self, timestep, states):
        update = {}
        intake = timestep * self.parameters['intake_rate']
        if states['external']['GLC'] >= intake:
            update = {
                'external': {'GLC': -2, 'MASS': 1},
                'internal': {'GLC': 2}}

        return update

class ToyDeriveVolume(Deriver):
    def __init__(self, initial_parameters={}):
        ports = {
            'compartment': ['MASS', 'DENSITY', 'VOLUME']}
        parameters = {}

        super(ToyDeriveVolume, self).__init__(ports, parameters)

    def ports_schema(self):
        return {
            port_id: {
                key: {
                    '_updater': 'set' if key == 'VOLUME' else 'accumulate',
                    '_default': 0.0,
                    '_emit': True}
                for key in keys}
            for port_id, keys in self.ports.items()}

    def next_update(self, timestep, states):
        volume = states['compartment']['MASS'] / states['compartment']['DENSITY']
        update = {
            'compartment': {'VOLUME': volume}}

        return update

class ToyDeath(Process):
    def __init__(self, initial_parameters={}):
        self.targets = initial_parameters.get('targets', [])
        ports = {
            'compartment': ['VOLUME'],
            'global': self.targets}
        super(ToyDeath, self).__init__(ports, {})

    def ports_schema(self):
        return {
            'compartment': {
                'VOLUME': {
                    '_default': 0.0,
                    '_emit': True}},
            'global': {
                target: {
                    '_default': None}
                for target in self.targets}}

    def next_update(self, timestep, states):
        volume = states['compartment']['VOLUME']
        update = {}

        if volume > 1.0:
            # kill the cell
            update = {
                'global': {
                    '_delete': [
                        (target,)
                        for target in self.targets]}}

        return update

class ToyCompartment(TreeCompartment):
    '''
    a toy compartment for testing

    '''
    def __init__(self, config):
        self.config = config

    def generate_processes(self, config):
        return {
            'metabolism': ToyMetabolism(
                {'mass_conversion_rate': 0.5}), # example of overriding default parameters
            'transport': ToyTransport(),
            'death': ToyDeath({'targets': [
                'metabolism',
                'transport']}),
            'external_volume': ToyDeriveVolume(),
            'internal_volume': ToyDeriveVolume()
        }

    def generate_topology(self, config):
        return{
            'metabolism': {
                'pool': ('cytoplasm',)},
            'transport': {
                'external': ('periplasm',),
                'internal': ('cytoplasm',)},
            'death': {
                'global': tuple(),
                'compartment': ('cytoplasm',)},
            'external_volume': {
                'compartment': ('periplasm',)},
            'internal_volume': {
                'compartment': ('cytoplasm',)}}


def test_compartment():
    toy_compartment = ToyCompartment({})
    settings = {
        'timestep': 1,
        'total_time': 10,
        'initial_state': {
            'periplasm': {
                'GLC': 20,
                'MASS': 100,
                'DENSITY': 10},
            'cytoplasm': {
                'GLC': 0,
                'MASS': 3,
                'DENSITY': 10}}}
    data = simulate_compartment_in_experiment(toy_compartment, settings)

if __name__ == '__main__':
    timeseries = test_compartment()
    print(timeseries)
