from __future__ import absolute_import, division, print_function

import copy
import uuid
import random
import datetime

import numpy as np
import logging as log

import pprint
pretty=pprint.PrettyPrinter(indent=2)
def pp(x):
    pretty.pprint(x)

from vivarium.core.process import (
    Process,
    divider_library)
from vivarium.utils.units import Quantity
from vivarium.utils.dict_utils import merge_dicts, deep_merge, deep_merge_check
from vivarium.core.emitter import get_emitter

# processes
from vivarium.processes.derive_globals import DeriveGlobals, AVOGADRO
from vivarium.processes.derive_counts import DeriveCounts
from vivarium.processes.derive_concentrations import DeriveConcentrations
from vivarium.processes.tree_mass import TreeMass

deriver_library = {
    'mmol_to_counts': DeriveCounts,
    'counts_to_mmol': DeriveConcentrations,
    'mass': TreeMass,
    'globals': DeriveGlobals,
}


def key_for_value(d, looking):
    found = None
    for key, value in d.items():
        if looking == value:
            found = key
            break
    return found


def get_in(d, path):
    if path:
        head = path[0]
        if head in d:
            return get_in(d[head], path[1:])
    else:
        return d


def assoc_in(d, path, value):
    if path:
        return dict(d, **{path[0]: assoc_in(d.get(path[0], {}), path[1:], value)})
    else:
        return value


def assoc_path(d, path, value):
    if path:
        head = path[0]
        if len(path) == 1:
            d[head] = value
        else:
            if head not in d:
                d[head] = {}
            assoc_path(d[head], path[1:], value)
    else:
        value


def update_in(d, path, f):
    if path:
        head = path[0]
        if len(path) == 1:
            d[head] = f(d.get(head, None))
        else:
            if not head in d:
                d[head] = {}
            update_in(d[head], path[1:], f)


def dissoc(d, removing):
    return {
        key: value
        for key, value in d.items()
        if not key in removing}
        

def select_keys(d, keys):
    return {
        key: d.get(key)
        for key in keys}


DEFAULT_TIMESTEP = 1.0

INFINITY = float('inf')
VERBOSE = False

class topologyError(Exception):
    pass

def npize(d):
    ''' Turn a dict into an ordered set of keys and values. '''

    ordered = [[key, value] for key, value in d.items()]
    keys = [key for key, _ in ordered]
    values = np.array([value for _, value in ordered], np.float64)

    return keys, values

## updater functions
# these function take in a variable key, the entire store's dict,
# the variable's current value, the variable's current update,
# and returns a new value, and other updates
def update_accumulate(current_value, new_value):
    return current_value + new_value

def update_set(current_value, new_value):
    return new_value

def update_merge(current_value, new_value):
    # merge dicts, with new_value replacing any shared keys with current_value
    update = current_value.copy()
    for k, v in current_value.items():
        new = new_value.get(k)
        if isinstance(new, dict):
            update[k] = deep_merge(dict(v), new)
        else:
            update[k] = new
    return update

updater_library = {
    'accumulate': update_accumulate,
    'set': update_set,
    'merge': update_merge}



def schema_for(port, keys, initial_state, default=0.0, updater='accumulate'):
    return {
        key: {
            '_default': initial_state.get(
                port, {}).get(key, default),
            '_updater': updater}
        for key in keys}

def always_true(x):
    return True

def identity(y):
    return y


class Store(object):
    schema_keys = set([
        '_default',
        '_updater',
        '_divider',
        '_value',
        '_properties',
        '_emit',
        # '_units',
    ])

    def __init__(self, config, outer=None, source=None):
        self.outer = outer
        self.inner = {}
        self.subschema = {}
        self.properties = {}
        self.default = None
        self.updater = None
        self.value = None
        self.units = None
        self.divider = None
        self.emit = False
        self.sources = {}

        self.apply_config(config, source)

    def check_default(self, new_default):
        if self.default is not None and new_default != self.default:
            print('_default schema conflict: {} and {}. selecting {}'.format(
                self.default, new_default, new_default))
        return new_default

    def check_value(self, new_value):
        if self.value is not None and new_value != self.value:
            raise Exception('_value schema conflict: {} and {}'.format(new_value, self.value))
        return new_value

    def apply_subschema_config(self, config, subschema_key):
        self.subschema = deep_merge(
            self.subschema,
            config.get(subschema_key, {}))
        return {
            key: value
            for key, value in config.items()
            if key != subschema_key}

    def apply_config(self, config, source=None):
        if '*' in config:
            config = self.apply_subschema_config(config, '*')

        if '_subschema' in config:
            if source:
                self.sources[source] = config['_subschema']
            config = self.apply_subschema_config(config, '_subschema')

        if self.schema_keys & config.keys():
            # self.units = config.get('_units', self.units)
            if '_default' in config:
                self.default = self.check_default(config.get('_default'))
            if '_value' in config:
                self.value = self.check_value(config.get('_value'))
                if isinstance(self.value, Quantity):
                    self.units = self.value.units

            self.updater = config.get('_updater', self.updater or 'accumulate')
            if isinstance(self.updater, str):
                self.updater = updater_library[self.updater]
            self.divider = config.get('_divider', self.divider)
            if isinstance(self.divider, str):
                self.divider = divider_library[self.divider]
            if isinstance(self.divider, dict) and isinstance(self.divider['divider'], str):
                self.divider['divider'] = divider_library[self.divider['divider']]

            self.properties = deep_merge(
                self.properties,
                config.get('_properties', {}))

            self.emit = config.get('_emit', self.emit)

            if source:
                self.sources[source] = config

        else:
            self.value = None

            for key, child in config.items():
                if not key in self.inner:
                    self.inner[key] = Store(child, outer=self, source=source)
                else:
                    self.inner[key].apply_config(child, source=source)

    def get_updater(self, update):
        updater = self.updater
        if '_updater' in update:
            updater = update['_updater']
            if isinstance(updater, str):
                updater = updater_library[updater]
        return updater

    def get_config(self, sources=False):
        config = {}
        if self.properties:
            config['_properties'] = self.properties
        if self.subschema:
            config['_subschema'] = self.subschema
        if sources and self.sources:
            config['_sources'] = self.sources

        if self.inner:
            child_config = {
                key: child.get_config(sources)
                for key, child in self.inner.items()}
            config.update(child_config)
        else:
            config.update({
                '_default': self.default,
                '_value': self.value})
            if self.updater:
                config['_updater'] = self.updater
            if self.divider:
                config['_divider'] = self.divider
            if self.units:
                config['_units'] = self.units
            if self.emit:
                config['_emit'] = self.emit

        return config

    def top(self):
        if self.outer:
            return self.outer.top()
        else:
            return self

    def path_for(self):
        if self.outer:
            key = key_for_value(self.outer.inner, self)
            above = self.outer.path_for()
            return above + (key,)
        else:
            return tuple()

    def get_value(self, condition=None, f=None):
        if self.inner:
            if condition is None:
                condition = always_true

            if f is None:
                f = identity

            return {
                key: f(child.get_value(condition, f))
                for key, child in self.inner.items()
                if condition(child)}
        else:
            if self.subschema:
                return {}
            # elif self.units:
            #     return self.value * self.units
            else:
                return self.value

    def get_path(self, path):
        if path:
            step = path[0]
            if step == '..':
                child = self.outer
            else:
                child = self.inner.get(step)

            if child:
                return child.get_path(path[1:])
            else:
                # TODO: more handling for bad paths?
                return None
        else:
            return self

    def get_paths(self, paths):
        return {
            key: self.get_path(path)
            for key, path in paths.items()}

    def get_values(self, paths):
        return {
            key: self.get_in(path)
            for key, path in paths.items()}

    def get_in(self, path):
        return self.get_path(path).get_value()

    def get_template(self, template):
        '''
        Pass in a template dict with None for each value you want to
        retrieve from the tree!
        '''

        state = {}
        for key, value in template.items():
            child = self.inner[key]
            if value is None:
                state[key] = child.get_value()
            else:
                state[key] = child.get_template(value)
        return state

    def emit_data(self):
        data = {}
        if self.inner:
            for key, child in self.inner.items():
                child_data = child.emit_data()
                if child_data is not None or child_data == 0:
                    data[key] = child_data
            return data
        else:
            if self.emit:
                if isinstance(self.value, Process):
                    return self.value.pull_data()
                else:
                    if self.units:
                        return self.value.to(self.units).magnitude
                    else:
                        return self.value

    def delete_path(self, path):
        if not path:
            self.inner = {}
            self.value = None
            return self
        else:
            target = self.get_path(path[:-1])
            remove = path[-1]
            if remove in target.inner:
                lost = target.inner[remove]
                del target.inner[remove]
                return lost

    def divide_value(self):
        if self.divider:
            # divider is either a function or a dict with topology
            if isinstance(self.divider, dict):
                divider = self.divider['divider']
                topology = self.divider['topology']
                state = self.outer.get_values(topology)
                return divider(self.value, state)
            else:
                return self.divider(self.value)
        elif self.inner:
            daughters = [{}, {}]
            for key, child in self.inner.items():
                division = child.divide_value()
                if division:
                    for daughter, divide in zip(daughters, division):
                        daughter[key] = divide
            return daughters

    def reduce(self, reducer, initial=None):
        value = initial
        for path, node in self.depth():
            value = reducer(value, path, node)
        return value

    def reduce_to(self, path, reducer, initial=None):
        value = self.reduce(reducer, initial)
        update = assoc_path({}, path, value)
        self.apply_update(update)

    def set_value(self, value):
        if self.inner or self.subschema:
            for child, child_value in value.items():
                if not child in self.inner:
                    if self.subschema:
                        self.inner[child] = Store(self.subschema, self)
                    else:
                        pass

                        # TODO: continue to ignore extra keys?
                        # print("setting value that doesn't exist in tree {} {}".format(
                        #     child, child_value))

                if child in self.inner:
                    self.inner[child].set_value(child_value)
        else:
            self.value = value

    def apply_defaults(self):
        if self.inner:
            for child in self.inner.values():
                child.apply_defaults()
        else:
            if self.value is None:
                self.value = self.default

    def apply_update(self, update):
        if self.inner or self.subschema:
            topology_updates = {}

            if '_delete' in update:
                # delete a list of paths
                for path in update['_delete']:
                    self.delete_path(path)

                update = dissoc(update, ['_delete'])

            if '_generate' in update:
                # generate a list of new compartments
                for generate in update['_generate']:
                    self.generate(
                        generate['path'],
                        generate['processes'],
                        generate['topology'],
                        generate['initial_state'])
                    assoc_path(
                        topology_updates,
                        generate['path'],
                        generate['topology'])
                self.apply_subschemas()
                self.apply_defaults()

                update = dissoc(update, '_generate')

            if '_divide' in update:
                # use dividers to find initial states for daughters
                divide = update['_divide']
                mother = divide['mother']
                daughters = divide['daughters']
                initial_state = self.inner[mother].get_value(
                    condition=lambda child: not(isinstance(child.value, Process)),
                    f=lambda child: copy.deepcopy(child))
                states = self.inner[mother].divide_value()

                for daughter, state in zip(daughters, states):
                    daughter_id = daughter['daughter']

                    # use initial state as default, merge in divided values
                    initial_state = deep_merge(
                        initial_state,
                        state)

                    self.generate(
                        daughter['path'],
                        daughter['processes'],
                        daughter['topology'],
                        daughter['initial_state'])
                    assoc_path(
                        topology_updates,
                        daughter['path'],
                        daughter['topology'])

                    self.apply_subschemas()
                    self.inner[daughter_id].set_value(initial_state)
                    self.apply_defaults()
                self.delete_path((mother,))

                update = dissoc(update, '_divide')

            for key, value in update.items():
                if key in self.inner:
                    child = self.inner[key]
                    inner_updates = child.apply_update(value)
                    if inner_updates:
                        topology_updates = deep_merge(
                            topology_updates,
                            {key: inner_updates})
                elif self.subschema:
                    self.inner[key] = Store(self.subschema, self)
                    self.inner[key].set_value(value)
                    self.inner[key].apply_defaults()

            return topology_updates

        else:
            if isinstance(update, dict) and '_reduce' in update:
                reduction = update['_reduce']
                top = self.get_path(reduction.get('from'))
                update = top.reduce(
                    reduction['reducer'],
                    initial=reduction['initial'])

            updater = self.updater
            if isinstance(update, dict) and self.schema_keys & update.keys():
                if '_updater' in update:
                    updater = self.get_updater(update)
                    update = update.get('_value', self.default)

            # if self.units:
            #     units_value = updater(self.value * self.units, update)
            #     self.value = units_value.to(self.units).magnitude
            # else:
            self.value = updater(self.value, update)

    def child_value(self, key):
        if key in self.inner:
            return self.inner[key].get_value()

    def state_for(self, path, keys):
        state = self.get_path(path)
        if state is None:
            return {}
        elif keys and keys[0] == '*':
            return state.get_value()
        else:
            return {
                key: state.child_value(key)
                for key in keys}

    def depth(self, path=()):
        base = [(path, self)]
        for key, child in self.inner.items():
            down = tuple(path + (key,))
            base += child.depth(down)
        return base

    def processes(self, path=()):
        return {
            path: state
            for path, state in self.depth()
            if state.value and isinstance(state.value, Process)}

    def establish_path(self, path, config, initial=None, source=None):
        if len(path) > 0:
            path_step = path[0]
            remaining = path[1:]

            if path_step == '..':
                if not self.outer:
                    raise Exception('outer does not exist for path: {}'.format(path))
                return self.outer.establish_path(
                    remaining, config,
                    initial=initial.get('..') if initial and isinstance(
                        initial, dict) else None,
                    source=source)

            else:
                if not path_step in self.inner:
                    self.inner[path_step] = Store({}, outer=self, source=source)
                return self.inner[path_step].establish_path(
                    remaining, config,
                    initial=initial.get(path_step) if initial and isinstance(
                        initial, dict) else None,
                    source=source)

        else:
            self.apply_config(config, source=source)
            if initial:
                self.value = initial
            return self

    def apply_subschema(self, subschema=None):
        if subschema is None:
            subschema = self.subschema
        for child_key, child in self.inner.items():
            child.apply_config(subschema)

    def apply_subschemas(self):
        if self.subschema:
            self.apply_subschema()
        for child in self.inner.values():
            child.apply_subschemas()

    def update_subschema(self, path, subschema):
        target = self.get_path(path)
        if target.subschema is None:
            target.subschema = subschema
        else:
            target.subschema = deep_merge(
                target.subschema,
                subschema)
        return target

    def generate_paths(self, processes, topology, initial_state):
        for key, subprocess in processes.items():
            subtopology = topology[key]
            if isinstance(subprocess, Process):
                process_state = Store({
                    '_value': subprocess,
                    '_updater': 'set'}, outer=self)
                self.inner[key] = process_state
                for port, targets in subprocess.ports_schema().items():
                    if port not in subtopology:
                        raise Exception('topology conflict: {} process does not have {} port'.format(key, port))
                    path = subtopology[port]
                    if path:
                        initial = get_in(initial_state, path)
                        for target, schema in targets.items():
                            source = self.path_for() + (key,)
                            if target == '*':
                                glob = self.establish_path(
                                    path, {
                                        '_subschema': schema},
                                    source=source)
                                glob.apply_subschema()
                                glob.apply_defaults()
                            else:
                                subpath = tuple(path) + (target,)
                                self.establish_path(
                                    subpath,
                                    schema,
                                    initial=initial.get(
                                        target) if initial and isinstance(
                                            initial, dict) else None,
                                    source=source)
            else:
                if not key in self.inner:
                    self.inner[key] = Store({}, outer=self)
                substate = initial_state.get(key, {})
                self.inner[key].generate_paths(
                    subprocess,
                    subtopology,
                    substate)

    def generate(self, path, processes, topology, initial_state):
        target = self.establish_path(path, {})
        target.generate_paths(processes, topology, initial_state)
        target.set_value(initial_state)
        target.apply_defaults()


def generate_derivers(processes, topology):
    deriver_processes = {}
    deriver_topology = {}
    for process_key, node in processes.items():
        subtopology = topology[process_key]
        if isinstance(node, Process):
            for deriver_key, config in node.derivers().items():
                if deriver_key not in deriver_processes:
                    # generate deriver process
                    deriver_config = config.get('config', {})
                    generate = config['deriver']
                    if isinstance(generate, str):
                        generate = deriver_library[generate]

                    deriver = generate(deriver_config)
                    deriver_processes[deriver_key] = deriver

                    # generate deriver topology
                    deriver_topology[deriver_key] = {}
                    for target, source in config.get('port_mapping', {}).items():
                        path = subtopology[source]
                        deriver_topology[deriver_key][target] = path
        else:
            subderivers = generate_derivers(node, subtopology)
            deriver_processes[process_key] = subderivers['processes']
            deriver_topology[process_key] = subderivers['topology']
    return {
        'processes': deriver_processes,
        'topology': deriver_topology}


class Compartment(object):
    def __init__(self, config):
        self.config = config

    def generate_processes(self, config):
        return {}

    def generate_topology(self, config):
        return {}

    def generate(self, config=None, path=tuple()):
        if config is None:
            config = {}

        processes = self.generate_processes(config)
        topology = self.generate_topology(config)

        # add derivers
        derivers = generate_derivers(processes, topology)
        processes = deep_merge(processes, derivers['processes'])
        topology = deep_merge(topology, derivers['topology'])

        return {
            'processes': assoc_in({}, path, processes),
            'topology': assoc_in({}, path, topology)}

    def get_parameters(self):
        processes = self.generate_processes({})
        return {
            process_id: process.parameters
            for process_id, process in processes.items()}


def generate_state(processes, topology, initial_state):
    state = Store({})
    state.generate_paths(processes, topology, initial_state)
    state.set_value(initial_state)
    state.apply_defaults()

    return state


def normalize_path(path):
    progress = []
    for step in path:
        if step == '..' and len(progress) > 0:
            progress = progress[:-1]
        else:
            progress.append(step)
    return progress


def timestamp(dt=None):
    if not dt:
        dt = datetime.datetime.now()

    return "%04d%02d%02d.%02d%02d%02d" % (
        dt.year, dt.month, dt.day,
        dt.hour, dt.minute, dt.second)


class Experiment(object):
    def __init__(self, config):
        self.config = config
        self.experiment_id = config.get('experiment_id', uuid.uuid1())
        self.description = config.get('description','')
        self.processes = config['processes']
        self.topology = config['topology']
        self.initial_state = config['initial_state']

        self.state = generate_state(
            self.processes,
            self.topology,
            self.initial_state)

        emitter_config = config.get('emitter', {})
        emitter_config['experiment_id'] = self.experiment_id
        self.emitter = get_emitter(emitter_config)

        self.local_time = 0.0

        # run the derivers
        self.send_updates([])

        # run emitter
        self.emit_configuration()
        self.emit_data()

        print('experiment {}'.format(self.experiment_id))

        print('\nPROCESSES:')
        pp(self.processes)

        print('\nTOPOLOGY:')
        pp(self.topology)

        # print('\nSTATE:')
        # pp(self.state.get_value())

        # print('\nCONFIG:')
        # pp(self.state.get_config())

    def emit_configuration(self):
        data = {
            'time_created': timestamp(),
            'experiment_id': self.experiment_id,
            'description': self.description,
            'processes': self.processes,
            'topology': self.topology,
            'state': self.state.get_config()}
        emit_config = {
            'table': 'configuration',
            'data': data}
        self.emitter.emit(emit_config)

    def absolute_update(self, path, new_update):
        absolute = {}
        for port, update in new_update.items():
            topology = get_in(self.topology, path + (port,))
            if topology is not None:
                state_path = path[:-1] + topology
                normal_path = normalize_path(state_path)
                absolute = assoc_in(absolute, normal_path, update)
        return absolute

    def process_update(self, path, state, interval):
        process = state.value
        process_topology = get_in(self.topology, path)
        ports = process.find_states(state.outer, process_topology)
        update = process.next_update(interval, ports)
        absolute = self.absolute_update(path, update)
        return absolute

    def apply_update(self, update):
        topology_updates = self.state.apply_update(update)
        if topology_updates:
            # print('topology updates for update {}: {}'.format(update, topology_updates))
            self.topology = deep_merge(self.topology, topology_updates)

    def run_derivers(self, derivers):
        for path, deriver in derivers.items():
            # timestep shouldn't influence derivers
            update = self.process_update(path, deriver, 0)
            self.apply_update(update)

    # def emit_paths(self, paths):
    #     emit_config = {
    #         'table': 'history',
    #         'data': data}
    #     self.emitter_emit(emit_config)

    def emit_data(self):
        data = self.state.emit_data()
        data.update({
            'time': self.local_time})
        emit_config = {
            'table': 'history',
            'data': data}
        self.emitter.emit(emit_config)

    def send_updates(self, updates, derivers=None):
        for update in updates:
            self.apply_update(update)
        if derivers is None:
            derivers = {
                path: state
                for path, state in self.state.depth()
                if state.value is not None and isinstance(state.value, Process) and state.value.is_deriver()}
        self.run_derivers(derivers)

    def update(self, timestep):
        ''' Run each process for the given time step and update the related states. '''

        time = 0

        def empty_front(t):
            return {
                'time': t,
                'update': {}}

        # keep track of which processes have simulated until when
        front = {}

        while time < timestep:
            full_step = INFINITY

            if VERBOSE:
                for state_id in self.states:
                    print('{}: {}'.format(time, self.states[state_id].to_dict()))

            processes = {}
            derivers = {}
            for path, state in self.state.depth():
                if state.value is not None and isinstance(state.value, Process):
                    if state.value.is_deriver():
                        derivers[path] = state
                    else:
                        processes[path] = state

            front = {
                path: process
                for path, process in front.items()
                if path in processes}

            for path, state in processes.items():
                if not path in front:
                    front[path] = empty_front(time)
                process_time = front[path]['time']

                if process_time <= time:
                    process = state.value
                    future = min(process_time + process.local_timestep(), timestep)
                    interval = future - process_time
                    update = self.process_update(path, state, interval)

                    if interval < full_step:
                        full_step = interval
                    front[path]['time'] = future
                    front[path]['update'] = update

            if full_step == INFINITY:
                # no processes ran, jump to next process
                next_event = timestep
                for process_name in front.keys():
                    if front[path]['time'] < next_event:
                        next_event = front[path]['time']
                time = next_event
            else:
                # at least one process ran, apply updates and continue
                future = time + full_step

                updates = []
                paths = []
                for path, advance in front.items():
                    if advance['time'] <= future:
                        new_update = advance['update']
                        new_update['_path'] = path
                        updates.append(new_update)
                        advance['update'] = {}
                        paths.append(path)

                self.send_updates(updates, derivers)
                # self.emit_paths(paths)
                self.emit_data()

                time = future

        for process_name, advance in front.items():
            assert advance['time'] == time == timestep
            assert len(advance['update']) == 0

        self.local_time += timestep

        # run emitters
        # self.emit_data()

    def update_interval(self, time, interval):
        while self.local_time < time:
            self.update(interval)


def test_recursive_store():
    environment_config = {
        'environment': {
            'temperature': {
                '_default': 0.0,
                '_updater': 'accumulate'},
            'fields': {
                (0, 1): {
                    'enzymeX': {
                        '_default': 0.0,
                        '_updater': 'set'},
                    'enzymeY': {
                        '_default': 0.0,
                        '_updater': 'set'}},
                (0, 2): {
                    'enzymeX': {
                        '_default': 0.0,
                        '_updater': 'set'},
                    'enzymeY': {
                        '_default': 0.0,
                        '_updater': 'set'}}},
            'agents': {
                '1': {
                    'location': {
                        '_default': (0, 0),
                        '_updater': 'set'},
                    'boundary': {
                        'external': {
                            '_default': 0.0,
                            '_updater': 'set'},
                        'internal': {
                            '_default': 0.0,
                            '_updater': 'set'}},
                    'transcripts': {
                        'flhDC': {
                            '_default': 0,
                            '_updater': 'accumulate'},
                        'fliA': {
                            '_default': 0,
                            '_updater': 'accumulate'}},
                    'proteins': {
                        'ribosome': {
                            '_default': 0,
                            '_updater': 'set'},
                        'flagella': {
                            '_default': 0,
                            '_updater': 'accumulate'}}},
                '2': {
                    'location': {
                        '_default': (0, 0),
                        '_updater': 'set'},
                    'boundary': {
                        'external': {
                            '_default': 0.0,
                            '_updater': 'set'},
                        'internal': {
                            '_default': 0.0,
                            '_updater': 'set'}},
                    'transcripts': {
                        'flhDC': {
                            '_default': 0,
                            '_updater': 'accumulate'},
                        'fliA': {
                            '_default': 0,
                            '_updater': 'accumulate'}},
                    'proteins': {
                        'ribosome': {
                            '_default': 0,
                            '_updater': 'set'},
                        'flagella': {
                            '_default': 0,
                            '_updater': 'accumulate'}}}}}}

    state = Store(environment_config)
    state.apply_update({})
    state.state_for(['environment'], ['temperature'])

def test_in():
    blank = {}
    path = ['where', 'are', 'we']
    assoc_path(blank, path, 5)
    print(blank)
    print(get_in(blank, path))
    update_in(blank, path, lambda x: x + 6)
    print(blank)


def test_timescales():
    class Slow(Process):
        def __init__(self):
            self.timestep = 3.0
            self.ports = {
                'state': ['base']}

        def ports_schema(self):
            return {
                'state': {
                    'base': {
                        '_default': 1.0}}}

        def local_timestep(self):
            return self.timestep

        def next_update(self, timestep, states):
            base = states['state']['base']
            next_base = timestep * base * 0.1

            return {
                'state': {'base': next_base}}

    class Fast(Process):
        def __init__(self):
            self.timestep = 0.1
            self.ports = {
                'state': ['base', 'motion']}

        def ports_schema(self):
            return {
                'state': {
                    'base': {
                        '_default': 1.0},
                    'motion': {
                        '_default': 0.0}}}

        def local_timestep(self):
            return self.timestep

        def next_update(self, timestep, states):
            base = states['state']['base']
            motion = timestep * base * 0.001

            return {
                'state': {'motion': motion}}

    processes = {
        'slow': Slow(),
        'fast': Fast()}

    states = {
        'state': {
            'base': 1.0,
            'motion': 0.0}}

    topology = {
        'slow': {'state': ('state',)},
        'fast': {'state': ('state',)}}

    emitter = {'type': 'null'}
    experiment = Experiment({
        'processes': processes,
        'topology': topology,
        'emitter': emitter,
        'initial_state': states})

    experiment.update(10.0)



if __name__ == '__main__':
    test_recursive_store()
    test_in()
    test_timescales()
