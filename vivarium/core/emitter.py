from __future__ import absolute_import, division, print_function

from pymongo import MongoClient
from confluent_kafka import Producer
import json
import copy

from vivarium.actor.actor import delivery_report
from vivarium.library.dict_utils import (
    merge_dicts, value_in_embedded_dict, get_path_list_from_dict, \
    get_value_from_path, make_path_dict)

HISTORY_INDEXES = [
    'time',
    'type',
    'simulation_id',
    'experiment_id']

CONFIGURATION_INDEXES = [
    'type',
    'simulation_id',
    'experiment_id']


def create_indexes(table, columns):
    '''Create all of the necessary indexes for the given table name.'''
    for column in columns:
        table.create_index(column)

def get_emitter(config):
    '''
    Get an emitter based on the provided config.

    config is a dict and requires three keys:
    * type: Type of emitter ('kafka' for a kafka emitter).
    * emitter: Any configuration the emitter type requires to initialize.
    * keys: A list of state keys to emit for each state label.
    '''

    if config is None:
        config = {'type': 'print'}
    emitter_type = config.get('type', 'print')

    if emitter_type == 'kafka':
        emitter = KafkaEmitter(config)
    elif emitter_type == 'database':
        emitter = DatabaseEmitter(config)
    elif emitter_type == 'null':
        emitter = NullEmitter(config)
    elif emitter_type == 'timeseries':
        emitter = TimeSeriesEmitter(config)
    else:
        emitter = Emitter(config)

    return emitter

def configure_emitter(config, processes, topology):
    emitter_config = config.get('emitter', {})
    emitter_config['keys'] = get_emitter_keys(processes, topology)
    emitter_config['experiment_id'] = config.get('experiment_id')
    emitter_config['simulation_id'] = config.get('simulation_id')
    return get_emitter(emitter_config)

def path_timeseries_from_data(data):
    embedded_timeseries = timeseries_from_data(data)
    times_vector = embedded_timeseries.pop('time')
    path_timeseries = make_path_dict(embedded_timeseries)
    path_timeseries['time'] = times_vector
    return path_timeseries

def timeseries_from_data(data):
    times_vector = list(data.keys())
    embedded_timeseries = {}
    for time, value in data.items():
        if isinstance(value, dict):
            embedded_timeseries = value_in_embedded_dict(value, embedded_timeseries)
        else:
            pass

    embedded_timeseries['time'] = times_vector
    return embedded_timeseries


class Emitter(object):
    '''
    Emit data to terminal
    '''
    def __init__(self, config):
        self.config = config

    def emit(self, data):
        print(data)

    def get_timeseries(self):
        raise Exception('emitter does not get timeseries')
        return {}

class NullEmitter(Emitter):
    '''
    Don't emit anything
    '''
    def emit(self, data):
        pass

class TimeSeriesEmitter(Emitter):

    def __init__(self, config):
        keys = config.get('keys', {})
        self.saved_data = {}

    def emit(self, data):
        # save history data
        if data['table'] == 'history':
            emit_data = data['data']
            time = emit_data.pop('time')
            self.saved_data[time] = emit_data

    def get_data(self):
        return self.saved_data

    def get_path_timeseries(self):
        return path_timeseries_from_data(self.saved_data)

    def get_timeseries(self):
        return timeseries_from_data(self.saved_data)


class KafkaEmitter(Emitter):
    '''
    Emit data to kafka

    example:
    config = {
        'host': 'localhost:9092',
        'topic': 'EMIT'}
    '''
    def __init__(self, config):
        self.config = config
        self.producer = Producer({
            'bootstrap.servers': self.config['host']})

    def emit(self, data):
        encoded = json.dumps(data, ensure_ascii=False).encode('utf-8')

        self.producer.produce(
            self.config['topic'],
            encoded,
            callback=delivery_report)

        self.producer.flush(timeout=0.1)


class DatabaseEmitter(Emitter):
    '''
    Emit data to a mongoDB database

	example:
	config = {
        'host': 'localhost:27017',
        'database': 'DB_NAME'}
    '''
    client = None

    def __init__(self, config):
        self.config = config
        self.experiment_id = config.get('experiment_id')

        # create singleton instance of mongo client
        if DatabaseEmitter.client is None:
            DatabaseEmitter.client = MongoClient(config['host'])

        self.db = getattr(self.client, config.get('database', 'simulations'))
        self.history = getattr(self.db, 'history')
        self.configuration = getattr(self.db, 'configuration')
        self.phylogeny = getattr(self.db, 'phylogeny')
        create_indexes(self.history, HISTORY_INDEXES)
        create_indexes(self.configuration, CONFIGURATION_INDEXES)
        create_indexes(self.phylogeny, CONFIGURATION_INDEXES)

    def emit(self, data_config):
        data = data_config['data']
        data.update({
            'experiment_id': self.experiment_id})

        table = getattr(self.db, data_config['table'])
        table.insert_one(data)
