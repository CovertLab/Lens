'''
====================
Colony Shape Deriver
====================
'''

from __future__ import absolute_import, division, print_function

import alphashape
import numpy as np
from pytest import approx
from shapely.geometry.collection import GeometryCollection
from shapely.geometry.linestring import LineString
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.point import Point
from shapely.geometry.polygon import Polygon

from vivarium.core.experiment import get_in
from vivarium.core.process import Deriver
from vivarium.library.units import units
from vivarium.core.registry import assert_no_divide


class Variables():
    AREA = 'surface_area'
    MASS = 'mass'
    CIRCUMFERENCE = 'circumference'
    MAJOR_AXIS = 'major_axis'
    MINOR_AXIS = 'minor_axis'



def major_minor_axes(shape):
    '''Calculate the lengths of the major and minor axes of a shape

    We assume that the major and minor axes are the dimensions of the
    minimum bounding rectangle of the shape. Note that this is different
    from using PCA to find the axes, especially for highly asymmetrical
    and concave shapes.

    Arguments:
        shape (shapely.polygon.Polygon): The shape to compute axes for.

    Returns:
        tuple: A tuple with the major axis first and the minor axis
        second.
    '''
    rect = shape.minimum_rotated_rectangle
    points = list(rect.exterior.coords)
    points = [np.array(point) for point in points]
    # Shapely returns polygon coordinates in the order in which they
    # would appear while traversing the boundary, so we know that any 3
    # consecutive points span the major and minor axes
    dimension_1 = abs(np.linalg.norm(points[1] - points[0]))
    dimension_2 = abs(np.linalg.norm(points[2] - points[1]))
    major = max(dimension_1, dimension_2)
    minor = min(dimension_1, dimension_2)
    return major, minor


def gen_agent_colony_map(agents, colony_shapes):
    '''Create a map from agent to the colony to which the agent belongs

    An agent is considered within a colony if its ``location``
    :term:`variable` intersects with the colony's interior or border. If
    an agent intersects with multiple colonies, we choose arbitrarily
    which colony the agent belongs to.

    Points containing any nan coordinate are considered to be outside of
    all shapes.

    .. note:: An agent may be part of no colony at all, in which case it
        will not be included in the returned map.

    Arguments:
        agents (dict): Dictionary of ``agents`` :term:`port` state whose
            keys are agent IDs and whose values are agent state
            dictionaries.
        colony_shapes (list): List of polygons that define the colonies.

    Returns:
        dict: Map from agent ID to index of colony in ``colony_shapes``.
    '''
    agent_colony_map = {}
    for agent_id, agent_state in agents.items():
        loc = agent_state['boundary']['location']
        if np.any(np.isnan(loc)):
            continue
        point = Point(*loc)
        for i, colony_shape in enumerate(colony_shapes):
            if colony_shape.intersects(point):
                agent_colony_map[agent_id] = i
                break
    return agent_colony_map


class ColonyShapeDeriver(Deriver):
    '''Derives colony shape metrics from cell locations
    '''

    name = 'colony_shape_deriver'
    defaults = {
        'alpha': 1.0,
        'bounds': [1, 1],
    }

    def ports_schema(self):
        return {
            'agents': {
                '*': {
                    'boundary': {
                        'location': {
                            '_default': [
                                0.5 * bound for bound in
                                self.parameters['bounds']
                            ],
                        },
                        'mass': {
                            '_default': 1339 * units.fg,
                        },
                    },
                },
            },
            'colony_global': {
                'surface_area': {
                    '_default': [],
                    '_updater': 'set',
                    '_divider': assert_no_divide,
                    '_emit': True,
                },
                'major_axis': {
                    '_default': [],
                    '_updater': 'set',
                    '_divider': assert_no_divide,
                    '_emit': True,
                },
                'minor_axis': {
                    '_default': [],
                    '_updater': 'set',
                    '_divider': assert_no_divide,
                    '_emit': True,
                },
                'mass': {
                    '_default': [],
                    '_updater': 'set',
                    '_divider': assert_no_divide,
                    '_emit': True,
                },
                'circumference': {
                    '_default': [],
                    '_updater': 'set',
                    '_divider': assert_no_divide,
                    '_emit': True,
                }
            },
        }

    def next_update(self, timestep, states):
        agents = states['agents']
        points = [
            agent['boundary']['location']
            for agent in agents.values()
        ]
        points = [tuple(point) for point in points if np.all(~np.isnan(point))]
        alpha_shape = alphashape.alphashape(
            set(points), self.parameters['alpha'])
        if isinstance(alpha_shape, Polygon):
            shapes = [alpha_shape]
        elif isinstance(alpha_shape, (Point, LineString)):
            # We need at least 3 cells to form a colony polygon
            shapes = []
        else:
            assert isinstance(
                alpha_shape, (MultiPolygon, GeometryCollection))
            shapes = list(alpha_shape)

        agent_colony_map = gen_agent_colony_map(agents, shapes)

        # Calculate colony surface areas
        areas = [shape.area for shape in shapes]

        # Calculate colony major and minor axes based on bounding
        # rectangles
        major_axes = []
        minor_axes = []
        for shape in shapes:
            if isinstance(shape, Polygon):
                major, minor = major_minor_axes(shape)
                major_axes.append(major)
                minor_axes.append(minor)
            else:
                major_axes.append(0)
                minor_axes.append(0)

        # Calculate colony circumference
        circumference = [shape.length for shape in shapes]

        # Calculate colony masses and cell surface areas
        mass = [0] * len(shapes)
        cell_area = [0] * len(shapes)
        for agent_id, agent_state in agents.items():
            if agent_id not in agent_colony_map:
                # We ignore agents not in any colony
                continue
            colony_index = agent_colony_map[agent_id]
            agent_mass = get_in(agent_state, ('boundary', 'mass'), 0)
            mass[colony_index] += agent_mass

        return {
            'colony_global': {
                'surface_area': areas,
                'major_axis': major_axes,
                'minor_axis': minor_axes,
                'circumference': circumference,
                'mass': mass,
            }
        }

ColonyShapeDeriver()

class TestDeriveColonyShape():

    def calc_shape_metrics(self, points, alpha=None):
        config = {}
        if alpha is not None:
            config['alpha'] = alpha
        deriver = ColonyShapeDeriver(config)
        states = {
            'agents': {
                str(i): {
                    'boundary': {
                        'location': list(point),
                        'mass': 1.0 * units.fg,
                    },
                }
                for i, point in enumerate(points)
            },
            'colony_global': {
                'surface_area': [],
                'major_axis': [],
                'minor_axis': [],
                'circumference': [],
            }
        }
        # Timestep does not matter
        update = deriver.next_update(-1, states)
        return update['colony_global']

    def flatten(self, lst):
        return [
            elem
            for sublist in lst
            for elem in sublist
        ]

    def test_convex(self):
        #    *
        #   / \
        #  * * *
        #   \ /
        #    *
        points = [
            (1, 2),
            (0, 1), (1, 1), (2, 1),
            (1, 0),
        ]
        metrics = self.calc_shape_metrics(points)
        assert metrics['surface_area'] == [2]
        assert metrics['major_axis'] == approx([np.sqrt(2)])
        assert metrics['minor_axis'] == approx([np.sqrt(2)])
        assert metrics['circumference'] == approx([4 * np.sqrt(2)])
        assert metrics['mass'] == [5 * units.fg]

    def test_concave(self):
        # *-*-*-*-*
        # |       |
        # * * *-*-*
        # |  /
        # * *
        # |  \
        # * * *-*-*
        # |       |
        # *-*-*-*-*
        points = (
            [(i, 4) for i in range(5)]
            + [(i, 3) for i in range(5)]
            + [(i, 2) for i in range(2)]
            + [(i, 1) for i in range(5)]
            + [(i, 0) for i in range(5)]
        )
        metrics = self.calc_shape_metrics(points)
        assert metrics['surface_area'] == [11]
        assert metrics['major_axis'] == [4]
        assert metrics['minor_axis'] == [4]
        assert metrics['circumference'] == approx([18 + 2 * np.sqrt(2)])
        assert metrics['mass'] == [22 * units.fg]

    def test_ignore_outliers_and_nan(self):
        #    *
        #   / \
        #  * * *            *
        #   \ /
        #    *
        points = [
            (1, 2),
            (0, 1), (1, 1), (2, 1), (10, 1),
            (1, 0),
            (np.nan, np.nan),
        ]
        metrics = self.calc_shape_metrics(points)
        assert metrics['surface_area'] == [2]
        assert metrics['major_axis'] == approx([np.sqrt(2)])
        assert metrics['minor_axis'] == approx([np.sqrt(2)])
        assert metrics['circumference'] == approx([4 * np.sqrt(2)])
        assert metrics['mass'] == [5 * units.fg]

    def test_colony_too_diffuse(self):
        #    *
        #
        #  *   *
        #
        #    *
        points = [
            (1, 2),
            (0, 1), (2, 1),
            (1, 0),
        ]
        metrics = self.calc_shape_metrics(points)
        expected_metrics = {
            'surface_area': [],
            'major_axis': [],
            'minor_axis': [],
            'circumference': [],
            'mass': [],
        }
        assert metrics == expected_metrics

    def test_single_cell(self):
        #
        #  *
        #
        points = [
            (1, 1),
        ]
        metrics = self.calc_shape_metrics(points)
        expected_metrics = {
            'surface_area': [],
            'major_axis': [],
            'minor_axis': [],
            'circumference': [],
            'mass': [],
        }
        assert metrics == expected_metrics

    def test_two_cells(self):
        #
        #  *   *
        #
        points = [
            (0, 1), (2, 1),
        ]
        metrics = self.calc_shape_metrics(points)
        expected_metrics = {
            'surface_area': [],
            'major_axis': [],
            'minor_axis': [],
            'circumference': [],
            'mass': [],
        }
        assert metrics == expected_metrics

    def test_no_cells(self):
        points = []
        metrics = self.calc_shape_metrics(points)
        expected_metrics = {
            'surface_area': [],
            'major_axis': [],
            'minor_axis': [],
            'circumference': [],
            'mass': [],
        }
        assert metrics == expected_metrics

    def test_find_multiple_colonies(self):
        #    *          *
        #   / \        / \
        #  * * *      * * *
        #   \ /        \ /
        #    *          *
        points = [
            (1, 2), (11, 2),
            (0, 1), (1, 1), (2, 1), (10, 1), (11, 1), (12, 1),
            (1, 0), (11, 0),
        ]
        metrics = self.calc_shape_metrics(points)
        assert metrics['surface_area'] == [2, 2]
        assert metrics['major_axis'] == approx([np.sqrt(2), np.sqrt(2)])
        assert metrics['minor_axis'] == approx([np.sqrt(2), np.sqrt(2)])
        assert metrics['circumference'] == approx([4 * np.sqrt(2)] * 2)
        assert metrics['mass'] == [5 * units.fg] * 2
