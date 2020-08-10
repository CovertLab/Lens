from __future__ import absolute_import, division, print_function

import os
import math
import random
import itertools

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
from matplotlib.colors import hsv_to_rgb
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import LineCollection
import numpy as np

from vivarium.library.dict_utils import get_value_from_path


DEFAULT_BOUNDS = [10, 10]

# constants
PI = math.pi

# colors for phylogeny initial agents
HUES = [hue/360 for hue in np.linspace(0,360,30)]
DEFAULT_HUE = HUES[0]
DEFAULT_SV = [100.0/100.0, 70.0/100.0]
BASELINE_TAG_COLOR = [220/360, 1.0, 0.2]  # HSV
FLOURESCENT_SV = [0.75, 1.0]  # SV for fluorescent colors

def check_plt_backend():
    # reset matplotlib backend for non-interactive plotting
    plt.close('all')
    if plt.get_backend() == 'TkAgg':
        matplotlib.use('Agg')

class LineWidthData(Line2D):
    def __init__(self, *args, **kwargs):
        _lw_data = kwargs.pop('linewidth', 1)
        super().__init__(*args, **kwargs)
        self._lw_data = _lw_data

    def _get_lw(self):
        if self.axes is not None:
            ppd = 72./self.axes.figure.dpi
            trans = self.axes.transData.transform
            return ((trans((1, self._lw_data))-trans((0, 0)))*ppd)[1]
        else:
            return 1

    def _set_lw(self, lw):
        self._lw_data = lw

    _linewidth = property(_get_lw, _set_lw)

def plot_agent(ax, data, color, agent_shape):
    # location, orientation, length
    x_center = data['boundary']['location'][0]
    y_center = data['boundary']['location'][1]
    theta = data['boundary']['angle'] / PI * 180 + 90 # rotate 90 degrees to match field
    length = data['boundary']['length']
    width = data['boundary']['width']

    # get color, convert to rgb
    rgb = hsv_to_rgb(color)

    if agent_shape == 'rectangle':
        # get bottom left position
        x_offset = (width / 2)
        y_offset = (length / 2)
        theta_rad = math.radians(theta)
        dx = x_offset * math.cos(theta_rad) - y_offset * math.sin(theta_rad)
        dy = x_offset * math.sin(theta_rad) + y_offset * math.cos(theta_rad)

        x = x_center - dx
        y = y_center - dy

        # Create a rectangle
        shape = patches.Rectangle(
            (x, y), width, length,
            angle=theta,
            linewidth=2,
            edgecolor='w',
            facecolor=rgb
        )
        ax.add_patch(shape)

    elif agent_shape == 'segment':
        membrane_width = 0.1
        membrane_color = [1, 1, 1]
        radius = width / 2

        # get the two ends
        length_offset = (length / 2) - radius
        theta_rad = math.radians(theta)
        dx = - length_offset * math.sin(theta_rad)
        dy = length_offset * math.cos(theta_rad)

        x1 = x_center - dx
        y1 = y_center - dy
        x2 = x_center + dx
        y2 = y_center + dy

        # segment plot
        membrane = LineWidthData(
            [x1, x2], [y1, y2],
            color=membrane_color,
            linewidth=width,
            solid_capstyle='round')
        line = LineWidthData(
            [x1, x2], [y1, y2],
            color=rgb,
            linewidth=width-membrane_width,
            solid_capstyle='round')
        ax.add_line(membrane)
        ax.add_line(line)


def plot_agents(
    ax, agents, agent_colors={}, agent_shape='segment', dead_color=None
):
    '''
    - ax: the axis for plot
    - agents: a dict with {agent_id: agent_data} and
        agent_data a dict with keys location, angle, length, width
    - agent_colors: dict with {agent_id: hsv color}
    - dead_color: List of 3 floats that define HSV color to use for dead
      cells. Dead cells only get treated differently if this is set.
    '''
    for agent_id, agent_data in agents.items():
        color = agent_colors.get(agent_id, [DEFAULT_HUE]+DEFAULT_SV)
        if dead_color:
            if agent_data['boundary']['dead']:
                color = dead_color
        plot_agent(ax, agent_data, color, agent_shape)

def mutate_color(baseline_hsv):
    mutation = 0.1
    new_hsv = [
        (n + np.random.uniform(-mutation, mutation))
        for n in baseline_hsv]
    # wrap hue around
    new_hsv[0] = new_hsv[0] % 1
    # reflect saturation and value
    if new_hsv[1] > 1:
        new_hsv[1] = 2 - new_hsv[1]
    if new_hsv[2] > 1:
        new_hsv[2] = 2 - new_hsv[2]
    return new_hsv

def color_phylogeny(ancestor_id, phylogeny, baseline_hsv, phylogeny_colors={}):
    """
    get colors for all descendants of the ancestor
    through recursive calls to each generation
    """
    phylogeny_colors.update({ancestor_id: baseline_hsv})
    daughter_ids = phylogeny.get(ancestor_id)
    if daughter_ids:
        for daughter_id in daughter_ids:
            daughter_color = mutate_color(baseline_hsv)
            color_phylogeny(daughter_id, phylogeny, daughter_color)
    return phylogeny_colors

def get_phylogeny_colors_from_names(agent_ids):
    '''Get agent colors using phlogeny saved in agent_ids
    This assumes the names use daughter_phylogeny_id() from meta_division
    '''

    # make phylogeny with {mother_id: [daughter_1_id, daughter_2_id]}
    phylogeny = {agent_id: [] for agent_id in agent_ids}
    for agent1, agent2 in itertools.combinations(agent_ids, 2):
        if agent1 == agent2[0:-1]:
            phylogeny[agent1].append(agent2)
        elif agent2 == agent1[0:-1]:
            phylogeny[agent2].append(agent1)

    # get initial ancestors
    daughters = list(phylogeny.values())
    daughters = set([item for sublist in daughters for item in sublist])
    mothers = set(list(phylogeny.keys()))
    ancestors = list(mothers - daughters)

    # agent colors based on phylogeny
    agent_colors = {agent_id: [] for agent_id in agent_ids}
    for agent_id in ancestors:
        hue = random.choice(HUES)  # select random initial hue
        initial_color = [hue] + DEFAULT_SV
        agent_colors.update(color_phylogeny(agent_id, phylogeny, initial_color))

    return agent_colors


def plot_snapshots(data, plot_config):
    '''Plot snapshots of the simulation over time

    The snapshots depict the agents and environmental molecule
    concentrations.

    Arguments:
        data (dict): A dictionary with the following keys:

            * **agents** (:py:class:`dict`): A mapping from times to
              dictionaries of agent data at that timepoint. Agent data
              dictionaries should have the same form as the hierarchy
              tree rooted at ``agents``.
            * **fields** (:py:class:`dict`): A mapping from times to
              dictionaries of environmental field data at that
              timepoint.  Field data dictionaries should have the same
              form as the hierarchy tree rooted at ``fields``.
            * **config** (:py:class:`dict`): The environmental
              configuration dictionary  with the following keys:

                * **bounds** (:py:class:`tuple`): The dimensions of the
                  environment.

        plot_config (dict): Accepts the following configuration options.
            Any options with a default is optional.

            * **n_snapshots** (:py:class:`int`): Number of snapshots to
              show per row (i.e. for each molecule). Defaults to 6.
            * **agent_shape** (:py:class:`str`): the shape of the agents.
              select from **rectangle**, **segment**
            * **phylogeny_names** (:py:class:`bool`): This selects agent
              colors based on phylogenies seved in their names using
              meta_division.py daughter_phylogeny_id()
            * **out_dir** (:py:class:`str`): Output directory, which is
              ``out`` by default.
            * **filename** (:py:class:`str`): Base name of output file.
              ``snapshots`` by default.
            * **skip_fields** (:py:class:`Iterable`): Keys of fields to
              exclude from the plot. This takes priority over
              ``include_fields``.
            * **include_fields** (:py:class:`Iterable`): Keys of fields
              to plot.
            * **field_label_size** (:py:class:`float`): Font size of the
              field label.
            * **dead_color** (:py:class:`list` of 3 :py:class:`float`s):
              Color for dead cells in HSV. Defaults to [0, 0, 0], which
              is black.
    '''
    check_plt_backend()

    n_snapshots = plot_config.get('n_snapshots', 6)
    out_dir = plot_config.get('out_dir', 'out')
    filename = plot_config.get('filename', 'snapshots')
    agent_shape = plot_config.get('agent_shape', 'segment')
    phylogeny_names = plot_config.get('phylogeny_names', True)
    skip_fields = plot_config.get('skip_fields', [])
    include_fields = plot_config.get('include_fields', None)
    field_label_size = plot_config.get('field_label_size', 20)
    dead_color = plot_config.get('dead_color', [0, 0, 0])

    # get data
    agents = data.get('agents', {})
    fields = data.get('fields', {})
    config = data.get('config', {})
    bounds = config.get('bounds', DEFAULT_BOUNDS)
    edge_length_x = bounds[0]
    edge_length_y = bounds[1]

    # time steps that will be used
    if agents and fields:
        assert set(list(agents.keys())) == set(list(fields.keys())), 'agent and field times are different'
        time_vec = list(agents.keys())
    elif agents:
        time_vec = list(agents.keys())
    elif fields:
        time_vec = list(fields.keys())
    else:
        raise Exception('No agents or field data')

    time_indices = np.round(np.linspace(0, len(time_vec) - 1, n_snapshots)).astype(int)
    snapshot_times = [time_vec[i] for i in time_indices]

    # get fields id and range
    if fields:
        if include_fields is None:
            field_ids = set(fields[time_vec[0]].keys())
        else:
            field_ids = set(include_fields)
        field_ids -= set(skip_fields)
        field_range = {}
        for field_id in field_ids:
            field_min = min([min(min(field_data[field_id])) for t, field_data in fields.items()])
            field_max = max([max(max(field_data[field_id])) for t, field_data in fields.items()])
            field_range[field_id] = [field_min, field_max]

    # get agent ids
    agent_ids = set()
    if agents:
        for time, time_data in agents.items():
            current_agents = list(time_data.keys())
            agent_ids.update(current_agents)
        agent_ids = list(agent_ids)

        # set agent colors
        if phylogeny_names:
            agent_colors = get_phylogeny_colors_from_names(agent_ids)
        else:
            agent_colors = {}
            for agent_id in agent_ids:
                hue = random.choice(HUES)  # select random initial hue
                color = [hue] + DEFAULT_SV
                agent_colors[agent_id] = color

    # make the figure
    n_rows = max(len(field_ids), 1)
    n_cols = n_snapshots + 1  # one column for the colorbar
    figsize = (12 * n_cols, 12 * n_rows)
    max_dpi = min([2**16 // dim for dim in figsize]) - 1
    fig = plt.figure(figsize=figsize, dpi=min(max_dpi, 100))
    grid = plt.GridSpec(n_rows, n_cols, wspace=0.2, hspace=0.2)
    plt.rcParams.update({'font.size': 36})

    # plot snapshot data in each subsequent column
    for col_idx, (time_idx, time) in enumerate(zip(time_indices, snapshot_times)):
        if field_ids:
            for row_idx, field_id in enumerate(field_ids):

                ax = init_axes(
                    fig, edge_length_x, edge_length_y, grid, row_idx,
                    col_idx, time, field_id, field_label_size,
                )
                ax.tick_params(
                    axis='both', which='both', bottom=False, top=False,
                    left=False, right=False,
                )

                # transpose field to align with agents
                field = np.transpose(np.array(fields[time][field_id])).tolist()
                vmin, vmax = field_range[field_id]
                im = plt.imshow(field,
                                origin='lower',
                                extent=[0, edge_length_x, 0, edge_length_y],
                                vmin=vmin,
                                vmax=vmax,
                                cmap='BuPu')
                if agents:
                    agents_now = agents[time]
                    plot_agents(
                        ax, agents_now,agent_colors, agent_shape,
                        dead_color,
                    )


                # colorbar in new column after final snapshot
                if col_idx == n_snapshots - 1:
                    cbar_col = col_idx + 1
                    ax = fig.add_subplot(grid[row_idx, cbar_col])
                    if row_idx == 0:
                        ax.set_title('Concentration (mmol/L)', y=1.08)
                    ax.axis('off')
                    if vmin == vmax:
                        continue
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("left", size="5%", pad=0.0)
                    fig.colorbar(im, cax=cax, format='%.6f')
                    ax.axis('off')
        else:
            row_idx = 0
            ax = init_axes(
                fig, bounds[0], bounds[1], grid, row_idx, col_idx,
                time, ""
            )
            if agents:
                agents_now = agents[time]
                plot_agents(ax, agents_now, agent_colors, agent_shape)

    fig_path = os.path.join(out_dir, filename)
    fig.subplots_adjust(wspace=0.7, hspace=0.1)
    fig.savefig(fig_path, bbox_inches='tight')
    plt.close(fig)

def get_fluorescent_color(baseline_hsv, tag_color, intensity):
    # move color towards bright fluoresence color when intensity = 1
    new_hsv = baseline_hsv[:]
    distance = [a - b for a, b in zip(tag_color, new_hsv)]

    # if hue distance > 180 degrees, go around in the other direction
    if distance[0] > 0.5:
        distance[0] = 1 - distance[0]
    elif distance[0] < -0.5:
        distance[0] = 1 + distance[0]

    new_hsv = [a + intensity * b for a, b in zip(new_hsv, distance)]
    new_hsv[0] = new_hsv[0] % 1

    return new_hsv

def plot_tags(data, plot_config):
    '''Plot snapshots of the simulation over time

    The snapshots depict the agents and the levels of tagged molecules
    in each agent by agent color intensity.

    Arguments:
        data (dict): A dictionary with the following keys:

            * **agents** (:py:class:`dict`): A mapping from times to
              dictionaries of agent data at that timepoint. Agent data
              dictionaries should have the same form as the hierarchy
              tree rooted at ``agents``.
            * **config** (:py:class:`dict`): The environmental
              configuration dictionary  with the following keys:

                * **bounds** (:py:class:`tuple`): The dimensions of the
                  environment.

        plot_config (dict): Accepts the following configuration options.
            Any options with a default is optional.

            * **n_snapshots** (:py:class:`int`): Number of snapshots to
              show per row (i.e. for each molecule). Defaults to 6.
            * **out_dir** (:py:class:`str`): Output directory, which is
              ``out`` by default.
            * **filename** (:py:class:`str`): Base name of output file.
              ``tags`` by default.
            * **tagged_molecules** (:py:class:`typing.Iterable`): The
              tagged molecules whose concentrations will be indicated by
              agent color. Each molecule should be specified as a
              :py:class:`tuple` of the path in the agent compartment
              to where the molecule's count can be found, with the last
              value being the molecule's count variable.
            * **background_color** (:py:class:`str`): use matplotlib colors,
              ``black`` by default
            * **tag_label_size** (:py:class:`float`): The font size for
              the tag name label
    '''
    check_plt_backend()

    n_snapshots = plot_config.get('n_snapshots', 6)
    out_dir = plot_config.get('out_dir', 'out')
    filename = plot_config.get('filename', 'tags')
    agent_shape = plot_config.get('agent_shape', 'segment')
    background_color = plot_config.get('background_color', 'black')
    tagged_molecules = plot_config['tagged_molecules']
    tag_path_name_map = plot_config.get('tag_path_name_map', {})
    tag_label_size = plot_config.get('tag_label_size', 20)

    if tagged_molecules == []:
        raise ValueError('At least one molecule must be tagged.')

    # get data
    agents = data['agents']
    config = data.get('config', {})
    bounds = config['bounds']
    edge_length_x, edge_length_y = bounds

    # time steps that will be used
    time_vec = list(agents.keys())
    time_indices = np.round(
        np.linspace(0, len(time_vec) - 1, n_snapshots)
    ).astype(int)
    snapshot_times = [time_vec[i] for i in time_indices]

    # get tag ids and range
    tag_ranges = {}
    tag_colors = {}

    for time, time_data in agents.items():
        for agent_id, agent_data in time_data.items():
            volume = agent_data.get('boundary', {}).get('volume', 0)
            for tag_id in tagged_molecules:
                count = get_value_from_path(agent_data, tag_id)
                conc = count / volume if volume else 0
                if tag_id in tag_ranges:
                    tag_ranges[tag_id] = [
                        min(tag_ranges[tag_id][0], conc),
                        max(tag_ranges[tag_id][1], conc)]
                else:
                    # add new tag
                    tag_ranges[tag_id] = [conc, conc]

                    # select random initial hue
                    hue = random.choice(HUES)
                    tag_color = [hue] + FLOURESCENT_SV
                    tag_colors[tag_id] = tag_color

    # make the figure
    n_rows = len(tagged_molecules)
    n_cols = n_snapshots + 1  # one column for the colorbar
    figsize = (12 * n_cols, 12 * n_rows)
    max_dpi = min([2**16 // dim for dim in figsize]) - 1
    fig = plt.figure(figsize=figsize, dpi=min(max_dpi, 100))
    grid = plt.GridSpec(n_rows, n_cols, wspace=0.2, hspace=0.2)
    plt.rcParams.update({'font.size': 36})

    # plot tags
    for row_idx, tag_id in enumerate(tag_ranges.keys()):
        used_agent_colors = []
        concentrations = []
        for col_idx, (time_idx, time) in enumerate(
            zip(time_indices, snapshot_times)
        ):
            tag_name = tag_path_name_map.get(tag_id, tag_id)
            ax = init_axes(
                fig, edge_length_x, edge_length_y, grid,
                row_idx, col_idx, time, tag_name, tag_label_size,
            )
            ax.tick_params(
                axis='both', which='both', bottom=False, top=False,
                left=False, right=False,
            )
            ax.set_facecolor(background_color)

            # update agent colors based on tag_level
            min_tag, max_tag = tag_ranges[tag_id]
            agent_tag_colors = {}
            for agent_id, agent_data in agents[time].items():
                agent_color = BASELINE_TAG_COLOR

                # get current tag concentration, and determine color
                count = get_value_from_path(agent_data, tag_id)
                volume = agent_data.get('boundary', {}).get('volume', 0)
                level = count / volume if volume else 0
                if min_tag != max_tag:
                    concentrations.append(level)
                    intensity = max((level - min_tag), 0)
                    intensity = min(intensity / (max_tag - min_tag), 1)
                    tag_color = tag_colors[tag_id]
                    agent_color = get_fluorescent_color(
                        BASELINE_TAG_COLOR, tag_color, intensity)
                    agent_rgb = matplotlib.colors.hsv_to_rgb(agent_color)
                    used_agent_colors.append(agent_rgb)

                agent_tag_colors[agent_id] = agent_color

            plot_agents(ax, agents[time], agent_tag_colors, agent_shape)

            # colorbar in new column after final snapshot
            if col_idx == n_snapshots - 1:
                cbar_col = col_idx + 1
                ax = fig.add_subplot(grid[row_idx, cbar_col])
                if row_idx == 0:
                    ax.set_title('Concentration (counts/fL)', y=1.08)
                ax.axis('off')
                if min_tag == max_tag:
                    continue
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("left", size="5%", pad=0.0)
                norm = matplotlib.colors.Normalize()
                # Sort colors and concentrations by concentration
                sorted_idx = np.argsort(concentrations)
                norm.autoscale(used_agent_colors)
                cmap = matplotlib.colors.ListedColormap(
                    np.array(used_agent_colors)[sorted_idx])
                mappable = matplotlib.cm.ScalarMappable(norm, cmap)
                mappable.set_array(np.array(concentrations)[sorted_idx])
                mappable.set_clim(min_tag, max_tag)
                fig.colorbar(mappable, cax=cax, format='%.6f')

    fig_path = os.path.join(out_dir, filename)
    fig.subplots_adjust(wspace=0.7, hspace=0.1)
    fig.savefig(fig_path, bbox_inches='tight')
    plt.close(fig)

def initialize_spatial_figure(bounds, fontsize=18):

    x_length = bounds[0]
    y_length = bounds[1]

    # set up figure
    n_ticks = 4
    plot_buffer = 0.02
    buffer = plot_buffer * min(bounds)
    min_edge = min(x_length, y_length)
    x_scale = x_length/min_edge
    y_scale = y_length/min_edge

    # make the figure
    fig = plt.figure(figsize=(8*x_scale, 8*y_scale))
    plt.rcParams.update({'font.size': fontsize, "font.family": "Times New Roman"})

    plt.xlim((0-buffer, x_length+buffer))
    plt.ylim((0-buffer, y_length+buffer))
    plt.xlabel(u'\u03bcm')
    plt.ylabel(u'\u03bcm')

    # specify the number of ticks for each edge
    [x_bins, y_bins] = [int(n_ticks * edge / min_edge) for edge in [x_length, y_length]]
    plt.locator_params(axis='y', nbins=y_bins)
    plt.locator_params(axis='x', nbins=x_bins)
    ax = plt.gca()

    return fig, ax

def get_agent_trajectories(agents, times):
    trajectories = {}
    for agent_id, series in agents.items():
        time_indices = series['boundary']['location']['time_index']
        series_times = [times[time_index] for time_index in time_indices]

        positions = series['boundary']['location']['value']
        angles = series['boundary']['angle']['value']
        series_values = [[x, y, theta] for ((x, y), theta) in zip(positions, angles)]

        trajectories[agent_id] = {
            'time': series_times,
            'value': series_values,
        }
    return trajectories

def get_agent_type_colors(agent_ids):
    """ get colors for each agent id by agent type
    Assumes that agents of the same type share the beginning
    of their name, followed by '_x' with x as a single number
    TODO -- make this more general for more digits and other comparisons"""
    agent_type_colors = {}
    agent_types = {}
    for agent1, agent2 in itertools.combinations(agent_ids, 2):
        if agent1[0:-2] == agent2[0:-2]:
            agent_type = agent1[0:-2]
            if agent_type not in agent_type_colors:
                color = plt.rcParams['axes.prop_cycle'].by_key()['color'][len(agent_type_colors)]
                agent_type_colors[agent_type] = color
            else:
                color = agent_type_colors[agent_type]
            agent_types[agent1] = agent_type
            agent_types[agent2] = agent_type
    for agent in agent_ids:
        if agent not in agent_types:
            color = plt.rcParams['axes.prop_cycle'].by_key()['color'][len(agent_type_colors)]
            agent_type_colors[agent] = color
            agent_types[agent] = agent

    return agent_types, agent_type_colors

def plot_agent_trajectory(agent_timeseries, config, out_dir='out', filename='trajectory'):
    check_plt_backend()

    # trajectory plot settings
    legend_fontsize = 18
    markersize = 25

    bounds = config.get('bounds', DEFAULT_BOUNDS)
    field = config.get('field')
    rotate_90 = config.get('rotate_90', False)

    # get agents
    times = np.array(agent_timeseries['time'])
    agents = agent_timeseries['agents']
    agent_types, agent_type_colors = get_agent_type_colors(list(agents.keys()))

    if rotate_90:
        field = rotate_field_90(field)
        for agent_id, series in agents.items():
            agents[agent_id] = rotate_agent_series_90(series, bounds)
        bounds = rotate_bounds_90(bounds)

    # get each agent's trajectory
    trajectories = get_agent_trajectories(agents, times)

    # initialize a spatial figure
    fig, ax = initialize_spatial_figure(bounds, legend_fontsize)

    # move x axis to top
    ax.tick_params(labelbottom=False,labeltop=True,bottom=False,top=True)
    ax.xaxis.set_label_coords(0.5, 1.12)

    if field is not None:
        field = np.transpose(field)
        shape = field.shape
        im = plt.imshow(field,
                        origin='lower',
                        extent=[0, shape[1], 0, shape[0]],
                        cmap='Greys')
        # colorbar for field concentrations
        cbar = plt.colorbar(im, pad=0.02, aspect=50, shrink=0.7)
        cbar.set_label('concentration', rotation=270, labelpad=20)

    for agent_id, trajectory_data in trajectories.items():
        agent_trajectory = trajectory_data['value']

        # convert trajectory to 2D array
        locations_array = np.array(agent_trajectory)
        x_coord = locations_array[:, 0]
        y_coord = locations_array[:, 1]

        # get agent type and color
        agent_type = agent_types[agent_id]
        agent_color = agent_type_colors[agent_type]

        # plot line
        ax.plot(x_coord, y_coord, linewidth=2, color=agent_color, label=agent_type)
        ax.plot(x_coord[0], y_coord[0],
                 color=(0.0, 0.8, 0.0), marker='.', markersize=markersize)  # starting point
        ax.plot(x_coord[-1], y_coord[-1],
                 color='r', marker='.', markersize=markersize)  # ending point

    # create legend for agent types
    agent_labels = [
        mlines.Line2D([], [], color=agent_color, linewidth=2, label=agent_type)
        for agent_type, agent_color in agent_type_colors.items()]
    agent_legend = plt.legend(
        title='agent type', handles=agent_labels, loc='upper center',
        bbox_to_anchor=(0.3, 0.0), ncol=2, prop={'size': legend_fontsize})
    ax.add_artist(agent_legend)

    # create a legend for start/end markers
    start = mlines.Line2D([], [],
            color=(0.0, 0.8, 0.0), marker='.', markersize=markersize, linestyle='None', label='start')
    end = mlines.Line2D([], [],
            color='r', marker='.', markersize=markersize, linestyle='None', label='end')
    marker_legend = plt.legend(
        title='trajectory', handles=[start, end], loc='upper center',
        bbox_to_anchor=(0.7, 0.0), ncol=2, prop={'size': legend_fontsize})
    ax.add_artist(marker_legend)

    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.7, hspace=0.1)
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close(fig)


def rotate_bounds_90(bounds):
    return [bounds[1], bounds[0]]

def rotate_field_90(field):
    return np.rot90(field, 3)  # rotate 3 times for 270

def rotate_agent_series_90(series, bounds):
    location_series = series['boundary']['location']
    angle_series = series['boundary']['angle']

    if isinstance(location_series, dict):
        # this ran with time_indexed_timeseries_from_data
        series['boundary']['location']['value'] = [[y, bounds[0] - x] for [x, y] in location_series['value']]
        series['boundary']['angle']['value'] = [theta + PI / 2 for theta in angle_series['value']]
    else:
        series['boundary']['location'] = [[y, bounds[0] - x] for [x, y] in location_series]
        series['boundary']['angle'] = [theta + PI / 2 for theta in angle_series]
    return series

def plot_temporal_trajectory(agent_timeseries, config, out_dir='out', filename='temporal'):
    check_plt_backend()

    bounds = config.get('bounds', DEFAULT_BOUNDS)
    field = config.get('field')
    rotate_90 = config.get('rotate_90', False)

    # get agents
    times = np.array(agent_timeseries['time'])
    agents = agent_timeseries['agents']

    if rotate_90:
        field = rotate_field_90(field)
        for agent_id, series in agents.items():
            agents[agent_id] = rotate_agent_series_90(series, bounds)
        bounds = rotate_bounds_90(bounds)

    # get each agent's trajectory
    trajectories = get_agent_trajectories(agents, times)

    # initialize a spatial figure
    fig, ax = initialize_spatial_figure(bounds)

    if field is not None:
        field = np.transpose(field)
        shape = field.shape
        im = plt.imshow(field,
                        origin='lower',
                        extent=[0, shape[1], 0, shape[0]],
                        cmap='Greys'
                        )

    for agent_id, trajectory_data in trajectories.items():
        agent_trajectory = trajectory_data['value']

        # convert trajectory to 2D array
        locations_array = np.array(agent_trajectory)
        x_coord = locations_array[:, 0]
        y_coord = locations_array[:, 1]

        # make multi-colored trajectory
        points = np.array([x_coord, y_coord]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=plt.get_cmap('cool'))
        lc.set_array(times)
        lc.set_linewidth(6)

        # plot line
        line = plt.gca().add_collection(lc)

    # color bar
    cbar = plt.colorbar(line, ticks=[times[0], times[-1]], aspect=90, shrink=0.4)
    cbar.set_label('time (s)', rotation=270)

    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.7, hspace=0.1)
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close(fig)


def plot_motility(timeseries, out_dir='out', filename='motility_analysis'):
    check_plt_backend()

    expected_velocity = 14.2  # um/s (Berg)
    expected_angle_between_runs = 68 # degrees (Berg)

    # time of motor behavior without chemotaxis
    expected_run_duration = 0.42  # s (Berg)
    expected_tumble_duration = 0.14  # s (Berg)

    times = timeseries['time']
    agents = timeseries['agents']

    motility_analysis = {
        agent_id: {
            'velocity': [],
            'angular_velocity': [],
            'angle_between_runs': [],
            'angle': [],
            'thrust': [],
            'torque': [],
            'run_duration': [],
            'tumble_duration': [],
            'run_time': [],
            'tumble_time': [],
        }
        for agent_id in list(agents.keys())}

    for agent_id, agent_data in agents.items():

        boundary_data = agent_data['boundary']
        cell_data = agent_data['cell']
        previous_time = times[0]
        previous_angle = boundary_data['angle'][0]
        previous_location = boundary_data['location'][0]
        previous_run_angle = boundary_data['angle'][0]
        previous_motor_state = cell_data['motor_state'][0]  # 1 for tumble, 0 for run
        run_duration = 0.0
        tumble_duration = 0.0
        dt = 0.0

        # go through each time point for this agent
        for time_idx, time in enumerate(times):
            motor_state = cell_data['motor_state'][time_idx]
            angle = boundary_data['angle'][time_idx]
            location = boundary_data['location'][time_idx]
            thrust = boundary_data['thrust'][time_idx]
            torque = boundary_data['torque'][time_idx]

            # get velocity
            if time != times[0]:
                dt = time - previous_time
                distance = (
                    (location[0] - previous_location[0]) ** 2 +
                    (location[1] - previous_location[1]) ** 2
                        ) ** 0.5
                velocity = distance / dt  # um/sec

                angle_change = ((angle - previous_angle) / PI * 180) % 360
                if angle_change > 180:
                    angle_change = 360 - angle_change
                angular_velocity = angle_change/ dt
            else:
                velocity = 0.0
                angular_velocity = 0.0

            # get angle change between runs
            angle_between_runs = None
            if motor_state == 0:  # run
                if previous_motor_state == 1:
                    angle_between_runs = angle - previous_run_angle
                previous_run_angle = angle

            # get run and tumble durations
            if motor_state == 0:  # run
                if previous_motor_state == 1:
                    # the run just started -- save the previous tumble time and reset to 0
                    motility_analysis[agent_id]['tumble_duration'].append(tumble_duration)
                    motility_analysis[agent_id]['tumble_time'].append(time)
                    tumble_duration = 0
                elif previous_motor_state == 0:
                    # the run is continuing
                    run_duration += dt
            elif motor_state == 1:
                if previous_motor_state == 0:
                    # the tumble just started -- save the previous run time and reset to 0
                    motility_analysis[agent_id]['run_duration'].append(run_duration)
                    motility_analysis[agent_id]['run_time'].append(time)
                    run_duration = 0
                elif previous_motor_state == 1:
                    # the tumble is continuing
                    tumble_duration += dt

            # save data
            motility_analysis[agent_id]['velocity'].append(velocity)
            motility_analysis[agent_id]['angular_velocity'].append(angular_velocity)
            motility_analysis[agent_id]['angle'].append(angle)
            motility_analysis[agent_id]['thrust'].append(thrust)
            motility_analysis[agent_id]['torque'].append(torque)
            motility_analysis[agent_id]['angle_between_runs'].append(angle_between_runs)

            # save previous location and time
            previous_location = location
            previous_angle = angle
            previous_time = time
            previous_motor_state = motor_state

    # plot results
    cols = 1
    rows = 7
    fig = plt.figure(figsize=(6 * cols, 1.2 * rows))
    plt.rcParams.update({'font.size': 12})

    # plot velocity
    ax1 = plt.subplot(rows, cols, 1)
    for agent_id, analysis in motility_analysis.items():
        velocity = analysis['velocity']
        mean_velocity = np.mean(velocity)
        ax1.plot(times, velocity, label=agent_id)
        ax1.axhline(y=mean_velocity, linestyle='dashed', label='mean_' + agent_id)
    ax1.axhline(y=expected_velocity, color='r', linestyle='dashed', label='expected mean')
    ax1.set_ylabel(u'velocity \n (\u03bcm/sec)')
    ax1.set_xlabel('time')
    # ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # plot angular velocity
    ax2 = plt.subplot(rows, cols, 2)
    for agent_id, analysis in motility_analysis.items():
        angular_velocity = analysis['angular_velocity']
        ax2.plot(times, angular_velocity, label=agent_id)
    ax2.set_ylabel(u'angular velocity \n (degrees/sec)')
    ax2.set_xlabel('time')
    # ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # plot thrust
    ax3 = plt.subplot(rows, cols, 3)
    for agent_id, analysis in motility_analysis.items():
        thrust = analysis['thrust']
        ax3.plot(times, thrust, label=agent_id)
    ax3.set_ylabel('thrust')

    # plot torque
    ax4 = plt.subplot(rows, cols, 4)
    for agent_id, analysis in motility_analysis.items():
        torque = analysis['torque']
        ax4.plot(times, torque, label=agent_id)
    ax4.set_ylabel('torque')

    # plot angles between runs
    ax5 = plt.subplot(rows, cols, 5)
    for agent_id, analysis in motility_analysis.items():
        # convert to degrees
        angle_between_runs = [
            (angle / PI * 180) % 360 if angle is not None else None
            for angle in analysis['angle_between_runs']]
        # pair with time
        run_angle_points = [
            [t, angle] if angle < 180 else [t, 360 - angle]
            for t, angle in dict(zip(times, angle_between_runs)).items()
            if angle is not None]

        plot_times = [point[0] for point in run_angle_points]
        plot_angles = [point[1] for point in run_angle_points]
        mean_angle_change = np.mean(plot_angles)
        ax5.scatter(plot_times, plot_angles, label=agent_id)
        ax5.axhline(y=mean_angle_change, linestyle='dashed') #, label='mean_' + agent_id)
    ax5.set_ylabel(u'degrees \n between runs')
    ax5.axhline(y=expected_angle_between_runs, color='r', linestyle='dashed', label='expected')
    ax5.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # plot run durations
    ax6 = plt.subplot(rows, cols, 6)
    for agent_id, analysis in motility_analysis.items():
        run_duration = analysis['run_duration']
        run_time = analysis['run_time']
        mean_run_duration = np.mean(run_duration)
        ax6.scatter(run_time, run_duration, label=agent_id)
        ax6.axhline(y=mean_run_duration, linestyle='dashed')
    ax6.set_ylabel('run \n duration \n (s)')
    ax6.axhline(y=expected_run_duration, color='r', linestyle='dashed', label='expected')
    ax6.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # plot tumble durations
    ax7 = plt.subplot(rows, cols, 7)
    for agent_id, analysis in motility_analysis.items():
        tumble_duration = analysis['tumble_duration']
        tumble_time = analysis['tumble_time']
        mean_tumble_duration = np.mean(tumble_duration)
        ax7.scatter(tumble_time, tumble_duration, label=agent_id)
        ax7.axhline(y=mean_tumble_duration, linestyle='dashed')
    ax7.set_ylabel('tumble \n duration \n (s)')
    ax7.axhline(y=expected_tumble_duration, color='r', linestyle='dashed', label='expected')
    ax7.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    fig_path = os.path.join(out_dir, filename)
    plt.subplots_adjust(wspace=0.7, hspace=0.4)
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close(fig)


def init_axes(
    fig, edge_length_x, edge_length_y, grid, row_idx, col_idx, time,
    molecule, ylabel_size=20
):
    ax = fig.add_subplot(grid[row_idx, col_idx])
    if row_idx == 0:
        plot_title = 'time: {:.4f} s'.format(float(time))
        plt.title(plot_title, y=1.08)
    if col_idx == 0:
        ax.set_ylabel(
            molecule, fontsize=ylabel_size, rotation='horizontal',
            horizontalalignment='right',
        )
    ax.set(xlim=[0, edge_length_x], ylim=[0, edge_length_y], aspect=1)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    return ax
