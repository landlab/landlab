# coding: utf8
# ! /usr/env/python
"""
"""

from itertools import chain

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from landlab import Component, RasterModelGrid
from landlab.plot import imshow_grid
from landlab.utils.return_array import return_array_at_node


def _flatten_structure(l):
    """
    Examples
    --------
    >>> from landlab.components.profiler.base_profiler import (
    ...     _flatten_structure)
    >>> struct = [[1, 2, 3, 4],
    ...           [[2, 3, 4, 5],
    ...            [3, 4, 5, 6]],
    ...           [4, 5, 6, 7]]
    >>> out = _flatten_structure(struct)
    >>> out
    [[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6], [4, 5, 6, 7]]
    >>> assert _flatten_structure(None) is None
    """
    if l is None:
        return None
    if isinstance(l[0], (np.number, int, float)):
        return [l]
    else:
        return list(chain(*map(_flatten_structure, l)))


def _flatten_color(l):
    """
    Examples
    --------
    >>> from landlab.components.profiler.base_profiler import _flatten_color
    >>> c = [[(1, 2, 3), (1, 2, 3)],
    ...           [[(1, 2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3)],
    ...            [(1, 2, 3), (1, 2, 3), (1, 2, 3)]],
    ...           [(1, 2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3)]]
    >>> out = _flatten_color(c)
    >>> out
    [[(1, 2, 3), (1, 2, 3)],
     [(1, 2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3)],
     [(1, 2, 3), (1, 2, 3), (1, 2, 3)],
     [(1, 2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3)]]
    >>> assert _flatten_color(None) is None
    """
    if l is None:
        return None
    if isinstance(l[0], tuple):
        return [l]
    else:
        return list(chain(*map(_flatten_color, l)))


def _recursive_max(jagged):
    """
    Examples
    --------
    from landlab.components.profiler.base_profiler import _recursive_max
    >>> struct = [[1, 2, 3, 4],
    ...           [[2, 3, 4, 5],
    ...            [3, 4, 5, 6]],
    ...           [4, 5, 6, 7]]
    >>> _recursive_max(struct)
    7
    >>> _recursive_max([100])
    100
    """
    return max(_recursive_max(j) if hasattr(j, "__iter__") else j for j in jagged)


def _recursive_min(jagged):
    """
    Examples
    --------
    from landlab.components.profiler.base_profiler import _recursive_min
    >>> struct = [[1, 2, 3, 4],
    ...           [[2, 3, 4, 5],
    ...            [3, 4, 5, 6]],
    ...           [4, 5, 6, 7]]
    >>> _recursive_min(struct)
    1
    >>> _recursive_min([100])
    100
    """
    return min(_recursive_min(j) if hasattr(j, "__iter__") else j for j in jagged)


class _BaseProfiler(Component):
    """Base class to handle profilers.

    Primarily exists to handle plotting.
    """

    _name = "_BaseProfiler"

    _input_var_names = ()

    _output_var_names = ()

    _var_units = {}

    _var_mapping = {}

    _var_doc = {}

    def __init__(self, grid, stopping_field):
        super(_BaseProfiler, self).__init__(grid)

    def run_one_step(self):
        """Calculate the profile datastructure and distances along it."""
        # calculate the profile IDs datastructure
        self._create_profile_structure()

        # calculate the distance along profile datastructure
        self._calculate_distances()

    def _calculate_distances(self):
        """Get distances along the profile_IDs datastructure."""
        self._distance_along_profile = []
        end_distances = {}

        # set the starting values for the beginnings of each netwrok.
        for network in self._profile_structure:
            starting_node = network[0][0]
            end_distances[starting_node] = 0

        # for each network
        for network in self._profile_structure:

            network_values = []
            # for each segment in the network.
            for segment in network:
                starting_node = segment[0]

                total_distance = end_distances[starting_node]

                profile_values = []
                profile_values.append(total_distance)

                # itterate up the profile
                for j in range(len(segment) - 1):
                    if isinstance(self._grid, RasterModelGrid):
                        total_distance += self._grid.length_of_d8[
                            self._link_to_flow_receiver[segment[j + 1]]
                        ]
                        profile_values.append(total_distance)
                    else:
                        total_distance += self._grid.length_of_link[
                            self._link_to_flow_receiver[segment[j + 1]]
                        ]
                        profile_values.append(total_distance)
                network_values.append(np.array(profile_values))
                end_distances[segment[-1]] = total_distance
            self._distance_along_profile.append(network_values)

    def plot_profiles(
        self,
        field="topographic__elevation",
        colors=None,
        xlabel="Distance Along Profile",
        ylabel="Plotted Quantity",
        title="Extracted Profiles",
    ):
        """
        Plot distance-upstream vs at at-node or size (nnodes,) quantity.

        Parameters
        ----------
        field : field name or nnode array
            Array of  the at-node-field to plot against distance upstream.
            Default value is the at-node field 'topographic__elevation'.
        colors : sequence of RGB tuples, optional
            Sequence of RGB tuples to use with each stream segment. Assumed
            to have the same structure as the profile datastructure.
        xlabel : str, optional
            X-axis label, default is "Distance Along Profile".
        ylabel : str, optional
            Y-axis label, default value is "Plotted Quantity".
        title : str, optional
            Plot title, default value is "Extracted Profiles".
        """
        quantity = return_array_at_node(self._grid, field)

        # flatten datastructure
        x_dist = _flatten_structure(self._distance_along_profile)
        node_ids = _flatten_structure(self._profile_structure)
        colors = _flatten_color(colors)

        # create segments the way that line collection likes them.
        segments = []
        qmin = []
        qmax = []
        for idx, nodes in enumerate(node_ids):
            segments.append(list(zip(x_dist[idx], quantity[nodes])))
            qmin.append(min(quantity[nodes]))
            qmax.append(max(quantity[nodes]))

        # We need to set the plot limits.
        fig, ax = plt.subplots()
        ax.set_xlim(_recursive_min(x_dist), _recursive_max(x_dist))
        ax.set_ylim(min(qmin), max(qmax))

        line_segments = LineCollection(segments, colors=colors)
        ax.add_collection(line_segments)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

    def plot_profiles_in_map_view(
        self, field="topographic__elevation", colors=None, **kwargs
    ):
        """
        Plot profile locations in map view.

        Parameters
        ----------
        field : field name or nnode array
            Array of  the at-node-field to plot as the 2D map values.
            Default value is the at-node field 'topographic__elevation'.
        colors : f..
        **kwargs : keyword arguments, arguments
            Additional parameters to pass to imshow_grid
        """
        # make imshow_grid background
        imshow_grid(self._grid, field, **kwargs)
        ax = plt.gca()

        segments = []

        # flatten datastructure
        node_ids = _flatten_structure(self._profile_structure)
        colors = _flatten_color(colors)

        # create segments the way that line collection likes them.
        segments = []
        for idx, nodes in enumerate(node_ids):
            segments.append(
                list(zip(self._grid.x_of_node[nodes], self._grid.y_of_node[nodes]))
            )

        line_segments = LineCollection(segments, colors=colors)
        ax.add_collection(line_segments)
