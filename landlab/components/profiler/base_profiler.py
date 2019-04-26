# coding: utf8
# ! /usr/env/python
"""
"""

from itertools import chain

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from landlab import Component
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


class _BaseProfiler(Component):
    """Base class to handle profilers.

    Primarily exists to handle plotting.
    """

    def __init__(self, grid, stopping_field):
        super(_BaseProfiler, self).__init__(grid)

    def plot_profiles(
        self,
        field="topographic__elevation",
        colors=None,
        xlabel="Distance Upstream",
        ylabel="Plotted Quantity",
        title="Channel Long Profile",
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
            X-axis label, default is "Distance Upstream".
        ylabel : str, optional
            Y-axis label, default value is "Plotted Quantity".
        title : str, optional
            Plot title, default value is "Channel Long Profile".
        """
        quantity = return_array_at_node(self._grid, field)

        # flatten datastructure
        x_dist = _flatten_structure(self._distance_along_profile)
        node_ids = _flatten_structure(self._profile_structure)
        colors = _flatten_color(colors)

        # create segments the way that line collection likes them.
        segments = []
        for idx, nodes in enumerate(node_ids):
            segments.append(list(zip(x_dist[idx], quantity[nodes])))

        # We need to set the plot limits.
        fig, ax = plt.subplots()
        ax.set_xlim(
            min(min(min(self._distance_along_profile))),
            max(max(max(self._distance_along_profile))),
        )
        ax.set_ylim(quantity.min(), quantity.max())

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
