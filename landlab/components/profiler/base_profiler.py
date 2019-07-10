# coding: utf8
# ! /usr/env/python
"""Base class for profile constructors."""

from abc import ABC, abstractmethod

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from landlab import Component
from landlab.plot import imshow_grid
from landlab.utils.return_array import return_array_at_node


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


class _BaseProfiler(Component, ABC):
    """Base class to handle profilers.

    Primarily exists to handle plotting.
    """

    _name = "_BaseProfiler"

    _input_var_names = ()

    _output_var_names = ()

    _var_units = {}

    _var_mapping = {}

    _var_doc = {}

    def __init__(self, grid):
        super(_BaseProfiler, self).__init__(grid)

    def run_one_step(self):
        """Calculate the profile datastructure and distances along it."""
        # calculate the profile IDs datastructure.
        self._create_profile_structure()

    @abstractmethod
    def _create_profile_structure(self):
        """Private class for creating profile structure.

        Expectation is that this will be overridden to create the following
        three private attributes:

        self._net_ids
        self._distance_along_profile

        are each lists of numpy arrays, one array per segment.

        self._colors

        is a list of RGBA tuples, one tuple per segment.

        The order of segments is expected to be consistent between each of the
        three datastructures.
        """
        ...  # pragma: no cover

    @property
    def distance_along_profile(self):
        """List of distances along profile for each segment."""
        return self._distance_along_profile

    @property
    def network_ids(self):
        """List of node ids for each segment."""
        return self._net_ids

    @property
    def colors(self):
        """List of colors for each segment."""
        return self._colors

    def plot_profiles(
        self,
        field="topographic__elevation",
        xlabel="Distance Along Profile",
        ylabel="Plotted Quantity",
        title="Extracted Profiles",
    ):
        """
        Plot distance-upstream vs at at-node or size (nnodes,) quantity.

        Parameters
        ----------
        field : field name or nnode array
            Array of the at-node-field to plot against distance upstream.
            Default value is the at-node field 'topographic__elevation'.
        xlabel : str, optional
            X-axis label, default is "Distance Along Profile".
        ylabel : str, optional
            Y-axis label, default value is "Plotted Quantity".
        title : str, optional
            Plot title, default value is "Extracted Profiles".
        """
        quantity = return_array_at_node(self._grid, field)

        # create segments the way that line collection likes them.
        segments = []
        qmin = []
        qmax = []
        for idx, nodes in enumerate(self._net_ids):
            segments.append(
                list(zip(self._distance_along_profile[idx], quantity[nodes]))
            )
            qmin.append(min(quantity[nodes]))
            qmax.append(max(quantity[nodes]))

        # We need to set the plot limits.
        ax = plt.gca()
        ax.set_xlim(
            _recursive_min(self._distance_along_profile),
            _recursive_max(self._distance_along_profile),
        )
        ax.set_ylim(min(qmin), max(qmax))

        line_segments = LineCollection(segments, colors=self._colors)
        ax.add_collection(line_segments)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

    def plot_profiles_in_map_view(self, field="topographic__elevation", **kwargs):
        """
        Plot profile locations in map view.

        Parameters
        ----------
        field : field name or nnode array
            Array of the at-node-field to plot as the 2D map values.
            Default value is the at-node field 'topographic__elevation'.
        """
        # make imshow_grid background
        imshow_grid(self._grid, field, **kwargs)
        ax = plt.gca()

        # create segments the way that line collection likes them.
        segments = []
        for idx, nodes in enumerate(self._net_ids):
            segments.append(
                list(zip(self._grid.x_of_node[nodes], self._grid.y_of_node[nodes]))
            )

        line_segments = LineCollection(segments, colors=self._colors)
        ax.add_collection(line_segments)
