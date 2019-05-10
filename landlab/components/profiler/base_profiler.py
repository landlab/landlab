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
    >>> c = [[(1, 2, 3), (1, 2, 3)],
    ...           [[(1, 2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3)],
    ...            [(1, 2, 3), (1, 2, 3), (1, 2, 3)]],
    ...           [(1, 2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3)]]
    >>> out = _flatten_structure(c)
    >>> out
    [(1, 2, 3), (1, 2, 3),
     (1, 2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3),
     (1, 2, 3), (1, 2, 3), (1, 2, 3),
     (1, 2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3)]
    """
    if l is None:
        return None
    if isinstance(l[0], (np.number, int, float)):
        return [l]
    else:
        return list(chain(*map(_flatten_structure, l)))


def _verify_structure_and_color(profile_structure, colors):
    """
    Color structure must either be
        a) None - then matplotlib default color cycle will be used
        b) One RGBA tuple. Then all segments will be that color.
        c) One RGBA tuple per highest level of the profile structure.
            Then all segments within that level (e.g. watershed for the
            ChannelProfiler) will be the same level.

    Examples
    --------
    >>> from landlab.components.profiler.base_profiler import (
    ...     _verify_structure_and_color)

    If both are None

    >>> ps, c = _verify_structure_and_color(None, None)
    >>> assert ps is None
    >>> assert c is None

    If color is a list of one tuple.

    >>> ps, c = _verify_structure_and_color(None, [(0, 1, 0, 1)])
    >>> assert ps is None
    >>> c
    [(0, 1, 0, 1)]

    If color is a tuple.

    >>> ps, c = _verify_structure_and_color(None, (0, 1, 0, 1))
    >>> assert ps is None
    >>> c
    (0, 1, 0, 1)

    If structure and color are perfectly parallel data structures.
    >>> struct = [[1, 2, 3, 4],
    ...           [[2, 3, 4, 5],
    ...            [3, 4, 5, 6]],
    ...           [4, 5, 6, 7]]
    >>> color = [(1, 1, 1, 1),
    ...           [(1, 1, 1, 0),
    ...            (1, 1, 0, 1)],
    ...           (0, 1, 0, 1)]
    >>> ps, c = _verify_structure_and_color(struct, color)
    >>> ps
    [[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6], [4, 5, 6, 7]]
    >>> c
    [(1, 1, 1, 1), (1, 1, 1, 0), (1, 1, 0, 1), (0, 1, 0, 1)]

    If color is parallel with the top level of hierarchy.

    >>> color = [(1, 1, 1, 1),
    ...          (1, 1, 1, 0),
    ...          (1, 1, 0, 1)]
    >>> ps, c = _verify_structure_and_color(struct, color)
    >>> c
    [(1, 1, 1, 1), (1, 1, 1, 0), (1, 1, 1, 0), (1, 1, 0, 1)]

    Some bad cases

    >>> import pytest

    If top level of color  has too few elements.

    >>> color = [(1, 1, 1, 1),
    ...          (1, 1, 1, 0)]
    >>> with pytest.raises(ValueError):
    ...     _verify_structure_and_color(struct, color)

    If a lower level of color has more than one, but fewer than correct
    elements.

    >>> struct = [[1, 2, 3, 4],
    ...           [[2, 3, 4, 5],
    ...            [3, 4, 5, 6],
    ...            [6, 7, 8, 9]],
    ...           [4, 5, 6, 7]]
    >>> color = [(1, 1, 1, 1),
    ...           [(1, 1, 1, 0),
    ...            (1, 1, 0, 1)],
    ...           (0, 1, 0, 1)]
    >>> with pytest.raises(ValueError):
    ...     _verify_structure_and_color(struct, color)

    """
    if (
        (colors is None)
        or (len(colors) == 1 and isinstance(colors[0], tuple))
        or (isinstance(colors, tuple))
    ):
        return (_flatten_structure(profile_structure), colors)

    else:
        new_colors = []
        new_struct = []

        if len(profile_structure) != len(colors):
            msg = (
                "Number of colors is different than the number of "
                "top level elements in the profile structure (watersheds "
                "if using the ChannelProfiler). The number must either be "
                "the same, or you must provide only one color."
            )
            raise ValueError(msg)

        for i in range(len(profile_structure)):
            p = _flatten_structure(profile_structure[i])
            c = _flatten_structure(colors[i])

            if len(c) == 1:
                c = c * len(p)
            if len(c) != len(p):
                msg = (
                    "Number of colors is different than the number of "
                    "segments in the profile structure (channel segments per "
                    "watershed if using the ChannelProfiler)."
                )
                raise ValueError(msg)
            new_struct.extend(p)
            new_colors.extend(c)

        return (new_struct, new_colors)


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

        The colors keyword argument can handle the following options.

            a) None - The matplotlib defaults will be used.
            b) One RGBA tuple. Then all segments of the profiles will be that
               color.
            c) One RGBA tuple per highest level of the profile data structure.
               Then all segments within that level (e.g. watershed for the
               ChannelProfiler) will be the same color.
            d) One RGBA tuple per segment in the profile data structure. The
               structure of the color keyword argument data structure must
               mirror that of the profile data structure.

        Parameters
        ----------
        field : field name or nnode array
            Array of  the at-node-field to plot against distance upstream.
            Default value is the at-node field 'topographic__elevation'.
        colors : sequence of RGBA tuples, optional
            Sequence of RGBA tuples to use with each stream segment. See above.
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
        node_ids, colors = _verify_structure_and_color(self._profile_structure, colors)

        # create segments the way that line collection likes them.
        segments = []
        qmin = []
        qmax = []
        for idx, nodes in enumerate(node_ids):
            segments.append(list(zip(x_dist[idx], quantity[nodes])))
            qmin.append(min(quantity[nodes]))
            qmax.append(max(quantity[nodes]))

        # We need to set the plot limits.
        ax = plt.gca()
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

        The colors keyword argument can handle the following options.

            a) None - The matplotlib defaults will be used.
            b) One RGBA tuple. Then all segments of the profiles will be that
               color.
            c) One RGBA tuple per highest level of the profile data structure.
               Then all segments within that level (e.g. watershed for the
               ChannelProfiler) will be the same color.
            d) One RGBA tuple per segment in the profile data structure. The
               structure of the color keyword argument data structure must
               mirror that of the profile data structure.

        Parameters
        ----------
        field : field name or nnode array
            Array of  the at-node-field to plot as the 2D map values.
            Default value is the at-node field 'topographic__elevation'.
        colors : sequence of RGB tuples, optional
            Sequence of RGB tuples to use with each stream segment. See above.
        **kwargs : keyword arguments, arguments
            Additional parameters to pass to imshow_grid
        """
        # make imshow_grid background
        imshow_grid(self._grid, field, **kwargs)
        ax = plt.gca()

        segments = []

        # flatten datastructure
        node_ids, colors = _verify_structure_and_color(self._profile_structure, colors)

        # create segments the way that line collection likes them.
        segments = []
        for idx, nodes in enumerate(node_ids):
            segments.append(
                list(zip(self._grid.x_of_node[nodes], self._grid.y_of_node[nodes]))
            )

        line_segments = LineCollection(segments, colors=colors)
        ax.add_collection(line_segments)
