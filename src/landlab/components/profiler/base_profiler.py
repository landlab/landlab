# ! /usr/env/python
"""Base class for profile constructors."""

from abc import ABC
from abc import abstractmethod

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
    >>> struct = [[1, 2, 3, 4], [[2, 3, 4, 5], [3, 4, 5, 6]], [4, 5, 6, 7]]
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
    >>> struct = [[1, 2, 3, 4], [[2, 3, 4, 5], [3, 4, 5, 6]], [4, 5, 6, 7]]
    >>> _recursive_min(struct)
    1
    >>> _recursive_min([100])
    100
    """
    return min(_recursive_min(j) if hasattr(j, "__iter__") else j for j in jagged)


class _BaseProfiler(ABC, Component):
    """Base class to handle profilers.

    Primarily exists to handle plotting.
    """

    _name = "_BaseProfiler"

    _unit_agnostic = True

    _info = {}

    def __init__(self, grid):
        super().__init__(grid)

    def run_one_step(self):
        """Calculate the profile data structure and distances along it."""
        # calculate the profile IDs data structure.
        self._create_profile_structure()

    @abstractmethod
    def _create_profile_structure(self):
        """Private class for creating profile structure.

        Expectation is that this will be overridden to create the following
        three private attributes:

        self._nodes
        self._distance_along_profile

        are each lists of numpy arrays, one array per segment.

        self._colors

        is a list of RGBA tuples, one tuple per segment.

        The order of segments is expected to be consistent between each of the
        three data structures.
        """
        ...  # pragma: no cover

    @property
    def distance_along_profile(self):
        """List of distances along profile for each segment.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import (
        ...     FastscapeEroder,
        ...     FlowAccumulator,
        ...     ChannelProfiler,
        ... )
        >>> mg = RasterModelGrid((10, 10), xy_spacing=10)
        >>> np.random.seed(42)
        >>> z = mg.add_zeros("topographic__elevation", at="node")
        >>> z[mg.core_nodes] += np.random.randn(mg.core_nodes.size)
        >>> fa = FlowAccumulator(mg)
        >>> sp = FastscapeEroder(mg, K_sp=0.0001)
        >>> dt = 1000
        >>> for i in range(200):
        ...     fa.run_one_step()
        ...     sp.run_one_step(dt=dt)
        ...     z[mg.core_nodes] += 0.001 * dt
        ...
        >>> profiler = ChannelProfiler(mg)
        >>> profiler.run_one_step()
        >>> profiler.distance_along_profile
        [array([  0.,  10.,  20.,  30.,  40.,  50.])]
        """
        return self._distance_along_profile

    @property
    def nodes(self):
        """List of node ids for each segment.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import (
        ...     FastscapeEroder,
        ...     FlowAccumulator,
        ...     ChannelProfiler,
        ... )
        >>> mg = RasterModelGrid((10, 10), xy_spacing=10)
        >>> np.random.seed(42)
        >>> z = mg.add_zeros("topographic__elevation", at="node")
        >>> z[mg.core_nodes] += np.random.randn(mg.core_nodes.size)
        >>> fa = FlowAccumulator(mg)
        >>> sp = FastscapeEroder(mg, K_sp=0.0001)
        >>> dt = 1000
        >>> for i in range(200):
        ...     fa.run_one_step()
        ...     sp.run_one_step(dt=dt)
        ...     z[mg.core_nodes] += 0.001 * dt
        ...
        >>> profiler = ChannelProfiler(mg)
        >>> profiler.run_one_step()
        >>> profiler.nodes
        [array([59, 58, 57, 56, 46, 45])]
        """
        return self._nodes

    @property
    def colors(self):
        """List of colors for each segment.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import (
        ...     FastscapeEroder,
        ...     FlowAccumulator,
        ...     ChannelProfiler,
        ... )
        >>> mg = RasterModelGrid((10, 10), xy_spacing=10)
        >>> np.random.seed(42)
        >>> z = mg.add_zeros("topographic__elevation", at="node")
        >>> z[mg.core_nodes] += np.random.randn(mg.core_nodes.size)
        >>> fa = FlowAccumulator(mg)
        >>> sp = FastscapeEroder(mg, K_sp=0.0001)
        >>> dt = 1000
        >>> for i in range(200):
        ...     fa.run_one_step()
        ...     sp.run_one_step(dt=dt)
        ...     z[mg.core_nodes] += 0.001 * dt
        ...
        >>> profiler = ChannelProfiler(mg)
        >>> profiler.run_one_step()
        >>> np.round(profiler.colors, decimals=2)
        array([[0.27, 0.  , 0.33, 1.  ]])
        """
        return self._colors

    def plot_profiles(
        self,
        field="topographic__elevation",
        xlabel="Distance Along Profile",
        ylabel="Plotted Quantity",
        title="Extracted Profiles",
        color=None,
    ):
        """Plot distance-upstream vs at at-node or size (nnodes,) quantity.

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
        color : RGBA tuple or color string
            Color to use in order to plot all profiles the same color. Default
            is None, and the colors assigned to each profile are used.
        """
        quantity = return_array_at_node(self._grid, field)

        # create segments the way that line collection likes them.
        segments = []
        qmin = []
        qmax = []
        for idx, nodes in enumerate(self._nodes):
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

        line_segments = LineCollection(segments)
        colors = color or self._colors
        line_segments.set_color(colors)
        ax.add_collection(line_segments)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

    def plot_profiles_in_map_view(
        self, field="topographic__elevation", endpoints_only=False, color=None, **kwds
    ):
        """Plot profile locations in map view.

        Parameters
        ----------
        field : field name or nnode array
            Array of the at-node-field to plot as the 2D map values.
            Default value is the at-node field 'topographic__elevation'.
        endpoints_only : boolean
            Boolean where False (default) indicates every node along the
            profile is plotted, or True indicating only segment endpoints are
            plotted.
        color : RGBA tuple or color string
            Color to use in order to plot all profiles the same color. Default
            is None, and the colors assigned to each profile are used.
        **kwds : dictionary
            Keyword arguments to pass to imshow_grid.
        """
        # make imshow_grid background
        imshow_grid(self._grid, field, **kwds)
        ax = plt.gca()

        # create segments the way that line collection likes them.
        segments = []
        for nodes in self._nodes:
            if endpoints_only:
                select_nodes = [nodes[0], nodes[-1]]
                segments.append(
                    list(
                        zip(
                            self._grid.x_of_node[select_nodes],
                            self._grid.y_of_node[select_nodes],
                        )
                    )
                )

            else:
                segments.append(
                    list(zip(self._grid.x_of_node[nodes], self._grid.y_of_node[nodes]))
                )

        line_segments = LineCollection(segments)
        colors = color or self._colors
        line_segments.set_color(colors)
        ax.add_collection(line_segments)
