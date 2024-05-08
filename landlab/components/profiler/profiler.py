# ! /usr/env/python
"""profiler.py component to create profiles with user-defined endpoints."""
from collections import OrderedDict

import numpy as np
from matplotlib import cm
from matplotlib import colors
from matplotlib import pyplot as plt

from landlab.components.profiler.base_profiler import _BaseProfiler


class Profiler(_BaseProfiler):
    """Extract and plot profiles set up using points within a grid.

    The profile is constructed from the first to final point in ``endpoints``.
    Endpoints are located at grid nodes. Two successive endpoints bound a
    profile segment. A profile with one segment is a straight line. The
    segments of a profile with multiple segments meet at endpoints. The grid
    nodes along the profile are sampled, including the segment endpoints. The
    extracted quantity of the node is retained. No interpolation is conducted
    even for profile traces that reside between nodes.

    The structure of the profile in a model grid is diagrammed below. The grid
    contains nine columns and nine rows. The profile is constructed from three
    endpoints that bound two segments. Here, ``o`` indicates a segment
    endpoint, ``.`` and ``*`` are sample nodes of the first and second segment,
    respectively. ``X`` are nodes not included in the profile. The first
    segment begins in the lower-left and continues horizontally and almost
    reaches the right boundary. The second segment is joined to the first in
    the lower-right of the grid and it continues diagonally to the upper-left.
    Segments have seven sample points each (nodes at endpoints are also
    sampled). The segments share the second endpoint. Segment and sample
    ordering is dictated by the ordering of endpoints. If the horizontal
    segment is the first segment, the endpoints used to construct this profile
    must be ordered: lower-left, lower-right, and then upper-left.::

        X X X X X X X X X
        X o X X X X X X X
        X X * X X X X X X
        X X X * X X X X X
        X X X X * X X X X
        X X X X X * X X X
        X X X X X X * X X
        X o . . . . . o X
        X X X X X X X X X

    The node IDs and distances along the profile are stored in a data structure
    called ``data_structure``. It is a dictionary with keys indicating the
    segment IDs that are enumerated along the profile.

    By default, a unique color will be assigned to each segment. To change the
    color, a user can change values stored in ``data_structure``. Additionally,
    a ``cmap`` keyword argument can provide some user control over the color at
    the instantiation of the component.

    The data structure of the example above will look as follows:

    .. code-block:: python

        {
            0: {
                "ids": [10, 11, 12, 13, 14, 15, 16],
                "distances": [0, 1, 2, 3, 4, 5, 6],
                "color": (0.27, 0, 0.33, 1),
            },
            1: {
                "ids": [16, 24, 32, 40, 48, 56, 64],
                "distances": [6, 7.41, 8.83, 10.24, 11.66, 13.07, 14.49],
                "color": (0.13, 0.57, 0.55, 1),
            },
        }


    Examples
    --------
    Create a model grid with the same dimensions as the diagram above.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components import Profiler
    >>> import numpy as np
    >>> mg = RasterModelGrid((10, 10), 10)
    >>> mg.at_node["topographic__elevation"] = mg.node_x * mg.node_y

    Create a profile with three endpoints. This profile is laid out the same as
    the diagram above.

    >>> endpoints = [10, 16, 64]
    >>> profiler = Profiler(mg, endpoints)
    >>> profiler.run_one_step()

    The keys of the data structure are the segment ids.

    >>> profiler.data_structure.keys()
    odict_keys([0, 1])

    The data structure contains data of segment samples. Below is the first
    segment.

    >>> profiler.data_structure[0]["ids"]
    array([10, 11, 12, 13, 14, 15, 16])
    >>> profiler.data_structure[0]["distances"]
    array([ 0., 10., 20., 30., 40., 50., 60.])
    >>> np.round(profiler.data_structure[0]["color"], decimals=2)
    array([0.27, 0.  , 0.33, 1.  ])

    Note that the first node of the second segment is the same as the final
    node of the first segment.

    >>> profiler.data_structure[1]["ids"]
    array([16, 26, 35, 45, 54, 64])

    Alternative to nodes, profiles can be instantiated with coordinates.

    >>> profiler = Profiler(mg, [(10, 10), (70, 10), (10, 70)])

    Endpoints can also be set with a combination of coordinates and nodes.

    >>> profiler = Profiler(mg, [(10, 10), 16, (10, 70)])

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    None Listed

    """

    _name = "Profiler"

    _unit_agnostic = True

    def __init__(self, grid, endpoints, cmap="viridis"):
        """Instantiate Profiler.

        Parameters
        ----------
        grid : RasterModelGrid
            A landlab RasterModelGrid.
        endpoints : list of node id integers or coordinate tuples
            The endpoints that bound segments of the profile. Endpoints can be
            node ids and/or tuples of coordinates (x, y, where these
            coordinates are the measurement from the grid lower-left). The list
            can be a mix of node ids and coordinate tuples. The profile begins
            with the first element of `endpoints` and continues in the order of
            this list.
        cmap : str
            A valid matplotlib cmap string. Default is "viridis".
        """
        super().__init__(grid)

        self._cmap = plt.colormaps[cmap]

        if not isinstance(endpoints, list) or len(endpoints) < 2:
            raise ValueError(
                "`endpoints` must be a list of at least 2 node IDs or a "
                "list of at least two tuples where each tuple contains the "
                "x, y coordinates of endpoints."
            )

        # Check if `endpoints` are within grid bounds while setting
        # `_end_nodes`.

        self._end_nodes = []
        for point in endpoints:
            node, _ = self._get_node_and_coords(point)
            self._end_nodes.append(node)

    @property
    def data_structure(self):
        """OrderedDict defining the profile.

        The node IDs and distances along the profile are stored in
        ``data_structure``. It is a dictionary with keys of the segment ID.
        The value of each key is itself a dictionary of the segment attributes.
        First, 'ids' contains a list of the node IDs of segment samples ordered
        from the start to the end of the segment. It includes the endpoints.
        Second, 'distances' contains a list of along-profile distances that
        mirrors the list in 'ids'. Finally, 'color' is an RGBA tuple indicating
        the color for the segment.
        """
        return self._data_struct

    def _create_profile_structure(self):
        """Create the data structure of the profile.

        The profile is processed by segment. Segments are bound by successive
        endpoints. The cumulative distance along the profile is accumulated by
        iteratively adding segment lengths.
        """
        self._data_struct = OrderedDict()
        grid = self._grid
        endnodes = self._end_nodes
        cum_dist = 0

        for i_endpt in range(len(endnodes) - 1):
            # Get the endpoints and samples of the segment.

            start_node, start_xy = self._get_node_and_coords(endnodes[i_endpt])
            end_node, end_xy = self._get_node_and_coords(endnodes[i_endpt + 1])

            sample_nodes = self._get_sample_nodes(start_node, end_node)

            # Calculate the along-profile distance of samples along the
            # segment.

            n_samples = len(sample_nodes)
            sample_distances = np.empty(n_samples, dtype=float)

            for i_sample, node in enumerate(sample_nodes):
                sample_xy = grid.xy_of_node[node]

                pt = self._project_point_onto_line(sample_xy, start_xy, end_xy)
                d = grid.calc_distances_of_nodes_to_point(pt, node_subset=start_node)
                sample_distances[i_sample] = d

            # Store the segment data.

            self._data_struct[i_endpt] = {
                "ids": np.array(sample_nodes),
                "distances": sample_distances + cum_dist,
            }

            cum_dist += max(sample_distances)

        self._assign_colors()
        self._create_flat_structures()

    def _assign_colors(self, color_mapping=None):
        """Assign a unique color for each segment.

        Parameters
        ----------
        color_mapping : str
            Color map name.
        """
        if color_mapping is None:
            segment_count = len(self._data_struct)
            norm = colors.Normalize(vmin=0, vmax=segment_count)
            mappable = cm.ScalarMappable(norm=norm, cmap=self._cmap)
            color_mapping = {
                segment_id: mappable.to_rgba(idx)
                for idx, segment_id in enumerate(self._data_struct)
            }

        for segment_id in self._data_struct:
            self._data_struct[segment_id]["color"] = color_mapping[segment_id]

    def _create_flat_structures(self):
        """Create expected flattened structures for ids, distances, and colors."""
        self._nodes = []
        self._distance_along_profile = []
        self._colors = []

        for segment_id in self._data_struct:
            self._nodes.append(self._data_struct[segment_id]["ids"])
            self._distance_along_profile.append(
                self._data_struct[segment_id]["distances"]
            )
            self._colors.append(self._data_struct[segment_id]["color"])

    def _get_node_and_coords(self, point):
        """Get the node and coordinates for a point.

        This method handles the option that endpoints can be a node or tuple.
        The grid methods called here verify if the point is within the grid.
        """
        if isinstance(point, (float, int, np.integer)):
            return point, self._grid.xy_of_node[point]
        elif isinstance(point, (tuple, list, np.ndarray)) and len(point) == 2:
            return self._grid.find_nearest_node(point), point
        else:
            raise TypeError(
                "each element of `endpoints` must be a number "
                "representing a node id or a tuple of node x, y "
                "coordinates"
            )

    def _get_sample_nodes(self, start_node, end_node):
        """Get the profile sample nodes using Bresenham's line algorithm.

        Parameters
        ----------
        start_node, end_node : integer
            The node id of a profile endpoint.

        Returns
        -------
        list of integers
            The node ids of the profile samples.

        Notes
        -----
        See: https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
        """
        # Get node column and row numbers to act as the coordinates.

        y0, x0 = np.argwhere(self._grid.nodes == start_node)[0]
        y1, x1 = np.argwhere(self._grid.nodes == end_node)[0]

        dx = x1 - x0
        dy = y1 - y0

        trace_is_steep = abs(dy) > abs(dx)

        if trace_is_steep:
            x0, y0 = y0, x0
            x1, y1 = y1, x1

        flipped_nodes = x0 > x1

        if flipped_nodes:
            x0, x1 = x1, x0
            y0, y1 = y1, y0

        dx = x1 - x0
        dy = y1 - y0

        error = int(dx / 2.0)

        if y0 < y1:
            y_step = 1
        else:
            y_step = -1

        # Iterate within the bounding box to identify the profile sample nodes.

        samples = []
        y = y0

        for x in range(x0, x1 + 1):
            if trace_is_steep:
                coord = (x, y)
            else:
                coord = (y, x)

            samples.append(self._grid.grid_coords_to_node_id(*coord))

            error -= abs(dy)
            if error < 0:
                y += y_step
                error += dx

        if flipped_nodes:
            samples.reverse()

        return samples

    def _project_point_onto_line(self, p, ep0, ep1):
        """Get the coordinates along a line nearest to a point.

        Parameters
        ----------
        p : tuple of floats
            The x, y coordinates of the point to project onto the line.
        ep0, ep1 : tuple of floats
            The endpoints of the line. Each endpoint is a tuple of x, y
            coordinates.

        Returns
        -------
        tuple
            The x, y coordinates along a line (bounded by `ep1` and `ep2`) that
            is nearest to `p`.
        """
        dx, dy = ep1[0] - ep0[0], ep1[1] - ep0[1]
        determinant = dx * dx + dy * dy
        coeff = (dy * (p[1] - ep0[1]) + dx * (p[0] - ep0[0])) / determinant

        return ep0[0] + coeff * dx, ep0[1] + coeff * dy

    def plot_profiles_in_map_view(
        self, field="topographic__elevation", endpoints_only=True, **kwds
    ):
        """Plot profile locations in map view.

        This method overrides the method in ``_BaseProfiler`` to set the
        default of ``endpoints_only`` to True.

        Parameters
        ----------
        field : field name or nnode array
            Array of the at-node-field to plot as the 2D map values.
            Default value is the at-node field 'topographic__elevation'.
        endpoints_only : boolean
            Boolean where False indicates every node along the profile is
            plotted, or True (default) indicating only segment endpoints are
            plotted.
        **kwds : dictionary
            Keyword arguments to pass to imshow_grid.
        """
        super().plot_profiles_in_map_view(field, endpoints_only=endpoints_only, **kwds)
