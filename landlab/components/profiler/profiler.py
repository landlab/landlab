# coding: utf8
# ! /usr/env/python
"""profiler.py component to create profiles with user-defined endpoints."""
from collections import OrderedDict

import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

from landlab.components.profiler.base_profiler import _BaseProfiler


class Profiler(_BaseProfiler):
    """Extract and plot profiles defined by points within a grid.

    The profile is constructed from the first to final point in `endpoints`.
    Endpoints are located at grid nodes. Two successive endpoints bound a
    profile segment. A profile with one segment is a straight line. The
    segments of a profile with multiple segments meet at endpoints. Every grid
    node along the profile and between segment endpoints is sampled, including
    the segment endpoints.

    The structure of the profile in a model grid is diagramed below. The
    diagram is explained in this list:
    -   The grid contains nine columns and nine rows.
    -   The profile is constructed from three endpoints that bound two
        segments.
    -   In the diagram below, 'o' indicates a segment endpoint, '.' and '*' are
        sample nodes of the first and second segment, respectively. 'X' are
        nodes not included in the profile.
    -   Segments have seven sample points each (nodes at endpoints are also
        sampled).
    -   The first segment begins in the lower-left and continues horizontally
        and almost reaches the right boundary. The second segment is joined to
        the first in the lower-right of the grid and it continues diagonally to
        the upper-left.
    -   Segment and sample ordering is dictated by the ordering of endpoints.
        The endpoints used to construct this profile must be of been ordered
        lower-left, lower-right, and then upper-left because the horizontal
        segment was designated as being first.

        X X X X X X X X X
        X o X X X X X X X
        X X * X X X X X X
        X X X * X X X X X
        X X X X * X X X X
        X X X X X * X X X
        X X X X X X * X X
        X o . . . . . o X
        X X X X X X X X X

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import Profiler
    >>> mg = RasterModelGrid((10, 10), 10)
    >>> mg.at_node['topographic__elevation'] = mg.node_x * mg.node_y

    Create a profile with three endpoints.
    >>> endpoints = [10, 16, 64]
    >>> profiler = Profiler(mg, endpoints)
    profiler.run_one_step()

    The keys of the network structure are the segment ids.
    >>> profiler.network_structure.keys()
    odict_keys([0, 1])

    The nodes of the first segment.
    >>> profiler.network_structure[0]['ids']
    [10, 11, 12, 13, 14, 15, 16]

    The first node of the second segment is the same as the final node of the
    first segment.
    >>> profiler.network_structure[1]['ids']
    [16, 24, 32, 40, 48, 56, 64]

    Alternative to nodes, profiles can be instantiated with coordinates.
    >>> profiler = Profiler(mg, [(10, 10), (70, 10), (10, 70)])

    Endpoints can also be set with a combination of coordinates and nodes.
    >>> profiler = Profiler(mg, [(10, 10), 16, (10, 70)])
    """
    _name = "Profiler"

    def __init__(self, grid, endpoints, cmap="viridis"):
        """Instantiate Profiler.

        Parameters
        ----------
        grid : RasterModelGrid
            A landlab RasterModelGrid.
        endpoints : list of integers or list of tuples
            The endpoints that bound segments of the profile. Endpoints can be
            node ids and tuples of coordinates (x, y). Both node ids and
            coordinate tuples can be used in the same endpoint list. The
            profile begins with the first element of `endpoints` and continues
            in the order of this list.
        cmap : str
            A valid matplotlib cmap string. Default is "viridis".
        """
        super(_BaseProfiler, self).__init__(grid)

        if not isinstance(endpoints, list) and len(endpoints) > 1:
            msg = ('`endpoints` must be a list of at least 2 nodes or a list '
                   'of at least two tuples where each tuple contains the x, y '
                   'coordinates of endpoints.')
            raise ValueError(msg)

        self._grid = grid
        self._endpoints = endpoints
        self._cmap = plt.get_cmap(cmap)

    @property
    def network_structure(self):
        """OrderedDict defining the profile.

        The node IDs and distances upstream of the channel network are stored
        in `network_structure`. It is a dictionary with keys of the segment
        ID.

        First, 'ids' contains a list of the segment node IDs ordered from
        the start to the end of the profile. It includes the endpoints. Second,
        'distances' contains a list of along-profile distances that mirrors the
        list in 'ids'. Finally, 'color' is an RGBA tuple indicating the color
        for the segment.
        """
        return self._net_struct

    def _create_profile_structure(self):
        """Create the network data structure of the profile.

        The profile is processed by segment. Segments are bound by successive
        endpoints. The cumulative distance along the profile is accumulated by
        iteratively adding segment lengths.
        """
        self._net_struct = OrderedDict()
        grid = self._grid
        endpts = self._endpoints
        cum_dist = 0

        for i_endpt in range(len(endpts) - 1):
            # Get the endpoints and samples of the segment.

            start_node, start_xy = self._get_node_and_coords(endpts[i_endpt])
            end_node, end_xy = self._get_node_and_coords(endpts[i_endpt + 1])

            sample_nodes = self._get_sample_nodes(start_node, end_node)

            # Calculate the along-profile distance of samples along the
            # segment.

            n_samples = len(sample_nodes)
            sample_distances = np.empty(n_samples, dtype=float)

            for i_sample, node in enumerate(sample_nodes):
                sample_xy = grid.xy_of_node[node]

                pt = self._project_point_onto_line(sample_xy, start_xy, end_xy)
                d = grid.calc_distances_of_nodes_to_point(
                        pt, node_subset=start_node)
                sample_distances[i_sample] = d

            # Store the segment data.

            self._net_struct[i_endpt] = {
                    'ids': sample_nodes,
                    'distances': sample_distances + cum_dist}

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
            segment_count = len(self._net_struct)
            norm = mpl.colors.Normalize(vmin=0, vmax=segment_count)
            mappable = cm.ScalarMappable(norm=norm, cmap=self._cmap)
            color_mapping = {
                segment_id: mappable.to_rgba(idx)
                for idx, segment_id in enumerate(self._net_struct)
            }

        for segment_id in self._net_struct:
            self._net_struct[segment_id]["color"] = color_mapping[
                segment_id
            ]

    def _create_flat_structures(self):
        """Create expected flattened structures for ids, distances, and colors.
        """
        self._net_ids = []
        self._distance_along_profile = []
        self._colors = []

        for segment_id in self._net_struct:
            self._net_ids.append(self._net_struct[segment_id]["ids"])
            self._distance_along_profile.append(
                    self._net_struct[segment_id]["distances"])
            self._colors.append(self._net_struct[segment_id]["color"])

    def _get_node_and_coords(self, point):
        """Get the node and coordinates for a point.

        The point might be a node or tuple. This method handles this option.
        """
        if isinstance(point, (float, int)):
            return point, self._grid.xy_of_node[point]
        elif isinstance(point, tuple):
            return self._grid.find_nearest_node(point), point
        else:
            raise TypeError('each element of `endpoints` must be a number '
                            'representing a node id or a tuple of node x, y '
                            'coordinates')

    def _get_sample_nodes(self, start_node, end_node):
        """Get the profile sample nodes using Bresenham's line algorithm.

        Parameters
        ----------
        start_node, end_node : integer
            The node id of the profile endpoints.

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
            # Transpose the profile trace.
            x0, y0 = y0, x0
            x1, y1 = y1, x1

        swapped_nodes = x0 > x1

        if swapped_nodes:
            # Swap the start and end nodes.
            x0, x1 = x1, x0
            y0, y1 = y1, y0

        dx = x1 - x0
        dy = y1 - y0

        error = int(dx / 2.)
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

        if swapped_nodes:
            # Reverse the list if the endpoint nodes were swapped.
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
