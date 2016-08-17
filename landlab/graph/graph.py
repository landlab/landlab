"""Define a graph of nodes-links-patches.

Examples
--------

>>> from landlab.graph import Graph

>>> node_x, node_y = [0, 0, 0, 1, 1, 1, 2, 2, 2], [0, 1, 2, 0, 1, 2, 0, 1, 2]
>>> graph = Graph((node_y, node_x))
>>> graph.x_of_node
array([ 0.,  1.,  2.,  0.,  1.,  2.,  0.,  1.,  2.])
>>> graph.y_of_node
array([ 0.,  0.,  0.,  1.,  1.,  1.,  2.,  2.,  2.])

>>> links = ((0, 1), (1, 2),
...          (0, 3), (1, 4), (2, 5),
...          (3, 4), (4, 5),
...          (3, 6), (4, 7), (5, 8),
...          (6, 7), (7, 8))
>>> graph = Graph((node_y, node_x), links=links)
>>> graph.nodes_at_link # doctest: +NORMALIZE_WHITESPACE
array([[0, 1], [1, 2],
       [0, 3], [1, 4], [2, 5],
       [3, 4], [4, 5],
       [3, 6], [4, 7], [5, 8],
       [6, 7], [7, 8]])
>>> graph.node_at_link_head
array([1, 2, 3, 4, 5, 4, 5, 6, 7, 8, 7, 8])
>>> graph.node_at_link_tail
array([0, 1, 0, 1, 2, 3, 4, 3, 4, 5, 6, 7])

>>> graph.links_at_node # doctest: +NORMALIZE_WHITESPACE
array([[ 0,  2, -1, -1], [ 1,  3,  0, -1], [ 4,  1, -1, -1],
       [ 5,  7,  2, -1], [ 6,  8,  5,  3], [ 9,  6,  4, -1],
       [10,  7, -1, -1], [11, 10,  8, -1], [11,  9, -1, -1]])

>>> graph.link_dirs_at_node # doctest: +NORMALIZE_WHITESPACE
array([[-1, -1,  0,  0], [-1, -1,  1,  0], [-1,  1,  0,  0],
       [-1, -1,  1,  0], [-1, -1,  1,  1], [-1,  1,  1,  0],
       [-1,  1,  0,  0], [-1,  1,  1,  0], [ 1,  1,  0,  0]])

>>> patches = ((5, 3, 0, 2), (6, 4, 1, 3), (10, 8, 5, 7), (11, 9, 6, 8))
>>> graph = Graph((node_y, node_x), links=links, patches=patches)
>>> graph.links_at_patch
array([[ 3,  5,  2,  0],
       [ 4,  6,  3,  1],
       [ 8, 10,  7,  5],
       [ 9, 11,  8,  6]])
>>> graph.nodes_at_patch
array([[4, 3, 0, 1],
       [5, 4, 1, 2],
       [7, 6, 3, 4],
       [8, 7, 4, 5]])
"""
from six.moves import range

import numpy as np

from ..core.utils import as_id_array, argsort_points_by_x_then_y
from ..utils.jaggedarray import flatten_jagged_array
from ..utils.decorators import store_result_in_grid, read_only_array
from .sort import sort_graph, reindex_by_xy, reorder_links_at_patch
from .object.at_node import get_links_at_node
from .object.at_patch import get_nodes_at_patch
from .quantity.of_link import (get_angle_of_link, get_length_of_link,
                               get_midpoint_of_link)
from .quantity.of_patch import get_centroid_of_patch, get_area_of_patch

from .sort.sort import reverse_one_to_many, reorient_link_dirs


def _parse_sorting_opt(sorting):
    SORTING_OPTS = ('xy', 'ccw', 'ne')

    as_dict = None

    if isinstance(sorting, bool):
        as_dict = dict([(opt, format) for opt in SORTING_OPTS])
    elif isinstance(sorting, dict):
        as_dict = dict(sorting.items())
        for opt in SORTING_OPTS:
            sorting.setdefault(opt, True)

    return as_dict


def find_perimeter_nodes(graph):
    """Find nodes on the perimeter of a graph.

    Uses a convex hull to locate the perimeter nodes of a graph.

    Parameters
    ----------
    graph : graph_like
        A Graph of nodes (just requires *xy_of_node*).

    Returns
    -------
    ndarray of int
        Identifiers of the perimeter nodes.
    """
    from scipy.spatial import ConvexHull

    hull = ConvexHull(graph.xy_of_node, qhull_options='Qt')
    return as_id_array(hull.vertices)


class Graph(object):

    """Define the connectivity of a graph of nodes, links, and patches."""

    def __init__(self, nodes, links=None, patches=None, sorting=True):
        """Define a graph of connected nodes.

        Parameters
        ----------
        nodes : tuple of array_like
            Coordinates of nodes as (*y*, *x*).
        links : array_like of tuple
            Tail node and head node for each link in the graph.
        patches : array_like of tuple
            Links that define each patch.
        sort : bool, optional
            Sort elements by their *x* and then *y* coordinate and use
            counter-clockwise element ordering when ordering one set
            of elements around another.
        """
        sorting = _parse_sorting_opt(sorting)
        if sorting is None:
            raise ValueError('bad argument for sorting keyword')
        if len(nodes[0]) != len(nodes[1]):
            raise ValueError('length mismatch in node coordinates')

        self._sorting = sorting

        nodes = [np.asarray(coord, dtype=float) for coord in nodes]

        if patches is not None:
            if len(patches) > 0:
                patches = flatten_jagged_array(patches, dtype=int)
            else:
                patches = None

        self._xy_of_node = np.stack((nodes[1], nodes[0])).T.copy()
        # self._y_of_node, self._x_of_node = nodes[0], nodes[1]
        self._nodes = np.arange(len(nodes[0]), dtype=int)

        self._create_nodes_at_link(links)
        self._create_links_at_patch(patches)

        not sorting['ne'] or reorient_link_dirs(self)
        not sorting['xy'] or reindex_by_xy(self)
        not sorting['ccw'] or reorder_links_at_patch(self)

        self._origin = (0., 0.)

    def _create_nodes_at_link(self, links):
        """Set up node-link data structures."""
        if links is not None:
            self._nodes_at_link = np.asarray(links, dtype=np.int)
            return self._nodes_at_link

    def _create_links_at_patch(self, patches):
        """Set up patch data structures."""
        from .matrix.at_patch import links_at_patch

        if patches is not None:
            self._links_at_patch = links_at_patch(
                patches, nodes_at_link=self.nodes_at_link)
            return self._links_at_patch

    @property
    def ndim(self):
        return 2

    @property
    def xy_of_node(self):
        """Get x and y-coordinates of node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
        >>> graph = Graph((node_y, node_x))
        >>> graph.xy_of_node[:, 0]
        array([ 0.,  1.,  2.,  0.,  1.,  2.])
        >>> graph.xy_of_node[:, 1]
        array([ 0.,  0.,  0.,  1.,  1.,  1.])
        """
        return self._xy_of_node

    @property
    def x_of_node(self):
        """Get x-coordinate of node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
        >>> graph = Graph((node_y, node_x))
        >>> graph.x_of_node
        array([ 0.,  1.,  2.,  0.,  1.,  2.])
        """
        return self._xy_of_node[:, 0]
        # return self._x_of_node

    @property
    def y_of_node(self):
        """Get y-coordinate of node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
        >>> graph = Graph((node_y, node_x))
        >>> graph.y_of_node
        array([ 0.,  0.,  0.,  1.,  1.,  1.])
        """
        return self._xy_of_node[:, 1]
        # return self._y_of_node

    @property
    def nodes(self):
        """Get identifier for each node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
        >>> graph = Graph((node_y, node_x))
        >>> graph.nodes
        array([0, 1, 2, 3, 4, 5])
        """
        return self._nodes

    @property
    @store_result_in_grid()
    def perimeter_nodes(self):
        return find_perimeter_nodes(self)

    @property
    def number_of_nodes(self):
        """Get total number of nodes.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
        >>> graph = Graph((node_y, node_x))
        >>> graph.number_of_nodes
        6
        """
        return self._nodes.size

    @property
    def nodes_at_link(self):
        """Get nodes at either end of links.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.nodes_at_link # doctest: +NORMALIZE_WHITESPACE
        array([[0, 1], [1, 2],
               [0, 3], [1, 4], [2, 5],
               [3, 4], [4, 5],
               [3, 6], [4, 7], [5, 8],
               [6, 7], [7, 8]])
        """
        return self._nodes_at_link

    @property
    def node_at_link_tail(self):
        """Get nodes at link tail.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.node_at_link_tail
        array([0, 1, 0, 1, 2, 3, 4, 3, 4, 5, 6, 7])
        """
        return self._nodes_at_link[:, 0]

    @property
    def node_at_link_head(self):
        """Get nodes at link head.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.node_at_link_head
        array([1, 2, 3, 4, 5, 4, 5, 6, 7, 8, 7, 8])
        """
        return self._nodes_at_link[:, 1]

    @property
    def number_of_links(self):
        """Get nodes at link head.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.number_of_links
        12
        """
        try:
            return len(self._nodes_at_link)
        except AttributeError:
            return 0

    @property
    def links_at_patch(self):
        """Get the links that define a patch.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches)
        >>> graph.links_at_patch
        array([[3, 5, 2, 0],
               [4, 6, 3, 1]])
        """
        return self._links_at_patch

    @property
    def nodes_at_patch(self):
        """Get the nodes that define a patch.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = ([0, 1, 2, 0, 1, 2, 0, 1, 2],
        ...                   [0, 0, 0, 1, 1, 1, 2, 2, 2])
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches)
        >>> graph.nodes_at_patch
        array([[4, 3, 0, 1],
               [5, 4, 1, 2]])
        """
        try:
            return self._nodes_at_patch
        except AttributeError:
            self._nodes_at_patch = get_nodes_at_patch(self)
            return self._nodes_at_patch

    @property
    def patches_at_node(self):
        """Get the patches that touch each node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = ([0, 1, 2, 0, 1, 2],
        ...                   [0, 0, 0, 1, 1, 1])
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5))
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches)
        >>> graph.patches_at_node # doctest: +NORMALIZE_WHITESPACE
        array([[ 0, -1], [ 0,  1], [ 1, -1],
               [ 0, -1], [ 0,  1], [ 1, -1]])
        """
        try:
            return self._patches_at_node
        except AttributeError:
            self._patches_at_node = reverse_one_to_many(self.nodes_at_patch)
            return self._patches_at_node

    @property
    @store_result_in_grid()
    @read_only_array
    def patches_at_link(self):
        """Get the patches on either side of each link.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = ([0, 1, 2, 0, 1, 2],
        ...                   [0, 0, 0, 1, 1, 1])
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5))
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches)
        >>> graph.patches_at_link # doctest: +NORMALIZE_WHITESPACE
        array([[ 0, -1], [ 1, -1],
               [ 0, -1], [ 0,  1], [ 1, -1],
               [ 0, -1], [ 1, -1]])
        """
        return reverse_one_to_many(self._links_at_patch)

    @property
    def number_of_patches(self):
        """Get the number of patches.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches)
        >>> graph.number_of_patches
        2
        """
        try:
            return len(self._links_at_patch)
        except AttributeError:
            return 0

    @property
    def links_at_node(self):
        """Get links touching a node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x = [0, 1, 2, 0, 1, 2, 0, 1, 2]
        >>> node_y = [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.links_at_node # doctest: +NORMALIZE_WHITESPACE
        array([[ 0,  2, -1, -1], [ 1,  3,  0, -1], [ 4,  1, -1, -1],
               [ 5,  7,  2, -1], [ 6,  8,  5,  3], [ 9,  6,  4, -1],
               [10,  7, -1, -1], [11, 10,  8, -1], [11,  9, -1, -1]])
        """
        try:
            return self._links_at_node
        except AttributeError:
            (self._links_at_node,
             self._link_dirs_at_node) = self._create_links_and_dirs_at_node()
            return self._links_at_node

    def _create_links_and_dirs_at_node(self):
        return get_links_at_node(self, sort=self._sorting['ccw'])

    @property
    def link_dirs_at_node(self):
        """Get directions of links touching a node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.link_dirs_at_node # doctest: +NORMALIZE_WHITESPACE
        array([[-1, -1,  0,  0], [-1, -1,  1,  0], [-1,  1,  0,  0],
               [-1, -1,  1,  0], [-1, -1,  1,  1], [-1,  1,  1,  0],
               [-1,  1,  0,  0], [-1,  1,  1,  0], [ 1,  1,  0,  0]])
        """
        try:
            return self._link_dirs_at_node
        except AttributeError:
            (self._links_at_node,
             self._link_dirs_at_node) = self._create_links_and_dirs_at_node()
            return self._link_dirs_at_node

    @property
    @store_result_in_grid()
    @read_only_array
    def angle_of_link(self):
        """Get the angle of each link.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.graph import Graph

        >>> node_x, node_y = ([0, 1, 2, 0, 1, 2],
        ...                   [0, 0, 0, 1, 1, 1])
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.angle_of_link * 180. / np.pi
        array([  0.,   0.,  90.,  90.,  90.,   0.,   0.])
        """
        return get_angle_of_link(self)

    @property
    @store_result_in_grid()
    @read_only_array
    def length_of_link(self):
        """Get the length of links.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.graph import UniformRectilinearGraph

        >>> graph = UniformRectilinearGraph((2, 3), spacing=(1, 2))
        >>> graph.length_of_link
        array([ 2.,  2.,  1.,  1.,  1.,  2.,  2.])
        """
        return get_length_of_link(self)

    @property
    @store_result_in_grid()
    @read_only_array
    def midpoint_of_link(self):
        """Get the middle of links.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.graph import UniformRectilinearGraph

        >>> graph = UniformRectilinearGraph((2, 3), spacing=(1, 2))
        >>> graph.midpoint_of_link # doctest: +NORMALIZE_WHITESPACE
        array([[ 1. ,  0. ], [ 3. ,  0. ],
               [ 0. ,  0.5], [ 2. ,  0.5], [ 4. ,  0.5],
               [ 1. ,  1. ], [ 3. ,  1. ]])
        """
        return get_midpoint_of_link(self)

    @property
    @store_result_in_grid()
    @read_only_array
    def xy_of_link(self):
        return get_midpoint_of_link(self)

    @property
    @store_result_in_grid()
    @read_only_array
    def xy_of_patch(self):
        return get_centroid_of_patch(self)

    @property
    @store_result_in_grid()
    @read_only_array
    def area_of_patch(self):
        return get_area_of_patch(self)
