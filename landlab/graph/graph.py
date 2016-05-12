"""Define a graph of nodes-links-patches.

Examples
--------

>>> from landlab.graph import Graph

>>> node_x, node_y = [0, 0, 0, 1, 1, 1, 2, 2, 2], [0, 1, 2, 0, 1, 2, 0, 1, 2]
>>> graph = Graph((node_y, node_x))
>>> graph.x_of_node
array([ 0.,  0.,  0.,  1.,  1.,  1.,  2.,  2.,  2.])
>>> graph.y_of_node
array([ 0.,  1.,  2.,  0.,  1.,  2.,  0.,  1.,  2.])

>>> links = ((0, 1), (1, 2),
...          (0, 3), (1, 4), (2, 5),
...          (3, 4), (4, 5),
...          (3, 6), (4, 7), (5, 8),
...          (6, 7), (7, 8))
>>> graph = Graph((node_y, node_x), links=links, rot_sort=True)
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
array([[ 2,  0, -1, -1], [ 3,  1,  0, -1], [ 4,  1, -1, -1],
       [ 7,  5,  2, -1], [ 8,  6,  3,  5], [ 9,  4,  6, -1],
       [10,  7, -1, -1], [11,  8, 10, -1], [ 9, 11, -1, -1]])

>>> graph.link_dirs_at_node # doctest: +NORMALIZE_WHITESPACE
array([[-1, -1,  0,  0], [-1, -1,  1,  0], [-1,  1,  0,  0],
       [-1, -1,  1,  0], [-1, -1,  1,  1], [-1,  1,  1,  0],
       [-1,  1,  0,  0], [-1,  1,  1,  0], [ 1,  1,  0,  0]])

>>> patches = ((5, 3, 0, 2), (6, 4, 1, 3), (10, 8, 5, 7), (11, 9, 6, 8))
>>> graph = Graph((node_y, node_x), links=links, patches=patches)
>>> graph.links_at_patch
array([[ 5,  3,  0,  2],
       [ 6,  4,  1,  3],
       [10,  8,  5,  7],
       [11,  9,  6,  8]])
>>> graph.nodes_at_patch
array([[4, 1, 0, 3],
       [5, 2, 1, 4],
       [7, 4, 3, 6],
       [8, 5, 4, 7]])
"""
from six.moves import range

import numpy as np

from ..core.utils import as_id_array, argsort_points_by_x_then_y
from ..utils.jaggedarray import flatten_jagged_array
from .sort import sort_graph, reindex_by_xy, reorder_links_at_patch
from .object.at_node import get_links_at_node
from .object.at_patch import get_nodes_at_patch
from .quantity.of_link import get_angle_of_link

from .sort.sort import reverse_one_to_many, reorient_link_dirs


class Graph(object):

    """Define the connectivity of a graph of nodes, links, and patches."""

    def __init__(self, nodes, links=None, patches=None, xy_sort=True,
                 rot_sort=True, sort_opts=None):
        """Define a graph of connected nodes.

        Parameters
        ----------
        nodes : tuple of array_like
            Coordinates of nodes as (*y*, *x*).
        links : array_like of tuple
            Tail node and head node for each link in the graph.
        patches : array_like of tuple
            Links that define each patch.
        xy_sort : bool, optional
            Sort elements by their *x* and then *y* coordinate.
        rot_sort : bool, optional
            Use counter-clockwise element ordering when ordering one set
            of elements around another.
        """
        self._xy_sort = xy_sort
        self._rot_sort = rot_sort

        nodes = [np.asarray(coord, dtype=float) for coord in nodes]
        # if links is not None:
        #     links = np.asarray(links, dtype=int)

        if len(nodes[0]) != len(nodes[1]):
            raise ValueError('length mismatch in node coordinates')

        # if xy_sort:
        #     nodes, links, patches = sort_graph(nodes, links=links,
        #                                        patches=patches)
        # else:
        if patches is not None:
            if len(patches) > 0:
                patches = flatten_jagged_array(patches, dtype=int)
            else:
                patches = None

        self._y_of_node = nodes[0]
        self._x_of_node = nodes[1]

        self._nodes = np.arange(len(self._x_of_node), dtype=int)

        if links is not None:
            self._setup_links(links)

        if patches is not None:
            self._setup_patches(patches)

        if xy_sort:
            reindex_by_xy(self)

        if rot_sort and patches is not None:
            reorder_links_at_patch(self)

        if links is not None:
            reorient_link_dirs(self)

    def _setup_links(self, links):
        """Set up node-link data structures."""
        self._nodes_at_link = np.asarray(links, dtype=np.int)

    def _setup_patches(self, patches):
        """Set up patch data structures."""
        from .sort.ext.remap_element import reorder_links_at_patch
        from .quantity.ext.of_link import calc_midpoint_of_link
        from .matrix.at_patch import links_at_patch

        # if self._rot_sort:
        #     xy_of_link = np.empty((self.number_of_links, 2), dtype=float)
        #     calc_midpoint_of_link(self.nodes_at_link, self.x_of_node,
        #                           self.y_of_node, xy_of_link)
        #     reorder_links_at_patch(patches[0], patches[1], xy_of_link)

        self._links_at_patch = links_at_patch(patches,
                                              nodes_at_link=self.nodes_at_link)

    def _reorder_links_at_node(self):
        from .cfuncs import _reorder_links_at_node

        outward_angle = self.angle_of_link[self.links_at_node]
        outward_angle[np.where(self.link_dirs_at_node == -1)] -= np.pi
        outward_angle[np.where(self.link_dirs_at_node == 0)] = 2 * np.pi

        sorted_links = as_id_array(np.argsort(outward_angle))

        _reorder_links_at_node(self._links_at_node, sorted_links)
        _reorder_links_at_node(self._link_dirs_at_node, sorted_links)

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
        return self._x_of_node

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
        return self._y_of_node

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
        return len(self._nodes_at_link)

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
        array([[0, 3, 5, 2],
               [1, 4, 6, 3]])
        """
        return self._links_at_patch

    @property
    def nodes_at_patch(self):
        """Get the nodes that define a patch.

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
        >>> graph.nodes_at_patch
        array([[1, 4, 3, 0],
               [2, 5, 4, 1]])
        """
        try:
            return self._nodes_at_patch
        except AttributeError:
            self._nodes_at_patch = get_nodes_at_patch(self)
            return self._nodes_at_patch

    @property
    def patches_at_node(self):
        try:
            return self._patches_at_node
        except AttributeError:
            self._patches_at_node = reverse_one_to_many(self._nodes_at_patch)
            return self._patches_at_node

    @property
    def patches_at_link(self):
        try:
            return self._patches_at_link
        except AttributeError:
            self._patches_at_link = reverse_one_to_many(self._links_at_patch)
            return self._patches_at_link

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
        return len(self._links_at_patch)

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
        >>> graph = Graph((node_y, node_x), links=links, rot_sort=True)
        >>> graph.links_at_node # doctest: +NORMALIZE_WHITESPACE
        array([[ 0,  2, -1, -1], [ 1,  3,  0, -1], [ 4,  1, -1, -1],
               [ 5,  7,  2, -1], [ 6,  8,  5,  3], [ 9,  6,  4, -1],
               [10,  7, -1, -1], [11, 10,  8, -1], [11,  9, -1, -1]])
        """
        try:
            return self._links_at_node
        except AttributeError:
            (self._links_at_node,
             self._link_dirs_at_node) = self.get_links_at_node()
            return self._links_at_node

    def get_links_at_node(self):
        return get_links_at_node(self, sort=self._rot_sort)

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
        >>> graph = Graph((node_y, node_x), links=links, rot_sort=True)
        >>> graph.link_dirs_at_node # doctest: +NORMALIZE_WHITESPACE
        array([[-1, -1,  0,  0], [-1, -1,  1,  0], [-1,  1,  0,  0],
               [-1, -1,  1,  0], [-1, -1,  1,  1], [-1,  1,  1,  0],
               [-1,  1,  0,  0], [-1,  1,  1,  0], [ 1,  1,  0,  0]])
        """
        try:
            return self._link_dirs_at_node
        except AttributeError:
            (self._links_at_node,
             self._link_dirs_at_node) = get_links_at_node(self,
                                                          sort=self._rot_sort)
            return self._link_dirs_at_node

    @property
    def angle_of_link(self):
        try:
            return self._angle_of_link
        except AttributeError:
            self._angle_of_link = get_angle_of_link(self)
            return self._angle_of_link
