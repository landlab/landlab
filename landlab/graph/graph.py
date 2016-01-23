"""Define a graph of nodes-links-patches.

Examples
--------

>>> from landlab.graph import Graph

>>> node_x, node_y = [0, 0, 0, 1, 1, 1, 2, 2, 2], [0, 1, 2, 0, 1, 2, 0, 1, 2]
>>> graph = Graph((node_y, node_x))
>>> graph.x_of_node
array([0, 0, 0, 1, 1, 1, 2, 2, 2])
>>> graph.y_of_node
array([0, 1, 2, 0, 1, 2, 0, 1, 2])

>>> links = ((0, 1), (1, 2),
...          (0, 3), (1, 4), (2, 5),
...          (3, 4), (4, 5),
...          (3, 6), (4, 7), (5, 8),
...          (6, 7), (7, 8))
>>> graph = Graph((node_y, node_x), links=links, ccw=True)
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
array([[0, 1, 3, 4],
       [1, 2, 4, 5],
       [3, 4, 6, 7],
       [4, 5, 7, 8]])
"""
from six.moves import range

import numpy as np

from ..core.utils import as_id_array, argsort_points_by_x_then_y
from ..utils.jaggedarray import flatten_jagged_array
from .sort import sort_graph


class Graph(object):

    """Define the connectivity of a graph of nodes, links, and patches."""

    def __init__(self, nodes, links=None, patches=None, sort=False, ccw=False):
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
            Sort elements by their *x* and then *y* coordinate.
        ccw : bool, optional
            Use counter-clockwise element ordering when ordering one set
            of elements around another.
        """
        self._sort = sort
        self._ccw = ccw

        nodes = [np.asarray(coord) for coord in nodes]
        if links is not None:
            links = np.asarray(links, dtype=int)

        if len(nodes[0]) != len(nodes[1]):
            raise ValueError('length mismatch in node coordinates')

        if sort:
            nodes, links, patches = sort_graph(nodes, links=links,
                                               patches=patches)
        else:
            if patches is not None:
                patches = flatten_jagged_array(patches, dtype=int)

        self._y_of_node = nodes[0]
        self._x_of_node = nodes[1]

        self._nodes = np.arange(len(self._x_of_node), dtype=int)

        if links is not None:
            self._setup_links(links)

        if patches is not None:
            self._setup_patches(patches)

    def _setup_links(self, links):
        """Set up node-link data structures."""
        self._nodes_at_link = np.asarray(links, dtype=np.int)

    def _setup_patches(self, patches):
        """Set up patch data structures."""
        from .cfuncs import reorder_links_at_patch

        x = (self.x_of_node[self.node_at_link_head] +
             self.x_of_node[self.node_at_link_tail]) * .5
        y = (self.y_of_node[self.node_at_link_head] +
             self.y_of_node[self.node_at_link_tail]) * .5
        xy_of_link = np.vstack((x, y)).T

        if self._ccw:
            reorder_links_at_patch(patches[0], patches[1], xy_of_link)
        self._links_at_patch = get_links_at_patch(patches)

    @property
    def x_of_node(self):
        """Get x-coordinate of node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
        >>> graph = Graph((node_y, node_x))
        >>> graph.x_of_node
        array([0, 1, 2, 0, 1, 2])
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
        array([0, 0, 0, 1, 1, 1])
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
        array([[0, 1, 3, 4],
               [1, 2, 4, 5]])
        """
        try:
            return self._nodes_at_patch
        except AttributeError:
            self._nodes_at_patch = get_nodes_at_patch(self)
            return self._nodes_at_patch

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
        >>> graph = Graph((node_y, node_x), links=links, ccw=True)
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
        return get_links_at_node(self, sort=self._ccw)

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
             self._link_dirs_at_node) = get_links_at_node(self, sort=True)
            return self._link_dirs_at_node

    @property
    def angle_of_link(self):
        try:
            return self._angle_of_link
        except AttributeError:
            self._angle_of_link = get_angle_of_link(self)
            return self._angle_of_link


def get_angle_of_link(graph):
    """Get angles of links in a graph.

    Parameters
    ----------
    graph : `Graph`
        A `Graph`.

    Returns
    -------
    ndarray of float
        Angle of each link as measured in radians from link tail to head.
    """
    y = graph.y_of_node[graph.nodes_at_link[:, ::-1]]
    x = graph.x_of_node[graph.nodes_at_link[:, ::-1]]
    angles = np.arctan2(np.diff(y), np.diff(x)).reshape((-1, )) + np.pi
    angles[angles == 2. * np.pi] = 0.
    return angles


def get_links_at_node(graph, sort=False):
    """Set up data structures for node-to-link connectivity.

    Parameters
    ----------
    graph : Graph
        A `Graph`.

    nodes_at_link: ndarray of int
        Nodes at either end of a link (tail node, then head node).
    number_of_nodes: int, optional
        The total number of nodes. If not given, use the largest node in
        *nodes_at_link*.

    Returns
    -------
    tuple of ndarray
        Tuple of *links_at_node* and *link_dirs_at_node*.
    """
    from .cfuncs import _setup_links_at_node

    node_count = np.bincount(graph.nodes_at_link.flat)
    number_of_nodes = graph.number_of_nodes

    max_node_count = np.max(node_count)

    link_dirs_at_node = np.full((number_of_nodes, max_node_count), 0,
                                dtype=int)
    links_at_node = np.full((number_of_nodes, max_node_count), -1, dtype=int)

    _setup_links_at_node(graph.nodes_at_link, links_at_node, link_dirs_at_node)

    if sort:
        sort_links_at_node_by_angle(links_at_node, link_dirs_at_node,
                                    graph.angle_of_link, inplace=True)

    return links_at_node, link_dirs_at_node


def sort_links_at_node_by_angle(links_at_node, link_dirs_at_node,
                                angle_of_link, inplace=True):
    """Sort links as spokes about a hub.

    Parameters
    ----------
    links_at_node : ndarray of int, shape `(n_nodes, max_links_per_node)`
        Links entering or leaving each node.
    link_dirs_at_node : ndarray of int, shape `(n_nodes, max_links_per_node)`
        Direction of links entering or leaving each node.
    angle_of_link : ndarray of float, shape `(n_links, )`
        Angle (in radians) of each link as measured from its head to tail.

    Returns
    -------
    tuple of (links_at_node, link_dirs_at_node)
        The sorted arrays. If `inplace` is `True`, these are the input
        arrays.
    """
    from .cfuncs import _reorder_links_at_node

    outward_angle = angle_of_link[links_at_node]

    links_entering = np.where(link_dirs_at_node == 1)
    outward_angle[links_entering] += np.pi
    outward_angle[outward_angle >= 2 * np.pi] -= 2 * np.pi

    outward_angle[np.where(link_dirs_at_node == 0)] = 4 * np.pi

    sorted_links = np.argsort(outward_angle)

    if not inplace:
        links_at_node = links_at_node.copy()
        link_dirs_at_node = link_dirs_at_node.copy()

    _reorder_links_at_node(links_at_node, sorted_links)
    _reorder_links_at_node(link_dirs_at_node, sorted_links)

    return links_at_node, link_dirs_at_node


def get_links_at_patch(patches):
    """Set up data structure that describes link-patch connectivity.

    Parameters
    ----------
    patches: iterable of iterables
        List of links for each patch.

    Returns
    -------
    ndarray
        Links for each patch.
    """
    from .cfuncs import _setup_links_at_patch

    max_links_per_patch = np.diff(patches[1]).max()
    n_patches = len(patches[1]) - 1
    links_at_patch = np.full((n_patches, max_links_per_patch), -1, dtype=int)

    _setup_links_at_patch(patches[0], patches[1], links_at_patch)

    return links_at_patch


def get_nodes_at_patch(graph):
    """Set up data structure that describes node-patch connectivity.

    Parameters
    ----------
    links_at_patch: ndarray
        Links that define each patch.
    nodes_at_link: ndarray
        Nodes that define each link.

    Returns
    -------
    ndarray
        Nodes that define each patch.
    """
    nodes_at_patch = np.full(graph.links_at_patch.shape, -1, dtype=int)

    for patch, links in enumerate(graph.links_at_patch):
        unique_nodes = np.unique(graph.nodes_at_link[links])
        nodes_at_patch[patch, :len(unique_nodes)] = unique_nodes

    return nodes_at_patch
