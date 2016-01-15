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
array([[ 0,  2, -1, -1], [ 0,  1,  3, -1], [ 1,  4, -1, -1],
       [ 2,  5,  7, -1], [ 3,  5,  6,  8], [ 4,  6,  9, -1],
       [ 7, 10, -1, -1], [ 8, 10, 11, -1], [ 9, 11, -1, -1]])

>>> graph.link_dirs_at_node # doctest: +NORMALIZE_WHITESPACE
array([[-1, -1,  0,  0], [ 1, -1, -1,  0], [ 1, -1,  0,  0],
       [ 1, -1, -1,  0], [ 1,  1, -1, -1], [ 1,  1, -1,  0],
       [ 1, -1,  0,  0], [ 1,  1, -1,  0], [ 1,  1,  0,  0]])

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


class Graph(object):

    """Define the connectivity of a graph of nodes, links, and patches."""

    def __init__(self, nodes, links=None, patches=None, sort=False):
        """Define a graph of connected nodes.

        Parameters
        ----------
        nodes : tuple of array_like
            Coordinates of nodes as (*y*, *x*).
        links : array_like of tuple
            Tail node and head node for each link in the graph.
        patches : array_like of tuple
            Links that define each patch.
        """
        nodes = [np.asarray(coord) for coord in nodes]
        if links is not None:
            links = np.asarray(links, dtype=int)

        if len(nodes[0]) != len(nodes[1]):
            raise ValueError('length mismatch in node coordinates')

        if sort:
            nodes, links, patches = sort_graph(nodes, links, patches)
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
        self._links_at_node, self._link_dirs_at_node = _setup_links_at_node(
            self._nodes_at_link, number_of_nodes=self.number_of_nodes)

    def _setup_patches(self, patches):
        """Set up patch data structures."""
        self._links_at_patch = _setup_links_at_patch(patches)

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
            self._nodes_at_patch = _setup_nodes_at_patch(self._links_at_patch,
                                                         self._nodes_at_link)
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
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1, 2, 2, 2]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.links_at_node # doctest: +NORMALIZE_WHITESPACE
        array([[ 0,  2, -1, -1], [ 0,  1,  3, -1], [ 1,  4, -1, -1],
               [ 2,  5,  7, -1], [ 3,  5,  6,  8], [ 4,  6,  9, -1],
               [ 7, 10, -1, -1], [ 8, 10, 11, -1], [ 9, 11, -1, -1]])
        """
        return self._links_at_node

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
        array([[-1, -1,  0,  0], [ 1, -1, -1,  0], [ 1, -1,  0,  0],
               [ 1, -1, -1,  0], [ 1,  1, -1, -1], [ 1,  1, -1,  0],
               [ 1, -1,  0,  0], [ 1,  1, -1,  0], [ 1,  1,  0,  0]])
        """
        return self._link_dirs_at_node

    def _setup_angle_of_link(self):
        y = self.y_of_node[self.nodes_at_link]
        x = self.x_of_node[self.nodes_at_link]
        return np.arctan2(np.diff(y), np.diff(x)).reshape((-1, ))

    @property
    def angle_of_link(self):
        try:
            return self._angle_of_link
        except AttributeError:
            self._angle_of_link = self._setup_angle_of_link()
            return self._angle_of_link


def _find_links_at_node(node, nodes_at_link):
    """Find links and link directions at a node.

    Examples
    --------
    >>> nodes_at_link = ((0, 1), (1, 2),
    ...                  (0, 3), (1, 4), (2, 5),
    ...                  (3, 4), (4, 5),
    ...                  (3, 6), (4, 7), (5, 8),
    ...                  (6, 7), (7, 8))

    The first node has only two links, both of which are directed outward.

    >>> (links, dirs) = _find_links_at_node(0, nodes_at_link)
    >>> links
    array([0, 2])
    >>> dirs
    array([-1, -1])

    The fourth node has two links entering and two link leaving.

    >>> (links, dirs) = _find_links_at_node(4, nodes_at_link)
    >>> links
    array([3, 5, 6, 8])
    >>> dirs
    array([ 1,  1, -1, -1])
    """
    from .cfuncs import _find_links_at_node

    nodes_at_link = np.asarray(nodes_at_link, dtype=int)
    nodes_at_link.shape = (-1, 2)

    links_at_node = np.full(4, -1, dtype=int)
    link_dirs_at_node = np.full(4, 0, dtype=int)

    n_links = _find_links_at_node(node, nodes_at_link, links_at_node,
                                  link_dirs_at_node)

    return links_at_node[:n_links], link_dirs_at_node[:n_links]


def _setup_links_at_node(nodes_at_link, number_of_nodes=None):
    """Set up data structures for node-to-link connectivity.

    Parameters
    ----------
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

    node_count = np.bincount(nodes_at_link.flat)
    number_of_nodes = number_of_nodes or len(node_count)

    max_node_count = np.max(node_count)

    link_dirs_at_node = np.full((number_of_nodes, max_node_count), 0,
                                dtype=int)
    links_at_node = np.full((number_of_nodes, max_node_count), -1, dtype=int)

    _setup_links_at_node(nodes_at_link, links_at_node, link_dirs_at_node)

    return links_at_node, link_dirs_at_node


def _setup_links_at_patch(patches):
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


def _setup_nodes_at_patch(links_at_patch, nodes_at_link):
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
    nodes_at_patch = np.full(links_at_patch.shape, -1, dtype=int)

    for patch, links in enumerate(links_at_patch):
        unique_nodes = np.unique(nodes_at_link[links])
        nodes_at_patch[patch, :len(unique_nodes)] = unique_nodes

    return nodes_at_patch
