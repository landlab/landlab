from six.moves import range

import numpy as np


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
    node_count = np.bincount(nodes_at_link.flat)
    number_of_nodes = number_of_nodes or len(node_count)

    max_node_count = np.max(node_count)

    link_dirs_at_node = np.full((number_of_nodes, max_node_count), 0,
                                dtype=int)
    links_at_node = np.full((number_of_nodes, max_node_count), -1, dtype=int)

    for node in range(number_of_nodes):
        (link, is_outward) = np.where(nodes_at_link == node)
        is_outward[is_outward == 1] = -1
        is_outward[is_outward == 0] = 1

        link_dirs_at_node[node, :len(is_outward)] = is_outward
        links_at_node[node, :len(link)] = link

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
    max_links_per_patch = np.max([len(links) for links in patches])
    links_at_patch = np.full((len(patches), max_links_per_patch), -1,
                             dtype=int)

    for patch, links in enumerate(patches):
        links_at_patch[patch, :len(links)] = links

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


class Graph(object):
    def __init__(self, nodes, links=None, patches=None):
        """Define a graph of connected nodes.

        Examples
        --------
        >>> from landlab.graph import Graph

        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
        >>> graph = Graph((node_y, node_x))
        >>> graph.x_of_node
        array([0, 1, 2, 0, 1, 2])
        >>> graph.y_of_node
        array([0, 0, 0, 1, 1, 1])

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
               [ 2,  5,  7, -1], [ 3,  5,  6,  8], [ 4,  6,  9, -1]])

        >>> graph.link_dirs_at_node # doctest: +NORMALIZE_WHITESPACE
        array([[ 1,  1,  0,  0], [-1,  1,  1,  0], [-1,  1,  0,  0],
               [-1,  1,  1,  0], [-1, -1,  1,  1], [-1, -1,  1,  0]])

        >>> patches = ((0, 3, 5, 2), (1, 4, 6, 3))
        >>> graph = Graph((node_y, node_x), links=links, patches=patches)
        >>> graph.links_at_patch
        array([[0, 3, 5, 2],
               [1, 4, 6, 3]])
        >>> graph.nodes_at_patch
        array([[0, 1, 3, 4],
               [1, 2, 4, 5]])
        """
        if len(nodes[0]) != len(nodes[1]):
            raise ValueError('length mismatch in node coordinates')

        self._y_of_node = np.asarray(nodes[0])
        self._x_of_node = np.asarray(nodes[1])

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
        self._nodes_at_patch = _setup_nodes_at_patch(self._links_at_patch,
                                                     self._nodes_at_link)
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
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
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
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
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
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
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
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
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
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
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
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
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
        return self._nodes_at_patch

    @property
    def number_of_patches(self):
        """Get the number of patches.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
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
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.links_at_node # doctest: +NORMALIZE_WHITESPACE
        array([[ 0,  2, -1, -1], [ 0,  1,  3, -1], [ 1,  4, -1, -1],
               [ 2,  5,  7, -1], [ 3,  5,  6,  8], [ 4,  6,  9, -1]])
        """
        return self._links_at_node

    @property
    def link_dirs_at_node(self):
        """Get directions of links touching a node.

        Examples
        --------
        >>> from landlab.graph import Graph
        >>> node_x, node_y = [0, 1, 2, 0, 1, 2], [0, 0, 0, 1, 1, 1]
        >>> links = ((0, 1), (1, 2),
        ...          (0, 3), (1, 4), (2, 5),
        ...          (3, 4), (4, 5),
        ...          (3, 6), (4, 7), (5, 8),
        ...          (6, 7), (7, 8))
        >>> graph = Graph((node_y, node_x), links=links)
        >>> graph.link_dirs_at_node # doctest: +NORMALIZE_WHITESPACE
        array([[ 1,  1,  0,  0], [-1,  1,  1,  0], [-1,  1,  0,  0],
               [-1,  1,  1,  0], [-1, -1,  1,  1], [-1, -1,  1,  0]])
        """
        return self._link_dirs_at_node
