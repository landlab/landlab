import numpy as np
from six.moves import range

from ...core.utils import as_id_array
from ...utils.jaggedarray import JaggedArray
from .status import CORE_NODE, CLOSED_BOUNDARY


def _split_link_ends(link_ends):
    """
    Examples
    --------
    >>> from landlab.grid.unstructured.links import _split_link_ends
    >>> _split_link_ends(((0, 1, 2), (3, 4, 5)))
    (array([0, 1, 2]), array([3, 4, 5]))
    >>> _split_link_ends([(0, 3), (1, 4), (2, 5)])
    (array([0, 1, 2]), array([3, 4, 5]))
    >>> _split_link_ends((0, 3))
    (array([0]), array([3]))
    """
    links = np.array(list(link_ends), ndmin=2, dtype=np.int)
    if len(links) != 2:
        links = links.transpose()

    if links.size == 0:
        return (np.array([], dtype=np.int), np.array([], dtype=np.int))
    else:
        return links[0], links[1]


def link_is_active(status_at_link_ends):
    """link_is_active((status0, status1))
    Check if a link is active.

    Links are *inactive* if they connect two boundary nodes or touch a
    closed boundary. Otherwise, the link is *active*.

    Parameters
    ----------
    status0, status1 : sequence of array-like
        Status at link start and end

    Returns
    -------
    ndarray, boolean :
        Boolean array that indicates if a link is active.
    """
    (status_at_link_start,
     status_at_link_end) = _split_link_ends(status_at_link_ends)

    return (((status_at_link_start == CORE_NODE) &
             ~ (status_at_link_end == CLOSED_BOUNDARY)) |
            ((status_at_link_end == CORE_NODE) &
             ~ (status_at_link_start == CLOSED_BOUNDARY)))


def find_active_links(node_status, node_at_link_ends):
    """find_active_links(node_status, (node0, node1))
    IDs of active links.

    Parameters
    ----------
    node_status : ndarray
        Status of nodes.
    node0, node1 : sequence of array-like
        Node ID at link start and end.

    Returns
    -------
    ndarray :
        Links IDs of active links.

    Examples
    --------
    >>> from landlab.grid.unstructured.links import find_active_links
    >>> links = [(0, 2), (1, 3), (0, 1), (1, 2), (0, 3)]
    >>> status = np.array([0, 0, 0, 0])
    >>> find_active_links(status, links)
    array([0, 1, 2, 3, 4])
    """
    node_at_link_start, node_at_link_end = _split_link_ends(node_at_link_ends)

    if len(node_at_link_end) != len(node_at_link_start):
        raise ValueError('Link arrays must be the same length')

    status_at_link_ends = (node_status[node_at_link_start],
                           node_status[node_at_link_end])

    (active_link_ids, ) = np.where(link_is_active(status_at_link_ends))

    return as_id_array(active_link_ids)


def in_link_count_per_node(node_at_link_ends, number_of_nodes=None):
    """in_link_count_per_node((node0, node1), number_of_nodes=None)
    Number of links entering nodes.

    Parameters
    ----------
    node0, node1 : sequence of array-like
        Node ID at link start and end.
    number_of_nodes : int, optional
        Number of nodes in the grid

    Returns
    -------
    ndarray :
        Number of links entering nodes.

    Examples
    --------
    >>> from landlab.grid.unstructured.links import in_link_count_per_node
    >>> link_ends = [(0, 3), (1, 4), (2, 5), (3, 6), (4, 7), (5, 8)]
    >>> in_link_count_per_node(zip(*link_ends))
    array([0, 0, 0, 1, 1, 1, 1, 1, 1])
    """
    node_at_link_start, node_at_link_end = _split_link_ends(node_at_link_ends)

    # if len(node_at_link_end) != len(node_at_link_start):
    #    raise ValueError('Link arrays must be the same length')
    return as_id_array(np.bincount(node_at_link_end, minlength=number_of_nodes))


def out_link_count_per_node(node_at_link_ends, number_of_nodes=None):
    """out_link_count_per_node((node0, node1), number_of_nodes=None)
    Number of links leaving nodes.

    Parameters
    ----------
    node0, node1 : sequence of array-like
        Node ID at link start and end.
    number_of_nodes : int, optional
        Number of nodes in the grid

    Returns
    -------
    ndarray :
        Number of links leaving nodes.

    Examples
    --------
    >>> from landlab.grid.unstructured.links import out_link_count_per_node
    >>> out_link_count_per_node(([0, 1, 2, 3, 4, 5], [3, 4, 5, 6, 7, 8]))
    array([1, 1, 1, 1, 1, 1])
    >>> out_link_count_per_node(([0, 1, 2, 3, 4, 5], [3, 4, 5, 6, 7, 8]),
    ...     number_of_nodes=9)
    array([1, 1, 1, 1, 1, 1, 0, 0, 0])
    """
    node_at_link_start, node_at_link_end = _split_link_ends(node_at_link_ends)
    if len(node_at_link_end) != len(node_at_link_start):
        raise ValueError('Link arrays must be the same length')
    return as_id_array(np.bincount(node_at_link_start,
                                   minlength=number_of_nodes))


def link_count_per_node(node_at_link_ends, number_of_nodes=None):
    """link_count_per_node((node0, node1), number_of_nodes=None)
    Number of links per node.

    Parameters
    ----------
    node0, node1 : sequence of array-like
        Node ID at link start and end.
    number_of_nodes : int, optional
        Number of nodes in the grid

    Returns
    -------
    ndarray :
        Number of links per nodes.

    Examples
    --------
    >>> from landlab.grid.unstructured.links import link_count_per_node
    >>> link_count_per_node(([0, 1, 2, 3, 4, 5], [3, 4, 5, 6, 7, 8]))
    array([1, 1, 1, 2, 2, 2, 1, 1, 1])
    """
    in_count = in_link_count_per_node(node_at_link_ends)
    out_count = out_link_count_per_node(node_at_link_ends)

    node_count = number_of_nodes or max(len(in_count), len(out_count))

    if len(in_count) < node_count:
        in_count = np.pad(in_count, (0, node_count - len(in_count)),
                          mode='constant')
    if len(out_count) < node_count:
        out_count = np.pad(out_count, (0, node_count - len(out_count)),
                           mode='constant')

    return in_count + out_count


def _sort_links_by_node(node_at_link_ends, link_ids=None, sortby=0):
    sorted_links = np.argsort(node_at_link_ends[sortby])

    if link_ids is not None:
        return np.array(link_ids, dtype=np.int)[sorted_links]
    else:
        return as_id_array(sorted_links)


def in_link_ids_at_node(node_at_link_ends, link_ids=None, number_of_nodes=None):
    """in_link_ids_at_node((node0, node1), number_of_nodes=None)
    Links entering nodes.

    Parameters
    ----------
    node0, node1 : sequence of array-like
        Node ID at link start and end.
    number_of_nodes : int, optional
        Number of nodes in the grid

    Returns
    -------
    tuple :
        Tuple of link id array and offset into link id array.

    Examples
    --------
    >>> from landlab.grid.unstructured.links import in_link_ids_at_node
    >>> (links, count) = in_link_ids_at_node(([0, 1, 2, 3, 4, 5],
    ...                                       [3, 4, 5, 6, 7, 8]))
    >>> links
    array([0, 1, 2, 3, 4, 5])
    >>> count
    array([0, 0, 0, 1, 1, 1, 1, 1, 1])


    >>> (links, count) = in_link_ids_at_node(([0, 1, 2, 3, 4, 5],
    ...                                       [3, 4, 5, 6, 7, 8]),
    ...                                      link_ids=range(1, 7))
    >>> links
    array([1, 2, 3, 4, 5, 6])
    >>> count
    array([0, 0, 0, 1, 1, 1, 1, 1, 1])
    """
    node_at_link_ends = _split_link_ends(node_at_link_ends)

    link_ids = _sort_links_by_node(node_at_link_ends, link_ids=link_ids,
                                   sortby=1)
    links_per_node = in_link_count_per_node(node_at_link_ends,
                                            number_of_nodes=number_of_nodes)
    return link_ids, links_per_node


def out_link_ids_at_node(node_at_link_ends, link_ids=None, number_of_nodes=None):
    """out_link_ids_at_node((node0, node1), number_of_nodes=None)
    Links leaving nodes.

    Parameters
    ----------
    node0, node1 : sequence of array-like
        Node ID at link start and end.
    number_of_nodes : int, optional
        Number of nodes in the grid

    Returns
    -------
    tuple :
        Tuple of link id array and offset into link id array.

    Examples
    --------
    >>> from landlab.grid.unstructured.links import out_link_ids_at_node
    >>> (links, count) = out_link_ids_at_node(
    ...     ([0, 1, 2, 3, 4, 5], [3, 4, 5, 6, 7, 8]), link_ids=range(-1, 5),
    ...     number_of_nodes=9)
    >>> links
    array([-1,  0,  1,  2,  3,  4])
    >>> count
    array([1, 1, 1, 1, 1, 1, 0, 0, 0])


    >>> (links, count) = out_link_ids_at_node(
    ...     ([0, 1, 2, 3, 4, 5], [3, 4, 5, 6, 7, 8]), number_of_nodes=9)
    >>> links
    array([0, 1, 2, 3, 4, 5])
    >>> count
    array([1, 1, 1, 1, 1, 1, 0, 0, 0])
    """
    node_at_link_ends = _split_link_ends(node_at_link_ends)

    link_ids = _sort_links_by_node(node_at_link_ends, link_ids=link_ids,
                                   sortby=0)
    links_per_node = out_link_count_per_node(node_at_link_ends,
                                             number_of_nodes=number_of_nodes)
    return link_ids, links_per_node


def link_ids_at_node(node_at_link_ends, number_of_nodes=None):
    """link_ids_at_node((node0, node1), number_of_nodes=None)
    Links entering and leaving nodes.

    Parameters
    ----------
    node0, node1 : sequence of array-like
        Node ID at link start and end.
    number_of_nodes : int, optional
        Number of nodes in the grid

    Returns
    -------
    tuple :
        Tuple of link id array and offset into link id array.

    Examples
    --------
    >>> from landlab.grid.unstructured.links import link_ids_at_node
    >>> (links, count) = link_ids_at_node(
    ...     ([0, 1, 2, 3, 4, 5], [3, 4, 5, 6, 7, 8]), number_of_nodes=9)
    >>> links
    array([0, 1, 2, 0, 3, 1, 4, 2, 5, 3, 4, 5])
    >>> count
    array([1, 1, 1, 2, 2, 2, 1, 1, 1])
    """
    links_per_node = link_count_per_node(node_at_link_ends,
                                         number_of_nodes=number_of_nodes)

    in_links = JaggedArray(
        *in_link_ids_at_node(node_at_link_ends,
                             number_of_nodes=number_of_nodes))
    out_links = JaggedArray(
        *out_link_ids_at_node(node_at_link_ends,
                              number_of_nodes=number_of_nodes))

    links = np.empty(in_links.size + out_links.size, dtype=int)

    offset = 0
    for node, link_count in enumerate(links_per_node):
        links[offset:offset + link_count] = np.concatenate(
            (in_links.row(node), out_links.row(node)))
        offset += link_count

    return links, links_per_node


class LinkGrid(object):
    """Create a grid of links that enter and leave nodes.
    __init__((node0, node1), number_of_nodes=None)

    Parameters
    ----------
    node0, node1 : sequence of array-like
        Node ID at link start and end.
    number_of_nodes : int, optional
        Number of nodes in the grid

    Returns
    -------
    LinkGrid :
        A newly-created grid

    Examples
    --------
    >>> from landlab.grid.unstructured.links import LinkGrid
    >>> lgrid = LinkGrid([(0, 1, 0, 2, 0), (2, 3, 1, 3, 3)], 4)
    >>> lgrid.number_of_links
    5
    >>> lgrid.number_of_nodes
    4
    >>> lgrid.number_of_in_links_at_node(0)
    0
    >>> lgrid.number_of_out_links_at_node(0)
    3
    >>> lgrid.out_link_at_node(0)
    array([0, 2, 4])
    >>> lgrid.nodes_at_link_id(1)
    array([1, 3])

    >>> lgrid = LinkGrid([(0, 1, 0, 2, 0), (2, 3, 1, 3, 3)], 4,
    ...                  link_ids=range(1, 6))
    >>> lgrid.nodes_at_link
    array([[0, 2],
           [1, 3],
           [0, 1],
           [2, 3],
           [0, 3]])
    >>> lgrid.out_link_at_node(0)
    array([1, 3, 5])
    >>> lgrid.nodes_at_link_id(1)
    array([0, 2])
    """

    def __init__(self, link_ends, number_of_nodes, link_ids=None,
                 node_status=None):
        """Create a grid of links that enter and leave nodes.
        __init__((node0, node1), number_of_nodes=None)

        Parameters
        ----------
        node0, node1 : sequence of array-like
            Node ID at link start and end.
        number_of_nodes : int, optional
            Number of nodes in the grid

        Returns
        -------
        LinkGrid :
            A newly-created grid

        Examples
        --------
        >>> from landlab.grid.unstructured.links import LinkGrid
        >>> lgrid = LinkGrid([(0, 1, 0, 2, 0), (2, 3, 1, 3, 3)], 4)
        >>> lgrid.number_of_links
        5
        >>> lgrid.number_of_nodes
        4
        >>> lgrid.number_of_in_links_at_node(0)
        0
        >>> lgrid.number_of_out_links_at_node(0)
        3
        >>> lgrid.out_link_at_node(0)
        array([0, 2, 4])
        >>> lgrid.nodes_at_link_id(1)
        array([1, 3])

        >>> lgrid = LinkGrid([(0, 1, 0, 2, 0), (2, 3, 1, 3, 3)], 4,
        ...                  link_ids=range(1, 6))
        >>> lgrid.nodes_at_link
        array([[0, 2],
               [1, 3],
               [0, 1],
               [2, 3],
               [0, 3]])
        >>> lgrid.out_link_at_node(0)
        array([1, 3, 5])
        >>> lgrid.nodes_at_link_id(1)
        array([0, 2])
        """
        link_ends = _split_link_ends(link_ends)

        self._in_link_at_node = JaggedArray(
            *in_link_ids_at_node(link_ends, link_ids=link_ids,
                                 number_of_nodes=number_of_nodes)

        )
        self._out_link_at_node = JaggedArray(
            *out_link_ids_at_node(link_ends, link_ids=link_ids,
                                  number_of_nodes=number_of_nodes)
        )
        self._link_ends = np.array(link_ends)
        if link_ids is not None:
            self._link_id_map = dict(zip(link_ids, range(len(link_ids))))
            self._link_ids = link_ids

        self._number_of_links = len(link_ends[0])
        self._number_of_nodes = number_of_nodes

        self._node_status = node_status

    @property
    def number_of_links(self):
        """Number of links in the grid.
        """
        return self._number_of_links

    @property
    def number_of_nodes(self):
        """Number of nodes in the grid.
        """
        return self._number_of_nodes

    def number_of_in_links_at_node(self, node):
        """Number of links entering a node.

        Parameters
        ----------
        node : int
            Node ID

        Returns
        -------
        int :
            Number of links entering the node.

        Examples
        --------
        >>> from landlab.grid.unstructured.links import LinkGrid
        >>> lgrid = LinkGrid([(0, 1, 0, 2), (2, 3, 1, 3)], 4)
        >>> [lgrid.number_of_in_links_at_node(node) for node in range(4)]
        [0, 1, 1, 2]
        """
        return self._in_link_at_node.length_of_row(node)

    def number_of_out_links_at_node(self, node):
        """Number of links leaving a node.

        Parameters
        ----------
        node : int
            Node ID

        Returns
        -------
        int :
            Number of links leaving the node.

        Examples
        --------
        >>> from landlab.grid.unstructured.links import LinkGrid
        >>> lgrid = LinkGrid([(0, 1, 0, 2), (2, 3, 1, 3)], 4)
        >>> [lgrid.number_of_out_links_at_node(node) for node in range(4)]
        [2, 1, 1, 0]
        """
        return self._out_link_at_node.length_of_row(node)

    def number_of_links_at_node(self, node):
        """Number of links entering and leaving a node.

        Parameters
        ----------
        node : int
            Node ID

        Returns
        -------
        int :
            Number of links entering and leaving the node.

        Examples
        --------
        >>> from landlab.grid.unstructured.links import LinkGrid
        >>> lgrid = LinkGrid([(0, 1, 0, 2), (2, 3, 1, 3)], 4)
        >>> [lgrid.number_of_links_at_node(node) for node in range(4)]
        [2, 2, 2, 2]
        """
        return (self.number_of_in_links_at_node(node) +
                self.number_of_out_links_at_node(node))

    @property
    def node_at_link_start(self):
        return self._link_ends[0]

    @property
    def node_at_link_end(self):
        return self._link_ends[1]

    @property
    def nodes_at_link(self):
        return self._link_ends.T

    @property
    def link_id(self):
        try:
            return self._link_ids
        except AttributeError:
            return np.arange(self.number_of_links)

    def nodes_at_link_id(self, link_id):
        try:
            return self.nodes_at_link[self._link_id_map[link_id]]
        except AttributeError:
            return self.nodes_at_link[link_id]

    def in_link_at_node(self, node):
        """Links entering a node.

        Parameters
        ----------
        node : int
            Node ID

        Returns
        -------
        ndarray :
            Links entering the node

        Examples
        --------
        >>> from landlab.grid.unstructured.links import LinkGrid
        >>> lgrid = LinkGrid([(0, 1, 0, 2), (2, 3, 1, 3)], 4)
        >>> len(lgrid.in_link_at_node(0)) == 0
        True
        >>> lgrid.in_link_at_node(3)
        array([1, 3])
        """
        return self._in_link_at_node.row(node)

    def out_link_at_node(self, node):
        """Links leaving a node.

        Parameters
        ----------
        node : int
            Node ID

        Returns
        -------
        ndarray :
            Links leaving the node

        Examples
        --------
        >>> from landlab.grid.unstructured.links import LinkGrid
        >>> lgrid = LinkGrid([(0, 1, 0, 2), (2, 3, 1, 3)], 4)
        >>> lgrid.out_link_at_node(0)
        array([0, 2])
        >>> len(lgrid.out_link_at_node(3)) == 0
        True
        """
        return self._out_link_at_node.row(node)

    def iter_nodes(self):
        """Iterate of the nodes of the grid.

        Returns
        -------
        ndarray :
            Links entering and leaving each node

        Examples
        --------
        >>> from landlab.grid.unstructured.links import LinkGrid
        >>> lgrid = LinkGrid([(0, 1, 0, 2), (2, 3, 1, 3)], 4)
        >>> for link in lgrid.iter_nodes(): link
        array([0, 2])
        array([2, 1])
        array([0, 3])
        array([1, 3])
        """
        for node in range(self.number_of_nodes):
            yield np.concatenate((
                self.in_link_at_node(node),
                self.out_link_at_node(node),
            ))

    @property
    def node_status_at_link_start(self):
        return self._node_status[self.node_at_link_start]

    @property
    def node_status_at_link_end(self):
        return self._node_status[self.node_at_link_end]
