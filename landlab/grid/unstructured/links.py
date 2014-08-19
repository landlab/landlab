import numpy as np


from .status import CORE_NODE, CLOSED_BOUNDARY


def _split_link_ends(link_ends):
    """
    Examples
    --------
    >>> _split_link_ends(((0, 1, 2), (3, 4, 5)))
    ((0, 1, 2), (3, 4, 5))
    >>> _split_link_ends([(0, 3), (1, 4), (2, 5)])
    ((0, 1, 2), (3, 4, 5))
    >>> _split_link_ends((0, 3))
    ((0,), (3,))
    """
    if len(link_ends) < 2:
        raise ValueError('Link array must be at least of length 2')
    elif len(link_ends) == 2:
        start, end = link_ends
    else:
        start, end = zip(*link_ends)

    try:
        if len(start) == len(end):
            return start, end
        else:
            raise ValueError('Link arrays must be the same length')
    except TypeError:
        return (start, ), (end, )


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
    """
    node_at_link_start, node_at_link_end = _split_link_ends(node_at_link_ends)

    if len(node_at_link_end) != len(node_at_link_start):
        raise ValueError('Link arrays must be the same length')

    status_at_link_ends = (node_status[node_at_link_start],
                           node_status[node_at_link_end])

    (active_link_ids, ) = np.where(link_is_active(status_at_link_ends))

    return active_link_ids


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
    >>> link_ends = [(0, 3), (1, 4), (2, 5), (3, 6), (4, 7), (5, 8)]
    >>> in_link_count_per_node(zip(*link_ends))
    array([0, 0, 0, 1, 1, 1, 1, 1, 1])
    """
    node_at_link_start, node_at_link_end = _split_link_ends(node_at_link_ends)

    if len(node_at_link_end) != len(node_at_link_start):
        raise ValueError('Link arrays must be the same length')

    return np.bincount(node_at_link_end, minlength=number_of_nodes)


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
    >>> out_link_count_per_node(([0, 1, 2, 3, 4, 5], [3, 4, 5, 6, 7, 8]))
    array([1, 1, 1, 1, 1, 1])
    >>> out_link_count_per_node(([0, 1, 2, 3, 4, 5], [3, 4, 5, 6, 7, 8]),
    ...     number_of_nodes=9)
    array([1, 1, 1, 1, 1, 1, 0, 0, 0])
    """
    node_at_link_start, node_at_link_end = _split_link_ends(node_at_link_ends)
    if len(node_at_link_end) != len(node_at_link_start):
        raise ValueError('Link arrays must be the same length')
    return np.bincount(node_at_link_start, minlength=number_of_nodes)


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


def in_link_ids_at_node(node_at_link_ends, number_of_nodes=None):
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
    >>> (links, offset) = in_link_ids_at_node(([0, 1, 2, 3, 4, 5],
    ...                                        [3, 4, 5, 6, 7, 8]))
    >>> links
    array([0, 1, 2, 3, 4, 5])
    >>> offset
    array([0, 0, 0, 0, 1, 2, 3, 4, 5, 6])

    Links entering node with id 0, 4 and 8

    >>> for ind in (0, 4, 8): links[offset[ind]:offset[ind +1]]
    array([], dtype=int64)
    array([1])
    array([5])
    """
    node_at_link_start, node_at_link_end = _split_link_ends(node_at_link_ends)

    link_ids = np.argsort(node_at_link_end)
    links_per_node = in_link_count_per_node(node_at_link_ends,
                                            number_of_nodes=number_of_nodes)
    offset = np.empty(len(links_per_node) + 1, dtype=int)
    np.cumsum(links_per_node, out=offset[1:])
    offset[0] = 0
    return link_ids, offset


def out_link_ids_at_node(node_at_link_ends, number_of_nodes=None):
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
    >>> (links, offset) = out_link_ids_at_node(
    ...     ([0, 1, 2, 3, 4, 5], [3, 4, 5, 6, 7, 8]), number_of_nodes=9)
    >>> links
    array([0, 1, 2, 3, 4, 5])
    >>> offset
    array([0, 1, 2, 3, 4, 5, 6, 6, 6, 6])

    Links leaving node with id 0, 4 and 8

    >>> for ind in (0, 4, 8): links[offset[ind]:offset[ind +1]]
    array([0])
    array([4])
    array([], dtype=int64)
    """
    node_at_link_start, node_at_link_end = _split_link_ends(node_at_link_ends)

    link_ids = np.argsort(node_at_link_start)
    links_per_node = out_link_count_per_node(node_at_link_ends,
                                             number_of_nodes=number_of_nodes)
    offset = np.empty(len(links_per_node) + 1, dtype=int)
    np.cumsum(links_per_node, out=offset[1:])
    offset[0] = 0
    return link_ids, offset


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
    >>> (links, offset) = link_ids_at_node(
    ...     ([0, 1, 2, 3, 4, 5], [3, 4, 5, 6, 7, 8]), number_of_nodes=9)
    >>> links
    array([0, 1, 2, 0, 3, 1, 4, 2, 5, 3, 4, 5])
    >>> offset
    array([ 0,  1,  2,  3,  5,  7,  9, 10, 11, 12])
    >>> len(offset)
    10

    Links entering and leaving node with id 0, 4 and 8

    >>> for ind in (0, 4, 8): links[offset[ind]:offset[ind +1]]
    array([0])
    array([1, 4])
    array([5])
    """
    links_per_node = link_count_per_node(node_at_link_ends,
                                         number_of_nodes=number_of_nodes)

    offset = np.empty(len(links_per_node) + 1, dtype=int)
    np.cumsum(links_per_node, out=offset[1:])
    offset[0] = 0

    (in_link_ids, in_offset) = in_link_ids_at_node(
        node_at_link_ends, number_of_nodes=number_of_nodes)
    (out_link_ids, out_offset) = out_link_ids_at_node(
        node_at_link_ends, number_of_nodes=number_of_nodes)

    in_link_count = in_link_count_per_node(node_at_link_ends,
                                           number_of_nodes=number_of_nodes)

    link_ids = np.zeros(offset[-1], dtype=int)
    for node in xrange(number_of_nodes):
        low = offset[node]
        middle = low + in_link_count[node]
        high = offset[node + 1]
        link_ids[low:middle] = in_link_ids[in_offset[node]:in_offset[node + 1]]
        link_ids[middle:high] = out_link_ids[out_offset[node]:out_offset[node + 1]]

    return link_ids, offset
