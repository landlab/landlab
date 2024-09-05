import numpy as np

from ...core.utils import as_id_array


def get_links(nodes_at_link, sort=True):
    """Get links and their directions at each node.

    Parameters
    ----------
    nodes_at_link : ndarray of int, shape `(n_links, 2)`
        Node identifier for each link tail and head.

    Returns
    -------
    tuple of ndarray of int
        Tuple of link identifiers for each node, link directions for each
        node and offsets into these arrays for each node.

    Examples
    --------
    ::
        0 -- 1
        |    | \
        2 -- 3 - 4

    >>> import numpy as np
    >>> nodes_at_link = np.array([[2, 3], [3, 4], [2, 0], [3, 1], [4, 1], [0, 1]])
    >>> (links_at_node, link_dirs_at_node, offset_to_node) = get_links(nodes_at_link)
    >>> links_at_node
    array([2, 5, 3, 4, 5, 0, 2, 0, 1, 3, 1, 4])
    >>> link_dirs_at_node
    array([ 1, -1,  1,  1,  1, -1, -1,  1, -1, -1,  1, -1], dtype=int8)
    >>> offset_to_node
    array([ 0,  2,  5,  7, 10, 12])
    """
    sorted = np.argsort(nodes_at_link.reshape((-1,)), kind="stable")

    links_at_node = sorted // 2
    link_dirs_at_node = sorted % 2
    link_dirs_at_node[link_dirs_at_node == 0] = -1

    (_, links_per_node) = np.unique(nodes_at_link, return_counts=True)

    offset_to_node = np.empty(len(links_per_node) + 1, dtype=int)
    offset_to_node[0] = 0
    links_per_node.cumsum(out=offset_to_node[1:])

    return (
        as_id_array(links_at_node),
        np.asarray(link_dirs_at_node, dtype=np.int8),
        offset_to_node,
    )
