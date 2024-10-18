import numpy as np

from landlab.core.utils import as_id_array


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
    from landlab.graph.object.ext.at_node import get_links_at_node

    node_count = np.bincount(graph.nodes_at_link.flat)
    number_of_nodes = graph.number_of_nodes

    max_node_count = np.max(node_count)

    link_dirs_at_node = np.full((number_of_nodes, max_node_count), 0, dtype=np.int8)
    links_at_node = np.full((number_of_nodes, max_node_count), -1, dtype=int)

    get_links_at_node(graph.nodes_at_link, links_at_node, link_dirs_at_node)

    if sort:
        sort_links_at_node_by_angle(
            links_at_node, link_dirs_at_node, graph.angle_of_link, inplace=True
        )

    return links_at_node, link_dirs_at_node


def sort_links_at_node_by_angle(
    links_at_node, link_dirs_at_node, angle_of_link, inplace=False
):
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

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.graph.object.at_node import sort_links_at_node_by_angle

    ::

        (1) - 1 -> (3)
         |          ^
         2          3
         V          |
        (0) - 0 -> (2)


    >>> links_at_node = [[2, 0], [1, 2], [3, 0], [3, 1]]
    >>> link_dirs_at_node = [[1, -1], [-1, -1], [-1, 1], [1, 1]]
    >>> angle_of_link = np.array([0.0, 0.0, -90.0, 90.0]) * np.pi / 180.0

    >>> out = sort_links_at_node_by_angle(
    ...     links_at_node, link_dirs_at_node, angle_of_link
    ... )

    The first item of the returned tuple is links at each node sorted
    counterclockwise by angle.

    >>> out[0]
    array([[0, 2], [2, 1], [3, 0],  [1, 3]])

    The second item is the direction of the link (entering or leaving).

    >>> out[1]
    array([[-1,  1], [-1, -1], [-1,  1], [ 1,  1]], dtype=int8)

    Because the input arrays are lists, not numpy arrays, the sort is not
    in-place.

    >>> out[0] is not links_at_node
    True

    >>> links_at_node = np.asarray([[2, 0], [1, 2], [3, 0], [3, 1]], dtype=int)
    >>> link_dirs_at_node = np.asarray(
    ...     [[1, -1], [-1, -1], [-1, 1], [1, 1]], dtype=np.int8
    ... )
    >>> angle_of_link = np.array([0.0, 0.0, -90.0, 90.0]) * np.pi / 180.0

    >>> _ = sort_links_at_node_by_angle(
    ...     links_at_node, link_dirs_at_node, angle_of_link, inplace=True
    ... )
    >>> links_at_node
    array([[0, 2], [2, 1], [3, 0],  [1, 3]])
    >>> link_dirs_at_node
    array([[-1,  1], [-1, -1], [-1,  1], [ 1,  1]], dtype=int8)
    """
    from landlab.graph.object.ext.at_node import reorder_rows_inplace

    out = (
        np.asarray(links_at_node, dtype=int),
        np.asarray(link_dirs_at_node, dtype=np.int8),
    )

    if inplace and (out[0] is not links_at_node or out[1] is not link_dirs_at_node):
        raise ValueError(
            "links_at_node and link_dirs_at_node arrays must be ndarray for in-place sort"
        )

    if not inplace:
        if out[0] is links_at_node:
            out[0] = links_at_node.copy()
        if out[1] is link_dirs_at_node:
            out[1] = link_dirs_at_node.copy()

    links_at_node, link_dirs_at_node = out

    outward_angle = angle_of_link[links_at_node]

    links_entering = np.where(link_dirs_at_node == 1)
    outward_angle[links_entering] += np.pi
    outward_angle[outward_angle >= 2 * np.pi] -= 2 * np.pi

    outward_angle[np.where(link_dirs_at_node == 0)] = 4 * np.pi

    sorted_links = as_id_array(np.argsort(outward_angle))

    reorder_rows_inplace(links_at_node, sorted_links)
    reorder_rows_inplace(link_dirs_at_node, sorted_links)

    return links_at_node, link_dirs_at_node
