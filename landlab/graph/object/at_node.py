import numpy as np
from ...utils.jaggedarray import unravel


# def get_links_at_node(graph):
#     (links_at_node,
#      _,
#      offset_to_node) = calc_links_at_node(graph.nodes_at_link)
# 
#     return unravel(links_at_node, offset_to_node)


def get_link_dirs_at_node(graph):
    (_,
     link_dirs_at_node,
     offset_to_node) = calc_links_at_node(graph.nodes_at_link)

    return unravel(link_dirs_at_node, offset_to_node)


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
    # from .cfuncs import _setup_links_at_node
    from .ext.at_node import get_links_at_node

    sorted_by_node = np.argsort(graph.nodes_at_link.flat)

    node_count = np.bincount(graph.nodes_at_link.flat)
    number_of_nodes = graph.number_of_nodes

    max_node_count = np.max(node_count)

    link_dirs_at_node = np.full((number_of_nodes, max_node_count), 0,
                                dtype=int)
    links_at_node = np.full((number_of_nodes, max_node_count), -1, dtype=int)

    get_links_at_node(graph.nodes_at_link, links_at_node, link_dirs_at_node)

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
    from .ext.at_node import reorder_links_at_node

    outward_angle = angle_of_link[links_at_node]

    links_entering = np.where(link_dirs_at_node == 1)
    outward_angle[links_entering] += np.pi
    outward_angle[outward_angle >= 2 * np.pi] -= 2 * np.pi

    outward_angle[np.where(link_dirs_at_node == 0)] = 4 * np.pi

    sorted_links = np.argsort(outward_angle)

    if not inplace:
        links_at_node = links_at_node.copy()
        link_dirs_at_node = link_dirs_at_node.copy()

    reorder_links_at_node(links_at_node, sorted_links)
    reorder_links_at_node(link_dirs_at_node, sorted_links)

    return links_at_node, link_dirs_at_node
