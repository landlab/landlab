import numpy as np


def get_angle_of_link(graph, out=None):
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
    if out is None:
        out = np.empty(graph.number_of_links, dtype=float)

    y = graph.y_of_node[graph.nodes_at_link[:, ::-1]]
    x = graph.x_of_node[graph.nodes_at_link[:, ::-1]]

    np.arctan2(np.diff(y).flat, np.diff(x).flat, out=out)
    out += np.pi
    out[out == 2. * np.pi] == 0.

    return out


def get_midpoint_of_link(graph, out=None):
    from ext.remap_element import calc_midpoint_of_link

    if out is None:
        out = np.empty((graph.number_of_links, 2), dtype=float)

    calc_midpoint_of_link(graph.nodes_at_link, graph.x_of_node,
                          graph.y_of_node, out)

    return out
