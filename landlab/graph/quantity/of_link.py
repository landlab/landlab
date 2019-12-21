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
    out : ndarray of float, optional
        Buffer to place output.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.graph import UniformRectilinearGraph
    >>> from landlab.graph.quantity.of_link import get_angle_of_link

    >>> graph = UniformRectilinearGraph((2, 2))
    >>> get_angle_of_link(graph) * 180. / np.pi
    array([  0.,  90.,  90.,   0.])

    >>> angles = np.empty(graph.number_of_links, dtype=float)
    >>> rtn = get_angle_of_link(graph, out=angles)
    >>> angles is rtn
    True
    """
    if out is None:
        out = np.empty(graph.number_of_links, dtype=float)

    y = graph.y_of_node[graph.nodes_at_link]
    x = graph.x_of_node[graph.nodes_at_link]

    np.arctan2(np.diff(y).flat, np.diff(x).flat, out=out)

    return np.mod(out, 2.0 * np.pi, out=out)


def get_midpoint_of_link(graph, out=None):
    """Get xy position of the midpoint of links.

    Parameters
    ----------
    graph : `Graph`
        A `Graph`.
    out : ndarray of float, optional
        Buffer to place output.

    Returns
    -------
    ndarray of float, shape `(n_links, 2)`
        Coordinates of link midpoints.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.graph import UniformRectilinearGraph
    >>> from landlab.graph.quantity.of_link import get_midpoint_of_link

    >>> graph = UniformRectilinearGraph((2, 2))
    >>> get_midpoint_of_link(graph)
    array([[ 0.5,  0. ],
           [ 0. ,  0.5],
           [ 1. ,  0.5],
           [ 0.5,  1. ]])

    >>> points = np.empty((graph.number_of_links, 2), dtype=float)
    >>> rtn = get_midpoint_of_link(graph, out=points)
    >>> points is rtn
    True
    """
    from .ext.of_link import calc_midpoint_of_link

    if out is None:
        out = np.empty((graph.number_of_links, 2), dtype=float)

    calc_midpoint_of_link(graph.nodes_at_link, graph.x_of_node, graph.y_of_node, out)

    return out


def get_length_of_link(graph):
    nodes_at_link = graph.nodes_at_link
    dx = graph.x_of_node[nodes_at_link[:, 0]] - graph.x_of_node[nodes_at_link[:, 1]]
    dy = graph.y_of_node[nodes_at_link[:, 0]] - graph.y_of_node[nodes_at_link[:, 1]]
    return np.sqrt(dx ** 2 + dy ** 2)
