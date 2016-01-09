import numpy as np

from ..core.utils import as_id_array, argsort_points_by_x_then_y
from ..utils.jaggedarray import flatten_jagged_array


def sort_graph(nodes, links=None, patches=None):
    """Sort elements of a graph by x, then y.

    Parameters
    ----------
    nodes : tuple of ndarray
        Coordinate of nodes as (y, x).
    links : ndarray of int, shape `(n_links, 2)`, optional
        Indices into *nodes* array of link tail then head.
    patches : array_like of array_like, optional
        Indices into *links* array of the links that define each patch.

    Returns
    -------
    (nodes, links, patches)
        Sorted nodes, links, and patches.

    Examples
    --------
    o---o---o
    |   | / |
    o---o---o

    >>> from landlab.graph.sort import sort_graph
    >>> import numpy as np
    >>> x = np.array([1., 2., 2., 0., 1., 0.])
    >>> y = np.array([0., 0., 1., 0., 1., 1.])

    Sort a graph with just points - no links or patches.

    >>> _ = sort_graph((y, x))
    >>> y
    array([ 0.,  0.,  0.,  1.,  1.,  1.])
    >>> x
    array([ 0.,  1.,  2.,  0.,  1.,  2.])

    Sort the points and links of a graph.

    >>> x = np.array([1., 2., 2., 0., 1., 0.])
    >>> y = np.array([0., 0., 1., 0., 1., 1.])
    >>> links = np.array([[3, 0], [0, 4], [4, 5], [5, 3], [0, 1], [1, 2],
    ...                   [2, 0], [2, 4]])
    >>> _ = sort_graph((y, x), links)
    >>> links # doctest: +NORMALIZE_WHITESPACE
    array([[0, 1], [1, 2], [3, 0], [1, 4], [5, 1], [2, 5], [4, 3], [5, 4]])

    Sort the points, links, and patches of a graph.

    >>> x = np.array([1., 2., 2., 0., 1., 0.])
    >>> y = np.array([0., 0., 1., 0., 1., 1.])
    >>> links = np.array([[3, 0], [0, 4], [4, 5], [5, 3], [0, 1], [1, 2],
    ...                   [2, 0], [2, 4]])
    >>> patches = (np.array([1, 6, 7, 4, 5, 6, 0, 1, 2, 3]),
    ...            np.array([0, 3, 6, 10]))
    >>> _ = sort_graph((y, x), links, patches)
    >>> patches[0]
    array([1, 5, 4, 0, 3, 6, 2, 3, 4, 7])
    >>> patches[1]
    array([ 0,  3,  7, 10])
    """
    from .cfuncs import _remap_nodes_at_link, _remap_links_at_patch

    if patches is not None and links is None:
        raise ValueError('graph that has patches must also have links')

    if links is not None:
        links = np.asarray(links)
    if patches is not None:
        if len(patches) == 2 and isinstance(patches[0], np.ndarray):
            links_at_patch, offset_to_patch = patches
        else:
            links_at_patch, offset_to_patch = flatten_jagged_array(patches,
                                                                   dtype=int)
    else:
        links_at_patch, offset_to_patch = (None, None)

    sorted_nodes = sort_nodes(nodes)

    if links is not None:
        _remap_nodes_at_link(links, np.argsort(sorted_nodes, kind='mergesort'))

        midpoint_of_link = np.empty((len(links), 2), dtype=float)
        sorted_links = sort_links(links, nodes,
                                  midpoint_of_link=midpoint_of_link) 

    if patches is not None:
        _remap_links_at_patch(links_at_patch,
                              np.argsort(sorted_links, kind='mergesort'))
        sort_patches(links_at_patch, offset_to_patch, midpoint_of_link)

    return nodes, links, (links_at_patch, offset_to_patch)


def sort_nodes(nodes):
    """Sort nodes based on their position.

    Parameters
    ----------
    nodes : tuple of ndarray of float
        Coordinates of nodes as (*y*, *x*).

    Returns
    -------
    ndarray of int
        Array of indices that sort the nodes.

    Examples
    --------
    >>> from landlab.graph.sort import sort_nodes
    >>> import numpy as np
    >>> x = np.array([0. , 1., 2.])
    >>> y = np.array([ .5, 0., 1.])
    >>> sort_nodes((y, x))
    array([1, 0, 2])
    >>> x
    array([ 1.,  0.,  2.])
    >>> y
    array([ 0. ,  0.5,  1. ])
    """
    from .cfuncs import _remap_nodes_at_link

    sorted_nodes = argsort_points_by_x_then_y((nodes[1], nodes[0]))
    nodes[0][:] = nodes[0][sorted_nodes]
    nodes[1][:] = nodes[1][sorted_nodes]

    return sorted_nodes


def sort_links(nodes_at_link, nodes, midpoint_of_link=None):
    """Sort links by their midpoint.

    Parameters
    ----------
    nodes_at_link : ndarray of int, shape `(n_links, 2)`
        Node for each link tail and head.
    nodes : tuple of ndarray of float
        Node coordinates.
    midpoint_of_link : ndarray of float, shape `(n_links, 2)`, optional
        Buffer to store the link midpoints that were used for sorting.

    Returns
    -------
    ndarray of int
        Array of indices that sort the links.

    Examples
    --------
    >>> from landlab.graph.sort import sort_nodes
    >>> import numpy as np
    >>> nodes = np.array([[0, 0, 0, 1, 1, 1],
    ...                   [0, 1, 2, 0, 1, 2]])
    >>> links = np.array([[0, 1], [0, 3], [1, 2], [1, 4], [2, 5], [3, 4],
    ...                   [4, 5]])
    >>> sort_links(links, nodes)
    array([0, 2, 1, 3, 4, 5, 6])
    >>> links
    array([[0, 1],
           [1, 2],
           [0, 3],
           [1, 4],
           [2, 5],
           [3, 4],
           [4, 5]])
    """
    y_of_node, x_of_node = nodes

    y_at_tail = y_of_node[nodes_at_link[:, 0]]
    x_at_tail = x_of_node[nodes_at_link[:, 0]]

    y_at_head = y_of_node[nodes_at_link[:, 1]]
    x_at_head = x_of_node[nodes_at_link[:, 1]]

    if midpoint_of_link is None:
        midpoint_of_link = np.empty((len(nodes_at_link), 2), dtype=float)

    midpoint_of_link[:, 0] = (x_at_tail + x_at_head) * .5
    midpoint_of_link[:, 1] = (y_at_tail + y_at_head) * .5

    sorted_links = argsort_points_by_x_then_y(midpoint_of_link)
    nodes_at_link[:] = nodes_at_link[sorted_links]

    return sorted_links


def sort_patches(links_at_patch, offset_to_patch, xy_of_link):
    """Sort patches by their centroid.

    Parameters
    ----------
    links_at_patch : ndarray of int
        Links that define patches.
    offset_to_patch : ndarray of int
        Offsets into *links_at_patch* for each patch.
    xy_of_link : ndarray of float, shape `(n_links, 2)`
        Midpoint coordinates for each link.

    Examples
    --------
    >>> from landlab.graph.sort import sort_patches
    >>> import numpy as np
    >>> links_at_patch = np.array([0, 1, 2, 3, 2, 4])
    >>> offset_to_patch = np.array([0, 3, 6])
    >>> xy_of_link = np.array([[0.0, 0.5], [0.5, 1.0], [.5, .5],
    ...                        [0.5, 0.0], [1.0, 0.5]])
    >>> sort_patches(links_at_patch, offset_to_patch, xy_of_link)
    array([1, 0])
    >>> links_at_patch
    array([3, 2, 4, 0, 1, 2])
    >>> offset_to_patch
    array([0, 3, 6])
    """
    from .cfuncs import _calc_center_of_patch, _resort_patches

    n_patches = len(offset_to_patch) - 1
    xy_at_patch = np.empty((n_patches, 2), dtype=float)

    _calc_center_of_patch(links_at_patch, offset_to_patch,
                          xy_of_link, xy_at_patch)

    sorted_patches = argsort_points_by_x_then_y(xy_at_patch)
    _resort_patches(links_at_patch, offset_to_patch, sorted_patches)

    return sorted_patches
