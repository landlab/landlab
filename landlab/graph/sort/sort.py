"""Sort the elements of a graph.

This module provides functions that sort the elements of a graph
structure.
"""

import numpy as np

from ...core.utils import argsort_points_by_x_then_y
from ...core.utils import as_id_array
from ...utils.jaggedarray import flatten_jagged_array
from ..quantity.ext.of_element import mean_of_children_at_parent
from .ext.argsort import sort_id_array


def remap(src, mapping, out=None, inplace=False):
    """Remap elements in an id array.

    Parameters
    ----------
    src : ndarray of int
        Initial array of ids.
    mapping : ndarray of int
        Mapping of ids.
    out : ndarray of int, optional
        Buffer into which to place output.
    inplace : bool, optional
        Mapping of values will include inplace.

    Returns
    -------
    ndarray of int
        Array of remapped values.

    Examples
    --------
    >>> from landlab.graph.sort.sort import remap
    >>> import numpy as np

    >>> src = np.array([1, 2, 3, 4])
    >>> mapping = np.array([-1, 10, 20, 30, 40])
    >>> remap(src, mapping)
    array([10, 20, 30, 40])
    """
    from .ext.remap_element import remap_graph_element

    if inplace:
        out = src
    else:
        if out is None:
            out = src.copy()
        else:
            out[:] = src[:]

    remap_graph_element(out.reshape((-1,)), mapping)

    return out


def reverse_one_to_one(ids, minlength=None):
    """Reverse a one-to-one mapping.

    Parameters
    ----------
    ids : ndarray of int, shape `(N, )`
        Array of identifier mapping.
    minlength : int, optional
        A minimum number of identifiers for the output array.
    Returns
    -------
    ndarray of int, shape `(n, )`
        Array of the reverse mapping.

    Examples
    --------
    >>> from landlab.graph.sort.sort import reverse_one_to_one
    >>> ids = np.array([-1, -1, 6, 3, -1, 2, 4, 1, -1, 5, 0], dtype=int)
    >>> reverse_one_to_one(ids)
    array([10,  7,  5,  3,  6,  9,  2])
    """
    from .ext.remap_element import reverse_one_to_one

    if minlength is None:
        minlength = ids.max() + 1
    out = np.full((minlength,), -1, dtype=int)

    reverse_one_to_one(ids, out)

    return out


def reverse_one_to_many(ids, min_counts=0):
    """Reverse a one-to-many mapping.

    Parameters
    ----------
    ids : ndarray of int, shape `(M, N)`
        Array of identifier mapping.

    Returns
    -------
    ndarray of int, shape `(m, n)`
        Array of the reverse mapping.

    Examples
    --------
    >>> from landlab.graph.sort.sort import reverse_one_to_many
    >>> ids = np.array([[1, 2, 3], [-1, -1, -1], [2, 3, -1]], dtype=int)
    >>> reverse_one_to_many(ids)
    array([[-1, -1],
           [ 0, -1],
           [ 0,  2],
           [ 0,  2]])
    """
    from .ext.remap_element import reverse_one_to_many

    counts = np.bincount(ids.reshape((-1,)) + 1)
    max_counts = np.max((np.max(counts[1:]), min_counts))

    out = np.full((ids.max() + 1, max_counts), -1, dtype=int)

    reverse_one_to_many(ids, out)

    return out


def reorder_links_at_patch(graph):
    from ..matrix.ext.matrix import roll_id_matrix_rows
    from ..object.ext.at_patch import get_rightmost_edge_at_patch
    from ..quantity.of_link import get_midpoint_of_link
    from ..quantity.of_patch import get_area_of_patch
    from .ext.remap_element import reverse_element_order

    if graph.number_of_patches == 0:
        return

    xy_of_link = get_midpoint_of_link(graph)

    shift = np.empty(graph.number_of_patches, dtype=int)

    get_rightmost_edge_at_patch(graph.links_at_patch, xy_of_link, shift)
    roll_id_matrix_rows(graph.links_at_patch, -shift)

    before = graph.links_at_patch.copy()
    area_before = get_area_of_patch(graph)

    negative_areas = as_id_array(np.where(get_area_of_patch(graph) < 0.0)[0])
    reverse_element_order(graph.links_at_patch, negative_areas)
    # reverse_element_order(graph._links_at_patch, negative_areas)

    # graph._nodes_at_patch = get_nodes_at_patch(graph)
    if "nodes_at_patch" in graph._ds:
        graph._ds = graph._ds.drop("nodes_at_patch")

    if np.any(get_area_of_patch(graph) < 0.0):
        raise ValueError(
            (graph.links_at_patch, before, get_area_of_patch(graph), area_before)
        )


def reorient_link_dirs(graph):
    from ..quantity.of_link import get_angle_of_link

    if graph.number_of_links == 0:
        return

    angles = get_angle_of_link(graph)
    links_to_swap = (angles < 7.0 * np.pi / 4.0) & (angles >= np.pi * 0.75)
    graph.nodes_at_link[links_to_swap, :] = graph.nodes_at_link[links_to_swap, ::-1]


def sort_links_at_patch(links_at_patch, nodes_at_link, xy_of_node):
    """Reorder links around a patch to be counterclockwise.

    Examples
    --------
    >>> from landlab.graph.sort.sort import sort_links_at_patch
    >>> import numpy as np
    >>> xy_of_node = np.array(
    ...     [
    ...         [0.0, 0.0],
    ...         [0.0, 1.0],
    ...         [1.0, 1.0],
    ...         [1.0, 0.0],
    ...     ]
    ... )
    >>> nodes_at_link = np.array([[0, 1], [1, 2], [2, 3], [3, 0]])
    >>> links_at_patch = np.array([[0, 1, 3, 2]])
    >>> sort_links_at_patch(links_at_patch, nodes_at_link, xy_of_node)
    >>> links_at_patch
    array([[2, 1, 0, 3]])

    >>> xy_of_node = np.array(
    ...     [
    ...         [0.0, 0.0],
    ...         [0.0, 1.0],
    ...         [1.0, 1.0],
    ...         [1.0, 0.0],
    ...         [2.0, 0.0],
    ...     ]
    ... )
    >>> nodes_at_link = np.array([[0, 1], [1, 2], [2, 3], [3, 0], [2, 4], [3, 4]])
    >>> links_at_patch = np.array([[0, 1, 3, 2], [2, 4, 5, -1]])
    >>> sort_links_at_patch(links_at_patch, nodes_at_link, xy_of_node)
    >>> links_at_patch
    array([[ 2,  1,  0,  3],
           [ 4,  2,  5, -1]])

    >>> xy_of_node = np.array(
    ...     [
    ...         [0.0, 0.0],
    ...         [0.0, 1.0],
    ...         [1.0, 1.0],
    ...         [1.0, 0.0],
    ...         [2.0, 1.0],
    ...     ]
    ... )
    >>> nodes_at_link = np.array([[0, 1], [1, 2], [2, 3], [3, 0], [2, 4], [3, 4]])
    >>> links_at_patch = np.array([[0, 1, 3, 2], [2, -1, 4, 5]])
    >>> sort_links_at_patch(links_at_patch, nodes_at_link, xy_of_node)
    >>> links_at_patch
    array([[ 2,  1,  0,  3],
           [ 4,  2,  5, -1]])
    """
    xy_of_link = np.mean(xy_of_node[nodes_at_link], axis=1)

    xy_of_patch = np.empty((len(links_at_patch), 2), dtype=float)

    mean_of_children_at_parent(links_at_patch, xy_of_link[:, 0], xy_of_patch[:, 0])
    mean_of_children_at_parent(links_at_patch, xy_of_link[:, 1], xy_of_patch[:, 1])

    sort_spokes_at_hub(links_at_patch, xy_of_patch, xy_of_link, inplace=True)


def reindex_by_xy(graph):
    sorted_nodes = reindex_nodes_by_xy(graph)

    if "nodes_at_link" in graph.ds:
        sorted_links = reindex_links_by_xy(graph)
    else:
        sorted_links = None

    if "links_at_patch" in graph.ds:
        sorted_patches = reindex_patches_by_xy(graph)
    else:
        sorted_patches = None

    return sorted_nodes, sorted_links, sorted_patches


def reindex_patches_by_xy(graph):
    from ..quantity.of_patch import get_centroid_of_patch

    if graph.number_of_patches == 1:
        return np.array([0], dtype=int)

    xy_at_patch = get_centroid_of_patch(graph)

    y = xy_at_patch[:, 1]
    y_min, y_max = y.min(), y.max()
    if np.isclose(y_min, y_max):
        y_max = y_min + 1.0

    sorted_patches = argsort_points_by_x_then_y(
        (xy_at_patch[:, 0], np.round((y - y_min) / (y_max - y_min), decimals=5))
    )

    graph.links_at_patch[:] = graph.links_at_patch[sorted_patches, :]

    if "nodes_at_patch" in graph._ds:
        graph._ds = graph.ds.drop("nodes_at_patch")

    return sorted_patches


def reindex_links_by_xy(graph):
    from ..quantity.of_link import get_midpoint_of_link
    from .ext.remap_element import remap_graph_element_ignore

    xy_of_link = get_midpoint_of_link(graph)

    sorted_links = argsort_points_by_x_then_y(xy_of_link)

    # graph._nodes_at_link[:] = graph._nodes_at_link[sorted_links, :]
    graph.nodes_at_link[:] = graph.nodes_at_link[sorted_links, :]

    # if hasattr(graph, '_links_at_patch'):
    if "links_at_patch" in graph.ds:
        remap_graph_element_ignore(
            graph.links_at_patch.reshape((-1,)),
            as_id_array(np.argsort(sorted_links, kind="stable")),
            -1,
        )

    return sorted_links


def reindex_nodes_by_xy(graph):
    from .ext.remap_element import remap_graph_element

    graph.y_of_node[:] = np.round(graph.y_of_node, decimals=6)

    sorted_nodes = argsort_points_by_x_then_y((graph.x_of_node, graph.y_of_node))

    graph.y_of_node[:] = graph.y_of_node[sorted_nodes]
    graph.x_of_node[:] = graph.x_of_node[sorted_nodes]

    if "nodes_at_link" in graph.ds:
        remap_graph_element(
            graph.nodes_at_link.reshape((-1,)),
            as_id_array(np.argsort(sorted_nodes, kind="stable")),
        )

    if "nodes_at_patch" in graph.ds:
        remap_graph_element(
            graph.nodes_at_patch.reshape((-1,)),
            as_id_array(np.argsort(sorted_nodes, kind="stable")),
        )

    return sorted_nodes


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
    >>> x = np.array([1.0, 2.0, 2.0, 0.0, 1.0, 0.0])
    >>> y = np.array([0.0, 0.0, 1.0, 0.0, 1.0, 1.0])

    Sort a graph with just points - no links or patches.

    >>> _ = sort_graph((y, x))
    >>> y
    array([0., 0., 0., 1., 1., 1.])
    >>> x
    array([0., 1., 2., 0., 1., 2.])

    Sort the points and links of a graph.

    >>> x = np.array([1.0, 2.0, 2.0, 0.0, 1.0, 0.0])
    >>> y = np.array([0.0, 0.0, 1.0, 0.0, 1.0, 1.0])
    >>> links = np.array(
    ...     [[3, 0], [0, 4], [4, 5], [5, 3], [0, 1], [1, 2], [2, 0], [2, 4]]
    ... )
    >>> _ = sort_graph((y, x), links)
    >>> links
    array([[0, 1], [1, 2], [3, 0], [1, 4], [5, 1], [2, 5], [4, 3], [5, 4]])

    Sort the points, links, and patches of a graph.

    >>> x = np.array([1.0, 2.0, 2.0, 0.0, 1.0, 0.0])
    >>> y = np.array([0.0, 0.0, 1.0, 0.0, 1.0, 1.0])
    >>> links = np.array(
    ...     [[3, 0], [0, 4], [4, 5], [5, 3], [0, 1], [1, 2], [2, 0], [2, 4]]
    ... )
    >>> patches = (np.array([1, 6, 7, 4, 5, 6, 0, 1, 2, 3]), np.array([0, 3, 6, 10]))
    >>> _ = sort_graph((y, x), links, patches)
    >>> patches[0]
    array([1, 5, 4, 0, 3, 6, 2, 3, 4, 7])
    >>> patches[1]
    array([ 0,  3,  7, 10])
    """
    from .ext.remap_element import remap_graph_element

    if patches is not None and links is None:
        raise ValueError("graph that has patches must also have links")

    if links is not None:
        links = as_id_array(links)

    if patches is not None:
        if len(patches) == 2 and isinstance(patches[0], np.ndarray):
            links_at_patch, offset_to_patch = patches
        else:
            links_at_patch, offset_to_patch = flatten_jagged_array(patches, dtype=int)
        links_at_patch, offset_to_patch = (
            as_id_array(links_at_patch),
            as_id_array(offset_to_patch),
        )
    else:
        links_at_patch, offset_to_patch = (None, None)

    sorted_nodes = sort_nodes(nodes)

    if links is not None:
        remap_graph_element(
            links.reshape((-1,)),
            as_id_array(np.argsort(sorted_nodes, kind="mergesort")),
        )
        midpoint_of_link = np.empty((len(links), 2), dtype=float)
        sorted_links = sort_links(links, nodes, midpoint_of_link=midpoint_of_link)

    if patches is not None:
        remap_graph_element(
            links_at_patch, as_id_array(np.argsort(sorted_links, kind="mergesort"))
        )
        sort_patches(links_at_patch, offset_to_patch, midpoint_of_link)

    if links_at_patch is None:
        return nodes, links, None
    else:
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
    >>> x = np.array([0.0, 1.0, 2.0])
    >>> y = np.array([0.5, 0.0, 1.0])
    >>> sort_nodes((y, x))
    array([1, 0, 2])
    >>> x
    array([1., 0., 2.])
    >>> y
    array([0. , 0.5, 1. ])
    """
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
    >>> nodes = np.array([[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 2]])
    >>> links = np.array([[0, 1], [0, 3], [1, 2], [1, 4], [2, 5], [3, 4], [4, 5]])
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
    from ..quantity.ext.of_link import calc_midpoint_of_link

    y_of_node, x_of_node = np.asfarray(nodes[0]), np.asfarray(nodes[1])

    if midpoint_of_link is None:
        midpoint_of_link = np.empty((len(nodes_at_link), 2), dtype=float)

    calc_midpoint_of_link(nodes_at_link, x_of_node, y_of_node, midpoint_of_link)

    sorted_links = argsort_points_by_x_then_y(midpoint_of_link)
    nodes_at_link[:] = nodes_at_link[sorted_links]
    midpoint_of_link[:] = midpoint_of_link[sorted_links]

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
    >>> xy_of_link = np.array(
    ...     [[0.0, 0.5], [0.5, 1.0], [0.5, 0.5], [0.5, 0.0], [1.0, 0.5]]
    ... )
    >>> sort_patches(links_at_patch, offset_to_patch, xy_of_link)
    array([1, 0])
    >>> links_at_patch
    array([3, 2, 4, 0, 1, 2])
    >>> offset_to_patch
    array([0, 3, 6])
    """
    from .ext.remap_element import calc_center_of_patch
    from .ext.remap_element import reorder_patches

    n_patches = len(offset_to_patch) - 1
    xy_at_patch = np.empty((n_patches, 2), dtype=float)

    calc_center_of_patch(links_at_patch, offset_to_patch, xy_of_link, xy_at_patch)

    sorted_patches = argsort_points_by_x_then_y(xy_at_patch)
    reorder_patches(links_at_patch, offset_to_patch, sorted_patches)

    return sorted_patches


def sort_spokes_at_hub_on_graph(graph, spoke=None, at="node", inplace=False):
    """Order spokes of a graph clockwise around spokes.

    Parameters
    ----------
    graph : Graph-like
        A landlab graph.
    spoke : str
        Name of graph elements that make the spokes.
    at : {'node', 'corner', 'link', 'face', 'patch', 'cell'}
        Namve of graph elements that make the hubs.

    Returns
    -------
    ndarray of int, shape `(n_spokes, n_hubs)`
        Spokes ordered around each hub.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.graph import UniformRectilinearGraph, TriGraph
    >>> from landlab.graph.sort.sort import sort_spokes_at_hub_on_graph

    >>> graph = UniformRectilinearGraph((3, 3))
    >>> sort_spokes_at_hub_on_graph(graph, "link", at="node")
    array([[ 0,  2, -1, -1],
           [ 1,  3,  0, -1],
           [ 4,  1, -1, -1],
           [ 5,  7,  2, -1],
           [ 6,  8,  5,  3],
           [ 9,  6,  4, -1],
           [10,  7, -1, -1],
           [11, 10,  8, -1],
           [11,  9, -1, -1]])

    >>> graph = TriGraph((3, 3), node_layout="hex", sort=True)
    >>> sort_spokes_at_hub_on_graph(graph, "patch", at="node")
    array([[ 0,  2, -1, -1, -1, -1],
           [ 1,  3,  0, -1, -1, -1],
           [ 4,  1, -1, -1, -1, -1],
           [ 5,  2, -1, -1, -1, -1],
           [ 6,  8,  5,  2,  0,  3],
           [ 7,  9,  6,  3,  1,  4],
           [ 7,  4, -1, -1, -1, -1],
           [ 5,  8, -1, -1, -1, -1],
           [ 8,  6,  9, -1, -1, -1],
           [ 9,  7, -1, -1, -1, -1]])
    """
    sorted_spokes = argsort_spokes_at_hub_on_graph(graph, spoke=spoke, at=at)
    if spoke == "patch":
        plural = "patches"
    else:
        plural = spoke + "s"
    spokes_at_hub = getattr(graph, f"{plural}_at_{at}")

    if inplace:
        out = spokes_at_hub
    else:
        out = np.empty_like(spokes_at_hub)

    return np.take(spokes_at_hub, sorted_spokes, out=out)


def argsort_spokes_at_hub_on_graph(graph, spoke=None, at="node"):
    """Order spokes clockwise around spokes.

    Parameters
    ----------
    graph : Graph-like
        A landlab graph.
    spoke : str
        Name of graph elements that make the spokes.
    at : str
        Namve of graph elements that make the hubs.

    Returns
    -------
    ndarray of int, shape `(n_spokes, n_hubs)`
        Spokes ordered around each hub.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.graph import UniformRectilinearGraph
    >>> from landlab.graph.sort.sort import argsort_spokes_at_hub_on_graph

    >>> graph = UniformRectilinearGraph((3, 3))
    >>> argsort_spokes_at_hub_on_graph(graph, "link", at="node")
    array([[ 0,  1,  2,  3],
           [ 4,  5,  6,  7],
           [ 9, 10,  8, 11],
           [12, 13, 15, 14],
           [16, 17, 18, 19],
           [21, 22, 23, 20],
           [24, 27, 25, 26],
           [28, 30, 31, 29],
           [34, 35, 32, 33]])
    """
    angles = calc_angle_of_spoke_on_graph(graph, spoke=spoke, at=at, badval=np.inf)
    angles[angles < 0] += np.pi

    n_hubs, n_spokes = angles.shape
    ordered_angles = np.argsort(angles, kind="stable")
    ordered_angles += np.arange(n_hubs).reshape((-1, 1)) * n_spokes

    return as_id_array(ordered_angles)


def calc_angle_of_spoke_on_graph(graph, spoke=None, at="node", badval=None):
    """Calculate angles spokes make with a hub.

    Parameters
    ----------
    graph : Graph-like
        A landlab graph.
    spoke : str
        Name of graph elements that make the spokes.
    at : str
        Namve of graph elements that make the hubs.
    badval : float or iterable of float, optional
        Value to insert for missing spokes. If an iterable, use items as
        bad values to use for different spokes.

    Returns
    -------
    ndarray of float, shape `(n_spokes, n_hubs)`
        Angle that spoke elements make with each hub element.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.graph import UniformRectilinearGraph
    >>> from landlab.graph.sort.sort import calc_angle_of_spoke_on_graph

    >>> graph = UniformRectilinearGraph((3, 3))
    >>> np.rad2deg(
    ...     calc_angle_of_spoke_on_graph(graph, "link", at="node", badval=np.nan)
    ... )
    array([[   0.,   90.,   nan,   nan],
           [   0.,   90.,  180.,   nan],
           [  nan,   90.,  180.,   nan],
           [   0.,   90.,   nan,  270.],
           [   0.,   90.,  180.,  270.],
           [  nan,   90.,  180.,  270.],
           [   0.,   nan,   nan,  270.],
           [   0.,   nan,  180.,  270.],
           [  nan,   nan,  180.,  270.]])
    """
    xy_of_hub = getattr(graph, f"xy_of_{at}")
    xy_of_spoke = getattr(graph, f"xy_of_{spoke}")
    if spoke == "patch":
        plural = "patches"
    else:
        plural = spoke + "s"
    spokes_at_hub = getattr(graph, f"{plural}_at_{at}")

    angle_of_spoke = calc_angle_of_spoke(
        spokes_at_hub, xy_of_hub, xy_of_spoke, badval=badval
    )

    return angle_of_spoke


def sort_spokes_at_hub(spokes_at_hub, xy_of_hub, xy_of_spokes, inplace=False):
    if inplace:
        out = spokes_at_hub
    else:
        out = np.empty_like(spokes_at_hub)

    dx = np.subtract(xy_of_spokes[:, 0][spokes_at_hub], xy_of_hub[:, 0, None])
    dy = np.subtract(xy_of_spokes[:, 1][spokes_at_hub], xy_of_hub[:, 1, None])

    angle_of_spoke_at_hub = np.arctan2(dy, dx, where=spokes_at_hub != -1)
    angle_of_spoke_at_hub[angle_of_spoke_at_hub < 0.0] += np.pi * 2.0

    sort_id_array(spokes_at_hub, angle_of_spoke_at_hub, out)

    return out


def argsort_spokes_at_hub(spokes_at_hub, xy_of_hub, xy_of_spokes):
    angles = calc_angle_of_spoke(spokes_at_hub, xy_of_hub, xy_of_spokes, badval=np.inf)
    angles[angles < 0] += np.pi

    n_hubs, n_spokes = angles.shape
    ordered_angles = np.argsort(angles, kind="stable")
    ordered_angles += np.arange(n_hubs).reshape((-1, 1)) * n_spokes

    return as_id_array(ordered_angles)


def calc_angle_of_spoke(spokes_at_hub, xy_of_hub, xy_of_spoke, badval=None):
    xy_of_spoke = xy_of_spoke[spokes_at_hub.flat]
    x_of_spoke = xy_of_spoke[:, 0].reshape(spokes_at_hub.shape)
    y_of_spoke = xy_of_spoke[:, 1].reshape(spokes_at_hub.shape)

    x_of_hub, y_of_hub = xy_of_hub[:, 0], xy_of_hub[:, 1]

    dx = (x_of_spoke.T - x_of_hub).T
    dy = (y_of_spoke.T - y_of_hub).T

    angle_of_spoke = np.arctan2(dy, dx)
    angle_of_spoke[angle_of_spoke < 0.0] += np.pi * 2.0

    if badval is not None:
        try:
            badval_for_column = enumerate(badval)
        except TypeError:
            angle_of_spoke[spokes_at_hub == -1] = badval
        else:
            for col, val in badval_for_column:
                angle_of_spoke[spokes_at_hub[:, col] == -1, col] = val

    return angle_of_spoke
