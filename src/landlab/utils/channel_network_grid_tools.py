from collections.abc import Sequence
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike
from numpy.typing import NDArray

from landlab import NetworkModelGrid
from landlab import RasterModelGrid
from landlab.components.flow_director.flow_director_steepest import FlowDirectorSteepest
from landlab.core.utils import require_id_array
from landlab.utils.geometry.planar import find_nearest_node

"""
A collection of tools for defining a channel network on a cellular-like ModelGrid
(e.g., RasterModelGrid, HexModelGrid) and mapping values (e.g., flow, shear stress)
between the cellular-like ModelGrid and NetworkModelGrid representations of the network.
"""


def get_link_nodes(nmgrid: NetworkModelGrid) -> NDArray[np.integer]:
    """Get the downstream (head) and upstream (tail) NetworkModelGrid link nodes.
    The nodes listed in the NetworkModelGrid nodes_at_link attribute
    may not be ordered as [head node, tail node]. Output from this function
    should be used for all channel_network_grid_tools functions that require a
    link_nodes input.

    Parameters
    ----------
    nmgrid : NetworkModelGrid

    Returns
    -------
    link_nodes : array_like
        For a nmgrid of L links, returns a Lx2 array_like, the ith row of the
        array is the [downstream node id, upstream node id] of the ith link
    """

    fd = FlowDirectorSteepest(nmgrid, "topographic__elevation")
    fd.run_one_step()

    return np.column_stack(
        (fd.downstream_node_at_link(), fd.upstream_node_at_link())
    ).astype(int, copy=False)


def _link_to_points_and_dist(
    point_0: tuple[float, float],
    point_1: tuple[float, float],
    number_of_points: int = 1000,
) -> tuple[float, float, float]:
    """Create a series of points between two points.
    Given two points defined by coordinates x0,y0 and x1,y1, define a series
    of points between them and the distance from point x0,y0 to each point.

    Parameters
    ----------
    point_0 : tuple of 2 floats
        Point 0 coordinates x and y
    point_1 : tuple of 2 floats
        Point 1 coordinates x and y
    number_of_points : int
        Number of points to create along the reach. The default is 1000.

    Returns
    -------
    X : array_like
        X coordinate of points
    Y : array_like
        Y coordinate of points
    dist : array_like
        Linear distance between points

    """
    x0 = point_0[0]
    y0 = point_0[1]
    x1 = point_1[0]
    y1 = point_1[1]
    X = np.linspace(x0, x1, number_of_points)
    Y = np.linspace(y0, y1, number_of_points)
    dist = np.hypot(X - x0, Y - y0)

    return X, Y, dist


def _dist_func(x0: float, x1: float, y0: float, y1: float) -> float:
    return np.hypot(x0 - x1, y0 - y1)


def extract_channel_nodes(grid: ModelGrid, Ct: float) -> NDArray[np.integer]:
    """Extract the channel nodes from a cellular-type ModelGrid DEM representation.
    Interpret which nodes of the DEM on a cellular-type ModelGrid represent
    the channel network. The channel network is all nodes that have a drainage
    area greater than or equal to the average drainage area at which channels
    initiate in the DEM (Ct, based on field or remote sensing evidence).

    Use Ct = average drainage area at which colluvial channels begin to get the
    entire channel network.

    Use Ct = the drainage area at which cascade channels typically begin to get
    the portion of the channel network where sediment transport is primarily via
    fluvial processes.

    Parameters
    ----------
    grid : ModelGrid
        A cellular-type ModelGrid with node field "drainage_area"
    Ct : float
        Channel threshold drainage area

    Returns
    -------
    cn : array_like of int
         Array of all node ids included in the channel network.

    """
    return np.flatnonzero(grid.at_node["drainage_area"] >= Ct)


def extract_terrace_nodes(
    grid: RasterModelGrid,
    terrace_width: int,
    acn: NDArray[np.integer],
    fcn: NDArray[np.integer],
) -> NDArray[np.integer]:
    """Determine which RasterModelGrid nodes coincide with a channel terrace.
    This function is specific to the RasterModelGrid. Terrarce nodes are assumed
    to be a fixed width (number of nodes) from the channel nodes.


    Parameters
    ----------
    grid : RasterModelGrid
    terrace_width : int
        Width of terrace in number of nodes. If provided as float, will be rounded
        to nearest int.
    acn : array_like
        Array of all node IDs included in the channel network.
    fcn : array_like
        Array of all node IDs included in the fluvial channel network.

    Raises
    ------
    ValueError
        Occurs if terrace width is less than 1.

    Returns
    -------
    terrace_nodes : array_like
        Array of all node IDs included in the terrace.

    """
    # round to int in case provided as float
    terrace_width = round(terrace_width)
    if terrace_width < 1:
        raise ValueError(f"terrace width must be 1 or greater ({terrace_width})")

    acn = np.asarray(acn, dtype=int)
    current_nodes = np.asarray(fcn, dtype=int)
    terrace_nodes = np.array([], dtype=int)

    for _ in range(terrace_width):
        adj_dn = grid.diagonal_adjacent_nodes_at_node[current_nodes].ravel()
        adj_n = grid.adjacent_nodes_at_node[current_nodes].ravel()

        neighbors = np.unique(np.concatenate((adj_n, adj_dn)))
        neighbors = neighbors[neighbors != -1]

        terrace_nodes = np.setdiff1d(neighbors, acn, assume_unique=True)

        current_nodes = terrace_nodes

    return terrace_nodes


def min_distance_to_network(
    grid: ModelGrid, acn: NDArray[np.integer], node_id: int
) -> tuple[float, int]:
    """Find the shortest distance (as the crow flies) to the channel network.
    Measured from a node of a cellular-type ModelGrid to the cellular-type
    ModelGrid channel nodes. Returns the distance and the closest channel node.

    Parameters
    ----------
    grid : cellular-type ModelGrid
    acn : list of int
        Array of all node ids included in the channel network.
    node_id : int
        ID of node from which the distance will be determined.

    Returns
    -------
    offset : float
        Distance between node and channel network.
    mdn : int
        ID of channel node that is closest node.

    """
    x0, y0 = grid.node_x[node_id], grid.node_y[node_id]
    x_acn, y_acn = grid.node_x[acn], grid.node_y[acn]

    dist = np.hypot(x_acn - x0, y_acn - y0)

    idx = np.argmin(dist)
    offset = dist[idx]
    mdn = acn[idx]

    return float(offset), int(mdn)


def choose_from_repeated(
    sorted_array: ArrayLike,
    choose: Literal["first", "last"] = "last",
) -> NDArray[np.bool_]:
    """Mark the first/last element of repeated values in a **sorted** 1-D array.

    Parameters
    ----------
    sorted_array : array_like
        Assumed sorted by the grouping key.
    choose : {'first','last'}, optional
        Whether to mark the first or last item of each run.

    Examples
    --------
    >>> array = [0, 0, 0, 2, 2, 5, 6, 6, 6, 6, 6]
    >>> is_last = choose_from_repeated(array, choose="last")
    >>> is_last.astype(int)
    array([0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1])
    """
    a = np.asarray(sorted_array).ravel()

    same_as_previous = np.zeros(a.size, dtype=bool)

    if a.size <= 1:
        return np.ones(a.size, dtype=bool)

    same_as_previous[1:] = a[1:] == a[:-1]
    if choose == "last":
        keep_mask = np.ones_like(same_as_previous)
        keep_mask[:-1] = ~same_as_previous[1:]
    elif choose == "first":
        keep_mask = ~same_as_previous
    else:
        raise ValueError(f"choose must be 'first' or 'last', got {choose!r}")

    return keep_mask


def choose_unique(
    values: ArrayLike,
    order_by: Sequence[ArrayLike] | None = None,
    choose: Literal["first", "last"] = "last",
) -> NDArray[np.intp]:
    """Find indices of unique values, selecting one representative if repeated.

    Examples
    --------
    >>> choose_unique([0, 1, 0, 0, 1], order_by=([10.0, 11.0, 12.0, 13.0, 14],))
    array([3, 4])

    >>> choose_unique([1, 0, 0, 1, 0], order_by=([10.0, 11.0, 12.0, 13.0, 14],))
    array([3, 4])
    """
    values = np.asarray(values).ravel()

    order_by = (
        () if order_by is None else tuple(np.asarray(key).ravel() for key in order_by)
    )

    if any(key.size != values.size for key in order_by):
        raise ValueError("All `order_by` arrays must match `values` length")

    sorted_rows = np.lexsort(order_by + (values,))

    is_last = choose_from_repeated(values[sorted_rows], choose=choose)

    return np.sort(sorted_rows[is_last])


def map_network_links_to_coincident_nodes(
    grid: ModelGrid,
    nmgrid: NetworkModelGrid,
    link_nodes: ArrayLike,
    remove_duplicates: bool = False,
) -> dict[str, NDArray]:
    """Map the links of a NetworkModelGrid to the nodes of a cellular-type ModelGrid.
    This function finds each cellular-type ModelGrid node that is coincident with
    a NetworkModelGrid link (nodes whose associated cell intersects the link).
    Each coincident cellular-type ModelGrid node is then recorded in a mapper dictionary
    (nmg_link_to_mg_coincident_nodes_mapper) in terms of its x and y coordinates,
    the link it is mapped to, and the downstream distance of the node on the link.
    The downstream distance of the node on the link is defined as the distance
    from the upstream end (tail) of the link to the first (most downstream) point
    within the node's cell.


    Parameters
    ----------
    grid : cellular-type ModelGrid
    nmgrid : NetworkModelGrid
    link_nodes : array_like
        Head and tail node of each link
    remove_duplicates : bool, optional
        If True, when two or more links are coincident with the same node, which
        can occur at stream junctions, the node is assigned to the link with the
        largest drainage area. If False, the node is assigned to each coincident
        link. The default is False.

    Returns
    -------
    nmg_link_to_coincident_nodes_mapper: dict
        Each key of the dictionary contains an array_like whose length is equal to the
        number of coincident nodes. Keys include link ID, coincident node ID,
        downstream distance of the coincident node, x coordinate of the coincident
        node, y coordinate of the coincident node and drainage area of the link.

    """

    # validate that link_nodes is correct format
    require_id_array(
        link_nodes,
        shape=("n_links", 2),
        max_id=nmgrid.number_of_nodes - 1,
        bad_id=None,
        name="link_nodes",
    )

    # for each link in the network model grid, map nodes of the other grid to
    # the link
    link_ids_list = []
    nodes_list = []
    xs_list = []
    ys_list = []
    downstream_dists_list = []
    link_drainage_areas_list = []
    for link_id, lknd in enumerate(link_nodes):
        # x and y of downstream (head) node of link
        x0 = nmgrid.x_of_node[lknd[0]]
        y0 = nmgrid.y_of_node[lknd[0]]
        # x and y of upstream (tail) node of link
        x1 = nmgrid.x_of_node[lknd[1]]
        y1 = nmgrid.y_of_node[lknd[1]]

        # get x and y coordinates and downstream distance from the upstream
        # node for 1000 points generated from downstream node to upstream node
        Xs, Ys, dists = _link_to_points_and_dist(
            (x0, y0), (x1, y1), number_of_points=1000
        )
        downstream_dists = dists.max() - dists  # convert to distance from tail node
        # find the node closest to each of the points
        nodes = find_nearest_node(
            np.array([grid.node_x, grid.node_y]).T, np.array([Xs, Ys]).T
        )
        # using the x and y coordinates of the first (most downstream) point
        # within the node's cell to represent the node location on the link
        mask = choose_from_repeated(nodes, choose="first")
        nodes = nodes[mask]

        link_ids_list.append((np.ones(len(nodes)) * link_id).astype(int))
        nodes_list.append(nodes)
        xs_list.append(grid.node_x[nodes])
        ys_list.append(grid.node_y[nodes])
        downstream_dists_list.append(downstream_dists[mask])
        link_drainage_areas_list.append(
            np.full(
                nodes.size,
                nmgrid.at_link["drainage_area"][link_id],
                dtype=float,
            )
        )

    nmg_link_to_coincident_nodes_mapper = {
        "link_id": np.concatenate(link_ids_list),
        "coincident_node": np.concatenate(nodes_list),
        "x": np.concatenate(xs_list),
        "y": np.concatenate(ys_list),
        "coincident_node_downstream_dist": np.concatenate(downstream_dists_list),
        "link_drainage_area": np.concatenate(link_drainage_areas_list),
    }

    if remove_duplicates:
        values = nmg_link_to_coincident_nodes_mapper["coincident_node"]
        area = nmg_link_to_coincident_nodes_mapper["link_drainage_area"]
        idx = choose_unique(values=values, order_by=[area], choose="last")
        idx.sort()
        for key in nmg_link_to_coincident_nodes_mapper.keys():

            nmg_link_to_coincident_nodes_mapper[key] = (
                nmg_link_to_coincident_nodes_mapper[key][idx]
            )

    return nmg_link_to_coincident_nodes_mapper
