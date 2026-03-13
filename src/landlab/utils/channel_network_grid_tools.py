from collections.abc import Sequence
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike
from numpy.typing import NDArray

from landlab.components.flow_director.flow_director_steepest import FlowDirectorSteepest

"""
A collection of tools for mapping values (e.g., flow, shear stress) between
network model grid and raster model grid representations of a channel network.
"""


def get_link_nodes(nmgrid):
    """Get the downstream (head) and upstream (tail) nodes at a link from
    flow director. The network model grid nodes_at_link attribute may not be
    ordered according to flow direction. Output from this function should be
    used for all channel_network_grid_tools functions that require a link_nodes
    input

    Parameters
    ----------
    nmgrid : network model grid

    Returns
    -------
    link_nodes : np.array
        for a nmgrid of n nodes, returns a nx2 np array, the ith row of the
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
):
    """Given two points defined by coordinates x0,y0 and x1,y1, define a series
    of points between them and the distance from point x0,y0 to each point.

    Parameters
    ----------
    point_0 : tuple of 2 floats
        point 0 coordinates x and y
    point_1 : tuple of 2 floats
        point 1 coordinates x and y
    number_of_points : int
        number of points to create along the reach. The default is 1000.

    Returns
    -------
    X : np array
        x coordinate of points
    Y : np array
        y coordinate of points
    dist : np array
        linear distance between points

    """
    x0 = point_0[0]
    y0 = point_0[1]
    x1 = point_1[0]
    y1 = point_1[1]
    X = np.linspace(x0, x1, number_of_points)
    Y = np.linspace(y0, y1, number_of_points)
    dist = np.hypot(X - x0, Y - y0)

    return X, Y, dist


def _dist_func(x0, x1, y0, y1):
    return np.hypot(x0 - x1, y0 - y1)


def extract_channel_nodes(grid, Ct):
    """interpret which nodes of the DEM represent the channel network as all nodes
    that have a drainage area >= to the average drainage area at which
    channels initiate in the DEM (Ct, based on field or remote sensing evidence).

    Use Ct = average drainage area at which colluvial channels to get the entire
    channel network.

    Use Ct = the drainage area at which cascade channels typically begin to get
    a channel network where sediment transport is primarily via fluvial processes

    Parameters
    ----------
    grid : raster model grid
        raster model grid with node field "drainage_area"
    Ct : float
        Channel threshold drainage area

    Returns
    -------
    cn : np array of int
         array of all node ids included in the channel network

    """
    return np.flatnonzero(grid.at_node["drainage_area"] >= Ct)


def extract_terrace_nodes(grid, terrace_width, acn, fcn):
    """Determine which raster model grid nodes coincide with channel terraces,
    which presently are assumed to be a fixed width (number of nodes) from
    the channel nodes


    Parameters
    ----------
    grid : raster model grid
    terrace_width : int
        Width of terrace in number of nodes. If provided as float, will be rounded
        to nearest int.
    acn : np array
        array of all node IDs included in the channel network
    fcn : np array
        array of all node IDs included in the fluvial channel network

    Raises
    ------
    ValueError
        Occurs if terrace width less than 1.

    Returns
    -------
    terrace_nodes : np array
        array of all node IDs included in the terrace

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


def min_distance_to_network(grid, acn, node_id):
    """Determine the shortest distance (as the crow flies) from a node to the
    channel network and the closest channel node

    Parameters
    ----------
    grid : raster model grid
    acn : list of int
        array of all node ids included in the channel network
    node_id : int
        ID of node from which the distance will be determined

    Returns
    -------
    offset : float
        distance between node and channel network
    mdn : int
        ID of channel node that is closest node

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


def map_nmg_links_to_rmg_coincident_nodes(
    grid, nmgrid, link_nodes, remove_duplicates=False
):
    """Map links of a network model grid to all coincident raster model grid
    nodes (nodes whose associated cell intersect the link). Each coincident raster model
    grid node is defined in terms of its x and y coordinates, the link it is mapped to
    and the downstream distance of the node on the link. The downstream distance
    of the node on the link is defined as the distance from the upstream end (tail) of
    the link to the farthest-downstream edge of the node's cell.


    Parameters
    ----------
    grid : raster model grid
    nmgrid : network model grid
    link_nodes : np array
        head and tail node of each link
    remove_duplicates : bool
        if True, when two or more links are coincident with the same node,
        the node is assigned to the link with the largest drainage area. If False,
        the node is assigned to each coincident link. The default is False.

    Returns
    -------

    nmg_link_to_rmg_coincident_nodes_mapper: dict
        each key of the dictionary contains an np.array whose length is equal to the
        number of coincident nodes. Keys include link ID, coincident node ID,
        downstream distance of the coincident node, x coordinate of the coicident
        node, y coordinate of the coincident node and drainage area of the link.

    """
    linkIDs_list = []
    nodes_list = []
    Xs_list = []
    Ys_list = []
    dists_list = []
    link_drainage_areas_list = []
    for linkID, lknd in enumerate(link_nodes):
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
        dists = dists.max() - dists  # convert to distance from tail node
        nodes = grid.find_nearest_node((Xs, Ys))
        # using the x and y coordinate of the first (most downstream) point
        # within the node's cell to represent the node location on the link
        mask = choose_from_repeated(nodes, choose="first")
        nodes = nodes[mask]

        linkIDs_list.append((np.ones(len(nodes)) * linkID).astype(int))
        nodes_list.append(nodes)
        Xs_list.append(grid.node_x[nodes])
        Ys_list.append(grid.node_y[nodes])
        dists_list.append(dists[mask])
        link_drainage_areas_list.append(
            (np.ones(len(nodes)) * nmgrid.at_link["drainage_area"][linkID]).astype(
                float
            )
        )

    nmg_link_to_rmg_coincident_nodes_mapper = {
        "linkID": np.concatenate(linkIDs_list),
        "coincident_node": np.concatenate(nodes_list),
        "x": np.concatenate(Xs_list),
        "y": np.concatenate(Ys_list),
        "coincident_node_downstream_dist": np.concatenate(dists_list),
        "link_drainage_area": np.concatenate(link_drainage_areas_list),
    }

    if remove_duplicates:
        values = nmg_link_to_rmg_coincident_nodes_mapper["coincident_node"]
        area = nmg_link_to_rmg_coincident_nodes_mapper["link_drainage_area"]
        idx = choose_unique(values=values, order_by=[area], choose="last")
        idx.sort()
        for key in nmg_link_to_rmg_coincident_nodes_mapper.keys():

            nmg_link_to_rmg_coincident_nodes_mapper[key] = (
                nmg_link_to_rmg_coincident_nodes_mapper[key][idx]
            )

    return nmg_link_to_rmg_coincident_nodes_mapper
