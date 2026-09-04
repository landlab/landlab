from collections.abc import Sequence
from typing import Literal

import numpy as np
import pandas as pd
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
    and downstream distance (distance from the upstream end (tail) of the link to the
    farthest-downstream edge of the node's cell)


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

    nmg_link_to_rmg_coincident_nodes_mapper: pandas dataframe
        each row of the dataframe lists the link ID, the coincident node ID,
        the downstream distance of the coincident node, the x and y coordinates
        of the coincident node and the drainage area of the link

    """
    Lxy = []  # list of all nodes and node attributes that coincide with the
    # network model grid links
    # loop through all links in network model grid to determine raster grid cells
    # coincident with each link and equivalent distance from upstream (tail) node
    for linkID, lknd in enumerate(link_nodes):  # for each link in network grid

        x0 = nmgrid.x_of_node[lknd[0]]  # x and y of downstream link node
        y0 = nmgrid.y_of_node[lknd[0]]
        x1 = nmgrid.x_of_node[lknd[1]]  # x and y of upstream link node
        y1 = nmgrid.y_of_node[lknd[1]]

        # x and y coordinates and downstream distance from the upstream (tail)
        # node for 1000 points generated from downstream node to upstream node
        X, Y, dist = _link_to_points_and_dist((x0, y0), (x1, y1), number_of_points=1000)
        dist = dist.max() - dist  # convert to distance from tail node
        nodelist = []  # list of nodes along link
        for i, y in enumerate(Y):
            x = X[i]
            node = grid.find_nearest_node((x, y))
            # if node not already in list, append - many points will be in same cell;
            # only need to list cell once
            if node not in nodelist:
                nodelist.append(node)
                xy = {
                    "linkID": linkID,
                    "coincident_node": node,
                    "x": grid.node_x[node],
                    "y": grid.node_y[node],
                    "coincident_node_downstream_dist": dist[i],
                    "link_drainage_area": nmgrid.at_link["drainage_area"][linkID],
                }
                Lxy.append(xy)
    df = pd.DataFrame(Lxy)

    # if remove_duplicates, remove duplicate node id from link with smaller
    # contributing area.
    if remove_duplicates:
        values = df["coincident_node"].to_numpy()
        area = df["link_drainage_area"].to_numpy()
        idx = choose_unique(values=values, order_by=[area], choose="last")
        idx.sort()
        df = df.iloc[idx].reset_index(drop=True)

    return df


def _remove_small_tribs(
    rmg_nodes_to_nmg_links_mapper,
    nmg_link_to_rmg_coincident_nodes_mapper,
    remove_small_trib_factor,
):
    """remove rmg channel nodes that do not have an equivalent nmg link from the mapper
    by looking for the rmg channel nodes that represent first order channels that flow into
    a mainstem and much higher order channel"""

    for link in np.unique(nmg_link_to_rmg_coincident_nodes_mapper["linkID"].values):
        if link in rmg_nodes_to_nmg_links_mapper["linkID"].values:
            # first get the coincident rmg node that represents the inlet to the link
            # as the node with shortest downstream distance from inlet
            mask1 = nmg_link_to_rmg_coincident_nodes_mapper["linkID"] == link
            min_dist = nmg_link_to_rmg_coincident_nodes_mapper[
                "coincident_node_downstream_dist"
            ][mask1].min()
            mask2 = (
                nmg_link_to_rmg_coincident_nodes_mapper[
                    "coincident_node_downstream_dist"
                ]
                == min_dist
            )
            inlet_coincident_node = nmg_link_to_rmg_coincident_nodes_mapper[
                "coincident_node"
            ][mask1][mask2].iloc[0]
            # now get the contributing area of the rmg channel node mapped to the link
            # inlet (inlet_CA).
            # if a small tributary node is also mapped to the link inlet, there may be
            # more than one contributing area associated with the inlet
            mask3 = (
                rmg_nodes_to_nmg_links_mapper["coincident_node"]
                == inlet_coincident_node
            )
            inlet_CA_ = rmg_nodes_to_nmg_links_mapper["node_drainage_area"][
                mask3
            ].values
            if (
                len(inlet_CA_) > 1
            ):  # if there is more than one, remove the contributing area that is much less than the link contributing area
                # Where "much less" is defined as being less than the contributing area to the link divided by the factor "remove_small_trib_factor"
                mask4 = (
                    inlet_CA_
                    > rmg_nodes_to_nmg_links_mapper["link_drainage_area"].iloc[0]
                    / remove_small_trib_factor
                )
                inlet_CA = inlet_CA_[mask4]
                # if one or more areas are NOT much less than the contributing area to the link, pick the smallest
                if len(inlet_CA) >= 1:
                    inlet_CA = inlet_CA.min()
                # Or if all areas are much less than the the contributing area to the link, pick the smallest
                elif len(inlet_CA) == 0:
                    inlet_CA = inlet_CA_.min()
            else:
                inlet_CA = inlet_CA_.min()

            # Any nodes that have a contributing area less than the inlet_CA are removed
            mask5 = (rmg_nodes_to_nmg_links_mapper["linkID"] == link) & (
                rmg_nodes_to_nmg_links_mapper["node_drainage_area"] < inlet_CA
            )
            rmg_nodes_to_nmg_links_mapper = rmg_nodes_to_nmg_links_mapper.drop(
                rmg_nodes_to_nmg_links_mapper.index[mask5].values
            )
    return rmg_nodes_to_nmg_links_mapper


def map_rmg_nodes_to_nmg_links(
    grid,
    nmg_link_to_rmg_coincident_nodes_mapper,
    rmg_nodes,
    remove_small_trib_factor=None,
):
    """Map the nodes representing the channel location in a DEM to the closest
    network model grid location. Network model grid location is described in
    terms of link id and distance down link, measured from the inlet node (tail)
    of the link.

    Parameters
    ----------
    grid : raster model grid
        needs to have node field "drainage_area"
    nmg_link_to_rmg_coincident_nodes_mapper : pandas dataframe
        each row of the dataframe lists the link ID, the coincident node ID,
        the downstream distance of the coincident node, the x and y coordinates
        of the coincident node and the drainage area of the link
    rmg_nodes : np.array
        an array of node ids to be mapped to the nmg links
    remove_small_tribs : None or int
        If int, channel nodes whose contributing area is much less than the contributing
        area of the closest link are not matched to the link. Where "much less"
        is defined as being less than the contributing area to the link divided by the factor
        "remove_small_trib_factor". If None, then the outlet node of small tributaries
        could be mapped to larger, main stem reaches of the channel network model.
        Default is remove_small_trib_factor = None

    Returns
    -------
    rmg_nodes_to_nmg_links_mapper : pandas dataframe
        each row of the dataframe lists the node ID, the link ID the node has been
        mapped too, the closest nmg-link-coincident node ID, the drainage area
        of the link and the drainage area of the node

    """

    def dist_between_nmg_and_rmg_nodes(row, xc, yc):
        """distance between channel node and link node"""
        return _dist_func(xc, row["x"], yc, row["y"])

    link_ = []
    for n in rmg_nodes:  # for each rmg node
        xc = grid.node_x[n]
        yc = grid.node_y[n]
        # compute the distance to all link coincident rmg nodes
        dist = nmg_link_to_rmg_coincident_nodes_mapper.apply(
            lambda row: dist_between_nmg_and_rmg_nodes(row, xc, yc), axis=1
        )
        # pick closest coincident node and corresponding link
        # if more than one (which can happen because the confluence between two
        # links overlay the same node), pick link with largest contributing area
        mask = dist == dist.min()
        dist_min_links = nmg_link_to_rmg_coincident_nodes_mapper[
            [
                "linkID",
                "coincident_node_downstream_dist",
                "coincident_node",
                "link_drainage_area",
            ]
        ][mask]
        link = dist_min_links[
            dist_min_links["link_drainage_area"]
            == dist_min_links["link_drainage_area"].max()
        ].head(1)
        link["node_drainage_area"] = grid.at_node["drainage_area"][
            n
        ]  # add node drainage area to attributes
        link_.append(link)

    rmg_nodes_to_nmg_links_mapper = pd.concat(link_)
    rmg_nodes_to_nmg_links_mapper["node"] = rmg_nodes
    # organize column order in mapper
    rmg_nodes_to_nmg_links_mapper = rmg_nodes_to_nmg_links_mapper[
        [
            "node",
            "linkID",
            "coincident_node",
            "coincident_node_downstream_dist",
            "link_drainage_area",
            "node_drainage_area",
        ]
    ].reset_index(drop=True)

    if (
        remove_small_trib_factor
    ):  # check for small tributary nodes assigned to link and remove them
        rmg_nodes_to_nmg_links_mapper = _remove_small_tribs(
            rmg_nodes_to_nmg_links_mapper,
            nmg_link_to_rmg_coincident_nodes_mapper,
            remove_small_trib_factor,
        )

    return rmg_nodes_to_nmg_links_mapper
