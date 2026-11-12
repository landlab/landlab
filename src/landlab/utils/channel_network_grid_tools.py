<<<<<<< HEAD
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from landlab.components.flow_director.flow_director_steepest import FlowDirectorSteepest


def _dist_func(x0, x1, y0, y1):
    return np.hypot(x0 - x1, y0 - y1)
=======
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
>>>>>>> master


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


<<<<<<< HEAD
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


def plot_nmgrids(nmgrid_1, nmgrid_2):
    """show the links and link ids of two network model grids in one plot"""

    def plot_nmgrid(nmgrid, line_color, alpha, fontsize, label):
        xnode = nmgrid.x_of_node
        xlink = nmgrid.midpoint_of_link[:, 0]
        ynode = nmgrid.y_of_node
        ylink = nmgrid.midpoint_of_link[:, 1]
        for link, val in enumerate(nmgrid.nodes_at_link):
            xv = xnode[val]
            yv = ynode[val]
            if link == 0:
                plt.plot(xv, yv, color=line_color, alpha=alpha, label=label)
            else:
                plt.plot(xv, yv, color=line_color, alpha=alpha, label="_nolegend_")
            plt.text(
                xlink[link],
                ylink[link],
                str(link),
                size=fontsize,
                color=line_color,
                alpha=alpha,
            )

    plt.figure(figsize=(5, 5))
    plot_nmgrid(nmgrid_1, line_color="red", alpha=1, fontsize=12, label="nmgrid_1")
    plot_nmgrid(nmgrid_2, line_color="green", alpha=0.37, fontsize=20, label="nmgrid_2")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()


def create_df_of_link_points(nmgrid, link_nodes, number_of_points):
    """convert the network model grid to a point representation, with each link
    of the grid represented by a series of number_of_points points. Each point
    is described by x and y coordinates and the link that the point represents.
    Note that this function differs from map_nmg_links_to_rmg_coincident_nodes,
    which converts the network model grid to it's node representation.
=======
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
>>>>>>> master


    Parameters
    ----------
<<<<<<< HEAD
    nmgrid : network model grid
    link_nodes : np.array
        for a nmgrid of n nodes, a nx2 np array, the ith row of the
        array is the [downstream node id, upstream node id] of the ith link
    number_of_points : int
        Each link is converted to a series of number_of_point points.

    Returns
    -------
    pandas dataframe that lists the link ID and x and y coordinates of each point

    """
    X_ = np.array([])
    Y_ = np.array([])
    link_ = np.array([])
    for linkID, lknd in enumerate(link_nodes):  # for each link in nmgrid 1 1
=======
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
    nodes. Each coincident raster model grid node is defined in terms of its
    x and y coordinates, the link it is mapped to and distance downstream from
    the upstream (tail) end of the link.


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
        each row of the dataframe lists the link ID, the coincident node ID, the
        x and y coordinates and the downstream distance of the coincident node
        and the drainage area of the link

    """
    Lxy = []  # list of all nodes and node attributes that coincide with the
    # network model grid links
    # loop through all links in network model grid to determine raster grid cells
    # coincident with each link and equivalent distance from upstream (tail) node
    for linkID, lknd in enumerate(link_nodes):  # for each link in network grid
>>>>>>> master

        x0 = nmgrid.x_of_node[lknd[0]]  # x and y of downstream link node
        y0 = nmgrid.y_of_node[lknd[0]]
        x1 = nmgrid.x_of_node[lknd[1]]  # x and y of upstream link node
        y1 = nmgrid.y_of_node[lknd[1]]

<<<<<<< HEAD
        # convert link to a series of points
        X, Y, dist = _link_to_points_and_dist((x0, y0), (x1, y1), number_of_points)

        X_ = np.concatenate((X_, X))
        Y_ = np.concatenate((Y_, Y))
        link_ = np.concatenate((link_, (np.ones(len(X)) * linkID).astype(int)))

    return pd.DataFrame(data=zip(link_, X_, Y_), columns=["linkID", "X", "Y"])


def map_nmg1_links_to_nmg2_links(
    nmgrid_1, nmgrid_2, number_of_points=11, plot_grids=False
):
    """given two slightly different network model grids of the same channel network,
    map each link from one network model grid (nmgrid_1) to the closest (based on
    the mean distance between links) link of the other network model grid (nmgrid_2).
    If two or more links of nmgrid_2 are equally close to a link of nmgrid_1, the
    link with the largest drainage area is mapped to the nmgrid_1 link

    Parameters
    ----------
    nmgrid_1 : network model grid
        grid that values will be mapped to
    nmgrid_2 : network model grid
        grid that values will be mapped from
    number_of_points : int
        Each link is converted to a series of number_of_point points. The relative
        distance of each link to the other is determined using these points.
        The default is 11. Below 11, mapping may not match expected.

    Returns
    -------
    link_mapper : dict
        Keys are the id of all links in nmgrid_1. Values are the link IDs of nmgrid_2
        that are mapped to each nmgrid_1 link.


    WARNING: In some situations this function may not map as expected. Set plot_grids
    to True and inspect results
    """

    warnings.warn(
        "In some situations this function may not map as expected. Set plot_grids to True and inspect results"
    )

    def distance_between_links(row, XY):
        return _dist_func(
            row["X"], XY[0], row["Y"], XY[1]
        )  # ((row['x']-XY[0])**2+(row['y']-XY[1])**2)**.5

    # convert the network model grid to a point representation, as described by
    # the link ID, x and y value of each point
    nmgrid_1_link_points = create_df_of_link_points(
        nmgrid_1, nmgrid_1.nodes_at_link, number_of_points
    )
    nmg1_linkIDs = nmgrid_1_link_points["linkID"].astype(int).values

    nmgrid_2_link_points = create_df_of_link_points(
        nmgrid_2, nmgrid_2.nodes_at_link, number_of_points
    )
    nmg2_linkIDs = nmgrid_2_link_points["linkID"].astype(int).values
    # for each point of each link of nmgrid_1, find the closest nmgrid_2 point
    # and link. nmgrid_2 link with highest number of points closest to the
    # nmgrid_1 link is mapped to the nmgrid_1 link.

    sublist1 = nmgrid_1_link_points[["X", "Y"]]  # get points that represent nmgrid_1
    sublist2 = nmgrid_2_link_points[["X", "Y"]]  # get points that represent nmgrid_2
    distance_matrix = cdist(
        sublist1, sublist2, metric="euclidean"
    )  # create the distance matrix, which lists the distance between all nmgrid_1 and nmgrid_2 points
    distance_matrix_nodiag = distance_matrix  # fill the diagonal values with inf
    np.fill_diagonal(distance_matrix_nodiag, np.inf)
    closest_point_indices = np.argmin(
        distance_matrix_nodiag, axis=1
    )  # find the minimum values
    linkID_array = np.tile(
        nmg2_linkIDs, (len(nmg1_linkIDs), 1)
    )  # create a matrix of the nmg 2 link ids
    nmg2_link_matrix = linkID_array[
        np.arange(len(nmg1_linkIDs)), closest_point_indices
    ]  # get the link id of the closest node

    # now count the number of times each nmgrid_2 point was closest to nmgrid_1 link
    link_mapper = {}
    for linkID_1 in nmg1_linkIDs:
        linkIDs_2 = nmg2_link_matrix[nmg1_linkIDs == linkID_1]
        count = np.bincount(linkIDs_2)
        # nmgrid_2 link with highest count is matched to nmgrid_1 link
        # if only one nmgrid_2 link has highest count, that is the link
        if (count == count.max()).sum() == 1:
            linkID_2 = np.argmax(count)
        else:  # if two or more nmgrid_2 links have the hightest count, select the
            # one that drains the largest area
            links_with_same_count = np.arange(len(count))[count == count.max()]
            DAs_ = nmgrid_2.at_link["drainage_area"][links_with_same_count]
            linkID_2 = links_with_same_count[DAs_ == DAs_.max()][0]  # to remove bracket
        link_mapper[linkID_1] = linkID_2

    if plot_grids:
        plot_nmgrids(nmgrid_1, nmgrid_2)

    return link_mapper
=======
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
                    "dist": dist[i],
                    "drainage_area": nmgrid.at_link["drainage_area"][linkID],
                }
                Lxy.append(xy)
    df = pd.DataFrame(Lxy)

    # if remove_duplicates, remove duplicate node id from link with smaller
    # contributing area.
    if remove_duplicates:
        values = df["coincident_node"].to_numpy()
        area = df["drainage_area"].to_numpy()
        idx = choose_unique(values=values, order_by=[area], choose="last")
        idx.sort()
        df = df.iloc[idx].reset_index(drop=True)

    return df
>>>>>>> master
