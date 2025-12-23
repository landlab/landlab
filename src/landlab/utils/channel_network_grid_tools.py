import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
from landlab.components.flow_director.flow_director_steepest import FlowDirectorSteepest
import warnings


def _dist_func(x0, x1, y0, y1):
    return np.hypot(x0 - x1, y0 - y1)


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


    Parameters
    ----------
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

        x0 = nmgrid.x_of_node[lknd[0]]  # x and y of downstream link node
        y0 = nmgrid.y_of_node[lknd[0]]
        x1 = nmgrid.x_of_node[lknd[1]]  # x and y of upstream link node
        y1 = nmgrid.y_of_node[lknd[1]]

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
