import numpy as np
import pandas as pd
from landlab.components.flow_director.flow_director_steepest import FlowDirectorSteepest


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


def map_nmg1_links_to_nmg2_links(nmgrid_1, nmgrid_2, number_of_points=11):
    """given two slightly different network model grids of the same channel network,
    map links from one network model grid (nmgrid_1) to the closest links of the
    other network model grid (nmgrid_2). If two or more links of nmgrid_2 are equally
    close to a link of nmgrid_1, the link with the largest drainage area is mapped
    to the nmgrid_1 link.


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
    """

    def distance_between_links(row, XY):
        return _dist_func(
            row["X"], XY[0], row["Y"], XY[1]
        )  # ((row['x']-XY[0])**2+(row['y']-XY[1])**2)**.5

    # get the head and tail nodes of each link
    linknodes_1 = get_link_nodes(nmgrid_1)
    linknodes_2 = get_link_nodes(nmgrid_2)

    # convert the network model grid to a point representation, as described by
    # the link ID, x and y value of each point
    nmgrid_1_link_points = create_df_of_link_points(
        nmgrid_1, linknodes_1, number_of_points
    )
    nmgrid_2_link_points = create_df_of_link_points(
        nmgrid_2, linknodes_2, number_of_points
    )

    # for each point of each link of nmgrid_1, find the closest nmgrid_2 point
    # and link. nmgrid_2 link with highest number of points mapped to nmgrid_1
    # link is mapped to the nmgrid_1 link.
    link_mapper = {}
    for linkID, lknd in enumerate(linknodes_1):  # for each link in nmgrid 1

        sublist = nmgrid_1_link_points[["X", "Y"]][
            nmgrid_1_link_points["linkID"] == linkID
        ]
        LinkL = []  # id of nmg2 link that is closest to nmg1 rmg node
        for j in range(len(sublist)):
            XY = [sublist.iloc[j]["X"], sublist.iloc[j]["Y"]]
            distances = nmgrid_2_link_points.apply(
                lambda row: distance_between_links(row, XY), axis=1
            )  # compute the distance from the nmgrid_1 point and all nmgrid_2 points
            offset = (
                distances.min()
            )  # find the minimum distance between the nmg1 point and all nmg2 points
            mdl = (
                nmgrid_2_link_points["linkID"][(distances == offset)]
                .values[0]
                .astype(int)
            )  # get the nmg2 link id with point at minimum distance from nmg1 point, if more than one, pick the first one
            LinkL.append(mdl)
        Links = np.array(LinkL)
        # number of times each nmgrid_2 point was closest to nmgrid_1 link
        count = np.bincount(Links)

        # nmgrid_2 link with highest count is matched to nmg1 link
        # if only one nmgrid_2 link has highest count, that is the link
        if (count == count.max()).sum() == 1:
            Link = np.argmax(count)
        else:  # if two or more nmgrid_2 links have the hightest count, select the
            # one that drains the largest area
            links_with_same_count = np.arange(len(count))[count == count.max()]
            DAs_ = nmgrid_2.at_link["drainage_area"][links_with_same_count]
            Link = links_with_same_count[DAs_ == DAs_.max()][0]  # to remove bracket
        link_mapper[linkID] = Link

    return link_mapper
