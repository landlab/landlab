from landlab.components.flow_director.flow_director_steepest import FlowDirectorSteepest
import numpy as np
import pandas as pd

"""
A collection of tools for mapping values (e.g., flow, shear stress) between a
network model grid and raster model grid representation of a channel network.
"""


def _flatten_lol(lol):
    """
    for list (l) in list of lists (lol) for item (i) in l append i"

    Parameters
    ----------
    lol : list of lists

    Returns
    -------
    the list of lists concatenated into a single list
    """
    return [i for l in lol for i in l]


def get_link_nodes(nmgrid):
    """get the downstream (head) and upstream (tail) nodes at a link from
    flow director. The network model grid nodes_at_link attribute may not be
    ordered according to flow direction. Output from this function should be
    used for all ChannelNetworkGridTool functions that require a link_nodes input

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
    upstream_node_id = []
    downstream_node_id = []
    for i in range(nmgrid.number_of_links):
        upstream_node_id.append(fd.upstream_node_at_link()[i])
        downstream_node_id.append(fd.downstream_node_at_link()[i])
    # create np array and transpose, each row is [downstream node id, upstream node id]
    link_nodes = np.array([downstream_node_id, upstream_node_id]).T
    return link_nodes


def _link_to_points_and_dist(x0, y0, x1, y1, number_of_points=1000):
    """Given two points defined by coordinates x0,y0 and x1,y1, create a series
    of points between them.
    If x0,y0 and x1,y1 are the head and tail nodes of a link, the distance is the
    downstream distance along the link.

    Parameters
    ----------
    x0 : float
        point 1 coordinate x
    y0 : float
        point 1 coordinate y
    x1 : float
        point 2 coordinate x
    y1 : float
        point 2 coordinate y
    number_of_points : int
        number of points to create along the reach. The default is 1000.

    Returns
    -------
    X : np array
        x coordinate of points
    Y : np array
        y coordinate of points
    dist : np array
        linear distance to point from the point 2

    """
    # create number_of_points points along domain of link
    X = np.linspace(x0, x1, number_of_points)
    Xs = np.abs(X - x0)  # change begin value to zero
    # determine distance from upstream node to each point
    if Xs.max() == 0:  # if a vertical link (x is constant)
        Y = np.linspace(y0, y1, number_of_points)  # y
        dist = np.abs(Y - y0)
    else:
        Y = y0 + (y1 - y0) / np.abs(x1 - x0) * (Xs)  # y
        dist = ((Y - y0) ** 2 + Xs**2) ** 0.5
    return X, Y, dist


def _dist_func(x0, x1, y0, y1):
    """returns linear distance between two points"""
    return ((x0 - x1) ** 2 + (y0 - y1) ** 2) ** 0.5


def extract_channel_nodes(grid, Ct):
    """interpret which nodes of the DEM represent the channel network as all nodes
    that have a drainage area >= to the average drainage area at which
    channels initiate in the DEM (Ct, based on field or remote sensing evidence).

    Ct = average drainage area at which colluvial channels to get the entire
    channel network.
    Ct = the drainage area at which cascade channels typically begin to get
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
    cn_mask = grid.at_node["drainage_area"] >= Ct
    cn = grid.nodes.flatten()[cn_mask]
    return cn


def extract_terrace_nodes(grid, terrace_width, acn, fcn):
    """determine which raster model grid nodes coincide with channel terraces,
    which presently are asssumed to be a fixed width (number of nodes) from
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
    TerraceNodes : np array
        array of all node IDs included in the terrace

    """
    terrace_width = np.round(terrace_width).astype(
        int
    )  # round to int in case provided as float

    if terrace_width < 1:  # check that at least 1
        msg = "terrace width must be 1 or greater"
        raise ValueError(msg)

    for i in range(terrace_width):
        if i == 0:
            # diagonal adjacent nodes to channel nodes
            AdjDN = np.ravel(grid.diagonal_adjacent_nodes_at_node[fcn])
            # adjacent nodes to channel nodes
            AdjN = np.ravel(grid.adjacent_nodes_at_node[fcn])
        elif i > 0:
            # diagonal adjacent nodes to channel nodes
            AdjDN = grid.diagonal_adjacent_nodes_at_node[TerraceNodes]
            # adjacent nodes to channel nodes
            AdjN = grid.adjacent_nodes_at_node[TerraceNodes]
        # all adjacent nodes to channel nodes
        AllNodes = np.concatenate((AdjN, AdjDN))
        # unique adjacent nodes
        AllNodes = np.unique(AllNodes)
        # unique adjacent nodes, excluding all channel nodes.
        TerraceNodes = AllNodes[np.isin(AllNodes, acn, invert=True)]
        # finally, remove any -1 nodes, which represent adjacent nodes outside
        # of the model grid
        TerraceNodes = TerraceNodes[~(TerraceNodes == -1)]
    return TerraceNodes


def min_distance_to_network(grid, acn, node_id):
    """determine the shortest distance (as the crow flies) from a node to the channel network and
    the closest channel node

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

    def distance_to_network(grid, row):
        """compute distance between nodes"""
        return _dist_func(
            row["x"], grid.node_x[node_id], row["y"], grid.node_y[node_id]
        )

    xyDF = pd.DataFrame(
        np.array([grid.node_x[acn], grid.node_y[acn]]).T, columns=["x", "y"]
    )
    xyDF.index = acn
    nmg_dist = xyDF.apply(lambda row: distance_to_network(grid, row), axis=1)
    offset = nmg_dist.min()  # minimum distancce
    mdn = xyDF[nmg_dist == offset].index.values  # find closest node
    if len(mdn) > 1:
        mdn = mdn[0]  # pick first in list if more than one
    return offset, mdn


def map_nmg_links_to_rmg_coincident_nodes(
    grid, nmgrid, link_nodes, remove_duplicates=False
):
    """maps the links of the network model grid to all coincident raster model grid
    nodes. Each coincident raster model grid node is defined in terms of it's
    x and y coordinates, the the link it is mapped to and distance downstream from
    the upstream end (tail end) of the link.


    Parameters
    ----------
    grid : raster model grid
    nmgrid : network model grid
    link_nodes : np array
        head and tail node of each link
    remove_duplicates : bool
        if True, when two or more links are coincident with the same node,
        the node is assigned to the link with the larges drainage area. If False,
        the node is assigned to each coincident link. The default is False.

    Returns
    -------

    nmg_link_to_rmg_coincident_nodes_mapper: pandas dataframe
        each row of the dataframe lists the link ID, the coincident node ID, the
        x and y coordinates and the downstream distance of the coincident node
        and the drainage area of the link

    """
    Lnodelist = []  # list of lists of all nodes that coincide with each link
    Ldistlist = (
        []
    )  # list of lists of the distance on the link (measured from upstream link node) for all nodes that coincide with each link
    LlinkIDlist = []
    Lxlist = []
    Lylist = []
    Lxy = []  # list of all nodes the coincide with the network links
    # loop through all links in network grid to determine raster grid cells that coincide with each link
    # and equivalent distance from upstream node on link
    for linkID, lknd in enumerate(link_nodes):  # for each link in network grid

        x0 = nmgrid.x_of_node[lknd[0]]  # x and y of downstream link node
        y0 = nmgrid.y_of_node[lknd[0]]
        x1 = nmgrid.x_of_node[lknd[1]]  # x and y of upstream link node
        y1 = nmgrid.y_of_node[lknd[1]]

        # x and y coordinates and downstream distance from the upstream node
        # for 1000 points generated from downstream node to upstream node
        X, Y, dist = _link_to_points_and_dist(x0, y0, x1, y1, number_of_points=1000)
        dist = dist.max() - dist  # convert to downstream distance
        nodelist = []  # list of nodes along link
        distlist = []  # list of distance along link corresponding to node
        linkIDlist = []
        xlist = []
        ylist = []
        for i, y in enumerate(Y):
            x = X[i]
            node = grid.find_nearest_node((x, y))  # change to grid.find_nearest_node
            if (
                node not in nodelist
            ):  # if node not already in list, append - many points will be in same cell; only need to list cell once
                nodelist.append(node)
                distlist.append(dist[i])
                linkIDlist.append(linkID)
                xlist.append(grid.node_x[node])
                ylist.append(grid.node_y[node])
                xy = {
                    "linkID": linkID,
                    "coincident_node": node,
                    "x": grid.node_x[node],
                    "y": grid.node_y[node],
                    "dist": dist[i],
                    "drainage_area": nmgrid.at_link["drainage_area"][linkID],
                }
                Lxy.append(xy)

        Lnodelist.append(nodelist)
        Ldistlist.append(distlist)
        LlinkIDlist.append(linkIDlist)
        Lxlist.append(xlist)
        Lylist.append(ylist)

    nmg_link_to_rmg_coincident_nodes_mapper = pd.DataFrame(Lxy)

    # if remove_duplicates, select link with largest mean contributing area.
    if remove_duplicates:
        for link in range(len(link_nodes)):
            for other_link in range(len(link_nodes)):
                if link != other_link:
                    link_coin_nodes = Lnodelist[link]
                    other_link_coin_nodes = Lnodelist[other_link]
                    link_a = nmgrid.at_link["drainage_area"][link]
                    other_link_a = nmgrid.at_link["drainage_area"][other_link]
                    dup = np.intersect1d(link_coin_nodes, other_link_coin_nodes)
                    # if contributing area of link is larger than contributing area
                    # of other_link, remove dupilcate nodes from other link
                    if len(dup) > 0:
                        # print('link {} and link {} have duplicates: {}'.format(link, other_link, dup))
                        if link_a >= other_link_a:
                            mask = ~np.isin(other_link_coin_nodes, dup)
                            Lnodelist[other_link] = list(
                                np.array(other_link_coin_nodes)[mask]
                            )
                            Ldistlist[other_link] = list(
                                np.array(Ldistlist[other_link])[mask]
                            )
                            LlinkIDlist[other_link] = list(
                                np.array(LlinkIDlist[other_link])[mask]
                            )
                            Lxlist[other_link] = list(
                                np.array(Lxlist[other_link])[mask]
                            )
                            Lylist[other_link] = list(
                                np.array(Lylist[other_link])[mask]
                            )
                        else:
                            mask = ~np.isin(link_coin_nodes, dup)
                            Lnodelist[link] = list(np.array(link_coin_nodes)[mask])
                            Ldistlist[link] = list(np.array(Ldistlist[link])[mask])
                            LlinkIDlist[link] = list(np.array(LlinkIDlist[link])[mask])
                            Lxlist[link] = list(np.array(Lxlist[link])[mask])
                            Lylist[link] = list(np.array(Lylist[link])[mask])
        LinkIDs = _flatten_lol(LlinkIDlist)
        nmg_link_to_rmg_coincident_nodes_mapper = pd.DataFrame(
            np.array(
                [
                    LinkIDs,
                    _flatten_lol(Lnodelist),
                    _flatten_lol(Lxlist),
                    _flatten_lol(Lylist),
                    _flatten_lol(Ldistlist),
                    nmgrid.at_link["drainage_area"][np.array(LinkIDs)],
                ]
            ).T,
            columns=["linkID", "coincident_node", "x", "y", "dist", "drainage_area"],
        )
        nmg_link_to_rmg_coincident_nodes_mapper["linkID"] = (
            nmg_link_to_rmg_coincident_nodes_mapper["linkID"].astype(int)
        )
        nmg_link_to_rmg_coincident_nodes_mapper["coincident_node"] = (
            nmg_link_to_rmg_coincident_nodes_mapper["coincident_node"].astype(int)
        )

    return nmg_link_to_rmg_coincident_nodes_mapper
