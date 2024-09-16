import numpy as np


def angle_finder(grid, dn, cn, rn):
    """Find the interior angle between two vectors on a grid.

    Parameters
    ----------
    grid : ModelGrid
        A landlab grid.
    dn : int or array of int
        Node or nodes at the end of the first vector.
    cn : int or array of int
        Node or nodes at the vertex between vectors.
    rn : int or array of int
        Node or nodes at the end of the second vector.

    Returns
    -------
    float or array of float
        Angle between vectors (in radians).

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.lateral_erosion.node_finder import angle_finder

    >>> grid = RasterModelGrid((3, 4))
    >>> np.rad2deg(angle_finder(grid, 8, 5, 0))
    90.0
    >>> np.rad2deg(angle_finder(grid, (8, 9, 10, 6), 5, 6))
    array([135.,  90.,  45.,   0.])
    """
    vertex = np.take(grid.x_of_node, cn), np.take(grid.y_of_node, cn)
    vec_1 = [
        np.take(grid.x_of_node, dn) - vertex[0],
        np.take(grid.y_of_node, dn) - vertex[1],
    ]
    vec_2 = [
        np.take(grid.x_of_node, rn) - vertex[0],
        np.take(grid.y_of_node, rn) - vertex[1],
    ]

    return np.arccos(
        (vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1])
        / np.sqrt((vec_1[0] ** 2 + vec_1[1] ** 2) * (vec_2[0] ** 2 + vec_2[1] ** 2))
    )


def forty_five_node(donor, i, receiver, neighbors, diag_neigh):
    radcurv_angle = 0.67

    lat_node = 0
    # In Landlab 2019: diagonal list goes [NE, NW, SW, SE]. Node list are ordered as [E,N,W,S]
    # if water flows SE-N OR if flow NE-S or E-NW or E-SW, erode west node
    if (
        donor == diag_neigh[0]
        and receiver == neighbors[3]
        or donor == diag_neigh[3]
        and receiver == neighbors[1]
        or donor == neighbors[0]
        and receiver == diag_neigh[2]
        or donor == neighbors[0]
        and receiver == diag_neigh[1]
    ):
        lat_node = neighbors[2]
    # if flow is from SW-N or NW-S or W-NE or W-SE, erode east node
    elif (
        donor == diag_neigh[1]
        and receiver == neighbors[3]
        or donor == diag_neigh[2]
        and receiver == neighbors[1]
        or donor == neighbors[2]
        and receiver == diag_neigh[3]
        or donor == neighbors[2]
        and receiver == diag_neigh[0]
    ):
        lat_node = neighbors[0]
    # if flow is from SE-W or SW-E or S-NE or S-NW, erode north node
    elif (
        donor == diag_neigh[3]
        and receiver == neighbors[2]
        or donor == diag_neigh[2]
        and receiver == neighbors[0]
        or donor == neighbors[3]
        and receiver == diag_neigh[0]
        or donor == neighbors[3]
        and receiver == diag_neigh[1]
    ):
        lat_node = neighbors[1]
    # if flow is from NE-W OR NW-E or N-SE or N-SW, erode south node
    elif (
        donor == diag_neigh[0]
        and receiver == neighbors[2]
        or donor == diag_neigh[1]
        and receiver == neighbors[0]
        or donor == neighbors[1]
        and receiver == diag_neigh[3]
        or donor == neighbors[1]
        and receiver == diag_neigh[2]
    ):
        lat_node = neighbors[3]
    return lat_node, radcurv_angle


def ninety_node(donor, i, receiver, link_list, neighbors, diag_neigh):
    # if flow is 90 degrees
    if donor in diag_neigh and receiver in diag_neigh:
        radcurv_angle = 1.37
        # if flow is NE-SE or NW-SW, erode south node
        if (
            donor == diag_neigh[0]
            and receiver == diag_neigh[3]
            or donor == diag_neigh[1]
            and receiver == diag_neigh[2]
        ):
            lat_node = neighbors[3]
        # if flow is SW-NW or SE-NE, erode north node
        elif (
            donor == diag_neigh[2]
            and receiver == diag_neigh[1]
            or donor == diag_neigh[3]
            and receiver == diag_neigh[0]
        ):
            lat_node = neighbors[1]
        # if flow is SW-SE or NW-NE, erode east node
        elif (
            donor == diag_neigh[2]
            and receiver == diag_neigh[3]
            or donor == diag_neigh[1]
            and receiver == diag_neigh[0]
        ):
            lat_node = neighbors[0]
        # if flow is SE-SW or NE-NW, erode west node
        elif (
            donor == diag_neigh[3]
            and receiver == diag_neigh[2]
            or donor == diag_neigh[0]
            and receiver == diag_neigh[1]
        ):
            lat_node = neighbors[2]
    elif donor not in diag_neigh and receiver not in diag_neigh:
        radcurv_angle = 1.37
        # if flow is from east, erode west node
        if donor == neighbors[0]:
            lat_node = neighbors[2]
        # if flow is from north, erode south node
        elif donor == neighbors[1]:
            lat_node = neighbors[3]
        # if flow is from west, erode east node
        elif donor == neighbors[2]:
            lat_node = neighbors[0]
        # if flow is from south, erode north node
        elif donor == neighbors[3]:
            lat_node = neighbors[1]
    return lat_node, radcurv_angle


def straight_node(donor, i, receiver, neighbors, diag_neigh):
    # ***FLOW LINK IS STRAIGHT, NORTH TO SOUTH***#
    if donor == neighbors[1] or donor == neighbors[3]:
        # print "flow is stright, N-S from ", donor, " to ", flowdirs[i]
        radcurv_angle = 0.23
        # neighbors are ordered E,N,W, S
        # if the west cell is boundary (neighbors=-1), erode from east node
        if neighbors[2] == -1:
            lat_node = neighbors[0]
        elif neighbors[0] == -1:
            lat_node = neighbors[2]
        else:
            # if could go either way, choose randomly. 0 goes East, 1 goes west
            ran_num = np.random.randint(0, 2)
            if ran_num == 0:
                lat_node = neighbors[0]
            if ran_num == 1:
                lat_node = neighbors[2]
    # ***FLOW LINK IS STRAIGHT, EAST-WEST**#
    elif donor == neighbors[0] or donor == neighbors[2]:
        radcurv_angle = 0.23
        #  Node list are ordered as [E,N,W,S]
        # if the north cell is boundary (neighbors=-1), erode from south node
        if neighbors[1] == -1:
            lat_node = neighbors[3]
        elif neighbors[3] == -1:
            lat_node = neighbors[1]
        else:
            # if could go either way, choose randomly. 0 goes south, 1 goes north
            ran_num = np.random.randint(0, 2)
            if ran_num == 0:
                lat_node = neighbors[1]
            if ran_num == 1:
                lat_node = neighbors[3]
    # if flow is straight across diagonal, choose node to erode at random
    elif donor in diag_neigh and receiver in diag_neigh:
        radcurv_angle = 0.23
        if receiver == diag_neigh[0]:
            poss_diag_nodes = neighbors[0 : 1 + 1]
        elif receiver == diag_neigh[1]:
            poss_diag_nodes = neighbors[1 : 2 + 1]
        elif receiver == diag_neigh[2]:
            poss_diag_nodes = neighbors[2 : 3 + 1]
        elif receiver == diag_neigh[3]:
            poss_diag_nodes = [neighbors[3], neighbors[0]]
        ran_num = np.random.randint(0, 2)
        if ran_num == 0:
            lat_node = poss_diag_nodes[0]
        if ran_num == 1:
            lat_node = poss_diag_nodes[1]
    return lat_node, radcurv_angle


def node_finder(grid, i, flowdirs, drain_area):
    """Find lateral neighbor node of the primary node for straight, 45 degree,
    and 90 degree channel segments.

    Parameters
    ----------
    grid : ModelGrid
        A Landlab grid object
    i : int
        node ID of primary node
    flowdirs : array
        Flow direction array
    drain_area : array
        drainage area array

    Returns
    -------
    lat_node : int
        node ID of lateral node
    radcurv_angle : float
        inverse radius of curvature of channel at lateral node
    """
    # receiver node of flow is flowdirs[i]
    receiver = flowdirs[i]

    # find indicies of where flowdirs=i to find donor nodes.
    # will donor nodes always equal the index of flowdir list?
    inflow = np.where(flowdirs == i)

    # if there are more than 1 donors, find the one with largest drainage area

    if len(inflow[0]) > 1:
        drin = drain_area[inflow]
        drmax = max(drin)
        maxinfl = inflow[0][np.where(drin == drmax)]
        # if donor nodes have same drainage area, choose one randomly
        if len(maxinfl) > 1:
            ran_num = np.random.randint(0, len(maxinfl))
            maxinfln = maxinfl[ran_num]
            donor = [maxinfln]
        else:
            donor = maxinfl
        # if inflow is empty, no donor
    elif len(inflow[0]) == 0:
        donor = i
    # else donor is the only inflow
    else:
        donor = inflow[0]
    # now we have chosen donor cell, next figure out if inflow/outflow lines are
    # straight, 45, or 90 degree angle. and figure out which node to erode
    link_list = grid.links_at_node[i]
    # this gives list of active neighbors for specified node
    # the order of this list is: [E,N,W,S]
    neighbors = grid.active_adjacent_nodes_at_node[i]
    # this gives list of all diagonal neighbors for specified node
    # the order of this list is: [NE,NW,SW,SE]
    diag_neigh = grid.diagonal_adjacent_nodes_at_node[i]
    angle_diff = np.rad2deg(angle_finder(grid, donor, i, receiver))

    if (donor == flowdirs[i]) or (donor == i):
        # this is a sink. no lateral ero
        radcurv_angle = 0.0
        lat_node = 0
    elif np.isclose(angle_diff, 0.0) or np.isclose(angle_diff, 180.0):
        [lat_node, radcurv_angle] = straight_node(
            donor, i, receiver, neighbors, diag_neigh
        )
    elif np.isclose(angle_diff, 45.0) or np.isclose(angle_diff, 135.0):
        [lat_node, radcurv_angle] = forty_five_node(
            donor, i, receiver, neighbors, diag_neigh
        )
    elif np.isclose(angle_diff, 90.0):
        [lat_node, radcurv_angle] = ninety_node(
            donor, i, receiver, link_list, neighbors, diag_neigh
        )
    else:
        lat_node = 0
        radcurv_angle = 0.0

    dx = grid.dx
    # INVERSE radius of curvature.
    radcurv_angle = radcurv_angle / dx
    return int(lat_node), radcurv_angle
