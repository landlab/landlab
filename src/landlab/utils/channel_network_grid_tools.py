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
