#! /usr/env/python

"""
flow_direction_mfd.py: calculate multiple-flow-direction flow directions.

Works on both a regular or irregular grid. Also calculates flow proportions.

KRB Jan 2017
"""

import numpy as np

from landlab.core.utils import as_id_array
from landlab.grid.base import BAD_INDEX_VALUE


def flow_directions_mfd(
    elev,
    neighbors_at_node,
    links_at_node,
    active_link_dir_at_node,
    link_slope,
    baselevel_nodes=None,
    partition_method="slope",
):
    """Find multiple-flow-direction flow directions on a grid.

    Finds and returns flow directions and proportions for a given elevation
    grid. Each node is assigned multiple flow directions, toward all of the N
    neighboring nodes that are lower than it. If none of the neighboring nodes
    are lower, it is assigned to itself. Flow proportions can be calculated as
    proportional to slope (default) or proportional to the square root of
    slope, which is the solution to a steady kinematic wave.

    Parameters
    ----------
    elev : array_like
        Elevations at nodes.
    neighbors_at_node : array_like (num nodes, max neighbors at node)
        For each node, the link IDs of active links.
    links_at_node : array_like (num nodes, max neighbors at node)

    link_dir_at_node: array_like (num nodes, max neighbors at node)

        IDs of the head node for each link.
    link_slope : array_like
        slope of each link, defined POSITIVE DOWNHILL (i.e., a negative value
        means the link runs uphill from the fromnode to the tonode).
    baselevel_nodes : array_like, optional
        IDs of open boundary (baselevel) nodes.
    partition_method: string, optional
        Method for partitioning flow. Options include 'slope' (default) and
        'square_root_of_slope'.

    Returns
    -------
    receivers : ndarray of size (num nodes, max neighbors at node)
        For each node, the IDs of the nodes that receive its flow. For nodes
        that do not direct flow to all neighbors, BAD_INDEX_VALUE is given as
        a placeholder. The ID of the node itself is given if no other receiver
        is assigned.
    proportions : ndarray of size (num nodes, max neighbors at node)
        For each receiver, the proportion of flow (between 0 and 1) is given.
        A proportion of zero indicates that the link does not have flow along
        it.
    slopes: ndarray of size (num nodes, max neighbors at node)
        For each node in the array ``recievers``, the slope value (positive
        downhill) in the direction of flow. If no flow occurs (value of
        ``recievers`` is -1), then this array is set to 0.
    steepest_slope : ndarray
        The slope value (positive downhill) in the direction of flow.
    steepest_receiver : ndarray
        For each node, the node ID of the node connected by the steepest link.
        BAD_INDEX_VALUE is given if no flow emmanates from the node.
    sink : ndarray
        IDs of nodes that are flow sinks (they are their own receivers)
    receiver_links : ndarray of size (num nodes, max neighbors at node)
        ID of links that leads from each node to its receiver, or
        BAD_INDEX_VALUE if no flow occurs on this link.
    steepest_link : ndarray
        For each node, the link ID of the steepest link.
        BAD_INDEX_VALUE is given if no flow emmanates from the node.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> import numpy as np
    >>> from landlab.components.flow_director.flow_direction_mfd import (
    ...     flow_directions_mfd,
    ... )
    >>> grid = RasterModelGrid((3, 3), xy_spacing=(1, 1))
    >>> elev = grid.add_field(
    ...     "topographic__elevation",
    ...     grid.node_x + grid.node_y,
    ...     at="node",
    ... )

    For the first example, we will not pass any diagonal elements to the flow
    direction algorithm.

    >>> neighbors_at_node = grid.adjacent_nodes_at_node
    >>> links_at_node = grid.links_at_node
    >>> active_link_dir_at_node = grid.active_link_dirs_at_node
    >>> link_slope = np.arctan(grid.calc_grad_at_link(elev))
    >>> slopes_to_neighbors_at_node = (
    ...     link_slope[links_at_node] * active_link_dir_at_node
    ... )
    >>> (
    ...     receivers,
    ...     proportions,
    ...     slopes,
    ...     steepest_slope,
    ...     steepest_receiver,
    ...     sink,
    ...     receiver_links,
    ...     steepest_link,
    ... ) = flow_directions_mfd(
    ...     elev,
    ...     neighbors_at_node,
    ...     links_at_node,
    ...     active_link_dir_at_node,
    ...     link_slope,
    ...     baselevel_nodes=None,
    ...     partition_method="slope",
    ... )
    >>> receivers
    array([[ 0, -1, -1, -1],
           [ 1, -1, -1, -1],
           [ 2, -1, -1, -1],
           [ 3, -1, -1, -1],
           [-1, -1,  3,  1],
           [-1, -1,  4, -1],
           [ 6, -1, -1, -1],
           [-1, -1, -1,  4],
           [ 8, -1, -1, -1]])
    >>> proportions
    array([[1. , 0. , 0. , 0. ],
           [1. , 0. , 0. , 0. ],
           [1. , 0. , 0. , 0. ],
           [1. , 0. , 0. , 0. ],
           [0. , 0. , 0.5, 0.5],
           [0. , 0. , 1. , 0. ],
           [1. , 0. , 0. , 0. ],
           [0. , 0. , 0. , 1. ],
           [1. , 0. , 0. , 0. ]])
    >>> proportions.sum(axis=-1)
    array([1., 1., 1., 1., 1., 1., 1., 1., 1.])

    In the second example, we will pass diagonal elements to the flow direction
    algorithm.

    >>> dal = grid.active_d8
    >>> neighbors_at_node = np.hstack(
    ...     (grid.adjacent_nodes_at_node, grid.diagonal_adjacent_nodes_at_node)
    ... )
    >>> links_at_node = grid.d8s_at_node
    >>> active_link_dir_at_node = grid.active_d8_dirs_at_node

    We need to create a list of diagonal links since it doesn't exist.

    >>> diag_links = np.sort(np.unique(grid.d8s_at_node[:, 4:]))
    >>> diag_links = diag_links[diag_links > 0]
    >>> diag_grads = np.zeros(diag_links.shape)
    >>> where_active_diag = dal >= diag_links.min()
    >>> active_diags_inds = dal[where_active_diag] - diag_links.min()
    >>> diag_grads = grid.calc_grad_at_diagonal(elev)
    >>> ortho_grads = grid.calc_grad_at_link(elev)
    >>> link_slope = np.hstack((np.arctan(ortho_grads), np.arctan(diag_grads)))
    >>> (
    ...     receivers,
    ...     proportions,
    ...     slopes,
    ...     steepest_slope,
    ...     steepest_receiver,
    ...     sink,
    ...     receiver_links,
    ...     steepest_link,
    ... ) = flow_directions_mfd(
    ...     elev,
    ...     neighbors_at_node,
    ...     links_at_node,
    ...     active_link_dir_at_node,
    ...     link_slope,
    ...     baselevel_nodes=None,
    ...     partition_method="slope",
    ... )
    >>> receivers
    array([[ 0, -1, -1, -1, -1, -1, -1, -1],
           [ 1, -1, -1, -1, -1, -1, -1, -1],
           [ 2, -1, -1, -1, -1, -1, -1, -1],
           [ 3, -1, -1, -1, -1, -1, -1, -1],
           [-1, -1,  3,  1, -1, -1,  0, -1],
           [-1, -1,  4, -1, -1, -1, -1, -1],
           [ 6, -1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1,  4, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1,  4, -1]])
    >>> proportions
    array([[1.        , 0.        , 0.        , 0.        , 0.        ,
            0.        , 0.        , 0.        ],
           [1.        , 0.        , 0.        , 0.        , 0.        ,
            0.        , 0.        , 0.        ],
           [1.        , 0.        , 0.        , 0.        , 0.        ,
            0.        , 0.        , 0.        ],
           [1.        , 0.        , 0.        , 0.        , 0.        ,
            0.        , 0.        , 0.        ],
           [0.        , 0.        , 0.31091174,  0.31091174, 0.        ,
            0.        , 0.37817653, 0.        ],
           [0.        , 0.        , 1.        , 0.        , 0.        ,
            0.        , 0.        , 0.        ],
           [1.        , 0.        , 0.        , 0.        , 0.        ,
            0.        , 0.        , 0.        ],
           [0.        , 0.        , 0.        , 1.        , 0.        ,
            0.        , 0.        , 0.        ],
           [0.        , 0.        , 0.        , 0.        , 0.        ,
            0.        , 1.        , 0.        ]])
    >>> slopes
    array([[0.        , 0.        , 0.        , 0.        , 0.        ,
            0.        , 0.        , 0.        ],
           [0.        , 0.        , 0.        , 0.        , 0.        ,
            0.        , 0.        , 0.        ],
           [0.        , 0.        , 0.        , 0.        , 0.        ,
            0.        , 0.        , 0.        ],
           [0.        , 0.        , 0.        , 0.        , 0.        ,
            0.        , 0.        , 0.        ],
           [0.        , 0.        , 0.78539816, 0.78539816, 0.        ,
            0.        , 0.95531662, 0.        ],
           [0.        , 0.        , 0.78539816, 0.        , 0.        ,
            0.        , 0.        , 0.        ],
           [0.        , 0.        , 0.        , 0.        , 0.        ,
            0.        , 0.        , 0.        ],
           [0.        , 0.        , 0.        , 0.78539816, 0.        ,
            0.        , 0.        , 0.        ],
           [0.        , 0.        , 0.        , 0.        , 0.        ,
            0.        , 0.95531662, 0.        ]])
    >>> proportions.sum(axis=-1)
    array([1., 1., 1., 1., 1., 1., 1., 1., 1.])
    """
    # Calculate the number of nodes.
    num_nodes = len(elev)

    # Create a node array
    node_id = np.arange(num_nodes)

    # Calculate the maximum number of neighbors at node.
    max_number_of_neighbors = neighbors_at_node.shape[1]

    # Make a copy of neighbors_at_node so we can change it into the receiver
    # array.
    receivers = neighbors_at_node.copy()

    # Construct the array of slope to neighbors at node. This also will adjust
    # for the slope convention based on the direction of the link.
    slopes_to_neighbors_at_node = link_slope[links_at_node] * active_link_dir_at_node

    # Make a copy so this can be changed based on where no flow occurs.
    receiver_links = links_at_node.copy()

    # some of these potential recievers may have already been assigned as
    # BAD_INDEX_VALUE because the link was inactive. Make a mask of these for
    # future use. Also find the close nodes.
    inactive_link_to_neighbor = active_link_dir_at_node == 0
    closed_nodes = np.sum(np.abs(active_link_dir_at_node), 1) == 0
    # Now calculate where flow occurs.
    # First, make an elevation array of potential receivers.
    potential_receiver_elev = elev[neighbors_at_node]

    # now make an array of the same shape (for direct comparison) of the source
    # node elevation.
    source_node_elev = elev[np.tile(node_id, (max_number_of_neighbors, 1)).T]

    # find where flow does not occur (source is lower that receiver)
    flow_does_not_occur = source_node_elev <= potential_receiver_elev

    # Where the source is lower, set receivers to BAD_INDEX_VALUE
    receivers[flow_does_not_occur] = BAD_INDEX_VALUE

    # Where the link is not active, set receivers to BAD_INDEX_VALUE
    receivers[inactive_link_to_neighbor] = BAD_INDEX_VALUE

    # Next, find where a node drains to itself
    drains_to_self = receivers.sum(1) == -1 * max_number_of_neighbors

    # Where this occurs, set the receiver ID in the first column of receivers
    # to the node ID.
    receivers[drains_to_self, 0] = node_id[drains_to_self]

    # Finally, set the first element of the closed nodes to themselves.
    receivers[closed_nodes, 0] = node_id[closed_nodes]

    # Next, calculate flow proportions.
    # Copy slope array and mask by where flow is not occuring and where the
    # link is inactive.
    flow_slopes = slopes_to_neighbors_at_node.copy()
    flow_slopes[flow_does_not_occur] = 0.0
    flow_slopes[inactive_link_to_neighbor] = 0.0

    if partition_method == "square_root_of_slope":
        values_for_partitioning = flow_slopes**0.5
    elif partition_method == "slope":
        values_for_partitioning = flow_slopes
    else:
        raise ValueError("Keyword argument to partition_method invalid.")

    # Calculate proportions by normalizing by rowsums.
    denom = np.tile(values_for_partitioning.sum(1), (max_number_of_neighbors, 1)).T
    denom[denom <= 0] = 1  # to prevent runtime errors
    proportions = values_for_partitioning / denom
    proportions[drains_to_self, 0] = 1
    proportions[drains_to_self, 1:] = 0

    # Might need to sort by proportions and rearrange to follow expectations
    # of no BAD_INDEX_VALUE value in first column. KRB NOT SURE

    # mask the receiver_links by where flow doesn't occur to return
    receiver_links[flow_does_not_occur] = BAD_INDEX_VALUE
    receiver_links[inactive_link_to_neighbor] = BAD_INDEX_VALUE

    # identify the steepest link so that the steepest receiver, link, and slope
    # can be returned.
    slope_sort = np.argsort(
        np.argsort(flow_slopes, axis=1, kind="stable"), axis=1, kind="stable"
    ) == (max_number_of_neighbors - 1)
    steepest_slope = flow_slopes[slope_sort]

    # identify the steepest link and steepest receiever.
    steepest_link = receiver_links[slope_sort]
    steepest_receiver = receivers[slope_sort]
    steepest_receiver[drains_to_self] = node_id[drains_to_self]

    # Optionally, handle baselevel nodes: they are their own receivers
    if baselevel_nodes is not None:
        receivers[baselevel_nodes, 0] = node_id[baselevel_nodes]
        receivers[baselevel_nodes, 1:] = -1
        proportions[baselevel_nodes, 0] = 1
        proportions[baselevel_nodes, 1:] = 0
        receiver_links[baselevel_nodes, :] = BAD_INDEX_VALUE
        steepest_slope[baselevel_nodes] = 0.0

    # The sink nodes are those that are their own receivers (this will normally
    # include boundary nodes as well as interior ones; "pits" would be sink
    # nodes that are also interior nodes).
    (sink,) = np.where(node_id == receivers[:, 0])
    sink = as_id_array(sink)

    slopes_to_neighbors_at_node[flow_does_not_occur] = 0
    slopes_to_neighbors_at_node[inactive_link_to_neighbor] = 0

    return (
        receivers,
        proportions,
        slopes_to_neighbors_at_node,
        steepest_slope,
        steepest_receiver,
        sink,
        receiver_links,
        steepest_link,
    )


if __name__ == "__main__":  # pragma: no cover
    import doctest

    doctest.testmod()
