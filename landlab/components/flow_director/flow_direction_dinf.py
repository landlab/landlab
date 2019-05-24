# def direct_dinf(grid, elevs='topographic_elevation', baselevel_nodes=None):


"""
flow_direction_dinf.py: calculate Dinfinity flow direction on raster grids.

Calculates flow direction and proportion on a raster grid by the Dinfinity
algorithm of Tarboton 1997.

KRB Feb 2017
"""

import numpy as np

from landlab import VoronoiDelaunayGrid  # for type tests
from landlab import BAD_INDEX_VALUE, CLOSED_BOUNDARY
from landlab.core.utils import as_id_array
from landlab.utils.return_array import return_array_at_node

UNDEFINED_INDEX = BAD_INDEX_VALUE


def flow_directions_dinf(grid, elevs="topographic__elevation", baselevel_nodes=None):
    """Find Dinfinity flow directions and proportions on a raster grid.

    Finds and returns flow directions and proportions for a given elevation
    grid by the D infinity method (Tarboton, 1997). Each node is assigned two
    flow directions, toward the two neighboring nodes that are on the steepest
    subtriangle. Partitioning of flow is done based on the aspect of the
    subtriangle.

    This method does not support irregular grids.

    Parameters
    ----------
    grid : ModelGrid
        A grid of type Voroni.
    elevs : field name at node or array of length node
        The surface to direct flow across.
    baselevel_nodes : array_like, optional
        IDs of open boundary (baselevel) nodes.

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
        UNDEFINED_INDEX if no flow occurs on this link.
    steepest_link : ndarray
        For each node, the link ID of the steepest link.
        BAD_INDEX_VALUE is given if no flow emmanates from the node.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.flow_director.flow_direction_dinf import(
    ...                                                  flow_directions_dinf)

    Dinfinity routes flow based on the relative proportion of flow along the
    triangular facets around a central raster node.

    >>> grid = RasterModelGrid((3,3), xy_spacing=(1, 1))
    >>> _ = grid.add_field('topographic__elevation',
    ...                     2.*grid.node_x+grid.node_y,
    ...                     at = 'node')
    >>> (receivers, proportions, slopes,
    ... steepest_slope, steepest_receiver,
    ... sink, receiver_links, steepest_link) = flow_directions_dinf(grid)
    >>> receivers
    array([[ 0, -1],
           [ 0, -1],
           [ 1, -1],
           [ 0, -1],
           [ 3,  0],
           [ 4,  1],
           [ 3, -1],
           [ 6,  3],
           [ 7,  4]])
    >>> proportions
    array([[ 1.        ,  0.        ],
           [ 1.        , -0.        ],
           [ 1.        , -0.        ],
           [ 1.        ,  0.        ],
           [ 0.40966553,  0.59033447],
           [ 0.40966553,  0.59033447],
           [ 1.        ,  0.        ],
           [ 0.40966553,  0.59033447],
           [ 0.40966553,  0.59033447]])

    This method also works if the elevations are passed as an array instead of
    the (implied) field name 'topographic__elevation'.

    >>> z = grid['node']['topographic__elevation']
    >>> (receivers, proportions, slopes,
    ... steepest_slope, steepest_receiver,
    ... sink, receiver_links, steepest_link) = flow_directions_dinf(grid, z)
    >>> receivers
    array([[ 0, -1],
           [ 0, -1],
           [ 1, -1],
           [ 0, -1],
           [ 3,  0],
           [ 4,  1],
           [ 3, -1],
           [ 6,  3],
           [ 7,  4]])
    >>> slopes
    array([[-1.        , -2.12132034],
           [ 2.        ,  0.70710678],
           [ 2.        ,  0.70710678],
           [ 1.        , -0.70710678],
           [ 2.        ,  2.12132034],
           [ 2.        ,  2.12132034],
           [ 1.        , -0.70710678],
           [ 2.        ,  2.12132034],
           [ 2.        ,  2.12132034]])
    >>> proportions
    array([[ 1.        ,  0.        ],
           [ 1.        , -0.        ],
           [ 1.        , -0.        ],
           [ 1.        ,  0.        ],
           [ 0.40966553,  0.59033447],
           [ 0.40966553,  0.59033447],
           [ 1.        ,  0.        ],
           [ 0.40966553,  0.59033447],
           [ 0.40966553,  0.59033447]])
    """
    # grid type testing
    if isinstance(grid, VoronoiDelaunayGrid):
        raise NotImplementedError(
            "Dinfinity is currently implemented for" " Raster grids only"
        )
    # get elevs
    elevs = return_array_at_node(grid, elevs)

    # find where there are closed nodes.
    closed_nodes = grid.status_at_node == CLOSED_BOUNDARY

    closed_elevation = np.max(elevs[~closed_nodes]) + 1000

    elevs[closed_nodes] = closed_elevation

    # Step 1, some basic set-up, gathering information about the grid.

    # Calculate the number of nodes.
    num_nodes = len(elevs)

    # Set the number of receivers and facets.
    num_receivers = 2
    num_facets = 8

    # Create a node array
    node_id = np.arange(num_nodes)

    # create an array of the triangle numbers
    tri_numbers = np.arange(num_facets)

    # Step 3, create some triangle datastructures because landlab (smartly)
    # makes it hard to deal with diagonals.

    # create list of triangle neighbors at node. Use orientation associated
    # with tarboton's 1997 algorithm, orthogonal link first, then diagonal.
    # has shape, (nnodes, 8 triangles, 2 neighbors)
    n_at_node = grid.adjacent_nodes_at_node
    dn_at_node = grid.diagonal_adjacent_nodes_at_node

    triangle_neighbors_at_node = np.stack(
        [
            np.vstack((n_at_node[:, 0], dn_at_node[:, 0])),
            np.vstack((n_at_node[:, 1], dn_at_node[:, 0])),
            np.vstack((n_at_node[:, 1], dn_at_node[:, 1])),
            np.vstack((n_at_node[:, 2], dn_at_node[:, 1])),
            np.vstack((n_at_node[:, 2], dn_at_node[:, 2])),
            np.vstack((n_at_node[:, 3], dn_at_node[:, 2])),
            np.vstack((n_at_node[:, 3], dn_at_node[:, 3])),
            np.vstack((n_at_node[:, 0], dn_at_node[:, 3])),
        ],
        axis=-1,
    )
    triangle_neighbors_at_node = triangle_neighbors_at_node.swapaxes(0, 1)

    # next create, triangle links at node
    l_at_node = grid.d8s_at_node[:, :4]
    dl_at_node = grid.d8s_at_node[:, 4:]
    triangle_links_at_node = np.stack(
        [
            np.vstack((l_at_node[:, 0], dl_at_node[:, 0])),
            np.vstack((l_at_node[:, 1], dl_at_node[:, 0])),
            np.vstack((l_at_node[:, 1], dl_at_node[:, 1])),
            np.vstack((l_at_node[:, 2], dl_at_node[:, 1])),
            np.vstack((l_at_node[:, 2], dl_at_node[:, 2])),
            np.vstack((l_at_node[:, 3], dl_at_node[:, 2])),
            np.vstack((l_at_node[:, 3], dl_at_node[:, 3])),
            np.vstack((l_at_node[:, 0], dl_at_node[:, 3])),
        ],
        axis=-1,
    )
    triangle_links_at_node = triangle_links_at_node.swapaxes(0, 1)

    # next create link directions and active link directions at node
    # link directions
    ld_at_node = grid.link_dirs_at_node
    dld_at_node = grid.diagonal_dirs_at_node
    triangle_link_dirs_at_node = np.stack(
        [
            np.vstack((ld_at_node[:, 0], dld_at_node[:, 0])),
            np.vstack((ld_at_node[:, 1], dld_at_node[:, 0])),
            np.vstack((ld_at_node[:, 1], dld_at_node[:, 1])),
            np.vstack((ld_at_node[:, 2], dld_at_node[:, 1])),
            np.vstack((ld_at_node[:, 2], dld_at_node[:, 2])),
            np.vstack((ld_at_node[:, 3], dld_at_node[:, 2])),
            np.vstack((ld_at_node[:, 3], dld_at_node[:, 3])),
            np.vstack((ld_at_node[:, 0], dld_at_node[:, 3])),
        ],
        axis=-1,
    )
    triangle_link_dirs_at_node = triangle_link_dirs_at_node.swapaxes(0, 1)

    # need to create a list of diagonal links since it doesn't exist.
    diag_links = np.sort(np.unique(grid.d8s_at_node[:, 4:]))
    diag_links = diag_links[diag_links > 0]

    # calculate graidents across diagonals and orthogonals
    diag_grads = grid._calculate_gradients_at_d8_links(elevs)
    ortho_grads = grid.calc_grad_at_link(elevs)

    # finally compile link slopes
    link_slope = np.hstack((ortho_grads, diag_grads))

    # Construct the array of slope to triangles at node. This also will adjust
    # for the slope convention based on the direction of the links.
    # this is a (nnodes, 2, 8) array
    slopes_to_triangles_at_node = (
        link_slope[triangle_links_at_node] * triangle_link_dirs_at_node
    )

    # Step 3: make arrays necessary for the specific tarboton algorithm.
    # create a arrays
    ac = np.array([0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0])
    af = np.array([1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0])

    # construct d1 and d2, we know these because we know where the orthogonal
    # links are
    diag_length = ((grid.dx) ** 2 + (grid.dy) ** 2) ** 0.5

    # for irregular grids, d1 and d2 will need to be matricies
    d1 = np.array(
        [grid.dx, grid.dy, grid.dy, grid.dx, grid.dx, grid.dy, grid.dy, grid.dy]
    )
    d2 = np.array(
        [grid.dx, grid.dx, grid.dy, grid.dy, grid.dx, grid.dx, grid.dy, grid.dy]
    )

    thresh = np.arctan(d2 / d1)

    # Step 4, Initialize receiver and proportion arrays
    receivers = UNDEFINED_INDEX * np.ones((num_nodes, num_receivers), dtype=int)
    receiver_closed = UNDEFINED_INDEX * np.ones((num_nodes, num_receivers), dtype=int)
    proportions = np.zeros((num_nodes, num_receivers), dtype=float)
    receiver_links = UNDEFINED_INDEX * np.ones((num_nodes, num_receivers), dtype=int)
    slopes_to_receivers = np.zeros((num_nodes, num_receivers), dtype=float)

    # Step  5  begin the algorithm in earnest

    # construct e0, e1, e2 for all triangles at all nodes.
    # will be (nnodes, nfacets=8 for raster or nfacets = max number of patches
    # for irregular grids.

    # e0 is origin point of the facet
    e0 = elevs[node_id]

    # e1 is the point on the orthogoanal edges
    e1 = elevs[triangle_neighbors_at_node[:, 0, :]]
    # e2 is the point on the diagonal edges
    e2 = elevs[triangle_neighbors_at_node[:, 1, :]]

    # modification of original algorithm to address Landlab boundary conditions.
    # first,
    # for e1 and e2, mask out where nodes do not exits (e.g.
    # triangle_neighbors_at_node == -1)
    e1[triangle_neighbors_at_node[:, 0, :] == -1] = np.nan
    e2[triangle_neighbors_at_node[:, 1, :] == -1] = np.nan

    # loop through and calculate s1 and s2
    # this will only loop nfacets times.
    s1 = np.empty_like(e1)
    s2 = np.empty_like(e2)

    for i in range(num_facets):
        s1[:, i] = (e0 - e1[:, i]) / d1[i]
        s2[:, i] = (e1[:, i] - e2[:, i]) / d2[i]

    # calculate r and s, the direction and magnitude
    r = np.arctan2(s2, s1)
    s = ((s1 ** 2) + (s2 ** 2)) ** 0.5

    r[np.isnan(r)] = 0
    # adjust r if it doesn't sit in the realm of (0, arctan(d2,d1))
    too_small = r < 0
    radj = r.copy()
    radj[too_small] = 0
    s[too_small] = s1[too_small]

    # to consider two big, we need to look by triangle.
    for i in range(num_facets):
        too_big = r[:, i] > thresh[i]
        radj[too_big, i] = thresh[i]
        s[too_big, i] = (e0[too_big] - e2[too_big, i]) / diag_length

    # calculate the geospatial version of r based on radj
    rg = np.empty_like(r)
    for i in range(num_facets):
        rg[:, i] = (af[i] * radj[:, i]) + (ac[i] * np.pi / 2.0)

    # set slopes that are nan to below zero
    # if there is a flat slope, it should be chosen over the closed or non-existant
    # triangles that are represented by the nan values.
    s[np.isnan(s)] = -999.0

    # sort slopes
    # we've set slopes going to closed or non-existant triangles to -999.0, so
    # we shouldn't ever choose these.
    steepest_sort = np.argsort(s)

    # determine the steepest triangle
    steepest_triangle = tri_numbers[steepest_sort[:, -1]]

    # initialize arrays for the steepest rg and steepest s
    steepest_rg = np.empty_like(node_id, dtype=float)
    steepest_s = np.empty_like(node_id, dtype=float)

    closed_triangle_neighbors = closed_nodes[triangle_neighbors_at_node]
    for n in node_id:
        steepest_rg[n] = rg[n, steepest_sort[n, -1]]
        receiver_closed[n] = closed_triangle_neighbors[n, :, steepest_sort[n, -1]]
        steepest_s[n] = s[n, steepest_sort[n, -1]]
        receivers[n, :] = triangle_neighbors_at_node[n, :, steepest_sort[n, -1]]
        receiver_links[n, :] = triangle_links_at_node[n, :, steepest_sort[n, -1]]
        slopes_to_receivers[n, :] = slopes_to_triangles_at_node[
            n, :, steepest_sort[n, -1]
        ]

    # construct the baseline for proportions
    rg_baseline = np.array([0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4]) * np.pi / 2.0

    # calculate alpha1 and alpha 2
    alpha2 = (steepest_rg - rg_baseline[steepest_triangle]) * af[steepest_triangle]
    alpha1 = thresh[steepest_triangle] - alpha2

    # calculate proportions from alpha
    proportions[:, 0] = (alpha1) / (alpha1 + alpha2)
    proportions[:, 1] = (alpha2) / (alpha1 + alpha2)

    # where proportions == 0, set reciever  to -1
    receivers[proportions == 0] = -1

    # END OF THE Tarboton algorithm, start of work to make this code mesh
    # with other landlab flow directing algorithms.

    # identify what drains to itself, and set proportion and id values based on
    # that.

    # if proportions is nan, drain to self
    drains_to_self = np.isnan(proportions[:, 0])

    # if all slopes are leading out or flat, drain to self
    drains_to_self[steepest_s <= 0] = True

    # if both receiver nodes are closed, drain to self
    drains_to_two_closed = receiver_closed.sum(axis=1) == num_receivers
    drains_to_self[drains_to_two_closed] = True

    # if drains to one closed receiver, check that the open receiver actually
    # gets flow. If so, route all to the open receiver. If the receiver getting
    # all the flow is closed, then drain to self.
    all_flow_to_closed = np.sum(receiver_closed * proportions, axis=1) == 1
    drains_to_self[all_flow_to_closed] = True

    drains_to_one_closed = receiver_closed.sum(axis=1) == 1
    fix_flow = drains_to_one_closed * (~all_flow_to_closed)
    first_column_has_closed = np.array(receiver_closed[:, 0] * fix_flow, dtype=bool)
    second_column_has_closed = np.array(receiver_closed[:, 1] * fix_flow, dtype=bool)

    # remove the link to the closed node
    receivers[first_column_has_closed, 0] = -1
    receivers[second_column_has_closed, 1] = -1

    # change the proportions
    proportions[first_column_has_closed, 0] = 0.0
    proportions[first_column_has_closed, 1] = 1.0

    proportions[second_column_has_closed, 0] = 1.0
    proportions[second_column_has_closed, 1] = 0.0

    # set properties of drains to self.
    receivers[drains_to_self, 0] = node_id[drains_to_self]
    receivers[drains_to_self, 1] = -1

    proportions[drains_to_self, 0] = 1.0
    proportions[drains_to_self, 1] = 0.0

    # set properties of closed
    receivers[closed_nodes, 0] = node_id[closed_nodes]
    receivers[closed_nodes, 1] = -1
    proportions[closed_nodes, 0] = 1.0
    proportions[closed_nodes, 1] = 0.0

    # mask the receiver_links by where flow doesn't occur to return
    receiver_links[drains_to_self, :] = UNDEFINED_INDEX

    # identify the steepest link so that the steepest receiver, link, and slope
    # can be returned.
    slope_sort = np.argsort(np.argsort(slopes_to_receivers, axis=1), axis=1) == (
        num_receivers - 1
    )
    steepest_slope = slopes_to_receivers[slope_sort]
    steepest_slope[drains_to_self] = 0.0

    # identify the steepest link and steepest receiever.
    steepest_link = receiver_links[slope_sort]
    steepest_link[drains_to_self] = UNDEFINED_INDEX

    steepest_receiver = receivers[slope_sort]
    steepest_receiver[drains_to_self] = node_id[drains_to_self]

    # Optionally, handle baselevel nodes: they are their own receivers
    if baselevel_nodes is not None:
        receivers[baselevel_nodes, 0] = node_id[baselevel_nodes]
        receivers[baselevel_nodes, 1:] = -1
        proportions[baselevel_nodes, 0] = 1
        proportions[baselevel_nodes, 1:] = 0
        receiver_links[baselevel_nodes, :] = UNDEFINED_INDEX
        steepest_slope[baselevel_nodes] = 0.0

    # ensure that if there is a -1, it is in the second column.
    order_reversed = receivers[:, 0] == -1

    receivers_out = receivers.copy()
    receivers_out[order_reversed, 1] = receivers[order_reversed, 0]
    receivers_out[order_reversed, 0] = receivers[order_reversed, 1]

    proportions_out = proportions.copy()
    proportions_out[order_reversed, 1] = proportions[order_reversed, 0]
    proportions_out[order_reversed, 0] = proportions[order_reversed, 1]

    slopes_to_receivers_out = slopes_to_receivers.copy()
    slopes_to_receivers_out[order_reversed, 1] = slopes_to_receivers[order_reversed, 0]
    slopes_to_receivers_out[order_reversed, 0] = slopes_to_receivers[order_reversed, 1]

    receiver_links_out = receiver_links.copy()
    receiver_links_out[order_reversed, 1] = receiver_links[order_reversed, 0]
    receiver_links_out[order_reversed, 0] = receiver_links[order_reversed, 1]

    # The sink nodes are those that are their own receivers (this will normally
    # include boundary nodes as well as interior ones; "pits" would be sink
    # nodes that are also interior nodes).
    (sink,) = np.where(node_id == receivers[:, 0])
    sink = as_id_array(sink)

    return (
        receivers_out,
        proportions_out,
        slopes_to_receivers_out,
        steepest_slope,
        steepest_receiver,
        sink,
        receiver_links_out,
        steepest_link,
    )


if __name__ == "__main__":  # pragma: no cover
    import doctest

    doctest.testmod()
