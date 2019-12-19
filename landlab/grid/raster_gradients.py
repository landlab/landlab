#! /usr/bin/env python
"""Calculate gradients on a raster grid.

Gradient calculators for raster grids
+++++++++++++++++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.raster_gradients.calc_grad_at_link
    ~landlab.grid.raster_gradients.calc_grad_at_active_link
    ~landlab.grid.raster_gradients.calc_grad_across_cell_faces
    ~landlab.grid.raster_gradients.calc_grad_across_cell_corners
    ~landlab.grid.raster_gradients.alculate_gradient_along_node_links

"""
from collections import deque

import numpy as np

from landlab.core.utils import make_optional_arg_into_id_array, radians_to_degrees
from landlab.grid import gradients
from landlab.grid.base import BAD_INDEX_VALUE, CLOSED_BOUNDARY
from landlab.utils.decorators import use_field_name_or_array


@use_field_name_or_array("node")
def calc_grad_at_link(grid, node_values, out=None):
    """Calculate gradients in node_values at links.

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.
    node_values : array_like or field name
        Values at nodes.
    out : ndarray, optional
        Buffer to hold result. If `None`, create a new array.

    Returns
    -------
    ndarray
        Gradients of the nodes values for each link.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((3, 3))
    >>> node_values = [0., 0., 0.,
    ...                1., 3., 1.,
    ...                2., 2., 2.]
    >>> grid.calc_grad_at_link(node_values)
    array([ 0.,  0.,  1.,  3.,  1.,  2., -2.,  1., -1.,  1.,  0.,  0.])

    >>> out = np.empty(grid.number_of_links, dtype=float)
    >>> rtn = grid.calc_grad_at_link(node_values, out=out)
    >>> rtn is out
    True
    >>> out
    array([ 0.,  0.,  1.,  3.,  1.,  2., -2.,  1., -1.,  1.,  0.,  0.])

    >>> grid = RasterModelGrid((3, 3), xy_spacing=(2, 1))
    >>> grid.calc_grad_at_link(node_values)
    array([ 0.,  0.,  1.,  3.,  1.,  1., -1.,  1., -1.,  1.,  0.,  0.])
    >>> _ = grid.add_field('node', 'elevation', node_values)
    >>> grid.calc_grad_at_link('elevation')
    array([ 0.,  0.,  1.,  3.,  1.,  1., -1.,  1., -1.,  1.,  0.,  0.])

    LLCATS: LINF GRAD
    """
    grads = gradients.calc_diff_at_link(grid, node_values, out=out)
    grads /= grid.length_of_link[: grid.number_of_links]

    #    n_vertical_links = (grid.shape[0] - 1) * grid.shape[1]
    #    diffs[:n_vertical_links] /= grid.dy
    #    diffs[n_vertical_links:] /= grid.dx

    return grads


@use_field_name_or_array("node")
def calc_grad_at_active_link(grid, node_values, out=None):
    """Calculate gradients over active links.

    .. deprecated:: 0.1
        Use :func:`calc_grad_across_cell_faces`
                or :func:`calc_grad_across_cell_corners` instead

    Calculates the gradient in quantity s at each active link in the grid.
    This is nearly identical to the method of the same name in ModelGrid,
    except that it uses a constant node spacing for link length to improve
    efficiency.

    Note that a negative gradient corresponds to a lower node in the
    direction of the link.

    Returns
    -------
    ndarray
        Gradients of the nodes values for each link.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((4, 5), xy_spacing=1.0)
    >>> u = [0., 1., 2., 3., 0.,
    ...      1., 2., 3., 2., 3.,
    ...      0., 1., 2., 1., 2.,
    ...      0., 0., 2., 2., 0.]
    >>> grad = grid.calc_grad_at_active_link(u)
    >>> grad # doctest: +NORMALIZE_WHITESPACE
    array([ 1.,  1., -1.,
            1.,  1., -1.,  1.,
           -1., -1., -1.,
            1.,  1., -1.,  1.,
           -1.,  0.,  1.])

    For greater speed, sending a pre-created numpy array as an argument
    avoids having to create a new one with each call:

    >>> grad = np.empty(grid.number_of_active_links)
    >>> rtn = grid.calc_grad_at_active_link(u, out=grad)
    >>> grad # doctest: +NORMALIZE_WHITESPACE
    array([ 1.,  1., -1.,
            1.,  1., -1.,  1.,
           -1., -1., -1.,
            1.,  1., -1.,  1.,
           -1.,  0.,  1.])
    >>> rtn is grad
    True

    >>> grid = RasterModelGrid((3, 3), xy_spacing=(2, 1))
    >>> node_values = [0., 0., 0.,
    ...                1., 3., 1.,
    ...                2., 2., 2.]
    >>> grid.calc_grad_at_active_link(node_values)
    array([ 3.,  1., -1., -1.])

    This function is *deprecated*. Instead, use ``calc_grad_at_link``.

    >>> grid = RasterModelGrid((3, 3), xy_spacing=(2, 1))
    >>> node_values = [0., 0., 0.,
    ...                1., 3., 1.,
    ...                2., 2., 2.]
    >>> grid.calc_grad_at_link(node_values)[grid.active_links]
    array([ 3.,  1., -1., -1.])

    LLCATS: LINF GRAD
    """
    if out is None:
        out = np.empty(len(grid.active_links), dtype=float)

    if len(out) != len(grid.active_links):
        raise ValueError("output buffer does not match that of the grid.")

    # grads = gradients.calc_diff_at_link(grid, node_values,
    #                                                  out=out)
    grads = gradients.calc_diff_at_link(grid, node_values)
    out[:] = grads[grid.active_links]
    out /= grid.length_of_link[grid.active_links]

    return out


@use_field_name_or_array("node")
def calc_grad_across_cell_faces(grid, node_values, *args, **kwds):
    """calc_grad_across_cell_faces(grid, node_values, [cell_ids], out=None)

    Get gradients across the faces of a cell.

    Calculate gradient of the value field provided by *node_values* across
    each of the faces of the cells of a grid. The returned gradients are
    ordered as right, top, left, and bottom.

    Note that the returned gradients are masked to exclude neighbor nodes which
    are closed. Beneath the mask is the value -1.

    Parameters
    ----------
    grid : RasterModelGrid
        Source grid.
    node_values : array_like or field name
        Quantity to take the gradient of defined at each node.
    cell_ids : array_like, optional
        If provided, cell ids to measure gradients. Otherwise, find gradients
        for all cells.
    out : array_like, optional
        Alternative output array in which to place the result.  Must
        be of the same shape and buffer length as the expected output.

    Returns
    -------
    (N, 4) Masked ndarray
        Gradients for each face of the cell.

    Examples
    --------
    Create a grid with two cells.

    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((3, 4))
    >>> x = np.array([0., 0., 0., 0.,
    ...               0., 0., 1., 1.,
    ...               3., 3., 3., 3.])

    A decrease in quantity across a face is a negative gradient.

    >>> grid.calc_grad_across_cell_faces(x) # doctest: +NORMALIZE_WHITESPACE
    masked_array(data =
     [[ 1.  3.  0.  0.]
     [ 0.  2. -1. -1.]],
                 mask =
     False,
           fill_value = 1e+20)

    >>> grid = RasterModelGrid((3, 4), xy_spacing=(1, 2))
    >>> grid.calc_grad_across_cell_faces(x) # doctest: +NORMALIZE_WHITESPACE
    masked_array(data =
     [[ 1.   1.5  0.   0. ]
     [ 0.   1.  -1.  -0.5]],
                  mask =
     False,
           fill_value = 1e+20)

    LLCATS: FINF GRAD
    """
    padded_node_values = np.empty(node_values.size + 1, dtype=float)
    padded_node_values[-1] = BAD_INDEX_VALUE
    padded_node_values[:-1] = node_values
    cell_ids = make_optional_arg_into_id_array(grid.number_of_cells, *args)
    node_ids = grid.node_at_cell[cell_ids]

    neighbors = grid.active_adjacent_nodes_at_node[node_ids]
    if BAD_INDEX_VALUE != -1:
        neighbors = np.where(neighbors == BAD_INDEX_VALUE, -1, neighbors)
    values_at_neighbors = padded_node_values[neighbors]
    masked_neighbor_values = np.ma.array(
        values_at_neighbors, mask=neighbors == BAD_INDEX_VALUE
    )
    values_at_nodes = node_values[node_ids].reshape(len(node_ids), 1)

    out = np.subtract(masked_neighbor_values, values_at_nodes, **kwds)

    out[:, (0, 2)] /= grid.dx
    out[:, (1, 3)] /= grid.dy

    return out


@use_field_name_or_array("node")
def calc_grad_across_cell_corners(grid, node_values, *args, **kwds):
    """calc_grad_across_cell_corners(grid, node_values, [cell_ids], out=None)

    Get gradients to diagonally opposite nodes.

    Calculate gradient of the value field provided by *node_values* to
    the values at diagonally opposite nodes. The returned gradients are
    ordered as upper-right, upper-left, lower-left and lower-right.

    Parameters
    ----------
    grid : RasterModelGrid
        Source grid.
    node_values : array_like or field name
        Quantity to take the gradient of defined at each node.
    cell_ids : array_like, optional
        If provided, cell ids to measure gradients. Otherwise, find gradients
        for all cells.
    out : array_like, optional
        Alternative output array in which to place the result.  Must
        be of the same shape and buffer length as the expected output.

    Returns
    -------
    (N, 4) ndarray
        Gradients to each diagonal node.

    Examples
    --------
    Create a grid with two cells.

    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((3, 4))
    >>> x = np.array([1., 0., 0., 1.,
    ...               0., 0., 1., 1.,
    ...               3., 3., 3., 3.])

    A decrease in quantity to a diagonal node is a negative gradient.

    >>> from math import sqrt
    >>> grid.calc_grad_across_cell_corners(x) * sqrt(2.)
    array([[ 3.,  3.,  1.,  0.],
           [ 2.,  2., -1.,  0.]])

    >>> grid = RasterModelGrid((3, 4), xy_spacing=(4, 3))
    >>> grid.calc_grad_across_cell_corners(x)
    array([[ 0.6,  0.6,  0.2,  0. ],
           [ 0.4,  0.4, -0.2,  0. ]])

    LLCATS: CNINF GRAD
    """
    cell_ids = make_optional_arg_into_id_array(grid.number_of_cells, *args)
    node_ids = grid.node_at_cell[cell_ids]

    values_at_diagonals = node_values[grid.diagonal_adjacent_nodes_at_node[node_ids]]
    values_at_nodes = node_values[node_ids].reshape(len(node_ids), 1)

    out = np.subtract(values_at_diagonals, values_at_nodes, **kwds)
    np.divide(out, np.sqrt(grid.dy ** 2.0 + grid.dx ** 2.0), out=out)

    return out


@use_field_name_or_array("node")
def calc_grad_along_node_links(grid, node_values, *args, **kwds):
    """calc_grad_along_node_links(grid, node_values, [cell_ids], out=None)

    Get gradients along links touching a node.

    Calculate gradient of the value field provided by *node_values* across
    each of the faces of the nodes of a grid. The returned gradients are
    ordered as right, top, left, and bottom. All returned values follow our
    standard sign convention, where a link pointing N or E and increasing in
    value is positive, a link pointing S or W and increasing in value is
    negative.

    Note that the returned gradients are masked to exclude neighbor nodes which
    are closed. Beneath the mask is the value -1.

    Parameters
    ----------
    grid : RasterModelGrid
        Source grid.
    node_values : array_like or field name
        Quantity to take the gradient of defined at each node.
    node_ids : array_like, optional
        If provided, node ids to measure gradients. Otherwise, find gradients
        for all nodes.
    out : array_like, optional
        Alternative output array in which to place the result.  Must
        be of the same shape and buffer length as the expected output.

    Returns
    -------
    (N, 4) Masked ndarray
        Gradients for each link of the node. Ordering is E,N,W,S.

    Examples
    --------
    Create a grid with nine nodes.

    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((3, 3))
    >>> x = np.array([0., 0., 0.,
    ...               0., 1., 2.,
    ...               2., 2., 2.])

    A decrease in quantity across a face is a negative gradient.

    >>> grid.calc_grad_along_node_links(x) # doctest: +NORMALIZE_WHITESPACE
    masked_array(data =
     [[-- -- -- --]
     [-- 1.0 -- --]
     [-- -- -- --]
     [1.0 -- -- --]
     [1.0 1.0 1.0 1.0]
     [-- -- 1.0 --]
     [-- -- -- --]
     [-- -- -- 1.0]
     [-- -- -- --]],
                 mask =
     [[ True  True  True  True]
     [ True False  True  True]
     [ True  True  True  True]
     [False  True  True  True]
     [False False False False]
     [ True  True False  True]
     [ True  True  True  True]
     [ True  True  True False]
     [ True  True  True  True]],
           fill_value = 1e+20)

    >>> grid = RasterModelGrid((3, 3), xy_spacing=(4, 2))
    >>> grid.calc_grad_along_node_links(x) # doctest: +NORMALIZE_WHITESPACE
    masked_array(data =
     [[-- -- -- --]
     [-- 0.5 -- --]
     [-- -- -- --]
     [0.25 -- -- --]
     [0.25 0.5 0.25 0.5]
     [-- -- 0.25 --]
     [-- -- -- --]
     [-- -- -- 0.5]
     [-- -- -- --]],
                 mask =
     [[ True  True  True  True]
     [ True False  True  True]
     [ True  True  True  True]
     [False  True  True  True]
     [False False False False]
     [ True  True False  True]
     [ True  True  True  True]
     [ True  True  True False]
     [ True  True  True  True]],
           fill_value = 1e+20)

    LLCATS: NINF LINF GRAD
    """
    padded_node_values = np.empty(node_values.size + 1, dtype=float)
    padded_node_values[-1] = BAD_INDEX_VALUE
    padded_node_values[:-1] = node_values
    node_ids = make_optional_arg_into_id_array(grid.number_of_nodes, *args)

    neighbors = grid.active_adjacent_nodes_at_node[node_ids]
    values_at_neighbors = padded_node_values[neighbors]
    masked_neighbor_values = np.ma.array(
        values_at_neighbors, mask=values_at_neighbors == BAD_INDEX_VALUE
    )
    values_at_nodes = node_values[node_ids].reshape(len(node_ids), 1)

    out = np.ma.empty_like(masked_neighbor_values, dtype=float)
    np.subtract(masked_neighbor_values[:, :2], values_at_nodes, out=out[:, :2], **kwds)
    np.subtract(values_at_nodes, masked_neighbor_values[:, 2:], out=out[:, 2:], **kwds)

    out[:, (0, 2)] /= grid.dx
    out[:, (1, 3)] /= grid.dy

    return out


def calc_unit_normals_at_cell_subtriangles(grid, elevs="topographic__elevation"):
    """Calculate unit normals on a cell.

    Calculate the eight unit normal vectors <a, b, c> to the eight
    subtriangles of a four-cornered (raster) cell.

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.
    elevs : str or ndarray, optional
        Field name or array of node values.

    Returns
    -------
    (n_ENE, n_NNE, n_NNW, n_WNW, n_WSW, n_SSW, n_SSE, n_ESE) :
        each a num-cells x length-3 array
        Len-8 tuple of the eight unit normal vectors <a, b, c> for the eight
        subtriangles in the cell. Order is from north of east, counter
        clockwise to south of east (East North East, North North East, North
        North West, West North West, West South West, South South West, South
        South East, East South East).

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> mg = RasterModelGrid((3, 3))
    >>> z = mg.node_x ** 2
    >>> eight_tris = mg.calc_unit_normals_at_cell_subtriangles(z)
    >>> type(eight_tris) is tuple
    True
    >>> len(eight_tris)
    8
    >>> eight_tris[0].shape == (mg.number_of_cells, 3)
    True
    >>> eight_tris # doctest: +NORMALIZE_WHITESPACE
    (array([[-0.9486833 ,  0.        ,  0.31622777]]),
     array([[-0.9486833 ,  0.        ,  0.31622777]]),
     array([[-0.70710678,  0.        ,  0.70710678]]),
     array([[-0.70710678,  0.        ,  0.70710678]]),
     array([[-0.70710678,  0.        ,  0.70710678]]),
     array([[-0.70710678,  0.        ,  0.70710678]]),
     array([[-0.9486833 ,  0.        ,  0.31622777]]),
     array([[-0.9486833 ,  0.        ,  0.31622777]]))

    LLCATS: CINF GRAD
    """

    # identify the grid neigbors at each location
    node_at_cell = grid.node_at_cell
    # calculate unit normals at all nodes.
    (
        n_ENE,
        n_NNE,
        n_NNW,
        n_WNW,
        n_WSW,
        n_SSW,
        n_SSE,
        n_ESE,
    ) = _calc_subtriangle_unit_normals_at_node(grid, elevs=elevs)

    # return only those at cell.
    return (
        n_ENE[node_at_cell, :],
        n_NNE[node_at_cell, :],
        n_NNW[node_at_cell, :],
        n_WNW[node_at_cell, :],
        n_WSW[node_at_cell, :],
        n_SSW[node_at_cell, :],
        n_SSE[node_at_cell, :],
        n_ESE[node_at_cell, :],
    )


def _calc_subtriangle_unit_normals_at_node(grid, elevs="topographic__elevation"):
    """Private Function: Calculate unit normals on subtriangles at all nodes.

    Calculate the eight unit normal vectors <a, b, c> to the eight
    subtriangles of a four-cornered (raster) cell. Unlike
    calc_unit_normals_at_node_subtriangles, this function also
    calculated unit normals at the degenerate part-cells around the
    boundary.

    On the grid boundaries where the cell is not fully defined, the unit normal
    is given as <nan, nan, nan>.

    Parameters
    ----------
    grid : RasterModelGrid
    A grid.
    elevs : str or ndarray, optional
    Field name or array of node values.

    Returns
    -------
    (n_ENE, n_NNE, n_NNW, n_WNW, n_WSW, n_SSW, n_SSE, n_ESE) :
    each a num-nodes x length-3 array
    Len-8 tuple of the eight unit normal vectors <a, b, c> for the eight
    subtriangles in the cell. Order is from north of east, counter
    clockwise to south of east (East North East, North North East, North
    North West, West North West, West South West, South South West, South
    South East, East South East).

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.raster_gradients import(
    ...       _calc_subtriangle_unit_normals_at_node)
    >>> mg = RasterModelGrid((3, 3))
    >>> z = mg.node_x ** 2
    >>> eight_tris = _calc_subtriangle_unit_normals_at_node(mg, z)
    >>> type(eight_tris) is tuple
    True
    >>> len(eight_tris)
    8
    >>> eight_tris[0].shape == (mg.number_of_nodes, 3)
    True
    >>> eight_tris[0] # doctest: +NORMALIZE_WHITESPACE
    array([[-0.70710678,  0.        ,  0.70710678],
               [-0.9486833 ,  0.        ,  0.31622777],
               [        nan,         nan,         nan],
               [-0.70710678,  0.        ,  0.70710678],
               [-0.9486833 ,  0.        ,  0.31622777],
               [        nan,         nan,         nan],
               [        nan,         nan,         nan],
               [        nan,         nan,         nan],
               [        nan,         nan,         nan]])

    LLCATS: CINF GRAD
    """

    try:
        z = grid.at_node[elevs]
    except TypeError:
        z = elevs
    #  cell has center node I
    # orthogonal neighbors P, R, T, V, counter clockwise from East
    # diagonal neihbors Q, S, U, W, counter clocwise from North East
    # There are 8 subtriangles that can be defined  with the following corners
    # (starting from the central node, and progressing counter-clockwise).
    # ENE: IPQ
    # NNE: IQR
    # NNW: IRS
    # WNW: IST
    # WSW: ITU
    # SSW: IUV
    # SSE: IVW
    # ESE: IWP

    # There are thus 8 vectors, IP, IQ, IR, IS, IT, IU, IV, IW

    # initialized difference matricies for cross product
    diff_xyz_IP = np.empty((grid.number_of_nodes, 3))  # East
    # ^this is the vector (xP-xI, yP-yI, zP-yI)
    diff_xyz_IQ = np.empty((grid.number_of_nodes, 3))  # Northeast
    diff_xyz_IR = np.empty((grid.number_of_nodes, 3))  # North
    diff_xyz_IS = np.empty((grid.number_of_nodes, 3))  # Northwest
    diff_xyz_IT = np.empty((grid.number_of_nodes, 3))  # West
    diff_xyz_IU = np.empty((grid.number_of_nodes, 3))  # Southwest
    diff_xyz_IV = np.empty((grid.number_of_nodes, 3))  # South
    diff_xyz_IW = np.empty((grid.number_of_nodes, 3))  # Southeast

    # identify the grid neigbors at each location
    node_at_cell = np.arange(grid.number_of_nodes)
    P = grid.adjacent_nodes_at_node[node_at_cell, 0]
    Q = grid.diagonal_adjacent_nodes_at_node[node_at_cell, 0]
    R = grid.adjacent_nodes_at_node[node_at_cell, 1]
    S = grid.diagonal_adjacent_nodes_at_node[node_at_cell, 1]
    T = grid.adjacent_nodes_at_node[node_at_cell, 2]
    U = grid.diagonal_adjacent_nodes_at_node[node_at_cell, 2]
    V = grid.adjacent_nodes_at_node[node_at_cell, 3]
    W = grid.diagonal_adjacent_nodes_at_node[node_at_cell, 3]

    # get x, y, z coordinates for each location
    x_I = grid.node_x[node_at_cell]
    y_I = grid.node_y[node_at_cell]
    z_I = z[node_at_cell]

    x_P = grid.node_x[P]
    y_P = grid.node_y[P]
    z_P = z[P]

    x_Q = grid.node_x[Q]
    y_Q = grid.node_y[Q]
    z_Q = z[Q]

    x_R = grid.node_x[R]
    y_R = grid.node_y[R]
    z_R = z[R]

    x_S = grid.node_x[S]
    y_S = grid.node_y[S]
    z_S = z[S]

    x_T = grid.node_x[T]
    y_T = grid.node_y[T]
    z_T = z[T]

    x_U = grid.node_x[U]
    y_U = grid.node_y[U]
    z_U = z[U]

    x_V = grid.node_x[V]
    y_V = grid.node_y[V]
    z_V = z[V]

    x_W = grid.node_x[W]
    y_W = grid.node_y[W]
    z_W = z[W]

    # calculate vectors by differencing
    diff_xyz_IP[:, 0] = x_P - x_I
    diff_xyz_IP[:, 1] = y_P - y_I
    diff_xyz_IP[:, 2] = z_P - z_I

    diff_xyz_IQ[:, 0] = x_Q - x_I
    diff_xyz_IQ[:, 1] = y_Q - y_I
    diff_xyz_IQ[:, 2] = z_Q - z_I

    diff_xyz_IR[:, 0] = x_R - x_I
    diff_xyz_IR[:, 1] = y_R - y_I
    diff_xyz_IR[:, 2] = z_R - z_I

    diff_xyz_IS[:, 0] = x_S - x_I
    diff_xyz_IS[:, 1] = y_S - y_I
    diff_xyz_IS[:, 2] = z_S - z_I

    diff_xyz_IT[:, 0] = x_T - x_I
    diff_xyz_IT[:, 1] = y_T - y_I
    diff_xyz_IT[:, 2] = z_T - z_I

    diff_xyz_IU[:, 0] = x_U - x_I
    diff_xyz_IU[:, 1] = y_U - y_I
    diff_xyz_IU[:, 2] = z_U - z_I

    diff_xyz_IV[:, 0] = x_V - x_I
    diff_xyz_IV[:, 1] = y_V - y_I
    diff_xyz_IV[:, 2] = z_V - z_I

    diff_xyz_IW[:, 0] = x_W - x_I
    diff_xyz_IW[:, 1] = y_W - y_I
    diff_xyz_IW[:, 2] = z_W - z_I

    # calculate cross product to get unit normal
    # cross product is orthogonal to both vectors, and is the normal
    # n = <a, b, c>, where plane is ax + by + cz = d
    nhat_ENE = np.cross(diff_xyz_IP, diff_xyz_IQ)  # <a, b, c>
    nhat_NNE = np.cross(diff_xyz_IQ, diff_xyz_IR)
    nhat_NNW = np.cross(diff_xyz_IR, diff_xyz_IS)
    nhat_WNW = np.cross(diff_xyz_IS, diff_xyz_IT)
    nhat_WSW = np.cross(diff_xyz_IT, diff_xyz_IU)
    nhat_SSW = np.cross(diff_xyz_IU, diff_xyz_IV)
    nhat_SSE = np.cross(diff_xyz_IV, diff_xyz_IW)
    nhat_ESE = np.cross(diff_xyz_IW, diff_xyz_IP)

    # now remove the bad subtriangles based on parts of the grid
    # make the bad subtriangle of length greater than one.
    bad = np.nan * np.ones((3,))

    # first, corners:
    corners = grid.nodes_at_corners_of_grid
    # lower left corner only has NNE and ENE
    nhat_NNW[corners[0], :] = bad
    nhat_WNW[corners[0], :] = bad
    nhat_WSW[corners[0], :] = bad
    nhat_SSW[corners[0], :] = bad
    nhat_SSE[corners[0], :] = bad
    nhat_ESE[corners[0], :] = bad

    # lower right corner only has NNW and WNW
    nhat_ENE[corners[1], :] = bad
    nhat_NNE[corners[1], :] = bad
    nhat_WSW[corners[1], :] = bad
    nhat_SSW[corners[1], :] = bad
    nhat_SSE[corners[1], :] = bad
    nhat_ESE[corners[1], :] = bad

    # upper left corner only has ESE and SSE
    nhat_ENE[corners[2], :] = bad
    nhat_NNE[corners[2], :] = bad
    nhat_NNW[corners[2], :] = bad
    nhat_WNW[corners[2], :] = bad
    nhat_WSW[corners[2], :] = bad
    nhat_SSW[corners[2], :] = bad

    # upper right corner only has WSW and SSW
    nhat_ENE[corners[3], :] = bad
    nhat_NNE[corners[3], :] = bad
    nhat_NNW[corners[3], :] = bad
    nhat_WNW[corners[3], :] = bad
    nhat_SSE[corners[3], :] = bad
    nhat_ESE[corners[3], :] = bad

    # next, sizes:
    # bottom row only has Norths
    bottom = grid.nodes_at_bottom_edge
    nhat_WSW[bottom, :] = bad
    nhat_SSW[bottom, :] = bad
    nhat_SSE[bottom, :] = bad
    nhat_ESE[bottom, :] = bad

    # left side only has Easts
    left = grid.nodes_at_left_edge
    nhat_NNW[left, :] = bad
    nhat_WNW[left, :] = bad
    nhat_WSW[left, :] = bad
    nhat_SSW[left, :] = bad

    # top row only has Souths
    top = grid.nodes_at_top_edge
    nhat_ENE[top, :] = bad
    nhat_NNE[top, :] = bad
    nhat_NNW[top, :] = bad
    nhat_WNW[top, :] = bad

    # right side only has Wests
    right = grid.nodes_at_right_edge
    nhat_ENE[right, :] = bad
    nhat_NNE[right, :] = bad
    nhat_SSE[right, :] = bad
    nhat_ESE[right, :] = bad

    # calculate magnitude of cross product so that the result is a unit normal
    nmag_ENE = np.sqrt(np.square(nhat_ENE).sum(axis=1))
    nmag_NNE = np.sqrt(np.square(nhat_NNE).sum(axis=1))
    nmag_NNW = np.sqrt(np.square(nhat_NNW).sum(axis=1))
    nmag_WNW = np.sqrt(np.square(nhat_WNW).sum(axis=1))
    nmag_WSW = np.sqrt(np.square(nhat_WSW).sum(axis=1))
    nmag_SSW = np.sqrt(np.square(nhat_SSW).sum(axis=1))
    nmag_SSE = np.sqrt(np.square(nhat_SSE).sum(axis=1))
    nmag_ESE = np.sqrt(np.square(nhat_ESE).sum(axis=1))

    # normalize the cross product with its magnitude so it is a unit normal
    # instead of a variable length normal.
    n_ENE = nhat_ENE / nmag_ENE.reshape(grid.number_of_nodes, 1)
    n_NNE = nhat_NNE / nmag_NNE.reshape(grid.number_of_nodes, 1)
    n_NNW = nhat_NNW / nmag_NNW.reshape(grid.number_of_nodes, 1)
    n_WNW = nhat_WNW / nmag_WNW.reshape(grid.number_of_nodes, 1)
    n_WSW = nhat_WSW / nmag_WSW.reshape(grid.number_of_nodes, 1)
    n_SSW = nhat_SSW / nmag_SSW.reshape(grid.number_of_nodes, 1)
    n_SSE = nhat_SSE / nmag_SSE.reshape(grid.number_of_nodes, 1)
    n_ESE = nhat_ESE / nmag_ESE.reshape(grid.number_of_nodes, 1)

    return (n_ENE, n_NNE, n_NNW, n_WNW, n_WSW, n_SSW, n_SSE, n_ESE)


def calc_slope_at_cell_subtriangles(
    grid, elevs="topographic__elevation", subtriangle_unit_normals=None
):
    """Calculate the slope (positive magnitude of gradient) at each of the
    eight cell subtriangles.

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.
    elevs : str or ndarray, optional
        Field name or array of node values.
    subtriangle_unit_normals : tuple of 8 (ncells, 3) arrays (optional)
        The unit normal vectors for the eight subtriangles of each cell,
        if already known. Order is from north of east, counter
        clockwise to south of east (East North East, North North East, North
        North West, West North West, West South West, South South West, South
        South East, East South East).

    Returns
    -------
    (s_ENE, s_NNE, s_NNW, s_WNW, s_WSW, s_SSW, s_SSE, s_ESE) :
        each a length num-cells array
        Len-8 tuple of the slopes (positive gradient magnitude) of each of the
        eight cell subtriangles, in radians. Order is from north of east,
        counter clockwise to south of east (East North East, North North East,
        North North West, West North West, West South West, South South West,
        South South East, East South East).

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> mg = RasterModelGrid((3, 3))
    >>> z = np.array([np.sqrt(3.), 0., 4./3.,
    ...               0., 0., 0.,
    ...               1., 0., 1./np.sqrt(3.)])
    >>> eight_tris = mg.calc_unit_normals_at_cell_subtriangles(z)
    >>> S = mg.calc_slope_at_cell_subtriangles(z, eight_tris)
    >>> S0 = mg.calc_slope_at_cell_subtriangles(z)
    >>> np.allclose(S, S0)
    True
    >>> type(S) is tuple
    True
    >>> len(S)
    8
    >>> len(S[0]) == mg.number_of_cells
    True
    >>> np.allclose(S[0], S[1])
    True
    >>> np.allclose(S[2], S[3])
    True
    >>> np.allclose(S[4], S[5])
    True
    >>> np.allclose(S[6], S[7])
    True
    >>> np.allclose(np.rad2deg(S[0])[0], 30.)
    True
    >>> np.allclose(np.rad2deg(S[2])[0], 45.)
    True
    >>> np.allclose(np.rad2deg(S[4])[0], 60.)
    True
    >>> np.allclose(np.cos(S[6])[0], 3./5.)
    True

    LLCATS: CINF GRAD
    """

    # calculate all subtriangle slopes
    (
        s_ENE,
        s_NNE,
        s_NNW,
        s_WNW,
        s_WSW,
        s_SSW,
        s_SSE,
        s_ESE,
    ) = _calc_subtriangle_slopes_at_node(
        grid, elevs=elevs, subtriangle_unit_normals=subtriangle_unit_normals
    )
    # return only those at cell
    if s_ENE.shape[0] == grid.number_of_nodes:
        node_at_cell = grid.node_at_cell
    else:
        node_at_cell = np.arange(grid.number_of_cells)

    return (
        s_ENE[node_at_cell],
        s_NNE[node_at_cell],
        s_NNW[node_at_cell],
        s_WNW[node_at_cell],
        s_WSW[node_at_cell],
        s_SSW[node_at_cell],
        s_SSE[node_at_cell],
        s_ESE[node_at_cell],
    )


def _calc_subtriangle_slopes_at_node(
    grid, elevs="topographic__elevation", subtriangle_unit_normals=None
):
    """Private Function: Calculate subtriangles slope at all nodes.

    Calculate the slope (positive magnitude of gradient) at each of the
    eight subtriangles, including those at not-full cells along the
    boundary.

    Those subtriangles that that don't exist because they are on the edge
    of the grid have slopes of NAN.

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.
    elevs : str or ndarray, optional
        Field name or array of node values.
    subtriangle_unit_normals : tuple of 8 (ncells, 3) or (nnodes, 3) arrays
        (optional)
        The unit normal vectors for the eight subtriangles of each cell or
        node,if already known. Order is from north of east, counter
        clockwise to south of east (East North East, North North East, North
        North West, West North West, West South West, South South West, South
        South East, East South East).

    Returns
    -------
    (s_ENE, s_NNE, s_NNW, s_WNW, s_WSW, s_SSW, s_SSE, s_ESE) :
        each a length num-cells array
        Len-8 tuple of the slopes (positive gradient magnitude) of each of the
        eight cell subtriangles, in radians. Order is from north of east,
        counter clockwise to south of east (East North East, North North East,
        North North West, West North West, West South West, South South West,
        South South East, East South East).

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.raster_gradients import(
    ...                         _calc_subtriangle_unit_normals_at_node,
    ...                         _calc_subtriangle_slopes_at_node)
    >>> mg = RasterModelGrid((3, 3))
    >>> z = np.array([np.sqrt(3.), 0., 4./3.,
    ...               0., 0., 0.,
    ...               1., 0., 1./np.sqrt(3.)])
    >>> eight_tris = _calc_subtriangle_unit_normals_at_node(mg, z)
    >>> S = _calc_subtriangle_slopes_at_node(mg, z, eight_tris)
    >>> S0 = _calc_subtriangle_slopes_at_node(mg, z)
    >>> np.allclose(S, S0, equal_nan=True)
    True
    >>> type(S) is tuple
    True
    >>> len(S)
    8
    >>> len(S[0]) == mg.number_of_nodes
    True
    >>> np.allclose(S[0][mg.core_nodes], S[1][mg.core_nodes])
    True
    >>> np.allclose(S[2][mg.core_nodes], S[3][mg.core_nodes])
    True
    >>> np.allclose(S[4][mg.core_nodes], S[5][mg.core_nodes])
    True
    >>> np.allclose(S[6][mg.core_nodes], S[7][mg.core_nodes])
    True
    >>> np.allclose(np.rad2deg(S[0][mg.core_nodes]), 30.)
    True
    >>> np.allclose(np.rad2deg(S[2][mg.core_nodes]), 45.)
    True
    >>> np.allclose(np.rad2deg(S[4])[mg.core_nodes], 60.)
    True
    >>> np.allclose(np.cos(S[6])[mg.core_nodes], 3./5.)
    True

    LLCATS: CINF GRAD
    """

    # verify that subtriangle_unit_normals is of the correct form.
    if subtriangle_unit_normals is not None:
        assert len(subtriangle_unit_normals) == 8
        assert subtriangle_unit_normals[0].shape[1] == 3
        assert subtriangle_unit_normals[1].shape[1] == 3
        assert subtriangle_unit_normals[2].shape[1] == 3
        assert subtriangle_unit_normals[3].shape[1] == 3
        assert subtriangle_unit_normals[4].shape[1] == 3
        assert subtriangle_unit_normals[5].shape[1] == 3
        assert subtriangle_unit_normals[6].shape[1] == 3
        assert subtriangle_unit_normals[7].shape[1] == 3
        (
            n_ENE,
            n_NNE,
            n_NNW,
            n_WNW,
            n_WSW,
            n_SSW,
            n_SSE,
            n_ESE,
        ) = subtriangle_unit_normals

        if subtriangle_unit_normals[7].shape[0] == grid.number_of_nodes:
            reshape_size = grid.number_of_nodes
        elif subtriangle_unit_normals[7].shape[0] == grid.number_of_cells:
            reshape_size = grid.number_of_cells
        else:
            ValueError("Subtriangles must be of lenght nnodes or ncells")
    else:
        n_ENE, n_NNE, n_NNW, n_WNW, n_WSW, n_SSW, n_SSE, n_ESE = _calc_subtriangle_unit_normals_at_node(
            grid, elevs
        )
        reshape_size = grid.number_of_nodes

    # combine z direction element of all eight so that the arccosine portion
    # only takes one function call.
    dotprod = np.empty((reshape_size, 8))
    dotprod[:, 0] = n_ENE[:, 2]  # by definition
    dotprod[:, 1] = n_NNE[:, 2]
    dotprod[:, 2] = n_NNW[:, 2]
    dotprod[:, 3] = n_WNW[:, 2]
    dotprod[:, 4] = n_WSW[:, 2]
    dotprod[:, 5] = n_SSW[:, 2]
    dotprod[:, 6] = n_SSE[:, 2]
    dotprod[:, 7] = n_ESE[:, 2]

    # take the inverse cosine of the z component to get the slope angle
    slopes_at_cell_subtriangles = np.arccos(dotprod)  #

    # split array into each subtriangle component.
    s_ENE = slopes_at_cell_subtriangles[:, 0].reshape(reshape_size)
    s_NNE = slopes_at_cell_subtriangles[:, 1].reshape(reshape_size)
    s_NNW = slopes_at_cell_subtriangles[:, 2].reshape(reshape_size)
    s_WNW = slopes_at_cell_subtriangles[:, 3].reshape(reshape_size)
    s_WSW = slopes_at_cell_subtriangles[:, 4].reshape(reshape_size)
    s_SSW = slopes_at_cell_subtriangles[:, 5].reshape(reshape_size)
    s_SSE = slopes_at_cell_subtriangles[:, 6].reshape(reshape_size)
    s_ESE = slopes_at_cell_subtriangles[:, 7].reshape(reshape_size)

    return (s_ENE, s_NNE, s_NNW, s_WNW, s_WSW, s_SSW, s_SSE, s_ESE)


def calc_aspect_at_cell_subtriangles(
    grid, elevs="topographic__elevation", subtriangle_unit_normals=None, unit="degrees"
):
    """Get tuple of arrays of aspect of each of the eight cell subtriangles.

    Aspect is returned as radians clockwise of north, unless input parameter
    units is set to 'degrees'.

    If subtriangle_unit_normals is provided the aspect will be calculated from
    these data.

    If it is not, it will be derived from elevation data at the nodes,
    which can either be a string referring to a grid field (default:
    'topographic__elevation'), or an nnodes-long numpy array of the
    values themselves.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    elevs : str or array (optional)
        Node field name or node array of elevations.
        If *subtriangle_unit_normals* is not provided, must be set, but unused
        otherwise.
    subtriangle_unit_normals : tuple of 8 (ncels, 3) arrays (optional)
        The unit normal vectors for the eight subtriangles of each cell,
        if already known. Order is from north of east, counter
        clockwise to south of east (East North East, North North East, North
        North West, West North West, West South West, South South West, South
        South East, East South East).
    unit : {'degrees', 'radians'}
        Controls the unit that the aspect is returned as.

    Returns
    -------
    (a_ENE, a_NNE, a_NNW, a_WNW, a_WSW, a_SSW, a_SSE, a_ESE) :
            each a length num-cells array
        Len-8 tuple of the aspect of each of the eight cell subtriangles.
        Aspect is returned as angle clockwise of north. Units are given as
        radians unless input parameter units is set to 'degrees'.
        Order is from north of east, counter clockwise to south of east (East
        North East, North North East, North North West, West North West, West
        South West, South South West, South South East, East South East).

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> mg = RasterModelGrid((3, 3))
    >>> z = np.array([1., 0., 1., 0., 0., 0., 1., 0., 1.])
    >>> eight_tris = mg.calc_unit_normals_at_cell_subtriangles(z)
    >>> A = mg.calc_aspect_at_cell_subtriangles(z, eight_tris)
    >>> A0 = mg.calc_aspect_at_cell_subtriangles(z)
    >>> np.allclose(A, A0)
    True
    >>> type(A) is tuple
    True
    >>> len(A)
    8
    >>> len(A[0]) == mg.number_of_cells
    True
    >>> A0  # doctest: +NORMALIZE_WHITESPACE
    (array([ 180.]), array([ 270.]), array([ 90.]), array([ 180.]),
     array([ 0.]), array([ 90.]), array([ 270.]), array([ 0.]))


    LLCATS: CINF SURF
    """

    # calculate all subtriangle slopes
    (
        angle_ENE,
        angle_NNE,
        angle_NNW,
        angle_WNW,
        angle_WSW,
        angle_SSW,
        angle_SSE,
        angle_ESE,
    ) = _calc_subtriangle_aspect_at_node(
        grid, elevs=elevs, subtriangle_unit_normals=subtriangle_unit_normals, unit=unit
    )
    # return only those at cell
    if angle_ESE.shape[0] == grid.number_of_nodes:
        node_at_cell = grid.node_at_cell
    else:
        node_at_cell = np.arange(grid.number_of_cells)

    if unit == "degrees" or unit == "radians":
        return (
            angle_ENE[node_at_cell],
            angle_NNE[node_at_cell],
            angle_NNW[node_at_cell],
            angle_WNW[node_at_cell],
            angle_WSW[node_at_cell],
            angle_SSW[node_at_cell],
            angle_SSE[node_at_cell],
            angle_ESE[node_at_cell],
        )
    else:
        raise TypeError("unit must be 'degrees' or 'radians'")


def _calc_subtriangle_aspect_at_node(
    grid, elevs="topographic__elevation", subtriangle_unit_normals=None, unit="degrees"
):
    """Private Function: Aspect of subtriangles at node.

    This function calculates the aspect of all subtriangles, including those
    that are at noded without cells (on the boundaries).

    Aspect is returned as radians clockwise of north, unless input parameter
    units is set to 'degrees'.

    If subtriangle_unit_normals is provided the aspect will be calculated from
    these data.

    If it is not, it will be derived from elevation data at the nodes,
    which can either be a string referring to a grid field (default:
    'topographic__elevation'), or an nnodes-long numpy array of the
    values themselves.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    elevs : str or array (optional)
        Node field name or node array of elevations.
        If *subtriangle_unit_normals* is not provided, must be set, but unused
        otherwise.
    subtriangle_unit_normals : tuple of 8 (ncells, 3) or (nnodes, 3) arrays
        (optional)
        The unit normal vectors for the eight subtriangles of each cell or
        node,if already known. Order is from north of east, counter
        clockwise to south of east (East North East, North North East, North
        North West, West North West, West South West, South South West, South
        South East, East South East).
    unit : {'degrees', 'radians'}
        Controls the unit that the aspect is returned as.

    Returns
    -------
    (a_ENE, a_NNE, a_NNW, a_WNW, a_WSW, a_SSW, a_SSE, a_ESE) :
            each a length num-cells array
        Len-8 tuple of the aspect of each of the eight cell subtriangles.
        Aspect is returned as angle clockwise of north. Units are given as
        radians unless input parameter units is set to 'degrees'.
        Order is from north of east, counter clockwise to south of east (East
        North East, North North East, North North West, West North West, West
        South West, South South West, South South East, East South East).

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.raster_gradients import (
    ...     _calc_subtriangle_unit_normals_at_node,
    ...     _calc_subtriangle_aspect_at_node)
    >>> mg = RasterModelGrid((3, 3))
    >>> z = np.array([1., 0., 1., 0., 0., 0., 1., 0., 1.])
    >>> eight_tris = _calc_subtriangle_unit_normals_at_node(mg, z)
    >>> A = _calc_subtriangle_aspect_at_node(mg, z, eight_tris)
    >>> A0 = _calc_subtriangle_aspect_at_node(mg, z)
    >>> np.allclose(A, A0, equal_nan=True)
    True
    >>> type(A) is tuple
    True
    >>> len(A)
    8
    >>> len(A[0]) == mg.number_of_nodes
    True
    >>> A0  # doctest: +NORMALIZE_WHITESPACE
    (array([  90.,  315.,   nan,   90.,  180.,   nan,   nan,   nan,   nan]),
     array([   0.,   90.,   nan,  135.,  270.,   nan,   nan,   nan,   nan]),
     array([  nan,   90.,    0.,   nan,   90.,  225.,   nan,   nan,   nan]),
     array([  nan,   45.,  270.,   nan,  180.,   90.,   nan,   nan,   nan]),
     array([  nan,   nan,   nan,   nan,    0.,   90.,   nan,  135.,  270.]),
     array([  nan,   nan,   nan,   nan,   90.,  315.,   nan,   90.,  180.]),
     array([  nan,   nan,   nan,   45.,  270.,   nan,  180.,   90.,   nan]),
     array([  nan,   nan,   nan,  270.,    0.,   nan,   90.,  225.,   nan]))

    LLCATS: CINF SURF
    """

    # verify that subtriangle_unit_normals is of the correct form.
    if subtriangle_unit_normals is not None:
        assert len(subtriangle_unit_normals) == 8
        assert subtriangle_unit_normals[0].shape[1] == 3
        assert subtriangle_unit_normals[1].shape[1] == 3
        assert subtriangle_unit_normals[2].shape[1] == 3
        assert subtriangle_unit_normals[3].shape[1] == 3
        assert subtriangle_unit_normals[4].shape[1] == 3
        assert subtriangle_unit_normals[5].shape[1] == 3
        assert subtriangle_unit_normals[6].shape[1] == 3
        assert subtriangle_unit_normals[7].shape[1] == 3

        if subtriangle_unit_normals[7].shape[0] == grid.number_of_nodes:
            reshape_size = grid.number_of_nodes
        elif subtriangle_unit_normals[7].shape[0] == grid.number_of_cells:
            reshape_size = grid.number_of_cells
        else:
            ValueError("Subtriangles must be of lenght nnodes or ncells")
        (
            n_ENE,
            n_NNE,
            n_NNW,
            n_WNW,
            n_WSW,
            n_SSW,
            n_SSE,
            n_ESE,
        ) = subtriangle_unit_normals

    # otherwise create it.
    else:
        n_ENE, n_NNE, n_NNW, n_WNW, n_WSW, n_SSW, n_SSE, n_ESE = _calc_subtriangle_unit_normals_at_node(
            grid, elevs
        )
        reshape_size = grid.number_of_nodes
    # calculate the aspect as an angle ccw from the x axis (math angle)
    angle_from_x_ccw_ENE = np.reshape(
        np.arctan2(n_ENE[:, 1], n_ENE[:, 0]), reshape_size
    )
    angle_from_x_ccw_NNE = np.reshape(
        np.arctan2(n_NNE[:, 1], n_NNE[:, 0]), reshape_size
    )
    angle_from_x_ccw_NNW = np.reshape(
        np.arctan2(n_NNW[:, 1], n_NNW[:, 0]), reshape_size
    )
    angle_from_x_ccw_WNW = np.reshape(
        np.arctan2(n_WNW[:, 1], n_WNW[:, 0]), reshape_size
    )
    angle_from_x_ccw_WSW = np.reshape(
        np.arctan2(n_WSW[:, 1], n_WSW[:, 0]), reshape_size
    )
    angle_from_x_ccw_SSW = np.reshape(
        np.arctan2(n_SSW[:, 1], n_SSW[:, 0]), reshape_size
    )
    angle_from_x_ccw_SSE = np.reshape(
        np.arctan2(n_SSE[:, 1], n_SSE[:, 0]), reshape_size
    )
    angle_from_x_ccw_ESE = np.reshape(
        np.arctan2(n_ESE[:, 1], n_ESE[:, 0]), reshape_size
    )
    # convert reference from math angle to angles clockwise from north
    # return as either  radians or degrees depending on unit.
    if unit == "degrees":
        return (
            radians_to_degrees(angle_from_x_ccw_ENE),
            radians_to_degrees(angle_from_x_ccw_NNE),
            radians_to_degrees(angle_from_x_ccw_NNW),
            radians_to_degrees(angle_from_x_ccw_WNW),
            radians_to_degrees(angle_from_x_ccw_WSW),
            radians_to_degrees(angle_from_x_ccw_SSW),
            radians_to_degrees(angle_from_x_ccw_SSE),
            radians_to_degrees(angle_from_x_ccw_ESE),
        )

    elif unit == "radians":
        angle_from_north_cw_ENE = (5.0 * np.pi / 2.0 - angle_from_x_ccw_ENE) % (
            2.0 * np.pi
        )
        angle_from_north_cw_NNE = (5.0 * np.pi / 2.0 - angle_from_x_ccw_NNE) % (
            2.0 * np.pi
        )
        angle_from_north_cw_NNW = (5.0 * np.pi / 2.0 - angle_from_x_ccw_NNW) % (
            2.0 * np.pi
        )
        angle_from_north_cw_WNW = (5.0 * np.pi / 2.0 - angle_from_x_ccw_WNW) % (
            2.0 * np.pi
        )
        angle_from_north_cw_WSW = (5.0 * np.pi / 2.0 - angle_from_x_ccw_WSW) % (
            2.0 * np.pi
        )
        angle_from_north_cw_SSW = (5.0 * np.pi / 2.0 - angle_from_x_ccw_SSW) % (
            2.0 * np.pi
        )
        angle_from_north_cw_SSE = (5.0 * np.pi / 2.0 - angle_from_x_ccw_SSE) % (
            2.0 * np.pi
        )
        angle_from_north_cw_ESE = (5.0 * np.pi / 2.0 - angle_from_x_ccw_ESE) % (
            2.0 * np.pi
        )

        return (
            angle_from_north_cw_ENE,
            angle_from_north_cw_NNE,
            angle_from_north_cw_NNW,
            angle_from_north_cw_WNW,
            angle_from_north_cw_WSW,
            angle_from_north_cw_SSW,
            angle_from_north_cw_SSE,
            angle_from_north_cw_ESE,
        )
    else:
        raise TypeError("unit must be 'degrees' or 'radians'")


def calc_unit_normals_at_patch_subtriangles(grid, elevs="topographic__elevation"):
    """Calculate unit normals on a patch.

    Calculate the four unit normal vectors <a, b, c> to the four possible
    subtriangles of a four-cornered (raster) patch.

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.
    elevs : str or ndarray, optional
        Field name or array of node values.

    Returns
    -------
    (n_TR, n_TL, n_BL, n_BR) : each a num-patches x length-3 array
        Len-4 tuple of the four unit normal vectors <a, b, c> for the four
        possible subtriangles in the patch. Order is (topright, topleft,
        bottomleft, bottomright).

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> mg = RasterModelGrid((4, 5))
    >>> z = mg.node_x ** 2
    >>> four_tris = mg.calc_unit_normals_at_patch_subtriangles(z)
    >>> type(four_tris) is tuple
    True
    >>> len(four_tris)
    4
    >>> np.allclose(four_tris[0], four_tris[1])
    True
    >>> np.allclose(four_tris[2], four_tris[3])
    True
    >>> np.allclose(four_tris[0], four_tris[2])
    True
    >>> np.allclose(np.square(four_tris[0]).sum(axis=1), 1.)
    True
    >>> four_tris[0]
    array([[-0.70710678,  0.        ,  0.70710678],
           [-0.9486833 ,  0.        ,  0.31622777],
           [-0.98058068,  0.        ,  0.19611614],
           [-0.98994949,  0.        ,  0.14142136],
           [-0.70710678,  0.        ,  0.70710678],
           [-0.9486833 ,  0.        ,  0.31622777],
           [-0.98058068,  0.        ,  0.19611614],
           [-0.98994949,  0.        ,  0.14142136],
           [-0.70710678,  0.        ,  0.70710678],
           [-0.9486833 ,  0.        ,  0.31622777],
           [-0.98058068,  0.        ,  0.19611614],
           [-0.98994949,  0.        ,  0.14142136]])

    LLCATS: PINF GRAD
    """
    try:
        z = grid.at_node[elevs]
    except TypeError:
        z = elevs
    # conceptualize patches as TWO sets of 3 nodes
    # the corners are PQRS, CC from NE
    diff_xyz_PQ = np.empty((grid.number_of_patches, 3))  # TOP
    # ^this is the vector (xQ-xP, yQ-yP, zQ-yP)
    diff_xyz_PS = np.empty((grid.number_of_patches, 3))  # RIGHT
    # we have RS and QR implicitly in PQ and PS - but store them too
    diff_xyz_RS = np.empty((grid.number_of_patches, 3))  # BOTTOM
    diff_xyz_QR = np.empty((grid.number_of_patches, 3))  # LEFT
    P = grid.nodes_at_patch[:, 0]
    Q = grid.nodes_at_patch[:, 1]
    R = grid.nodes_at_patch[:, 2]
    S = grid.nodes_at_patch[:, 3]
    x_P = grid.node_x[P]
    y_P = grid.node_y[P]
    z_P = z[P]
    x_Q = grid.node_x[Q]
    y_Q = grid.node_y[Q]
    z_Q = z[Q]
    x_R = grid.node_x[R]
    y_R = grid.node_y[R]
    z_R = z[R]
    x_S = grid.node_x[S]
    y_S = grid.node_y[S]
    z_S = z[S]
    diff_xyz_PQ[:, 0] = x_Q - x_P
    diff_xyz_PQ[:, 1] = y_Q - y_P
    diff_xyz_PQ[:, 2] = z_Q - z_P
    diff_xyz_PS[:, 0] = x_S - x_P
    diff_xyz_PS[:, 1] = y_S - y_P
    diff_xyz_PS[:, 2] = z_S - z_P
    diff_xyz_RS[:, 0] = x_S - x_R
    diff_xyz_RS[:, 1] = y_S - y_R
    diff_xyz_RS[:, 2] = z_S - z_R
    diff_xyz_QR[:, 0] = x_R - x_Q
    diff_xyz_QR[:, 1] = y_R - y_Q
    diff_xyz_QR[:, 2] = z_R - z_Q
    # make the other ones
    # cross product is orthogonal to both vectors, and is the normal
    # n = <a, b, c>, where plane is ax + by + cz = d
    nhat_topleft = np.cross(diff_xyz_PQ, diff_xyz_QR)  # <a, b, c>
    nhat_bottomright = np.cross(diff_xyz_PS, diff_xyz_RS)
    nhat_topright = np.cross(diff_xyz_PQ, diff_xyz_PS)
    nhat_bottomleft = np.cross(diff_xyz_QR, diff_xyz_RS)
    nmag_topleft = np.sqrt(np.square(nhat_topleft).sum(axis=1))
    nmag_bottomright = np.sqrt(np.square(nhat_bottomright).sum(axis=1))
    nmag_topright = np.sqrt(np.square(nhat_topright).sum(axis=1))
    nmag_bottomleft = np.sqrt(np.square(nhat_bottomleft).sum(axis=1))
    n_TR = nhat_topright / nmag_topright.reshape(grid.number_of_patches, 1)
    n_TL = nhat_topleft / nmag_topleft.reshape(grid.number_of_patches, 1)
    n_BL = nhat_bottomleft / nmag_bottomleft.reshape(grid.number_of_patches, 1)
    n_BR = nhat_bottomright / nmag_bottomright.reshape(grid.number_of_patches, 1)

    return (n_TR, n_TL, n_BL, n_BR)


def calc_slope_at_patch(
    grid,
    elevs="topographic__elevation",
    ignore_closed_nodes=True,
    subtriangle_unit_normals=None,
):
    """Calculate the slope (positive magnitude of gradient) at raster patches.

    Returns the mean of the slopes of the four possible patch subtriangles.

    If ignore_closed_nodes is True, closed nodes do not affect slope
    calculations. If more than one closed node is present in a patch, the
    patch slope is set to zero.

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.
    elevs : str or ndarray, optional
        Field name or array of node values.
    ignore_closed_nodes : bool
        If True, do not incorporate values at closed nodes into the calc.
    subtriangle_unit_normals : tuple of 4 (npatches, 3) arrays (optional)
        The unit normal vectors for the four subtriangles of each patch,
        if already known. Order is TR, TL, BL, BR.

    Returns
    -------
    slopes_at_patch : n_patches-long array
        The slope (positive gradient magnitude) of each patch, in radians.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> mg = RasterModelGrid((4, 5))
    >>> z = mg.node_x
    >>> S = mg.calc_slope_at_patch(elevs=z)
    >>> S.size == mg.number_of_patches
    True
    >>> np.allclose(S, np.pi/4.)
    True
    >>> z = mg.node_y**2
    >>> mg.calc_slope_at_patch(elevs=z).reshape((3, 4))
    array([[ 0.78539816,  0.78539816,  0.78539816,  0.78539816],
           [ 1.24904577,  1.24904577,  1.24904577,  1.24904577],
           [ 1.37340077,  1.37340077,  1.37340077,  1.37340077]])

    >>> from landlab import CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY
    >>> z = mg.node_x.copy()
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> mg.status_at_node[11] = CLOSED_BOUNDARY
    >>> mg.status_at_node[9] = FIXED_VALUE_BOUNDARY
    >>> z[11] = 100.  # this should get ignored now
    >>> z[9] = 2.  # this should be felt by patch 7 only
    >>> mg.calc_slope_at_patch(elevs=z, ignore_closed_nodes=True).reshape(
    ...     (3, 4)) * 4./np.pi
    array([[ 0.,  0.,  0.,  0.],
           [ 0.,  1.,  1.,  1.],
           [ 0.,  0.,  0.,  0.]])

    LLCATS: PINF GRAD
    """
    if subtriangle_unit_normals is not None:
        assert len(subtriangle_unit_normals) == 4
        assert subtriangle_unit_normals[0].shape[1] == 3
        assert subtriangle_unit_normals[1].shape[1] == 3
        assert subtriangle_unit_normals[2].shape[1] == 3
        assert subtriangle_unit_normals[3].shape[1] == 3
        n_TR, n_TL, n_BL, n_BR = subtriangle_unit_normals
    else:
        n_TR, n_TL, n_BL, n_BR = grid.calc_unit_normals_at_patch_subtriangles(elevs)
    dotprod_TL = n_TL[:, 2]  # by definition
    dotprod_BR = n_BR[:, 2]
    dotprod_TR = n_TR[:, 2]
    dotprod_BL = n_BL[:, 2]
    slopes_at_patch_TL = np.arccos(dotprod_TL)  # 1 node order
    slopes_at_patch_BR = np.arccos(dotprod_BR)  # 3
    slopes_at_patch_TR = np.arccos(dotprod_TR)  # 0
    slopes_at_patch_BL = np.arccos(dotprod_BL)  # 2
    if ignore_closed_nodes:
        badnodes = grid.status_at_node[grid.nodes_at_patch] == CLOSED_BOUNDARY
        tot_bad = badnodes.sum(axis=1)
        tot_tris = 4.0 - 3.0 * (tot_bad > 0)  # 4 where all good, 1 where not
        # now shut down the bad tris. Remember, one bad node => 3 bad tris.
        # anywhere where badnodes > 1 will have zero from summing, so div by 1
        # assert np.all(np.logical_or(np.isclose(tot_tris, 4.),
        #                             np.isclose(tot_tris, 1.)))
        corners_rot = deque(
            [
                slopes_at_patch_BR,
                slopes_at_patch_TR,
                slopes_at_patch_TL,
                slopes_at_patch_BL,
            ]
        )
        # note initial offset so we are centered around TR on first slice
        for i in range(4):
            for j in range(3):
                (corners_rot[j])[badnodes[:, i]] = 0.0
            corners_rot.rotate(-1)
    else:
        tot_tris = 4.0
    mean_slope_at_patch = (
        slopes_at_patch_TR
        + slopes_at_patch_TL
        + slopes_at_patch_BL
        + slopes_at_patch_BR
    ) / tot_tris

    return mean_slope_at_patch


def calc_grad_at_patch(
    grid,
    elevs="topographic__elevation",
    ignore_closed_nodes=True,
    subtriangle_unit_normals=None,
    slope_magnitude=None,
):
    """Calculate the components of the gradient of each raster patch.

    Returns the mean gradient of the four possible patch subtriangles,
    in radians.

    If ignore_closed_nodes is True, closed nodes do not affect gradient
    calculations. If more than one closed node is present in a patch, the
    patch gradients in both x and y directions are set to zero.

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.
    elevs : str or ndarray, optional
        Field name or array of node values.
    ignore_closed_nodes : bool
        If True, do not incorporate values at closed nodes into the calc.
    subtriangle_unit_normals : tuple of 4 (npatches, 3) arrays (optional)
        The unit normal vectors for the four subtriangles of each patch,
        if already known. Order is TR, TL, BL, BR.
    slope_magnitude : array with size num_patches (optional)
        The mean slope of each patch, if already known. Units must be the
        same as provided here!

    Returns
    -------
    gradient_tuple : (x_component_at_patch, y_component_at_patch)
        Len-2 tuple of arrays giving components of gradient in the x and y
        directions, in the units of *radians*.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> mg = RasterModelGrid((4, 5))
    >>> z = mg.node_y
    >>> (x_grad, y_grad) = mg.calc_grad_at_patch(elevs=z)
    >>> np.allclose(y_grad, np.pi/4.)
    True
    >>> np.allclose(x_grad, 0.)
    True

    >>> from landlab import CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY
    >>> z = mg.node_x.copy()
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> mg.status_at_node[11] = CLOSED_BOUNDARY
    >>> mg.status_at_node[[9, 2]] = FIXED_VALUE_BOUNDARY
    >>> z[11] = 100.  # this should get ignored now
    >>> z[9] = 2.  # this should be felt by patch 7 only
    >>> z[2] = 1.  # should be felt by patches 1 and 2
    >>> xgrad, ygrad = mg.calc_grad_at_patch(
    ...     elevs=z, ignore_closed_nodes=True)
    >>> (xgrad.reshape((3, 4)) * 4./np.pi)[1, 1:]
    array([ 1.,  1., -1.])
    >>> np.allclose(ygrad[1:3], xgrad[1:3])
    True

    LLCATS: PINF GRAD
    """
    if subtriangle_unit_normals is not None:
        assert len(subtriangle_unit_normals) == 4
        assert subtriangle_unit_normals[0].shape[1] == 3
        assert subtriangle_unit_normals[1].shape[1] == 3
        assert subtriangle_unit_normals[2].shape[1] == 3
        assert subtriangle_unit_normals[3].shape[1] == 3
        n_TR, n_TL, n_BL, n_BR = subtriangle_unit_normals
    else:
        n_TR, n_TL, n_BL, n_BR = grid.calc_unit_normals_at_patch_subtriangles(elevs)
    if slope_magnitude is not None:
        assert slope_magnitude.size == grid.number_of_patches
        slopes_at_patch = slope_magnitude
    else:
        slopes_at_patch = grid.calc_slope_at_patch(
            elevs=elevs,
            ignore_closed_nodes=ignore_closed_nodes,
            subtriangle_unit_normals=(n_TR, n_TL, n_BL, n_BR),
        )

    if ignore_closed_nodes:
        badnodes = grid.status_at_node[grid.nodes_at_patch] == CLOSED_BOUNDARY
        corners_rot = deque([n_BR, n_TR, n_TL, n_BL])
        # note initial offset so we are centered around TR on first slice
        for i in range(4):
            for j in range(3):
                (corners_rot[j])[badnodes[:, i], :] = 0.0
            corners_rot.rotate(-1)

    n_sum_x = n_TR[:, 0] + n_TL[:, 0] + n_BL[:, 0] + n_BR[:, 0]
    n_sum_y = n_TR[:, 1] + n_TL[:, 1] + n_BL[:, 1] + n_BR[:, 1]
    theta_sum = np.arctan2(-n_sum_y, -n_sum_x)
    x_slope_patches = np.cos(theta_sum) * slopes_at_patch
    y_slope_patches = np.sin(theta_sum) * slopes_at_patch

    return (x_slope_patches, y_slope_patches)


def calc_slope_at_node(
    grid,
    elevs="topographic__elevation",
    method="patch_mean",
    ignore_closed_nodes=True,
    return_components=False,
):
    """Array of slopes at nodes, averaged over neighboring patches.

    Produces a value for node slope (i.e., mean gradient magnitude)
    at each node in a manner analogous to a GIS-style slope map.
    If method=='patch_mean', it averages the gradient on each of the
    patches surrounding the node; if method=='Horn', it returns the
    resolved slope direction. Directional information can still be
    returned through use of the return_components keyword.
    All values are returned in radians, including the components;
    take the tan to recover the rise/run.

    Note that under these definitions, it is not always true that::

        mag, cmp = mg.calc_slope_at_node(z)
        mag**2 == cmp[0]**2 + cmp[1]**2  # only if method=='Horn'

    If ignore_closed_nodes is False, all proximal elevation values will be used
    in the calculation. If True, only unclosed nodes are used.

    This is a verion of this code specialized for a raster. It subdivides
    the four square patches around each node into subtriangles,
    in order to ensure more correct solutions that incorporate equally
    weighted information from all surrounding nodes on rough surfaces.

    Parameters
    ----------
    elevs : str or ndarray, optional
        Field name or array of node values.
    method : {'patch_mean', 'Horn'}
        Controls the slope algorithm. Current options are 'patch_mean',
        which takes the mean slope of each pf the four neighboring
        square patches, and 'Horn', which is the standard ArcGIS slope
        algorithm. These produce very similar solutions; the Horn method
        gives a vector mean and the patch_mean gives a scalar mean.
    ignore_closed_nodes : bool
        If True, do not incorporate values at closed nodes into the calc.
    return_components : bool
        If True, return a tuple, (array_of_magnitude,
        (array_of_slope_x_radians, array_of_slope_y_radians)).
        If false, return an array of floats of the slope magnitude.

    Returns
    -------
    float array or length-2 tuple of float arrays
        If return_components, returns (array_of_magnitude,
        (array_of_slope_x_radians, array_of_slope_y_radians)).
        If not return_components, returns an array of slope magnitudes.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RadialModelGrid, RasterModelGrid
    >>> mg = RasterModelGrid((5, 5))
    >>> z = mg.node_x
    >>> slopes = mg.calc_slope_at_node(elevs=z)
    >>> np.allclose(slopes, np.pi / 4.)
    True
    >>> mg = RasterModelGrid((4, 5), xy_spacing=2.)
    >>> z = - mg.node_y
    >>> slope_mag, cmp = mg.calc_slope_at_node(elevs=z,
    ...                                        return_components=True)
    >>> np.allclose(slope_mag, np.pi / 4.)
    True
    >>> np.allclose(cmp[0], 0.)
    True
    >>> np.allclose(cmp[1], - np.pi / 4.)
    True
    >>> mg = RasterModelGrid((4, 4))
    >>> z = mg.node_x ** 2 + mg.node_y ** 2
    >>> slopes, cmp = mg.calc_slope_at_node(z, return_components=True)
    >>> slopes
    array([ 0.95531662,  1.10991779,  1.32082849,  1.37713803,  1.10991779,
            1.20591837,  1.3454815 ,  1.38904403,  1.32082849,  1.3454815 ,
            1.39288142,  1.41562833,  1.37713803,  1.38904403,  1.41562833,
            1.43030663])
    >>> np.allclose(cmp[0].reshape((4, 4))[:, 0],
    ...             cmp[1].reshape((4, 4))[0, :])  # test radial symmetry
    True

    LLCATS: NINF GRAD SURF
    """
    if method not in ("patch_mean", "Horn"):
        raise ValueError("method name not understood")
    try:
        patches_at_node = grid.patches_at_node()
    except TypeError:  # was a property, not a fn (=> new style)
        if not ignore_closed_nodes:
            patches_at_node = np.ma.masked_where(
                grid.patches_at_node == -1, grid.patches_at_node, copy=False
            )
        else:
            patches_at_node = np.ma.masked_where(
                np.logical_not(grid.patches_present_at_node),
                grid.patches_at_node,
                copy=False,
            )
    # now, we also want to mask any "closed" patches (any node closed)
    closed_patches = (grid.status_at_node[grid.nodes_at_patch] == CLOSED_BOUNDARY).sum(
        axis=1
    ) > 0
    closed_patch_mask = np.logical_or(
        patches_at_node.mask, closed_patches[patches_at_node.data]
    )

    if method == "patch_mean":
        n_TR, n_TL, n_BL, n_BR = grid.calc_unit_normals_at_patch_subtriangles(elevs)

        mean_slope_at_patches = grid.calc_slope_at_patch(
            elevs=elevs,
            ignore_closed_nodes=ignore_closed_nodes,
            subtriangle_unit_normals=(n_TR, n_TL, n_BL, n_BR),
        )

        # now CAREFUL - patches_at_node is MASKED
        slopes_at_node_unmasked = mean_slope_at_patches[patches_at_node]
        slopes_at_node_masked = np.ma.array(
            slopes_at_node_unmasked, mask=closed_patch_mask
        )
        slope_mag = np.mean(slopes_at_node_masked, axis=1).data
        if return_components:
            (x_slope_patches, y_slope_patches) = grid.calc_grad_at_patch(
                elevs=elevs,
                ignore_closed_nodes=ignore_closed_nodes,
                subtriangle_unit_normals=(n_TR, n_TL, n_BL, n_BR),
                slope_magnitude=mean_slope_at_patches,
            )
            x_slope_unmasked = x_slope_patches[patches_at_node]
            x_slope_masked = np.ma.array(x_slope_unmasked, mask=closed_patch_mask)
            x_slope = np.mean(x_slope_masked, axis=1).data
            y_slope_unmasked = y_slope_patches[patches_at_node]
            y_slope_masked = np.ma.array(y_slope_unmasked, mask=closed_patch_mask)
            y_slope = np.mean(y_slope_masked, axis=1).data
            mean_grad_x = x_slope
            mean_grad_y = y_slope
    elif method == "Horn":
        z = np.empty(grid.number_of_nodes + 1, dtype=float)
        mean_grad_x = grid.empty(at="node", dtype=float)
        mean_grad_y = grid.empty(at="node", dtype=float)
        z[-1] = 0.0
        try:
            z[:-1] = grid.at_node[elevs]
        except TypeError:
            z[:-1] = elevs
        # proof code for bad indexing:
        diags = grid.diagonal_neighbors_at_node.copy()  # LL order
        orthos = grid.adjacent_nodes_at_node.copy()
        # these have closed node neighbors...
        for dirs in (diags, orthos):
            dirs[dirs == BAD_INDEX_VALUE] = -1  # indexing to work
        # now make an array like patches_at_node to store the interim calcs
        patch_slopes_x = np.ma.zeros(patches_at_node.shape, dtype=float)
        patch_slopes_y = np.ma.zeros(patches_at_node.shape, dtype=float)
        diff_E = z[orthos[:, 0]] - z[:-1]
        diff_W = z[:-1] - z[orthos[:, 2]]
        diff_N = z[orthos[:, 1]] - z[:-1]
        diff_S = z[:-1] - z[orthos[:, 3]]
        patch_slopes_x[:, 0] = z[diags[:, 0]] - z[orthos[:, 1]] + diff_E
        patch_slopes_x[:, 1] = z[orthos[:, 1]] - z[diags[:, 1]] + diff_W
        patch_slopes_x[:, 2] = z[orthos[:, 3]] - z[diags[:, 2]] + diff_W
        patch_slopes_x[:, 3] = z[diags[:, 3]] - z[orthos[:, 3]] + diff_E
        patch_slopes_y[:, 0] = z[diags[:, 0]] - z[orthos[:, 0]] + diff_N
        patch_slopes_y[:, 1] = z[diags[:, 1]] - z[orthos[:, 2]] + diff_N
        patch_slopes_y[:, 2] = z[orthos[:, 2]] - z[diags[:, 2]] + diff_S
        patch_slopes_y[:, 3] = z[orthos[:, 0]] - z[diags[:, 3]] + diff_S
        patch_slopes_x /= 2.0 * grid.dx
        patch_slopes_y /= 2.0 * grid.dy
        patch_slopes_x.mask = closed_patch_mask
        patch_slopes_y.mask = closed_patch_mask
        mean_grad_x = patch_slopes_x.mean(axis=1).data
        mean_grad_y = patch_slopes_y.mean(axis=1).data
        slope_mag = np.arctan(np.sqrt(np.square(mean_grad_x) + np.square(mean_grad_y)))
        if return_components:
            mean_grad_x = np.arctan(mean_grad_x)
            mean_grad_y = np.arctan(mean_grad_y)

    if return_components:
        return slope_mag, (mean_grad_x, mean_grad_y)

    else:
        return slope_mag
