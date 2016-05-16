#! /usr/bin/env python
"""Calculate gradients on a raster grid.

Gradient calculators for raster grids
++++++++++++++++++++++++++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.raster_gradients.calc_grad_at_link
    ~landlab.grid.raster_gradients.calc_grad_at_active_link
    ~landlab.grid.raster_gradients.calc_grad_across_cell_faces
    ~landlab.grid.raster_gradients.calc_grad_across_cell_corners
    ~landlab.grid.raster_gradients.alculate_gradient_along_node_links

"""
import numpy as np

from landlab.core.utils import make_optional_arg_into_id_array
from landlab.grid import gradients
from landlab.grid.base import BAD_INDEX_VALUE, CLOSED_BOUNDARY
from landlab.utils.decorators import use_field_name_or_array


@use_field_name_or_array('node')
def calc_grad_at_link(grid, node_values, out=None):
    """Calculate gradients in node_values at links.

    Construction::

        calc_grad_at_link(grid, node_values, out=None)

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

    >>> grid = RasterModelGrid((3, 3), spacing=(1, 2))
    >>> grid.calc_grad_at_link(node_values)
    array([ 0.,  0.,  1.,  3.,  1.,  1., -1.,  1., -1.,  1.,  0.,  0.])
    >>> _ = grid.add_field('node', 'elevation', node_values)
    >>> grid.calc_grad_at_link('elevation')
    array([ 0.,  0.,  1.,  3.,  1.,  1., -1.,  1., -1.,  1.,  0.,  0.])
    """
    grads = gradients.calc_diff_at_link(grid, node_values, out=out)
    grads /= grid.length_of_link

#    n_vertical_links = (grid.shape[0] - 1) * grid.shape[1]
#    diffs[:n_vertical_links] /= grid.dy
#    diffs[n_vertical_links:] /= grid.dx

    return grads


@use_field_name_or_array('node')
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
    >>> grid = RasterModelGrid(4, 5, 1.0)
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

    >>> grid = RasterModelGrid((3, 3), spacing=(1, 2))
    >>> node_values = [0., 0., 0.,
    ...                1., 3., 1.,
    ...                2., 2., 2.]
    >>> grid.calc_grad_at_active_link(node_values)
    array([ 3.,  1., -1., -1.])

    This function is *deprecated*. Instead, use ``calc_grad_at_link``.

    >>> grid = RasterModelGrid((3, 3), spacing=(1, 2))
    >>> node_values = [0., 0., 0.,
    ...                1., 3., 1.,
    ...                2., 2., 2.]
    >>> grid.calc_grad_at_link(node_values)[grid.active_links]
    array([ 3.,  1., -1., -1.])
    """
    if out is None:
        out = grid.empty(at='active_link')

    if len(out) != grid.number_of_active_links:
        raise ValueError('output buffer does not match that of the grid.')

    # grads = gradients.calculate_diff_at_active_links(grid, node_values,
    #                                                  out=out)
    grads = gradients.calc_diff_at_link(grid, node_values)
    out[:] = grads[grid.active_links]
    out /= grid.length_of_link[grid.active_links]

    return out


@use_field_name_or_array('node')
def calc_grad_across_cell_faces(grid, node_values, *args, **kwds):
    """Get gradients across the faces of a cell.

    Calculate gradient of the value field provided by *node_values* across
    each of the faces of the cells of a grid. The returned gradients are
    ordered as right, top, left, and bottom.

    Note that the returned gradients are masked to exclude neighbor nodes which
    are closed. Beneath the mask is the value -1.

    Construction::

        calc_grad_across_cell_faces(grid, node_values, [cell_ids],
                                             out=None)

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
    >>> grid = RasterModelGrid(3, 4)
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

    >>> grid = RasterModelGrid((3, 4), spacing=(2, 1))
    >>> grid.calc_grad_across_cell_faces(x) # doctest: +NORMALIZE_WHITESPACE
    masked_array(data =
     [[ 1.   1.5  0.   0. ]
     [ 0.   1.  -1.  -0.5]],
                  mask =
     False,
           fill_value = 1e+20)
    """
    padded_node_values = np.empty(node_values.size + 1, dtype=float)
    padded_node_values[-1] = BAD_INDEX_VALUE
    padded_node_values[:-1] = node_values
    cell_ids = make_optional_arg_into_id_array(grid.number_of_cells, *args)
    node_ids = grid.node_at_cell[cell_ids]

    neighbors = grid.active_neighbors_at_node(node_ids)
    if BAD_INDEX_VALUE != -1:
        neighbors = np.where(neighbors == BAD_INDEX_VALUE, -1, neighbors)
    values_at_neighbors = padded_node_values[neighbors]
    masked_neighbor_values = np.ma.array(
        values_at_neighbors, mask=neighbors == BAD_INDEX_VALUE)
    values_at_nodes = node_values[node_ids].reshape(len(node_ids), 1)

    out = np.subtract(masked_neighbor_values, values_at_nodes, **kwds)

    out[:, (0, 2)] /= grid.dx
    out[:, (1, 3)] /= grid.dy

    return out


@use_field_name_or_array('node')
def calc_grad_across_cell_corners(grid, node_values, *args, **kwds):
    """Get gradients to diagonally opposite nodes.

    Calculate gradient of the value field provided by *node_values* to
    the values at diagonally opposite nodes. The returned gradients are
    ordered as upper-right, upper-left, lower-left and lower-right.

    Construction::

        calc_grad_across_cell_corners(grid, node_values, [cell_ids],
                                               out=None)

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
    >>> grid = RasterModelGrid(3, 4)
    >>> x = np.array([1., 0., 0., 1.,
    ...               0., 0., 1., 1.,
    ...               3., 3., 3., 3.])

    A decrease in quantity to a diagonal node is a negative gradient.

    >>> from math import sqrt
    >>> grid.calc_grad_across_cell_corners(x) * sqrt(2.)
    array([[ 3.,  3.,  1.,  0.],
           [ 2.,  2., -1.,  0.]])

    >>> grid = RasterModelGrid((3, 4), spacing=(3, 4))
    >>> grid.calc_grad_across_cell_corners(x)
    array([[ 0.6,  0.6,  0.2,  0. ],
           [ 0.4,  0.4, -0.2,  0. ]])
    """
    cell_ids = make_optional_arg_into_id_array(grid.number_of_cells, *args)
    node_ids = grid.node_at_cell[cell_ids]

    values_at_diagonals = node_values[grid._get_diagonal_list(node_ids)]
    values_at_nodes = node_values[node_ids].reshape(len(node_ids), 1)

    out = np.subtract(values_at_diagonals, values_at_nodes, **kwds)
    np.divide(out, np.sqrt(grid.dy ** 2. + grid.dx ** 2.), out=out)

    return out


@use_field_name_or_array('node')
def calc_grad_along_node_links(grid, node_values, *args, **kwds):
    """Get gradients along links touching a node.

    Calculate gradient of the value field provided by *node_values* across
    each of the faces of the nodes of a grid. The returned gradients are
    ordered as right, top, left, and bottom. All returned values follow our
    standard sign convention, where a link pointing N or E and increasing in
    value is positive, a link pointing S or W and increasing in value is
    negative.

    Note that the returned gradients are masked to exclude neighbor nodes which
    are closed. Beneath the mask is the value -1.

    Construction::

        calc_grad_along_node_links(grid, node_values, [cell_ids],
                                            out=None)

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
    >>> grid = RasterModelGrid(3, 3)
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

    >>> grid = RasterModelGrid((3, 3), spacing=(2, 4))
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
    """
    padded_node_values = np.empty(node_values.size + 1, dtype=float)
    padded_node_values[-1] = BAD_INDEX_VALUE
    padded_node_values[:-1] = node_values
    node_ids = make_optional_arg_into_id_array(grid.number_of_nodes, *args)

    neighbors = grid.active_neighbors_at_node(node_ids, bad_index=-1)
    values_at_neighbors = padded_node_values[neighbors]
    masked_neighbor_values = np.ma.array(
        values_at_neighbors, mask=values_at_neighbors == BAD_INDEX_VALUE)
    values_at_nodes = node_values[node_ids].reshape(len(node_ids), 1)

    out = np.ma.empty_like(masked_neighbor_values, dtype=float)
    np.subtract(masked_neighbor_values[:, :2],
                values_at_nodes, out=out[:, :2], **kwds)
    np.subtract(values_at_nodes, masked_neighbor_values[:, 2:],
                out=out[:, 2:], **kwds)

    out[:, (0, 2)] /= grid.dx
    out[:, (1, 3)] /= grid.dy

    return out


def calc_unit_normals_at_patch_subtriangles(grid,
                                            elevs='topographic__elevation'):
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
    n_TR = nhat_topright/nmag_topright.reshape(grid.number_of_patches, 1)
    n_TL = nhat_topleft/nmag_topleft.reshape(grid.number_of_patches, 1)
    n_BL = nhat_bottomleft/nmag_bottomleft.reshape(
        grid.number_of_patches, 1)
    n_BR = nhat_bottomright/nmag_bottomright.reshape(
        grid.number_of_patches, 1)

    return (n_TR, n_TL, n_BL, n_BR)


def calc_slope_at_patch(grid, elevs='topographic__elevation',
                        subtriangle_unit_normals=None):
    """Calculate the slope (positive magnitude of gradient) at raster patches.

    Returns the mean of the slopes of the four possible patch subtriangles.

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.
    elevs : str or ndarray, optional
        Field name or array of node values.
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
    """
    if subtriangle_unit_normals is not None:
        assert len(subtriangle_unit_normals) == 4
        assert subtriangle_unit_normals[0].shape[1] == 3
        assert subtriangle_unit_normals[1].shape[1] == 3
        assert subtriangle_unit_normals[2].shape[1] == 3
        assert subtriangle_unit_normals[3].shape[1] == 3
        n_TR, n_TL, n_BL, n_BR = subtriangle_unit_normals
    else:
        n_TR, n_TL, n_BL, n_BR = (
            grid.calc_unit_normals_at_patch_subtriangles(elevs))
    dotprod_TL = n_TL[:, 2]  # by definition
    dotprod_BR = n_BR[:, 2]
    dotprod_TR = n_TR[:, 2]
    dotprod_BL = n_BL[:, 2]
    slopes_at_patch_TL = np.arccos(dotprod_TL)
    slopes_at_patch_BR = np.arccos(dotprod_BR)
    slopes_at_patch_TR = np.arccos(dotprod_TR)
    slopes_at_patch_BL = np.arccos(dotprod_BL)
    mean_slope_at_patch = (slopes_at_patch_TR + slopes_at_patch_TL +
                           slopes_at_patch_BL + slopes_at_patch_BR) / 4.

    return mean_slope_at_patch


def calc_grad_at_patch(grid, elevs='topographic__elevation',
                       subtriangle_unit_normals=None,
                       slope_magnitude=None):
    """Calculate the components of the gradient of each raster patch.

    Returns the mean gradient of the four possible patch subtriangles,
    in radians.

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.
    elevs : str or ndarray, optional
        Field name or array of node values.
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
    """
    if subtriangle_unit_normals is not None:
        assert len(subtriangle_unit_normals) == 4
        assert subtriangle_unit_normals[0].shape[1] == 3
        assert subtriangle_unit_normals[1].shape[1] == 3
        assert subtriangle_unit_normals[2].shape[1] == 3
        assert subtriangle_unit_normals[3].shape[1] == 3
        n_TR, n_TL, n_BL, n_BR = subtriangle_unit_normals
    else:
        n_TR, n_TL, n_BL, n_BR = \
            grid.calc_unit_normals_at_patch_subtriangles(elevs)
    if slope_magnitude is not None:
        assert slope_magnitude.size == grid.number_of_patches
        slopes_at_patch = slope_magnitude
    else:
        slopes_at_patch = grid.calc_slope_at_patch(
            elevs=elevs, subtriangle_unit_normals=(n_TR, n_TL, n_BL, n_BR))

    n_sum_x = n_TR[:, 0] + n_TL[:, 0] + n_BL[:, 0] + n_BR[:, 0]
    n_sum_y = n_TR[:, 1] + n_TL[:, 1] + n_BL[:, 1] + n_BR[:, 1]
    theta_sum = np.arctan2(-n_sum_y, -n_sum_x)
    x_slope_patches = np.cos(theta_sum)*slopes_at_patch
    y_slope_patches = np.sin(theta_sum)*slopes_at_patch

    return (x_slope_patches, y_slope_patches)


def calc_slope_at_node(grid, elevs='topographic__elevation',
                       method='patch_mean', return_components=False):
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
    >>> mg = RasterModelGrid((5, 5), 1.)
    >>> z = mg.node_x
    >>> slopes = mg.calc_slope_at_node(elevs=z)
    >>> np.allclose(slopes, np.pi / 4.)
    True
    >>> mg = RasterModelGrid((4, 5), 2.)
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
    """
    if method not in ('patch_mean', 'Horn'):
        raise ValueError('method name not understood')
    try:
        patches_at_node = grid.patches_at_node()
    except TypeError:  # was a property, not a fn (=> new style)
        patches_at_node = np.ma.masked_where(
            grid.patches_at_node == -1, grid.patches_at_node,
            copy=False)
    # now, we also want to mask any "closed" patches (any node closed)
    closed_patches = (grid.status_at_node[grid.nodes_at_patch] ==
                      CLOSED_BOUNDARY).sum(axis=1) > 0
    closed_patch_mask = np.logical_or(
        patches_at_node.mask, closed_patches[patches_at_node.data])

    if method == 'patch_mean':
        n_TR, n_TL, n_BL, n_BR = \
            grid.calc_unit_normals_at_patch_subtriangles(elevs)

        mean_slope_at_patches = grid.calc_slope_at_patch(
            elevs=elevs, subtriangle_unit_normals=(n_TR, n_TL, n_BL, n_BR))

        # now CAREFUL - patches_at_node is MASKED
        slopes_at_node_unmasked = mean_slope_at_patches[patches_at_node]
        slopes_at_node_masked = np.ma.array(slopes_at_node_unmasked,
                                            mask=closed_patch_mask)
        slope_mag = np.mean(slopes_at_node_masked, axis=1).data
        if return_components:
            (x_slope_patches, y_slope_patches) = grid.calc_grad_at_patch(
                elevs=elevs, subtriangle_unit_normals=(
                    n_TR, n_TL, n_BL, n_BR),
                slope_magnitude=mean_slope_at_patches)
            x_slope_unmasked = x_slope_patches[patches_at_node]
            x_slope_masked = np.ma.array(x_slope_unmasked,
                                         mask=closed_patch_mask)
            x_slope = np.mean(x_slope_masked, axis=1).data
            y_slope_unmasked = y_slope_patches[patches_at_node]
            y_slope_masked = np.ma.array(y_slope_unmasked,
                                         mask=closed_patch_mask)
            y_slope = np.mean(y_slope_masked, axis=1).data
            mean_grad_x = x_slope
            mean_grad_y = y_slope
    elif method == 'Horn':
        z = np.empty(grid.number_of_nodes + 1, dtype=float)
        mean_grad_x = grid.empty(at='node', dtype=float)
        mean_grad_y = grid.empty(at='node', dtype=float)
        z[-1] = 0.
        try:
            z[:-1] = grid.at_node[elevs]
        except TypeError:
            z[:-1] = elevs
        # proof code for bad indexing:
        diags = grid.diagonal_neighbors_at_node.copy()  # LL order
        orthos = grid.neighbors_at_node.copy()
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
        patch_slopes_x /= (2. * grid.dx)
        patch_slopes_y /= (2. * grid.dy)
        patch_slopes_x.mask = closed_patch_mask
        patch_slopes_y.mask = closed_patch_mask
        mean_grad_x = patch_slopes_x.mean(axis=1).data
        mean_grad_y = patch_slopes_y.mean(axis=1).data
        slope_mag = np.arctan(np.sqrt(np.square(mean_grad_x) +
                                      np.square(mean_grad_y)))
        if return_components:
            mean_grad_x = np.arctan(mean_grad_x)
            mean_grad_y = np.arctan(mean_grad_y)

    if return_components:
        return slope_mag, (mean_grad_x, mean_grad_y)

    else:
        return slope_mag
