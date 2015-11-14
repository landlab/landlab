#! /usr/bin/env python
"""Calculate gradients on a raster grid."""
import numpy as np

from landlab.core.utils import make_optional_arg_into_id_array
from landlab.grid import gradients
from landlab.grid.base import BAD_INDEX_VALUE
from landlab.utils.decorators import use_field_name_or_array


@use_field_name_or_array('node')
def calculate_gradients_at_links(grid, node_values, out=None):
    """Calculate gradients over links.

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.
    node_values : array_like
        Values at nodes.
    out : ndarray, optional
        Buffer to hold result. If `None`, create a new array.

    Returns
    -------
    ndarray
        Gradients of the nodes values for each link.

    .. deprecated:: 0.1
        Use :func:`calculate_gradient_across_cell_faces`
                or :func:`calculate_gradient_across_cell_corners` instead
    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((3, 3))
    >>> node_values = [0., 0., 0.,
    ...                1., 3., 1.,
    ...                2., 2., 2.]
    >>> grid.calculate_gradients_at_links(node_values)
    array([ 1.,  3.,  1.,  1., -1.,  1.,  0.,  0.,  2., -2.,  0.,  0.])

    >>> out = np.empty(grid.number_of_links, dtype=float)
    >>> rtn = grid.calculate_gradients_at_links(node_values, out=out)
    >>> rtn is out
    True
    >>> out
    array([ 1.,  3.,  1.,  1., -1.,  1.,  0.,  0.,  2., -2.,  0.,  0.])

    >>> grid = RasterModelGrid((3, 3), spacing=(1, 2))
    >>> grid.calculate_gradients_at_links(node_values)
    array([ 1.,  3.,  1.,  1., -1.,  1.,  0.,  0.,  1., -1.,  0.,  0.])

    >>> _ = grid.add_field('node', 'elevation', node_values)
    >>> grid.calculate_gradients_at_links('elevation')
    array([ 1.,  3.,  1.,  1., -1.,  1.,  0.,  0.,  1., -1.,  0.,  0.])
    """
    diffs = gradients.calculate_diff_at_links(grid, node_values, out=out)

    n_vertical_links = (grid.shape[0] - 1) * grid.shape[1]
    diffs[:n_vertical_links] /= grid.dy
    diffs[n_vertical_links:] /= grid.dx

    return diffs


@use_field_name_or_array('node')
def calculate_gradients_at_active_links(grid, node_values, out=None):
    """Calculate gradients over active links.

    .. deprecated:: 0.1
        Use :func:`calculate_gradient_across_cell_faces`
                or :func:`calculate_gradient_across_cell_corners` instead

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
    >>> grad = grid.calculate_gradients_at_active_links(u)
    >>> grad
    array([ 1.,  1., -1., -1., -1., -1., -1.,  0.,  1.,  1.,  1., -1.,  1.,
            1.,  1., -1.,  1.])

    For greater speed, sending a pre-created numpy array as an argument
    avoids having to create a new one with each call:

    >>> grad = np.empty(grid.number_of_active_links)
    >>> rtn = grid.calculate_gradients_at_active_links(u, out=grad)
    >>> grad
    array([ 1.,  1., -1., -1., -1., -1., -1.,   0.,  1.,  1.,  1.,
           -1.,  1.,  1.,  1., -1.,  1.])
    >>> rtn is grad
    True

    >>> grid = RasterModelGrid((3, 3), spacing=(1, 2))
    >>> node_values = [0., 0., 0.,
    ...                1., 3., 1.,
    ...                2., 2., 2.]
    >>> grid.calculate_gradients_at_active_links(node_values)
    array([ 3., -1.,  1., -1.])
    """
    diffs = gradients.calculate_diff_at_active_links(grid, node_values,
                                                     out=out)

    n_vertical_links = (grid.shape[0] - 1) * grid.shape[1]

    vertical_active_links = np.where(grid.active_links < n_vertical_links)
    horizontal_active_links = np.where(grid.active_links >= n_vertical_links)

    diffs[vertical_active_links] /= grid.dy
    diffs[horizontal_active_links] /= grid.dx

    return diffs


@use_field_name_or_array('node')
def calculate_gradient_across_cell_faces(grid, node_values, *args, **kwds):
    """calculate_gradient_across_cell_faces(grid, node_values, [cell_ids], out=None)
    Get gradients across the faces of a cell.

    Calculate gradient of the value field provided by *node_values* across
    each of the faces of the cells of a grid. The returned gradients are
    ordered as right, top, left, and bottom.

    Note that the returned gradients are masked to exclude neighbor nodes which
    are closed. Beneath the mask is the value numpy.iinfo(numpy.int32).max.

    Parameters
    ----------
    grid : RasterModelGrid
        Source grid.
    node_values : array_like
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

    >>> grid.calculate_gradient_across_cell_faces(x)
    masked_array(data =
     [[ 1.  3.  0.  0.]
     [ 0.  2. -1. -1.]],
                 mask =
     False,
           fill_value = 1e+20)

    >>> grid = RasterModelGrid((3, 4), spacing=(2, 1))
    >>> grid.calculate_gradient_across_cell_faces(x)
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

    neighbors = grid.get_active_neighbors_at_node(node_ids)
    if BAD_INDEX_VALUE != -1:
        neighbors = np.where(neighbors == BAD_INDEX_VALUE, -1, neighbors)
    values_at_neighbors = padded_node_values[neighbors]
    masked_neighbor_values = np.ma.array(
        values_at_neighbors, mask=values_at_neighbors == BAD_INDEX_VALUE)
    values_at_nodes = node_values[node_ids].reshape(len(node_ids), 1)

    out = np.subtract(masked_neighbor_values, values_at_nodes, **kwds)

    out[:, (0, 2)] /= grid.dx
    out[:, (1, 3)] /= grid.dy

    return out


@use_field_name_or_array('node')
def calculate_gradient_across_cell_corners(grid, node_values, *args, **kwds):
    """calculate_gradient_across_cell_corners(grid, node_values, [cell_ids], out=None)
    Get gradients to diagonally opposite nodes.

    Calculate gradient of the value field provided by *node_values* to
    the values at diagonally opposite nodes. The returned gradients are
    ordered as upper-right, upper-left, lower-left and lower-right.

    Parameters
    ----------
    grid : RasterModelGrid
        Source grid.
    node_values : array_like
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
    >>> grid.calculate_gradient_across_cell_corners(x) * sqrt(2.)
    array([[ 3.,  3.,  1.,  0.],
           [ 2.,  2., -1.,  0.]])

    >>> grid = RasterModelGrid((3, 4), spacing=(3, 4))
    >>> grid.calculate_gradient_across_cell_corners(x)
    array([[ 0.6,  0.6,  0.2,  0. ],
           [ 0.4,  0.4, -0.2,  0. ]])
    """
    cell_ids = make_optional_arg_into_id_array(grid.number_of_cells, *args)
    node_ids = grid.node_at_cell[cell_ids]

    values_at_diagonals = node_values[grid.get_diagonal_list(node_ids)]
    values_at_nodes = node_values[node_ids].reshape(len(node_ids), 1)

    out = np.subtract(values_at_diagonals, values_at_nodes, **kwds)
    np.divide(out, np.sqrt(grid.dy ** 2. + grid.dx ** 2.), out=out)

    return out


@use_field_name_or_array('node')
def calculate_gradient_along_node_links(grid, node_values, *args, **kwds):
    """calculate_gradient_along_node_links(grid, node_values, [cell_ids], out=None)
    Get gradients along links touching a node.

    Calculate gradient of the value field provided by *node_values* across
    each of the faces of the nodes of a grid. The returned gradients are
    ordered as right, top, left, and bottom. All returned values follow our
    standard sign convention, where a link pointing N or E and increasing in
    value is positive, a link pointing S or W and increasing in value is
    negative.

    Note that the returned gradients are masked to exclude neighbor nodes which
    are closed. Beneath the mask is the value numpy.iinfo(numpy.int32).max.

    Parameters
    ----------
    grid : RasterModelGrid
        Source grid.
    node_values : array_like
        Quantity to take the gradient of defined at each node.
    node_ids : array_like, optional
        If provided, node ids to measure gradients. Otherwise, find gradients
        for all nodes.
    out : array_like, optional
        Alternative output array in which to place the result.  Must
        be of the same shape and buffer length as the expected output.

    Returns
    -------
    (N, 4) ndarray
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

    >>> grid.calculate_gradient_along_node_links(x)
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
    >>> grid.calculate_gradient_along_node_links(x)
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

    neighbors = grid.get_active_neighbors_at_node(node_ids, bad_index=-1)
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
