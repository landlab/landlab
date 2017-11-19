#! /usr/bin/env python
"""Calculate steepest descent on a raster grid.

Steepest descent functions for raster grids
+++++++++++++++++++++++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.raster_steepest_descent._calc_steepest_descent_across_cell_corners
    ~landlab.grid.raster_steepest_descent._calc_steepest_descent_across_cell_faces

"""
import numpy as np

from landlab.core.utils import make_optional_arg_into_id_array
from landlab.utils.decorators import deprecated
from landlab.grid.raster_gradients import (
    calc_grad_across_cell_faces,
    calc_grad_across_cell_corners,
    calc_grad_along_node_links)


_VALID_ROUTING_METHODS = set(['d8', 'd4'])


def _assert_valid_routing_method(method):
    """Check if name is a valid routing method.

    Parameters
    ----------
    method : str
        Name of a routing method.

    Raises
    ------
    ValueError
        If the method name is invalid.
    """
    if method not in _VALID_ROUTING_METHODS:
        raise ValueError(
            '%s: routing method not understood. should be one of %s' %
            (method, ', '.join(_VALID_ROUTING_METHODS)))


def _calc_steepest_descent_across_cell_corners(grid, node_values, *args,
                                              **kwds):
    """Get steepest gradient to the diagonals of a cell.

    Calculate the gradients in *node_values* measure to the diagonals of cells
    IDs, *cell_ids*. Slopes upward from the cell are reported as positive.
    If *cell_ids* is not given, calculate gradients for all cells.

    Use the *return_node* keyword to return a tuple, with the first element
    being the gradients and the second the node id of the node in the direction
    of the minimum gradient, i.e., the steepest descent. Note the gradient
    value returned is probably thus negative.

    Parameters
    ----------
    grid : RasterModelGrid
        Input grid.
    node_values : array_like
        Values to take gradient of.
    cell_ids : array_like, optional
        IDs of grid cells to measure gradients.
    return_node: boolean, optional
        If `True`, return node IDs of the node that has the steepest descent.
    out : ndarray, optional
        Alternative output array in which to place the result.  Must
        be of the same shape and buffer length as the expected output.

    Returns
    -------
    ndarray :
        Calculated gradients to lowest node across cell faces.

    Examples
    --------
    Create a rectilinear grid that is 3 nodes by 3 nodes and so has one cell
    centered around node 4.

    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid(3, 3)
    >>> values_at_nodes = np.arange(9.)

    Calculate gradients to cell diagonals and choose the gradient to the
    lowest node.

    >>> from math import sqrt
    >>> grid._calc_steepest_descent_across_cell_corners(
    ...     values_at_nodes) * sqrt(2.)
    array([-4.])

    The steepest gradient is to node with id 0.

    >>> (_, ind) = grid._calc_steepest_descent_across_cell_corners(
    ...                values_at_nodes, return_node=True)
    >>> ind
    array([0])

    >>> grid = RasterModelGrid(3, 3)
    >>> node_values = grid.zeros()
    >>> node_values[0] = -1
    >>> grid._calc_steepest_descent_across_cell_corners(node_values, 0)
    array([-0.70710678])

    Get both the maximum gradient and the node to which the gradient is
    measured.

    >>> grid._calc_steepest_descent_across_cell_corners(node_values, 0,
    ...     return_node=True)
    (array([-0.70710678]), array([0]))

    LLCATS: NINF CNINF GRAD
    """
    return_node = kwds.pop('return_node', False)

    cell_ids = make_optional_arg_into_id_array(grid.number_of_cells, *args)

    grads = calc_grad_across_cell_corners(grid, node_values, cell_ids)

    if return_node:
        ind = np.argmin(grads, axis=1)
        node_ids = grid.diagonal_cells[grid.node_at_cell[cell_ids], ind]
        if 'out' not in kwds:
            out = np.empty(len(cell_ids), dtype=grads.dtype)
        out[:] = grads[range(len(cell_ids)), ind]
        return (out, node_ids)
    else:
        return grads.min(axis=1, **kwds)

@deprecated(use='FlowRouter', version=1.0)
def _calc_steepest_descent_across_cell_faces(grid, node_values, *args, **kwds):
    """Get steepest gradient across the faces of a cell.

    This method calculates the gradients in *node_values* across all four
    faces of the cell or cells with ID *cell_ids*. Slopes upward from the
    cell are reported as positive. If *cell_ids* is not given, calculate
    gradients for all cells.

    Use the *return_node* keyword to return a tuple, with the first element
    being the gradients and the second the node id of the node in the direction
    of the minimum gradient, i.e., the steepest descent. Note the gradient
    value returned is probably thus negative.

    Parameters
    ----------
    grid : RasterModelGrid
        Input grid.
    node_values : array_like
        Values to take gradient of.
    cell_ids : array_like, optional
        IDs of grid cells to measure gradients.
    return_node: boolean, optional
        Return node IDs of the node that has the steepest descent.
    out : ndarray, optional
        Alternative output array in which to place the result.  Must
        be of the same shape and buffer length as the expected output.

    Returns
    -------
    ndarray :
        Calculated gradients to lowest node across cell faces.

    Convention: gradients positive UP

    Examples
    --------
    Create a rectilinear grid that is 3 nodes by 3 nodes and so has one cell
    centered around node 4.

    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid(3, 3)
    >>> values_at_nodes = np.arange(9.)

    Calculate gradients across each cell face and choose the gradient to the
    lowest node.

    >>> grid._calc_steepest_descent_across_cell_faces(values_at_nodes)
    masked_array(data = [-3.],
                 mask = False,
           fill_value = 1e+20)
    <BLANKLINE>

    The steepest gradient is to node with id 1.

    >>> (_, ind) = grid._calc_steepest_descent_across_cell_faces(
    ...     values_at_nodes, return_node=True)
    >>> ind
    array([1])

    >>> grid = RasterModelGrid(3, 3)
    >>> node_values = grid.zeros()
    >>> node_values[1] = -1
    >>> grid._calc_steepest_descent_across_cell_faces(node_values, 0)
    masked_array(data = [-1.],
                 mask = False,
           fill_value = 1e+20)
    <BLANKLINE>

    Get both the maximum gradient and the node to which the gradient is
    measured.

    >>> grid._calc_steepest_descent_across_cell_faces(node_values, 0,
    ...     return_node=True)
    (array([-1.]), array([1]))

    LLCATS: NINF FINF GRAD
    """
    return_node = kwds.pop('return_node', False)

    cell_ids = make_optional_arg_into_id_array(grid.number_of_cells, *args)

    grads = calc_grad_across_cell_faces(grid, node_values, cell_ids)

    if return_node:
        ind = np.argmin(grads, axis=1)
        node_ids = grid.active_neighbors_at_node[grid.node_at_cell[
            cell_ids], ind]
        # node_ids = grid.neighbor_nodes[grid.node_at_cell[cell_ids], ind]
        if 'out' not in kwds:
            out = np.empty(len(cell_ids), dtype=grads.dtype)
        out[:] = grads[range(len(cell_ids)), ind]
        return (out, node_ids)
        # return (out, 3 - ind)
    else:
        return grads.min(axis=1, **kwds)
