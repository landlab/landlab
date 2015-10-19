#! /usr/bin/env python
import numpy as np

from landlab.grid import gradients
from landlab.grid.base import BAD_INDEX_VALUE


def _make_optional_arg_into_array(number_of_elements, *args):
    """Transform an optional argument into an array.

    Parameters
    ----------
    number_of_elements : int
        Number of elements in the grid.
    array : array_like
        Iterable to convert to an array.

    Returns
    -------
    ndarray
        Input array converted to a numpy array, or a newly-created numpy
        array.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.raster_funcs import _make_optional_arg_into_array
    >>> _make_optional_arg_into_array(4)
    array([0, 1, 2, 3])
    >>> _make_optional_arg_into_array(4, [0, 0, 0, 0])
    [0, 0, 0, 0]
    >>> _make_optional_arg_into_array(4, (1, 1, 1, 1))
    [1, 1, 1, 1]
    >>> _make_optional_arg_into_array(4, np.ones(4))
    array([ 1.,  1.,  1.,  1.])
    """
    if len(args) >= 2:
        raise ValueError('Number of arguments must be 0 or 1.')

    if len(args) == 0:
        ids = np.arange(number_of_elements)
    else:
        ids = args[0]
        if not isinstance(ids, list) and not isinstance(ids, np.ndarray):
            try:
                ids = list(ids)
            except TypeError:
                ids = [ids]
    return ids


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

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.raster_gradients import (
    ...     calculate_gradients_at_links)
    >>> grid = RasterModelGrid((3, 3))
    >>> node_values = [0., 0., 0., 1., 3., 1., 2., 2., 2.]
    >>> calculate_gradients_at_links(grid, node_values)
    array([ 1.,  3.,  1.,  1., -1.,  1.,  0.,  0.,  2., -2.,  0.,  0.])

    >>> out = np.empty(grid.number_of_links, dtype=float)
    >>> rtn = calculate_gradients_at_links(grid, node_values, out=out)
    >>> rtn is out
    True
    >>> out
    array([ 1.,  3.,  1.,  1., -1.,  1.,  0.,  0.,  2., -2.,  0.,  0.])

    .. deprecated:: 0.1
        Use :func:`calculate_gradient_across_cell_faces`
                or :func:`calculate_gradient_across_cell_corners` instead
    """
    diffs = gradients.calculate_diff_at_links(grid, node_values, out=out)
    return np.divide(diffs, grid.node_spacing, out=diffs)


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
    >>> rmg = RasterModelGrid(4, 5, 1.0)
    >>> u = [0., 1., 2., 3., 0.,
    ...      1., 2., 3., 2., 3.,
    ...      0., 1., 2., 1., 2.,
    ...      0., 0., 2., 2., 0.]
    >>> grad = rmg.calculate_gradients_at_active_links(u)
    >>> grad
    array([ 1.,  1., -1., -1., -1., -1., -1.,  0.,  1.,  1.,  1., -1.,  1.,
            1.,  1., -1.,  1.])

    For greater speed, sending a pre-created numpy array as an argument
    avoids having to create a new one with each call:

    >>> grad = np.empty(rmg.number_of_active_links)
    >>> rtn = rmg.calculate_gradients_at_active_links(u, out=grad)
    >>> grad
    array([ 1.,  1., -1., -1., -1., -1., -1.,   0.,  1.,  1.,  1.,
           -1.,  1.,  1.,  1., -1.,  1.])
    >>> rtn is grad
    True
    """
    diffs = gradients.calculate_diff_at_active_links(grid, node_values,
                                                     out=out)
    return np.divide(diffs, grid.node_spacing, out=diffs)


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
    node_values : array_link
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
    >>> from landlab.grid.raster_funcs import (
    ...     calculate_gradient_across_cell_faces)
    >>> grid = RasterModelGrid(3, 4)
    >>> x = np.array([0., 0., 0., 0., 0., 0., 1., 1., 3., 3., 3., 3.])

    A decrease in quantity across a face is a negative gradient.

    >>> calculate_gradient_across_cell_faces(grid, x)
    masked_array(data =
     [[ 1.  3.  0.  0.]
     [ 0.  2. -1. -1.]],
                 mask =
     False,
           fill_value = 1e+20)
    <BLANKLINE>
    """
    padded_node_values = np.empty(node_values.size + 1, dtype=float)
    padded_node_values[-1] = BAD_INDEX_VALUE
    padded_node_values[:-1] = node_values
    cell_ids = _make_optional_arg_into_array(grid.number_of_cells, *args)
    node_ids = grid.node_at_cell[cell_ids]

    neighbors = grid.get_neighbor_list(node_ids)
    if BAD_INDEX_VALUE != -1:
        neighbors = np.where(neighbors == BAD_INDEX_VALUE, -1, neighbors)
    values_at_neighbors = padded_node_values[neighbors]
    masked_neighbor_values = np.ma.array(
        values_at_neighbors, mask=values_at_neighbors == BAD_INDEX_VALUE)
    values_at_nodes = node_values[node_ids].reshape(len(node_ids), 1)

    out = np.subtract(masked_neighbor_values, values_at_nodes, **kwds)
    out = np.multiply(out, 1. / grid.node_spacing, out=out)

    return out
