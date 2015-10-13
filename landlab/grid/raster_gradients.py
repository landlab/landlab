#! /usr/bin/env python
import numpy as np

from landlab.grid import gradients


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
