#! /usr/bin/env python
import numpy as np

from landlab.grid import gradients


def calculate_gradients_at_links(grid, node_values, out=None):
    """Calculate gradients over links.

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

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> rmg = RasterModelGrid(4, 5, 1.0)
    >>> u = [0., 1., 2., 3., 0.,
    ...     1., 2., 3., 2., 3.,
    ...     0., 1., 2., 1., 2.,
    ...     0., 0., 2., 2., 0.]
    >>> u = np.array(u)
    >>> u
    array([ 0.,  1.,  2.,  3.,  0.,  1.,  2.,  3.,  2.,  3.,  0.,  1.,  2.,
            1.,  2.,  0.,  0.,  2.,  2.,  0.])
    >>> grad = rmg.calculate_gradients_at_active_links(u)
    >>> grad
    array([ 1.,  1., -1., -1., -1., -1., -1.,  0.,  1.,  1.,  1., -1.,  1.,
            1.,  1., -1.,  1.])

    For greater speed, sending a pre-created numpy array as an argument
    avoids having to create a new one with each call:

    >>> grad = np.zeros(rmg.number_of_active_links)
    >>> u = u*10
    >>> grad = rmg.calculate_gradients_at_active_links(u, grad)
    >>> grad
    array([ 10.,  10., -10., -10., -10., -10., -10.,   0.,  10.,  10.,  10.,
           -10.,  10.,  10.,  10., -10.,  10.])
    """
    diffs = gradients.calculate_diff_at_active_links(grid, node_values,
                                                     out=out)
    return np.divide(diffs, grid.node_spacing, out=diffs)
