#! /usr/bin/env python
"""Calculate gradients of quantities over links."""
import numpy as np
from landlab.utils.decorators import use_field_name_or_array


@use_field_name_or_array('node')
def calculate_gradients_at_active_links(grid, node_values, out=None):
    """Calculate gradients of node values over active links.

    Calculates the gradient in *quantity* node values at each active link in
    the grid.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    node_values : ndarray or field name
        Values at grid nodes.
    out : ndarray, optional
        Buffer to hold the result.

    Returns
    -------
    ndarray
        Gradients across active links.
    """
    if out is None:
        out = grid.empty(centering='active_link')
    return np.divide(node_values[grid.activelink_tonode] -
                     node_values[grid.activelink_fromnode],
                     grid.link_length[grid.active_links], out=out)


@use_field_name_or_array('node')
def calculate_gradients_at_links(grid, node_values, out=None):
    """Calculate gradients of node values over links.

    Calculates the gradient in *quantity* node_values at each link in the grid.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    node_values : ndarray or field name
        Values at grid nodes.
    out : ndarray, optional
        Buffer to hold the result.

    Returns
    -------
    ndarray
        Gradients across links.
    """
    if out is None:
        out = grid.empty(centering='link')
    return np.divide(node_values[grid.node_at_link_head] -
                     node_values[grid.node_at_link_tail],
                     grid.link_length, out=out)


@use_field_name_or_array('node')
def calculate_gradients_at_faces(grid, node_values, out=None):
    """Calculate gradients of node values over faces.

    Calculate and return gradient in *node_values* at each face in the grid.
    Gradients are calculated from the nodes at either end of the link that
    crosses each face.
    
    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    node_values : ndarray or field name
        Values at grid nodes.
    out : ndarray, optional
        Buffer to hold the result.

    Returns
    -------
    ndarray (x number of faces)
        Gradients across faces.
    
    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> rg = RasterModelGrid(3, 4, 10.0)
    >>> z = rg.add_zeros('node', 'topographic__elevation')
    >>> z[5] = 50.0
    >>> z[6] = 36.0
    >>> calculate_gradients_at_faces(rg, z)  # there are 7 faces
    array([ 5. ,  3.6,  5. , -1.4, -3.6, -5. , -3.6])

    >>> from landlab import HexModelGrid
    >>> hg = HexModelGrid(3, 3, 10.0)
    >>> z = rg.add_zeros('node', 'topographic__elevation')
    >>> z[4] = 50.0
    >>> z[5] = 36.0
    >>> calculate_gradients_at_faces(hg, z)  # there are 11 faces
    array([ 5. ,  5. ,  3.6,  3.6,  5. , -1.4, -3.6, -5. , -5. , -3.6, -3.6])
    """
    if out is None:
        out = grid.empty(centering='face')
    laf = grid.link_at_face
    return np.divide(node_values[grid.node_at_link_head[laf]] -
                     node_values[grid.node_at_link_tail[laf]],
                     grid.link_length[laf], out=out)


@use_field_name_or_array('node')
def calculate_diff_at_links(grid, node_values, out=None):
    """Calculate differences of node values over links.

    Calculates the difference in quantity *node_values* at each link in the
    grid.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    node_values : ndarray or field name
        Values at grid nodes.
    out : ndarray, optional
        Buffer to hold the result.

    Returns
    -------
    ndarray
        Differences across links.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> rmg = RasterModelGrid((3, 3))
    >>> z = np.zeros(9)
    >>> z[4] = 1.
    >>> rmg.calculate_diff_at_links(z)
    array([ 0.,  0.,  0.,  1.,  0.,  1., -1.,  0., -1.,  0.,  0.,  0.])
    """
    if out is None:
        out = grid.empty(centering='link')
    node_values = np.asarray(node_values)
    return np.subtract(node_values[grid.node_at_link_head],
                       node_values[grid.node_at_link_tail], out=out)


@use_field_name_or_array('node')
def calculate_diff_at_active_links(grid, node_values, out=None):
    """Calculate differences of node values over active links.

    Calculates the difference in quantity *node_values* at each active link
    in the grid.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    node_values : ndarray or field name
        Values at grid nodes.
    out : ndarray, optional
        Buffer to hold the result.

    Returns
    -------
    ndarray
        Differences across active links.
    """
    if out is None:
        out = grid.empty(centering='active_link')
    node_values = np.asarray(node_values)
    return np.subtract(node_values[grid.activelink_tonode],
                       node_values[grid.activelink_fromnode], out=out)
