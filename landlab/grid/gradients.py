#! /usr/bin/env python
"""Calculate gradients of quantities over links."""


def calculate_gradients_at_active_links(grid, node_values, out=None):
    """Calculate gradients of node values over active links.

    Calculates the gradient in *quantity* node values at each active link in
    the grid.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    node_values : ndarray
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


def calculate_gradients_at_links(grid, node_values, out=None):
    """Calculate gradients of node values over links.

    Calculates the gradient in *quantity* node_values at each link in the grid.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    node_values : ndarray
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
