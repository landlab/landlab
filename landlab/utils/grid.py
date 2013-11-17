
import numpy as np


def calculate_gradients_at_active_links(grid, node_values, out=None):
    """
    Calculates the gradient in *quantity* node_values at each active link in
    the grid.
    """
    if out is None:
        out = grid.empty(centering='active_link')
    return np.divide(node_values[grid.activelink_tonode] -
                     node_values[grid.activelink_fromnode],
                     grid.link_length[grid.active_links], out=out)


def calculate_diff_at_active_links(grid, node_values, out=None):
    """
    Calculates the difference in quantity *node_values* at each active link
    in the grid.
    """
    if out is None:
        out = grid.empty(centering='active_link')
    return np.subtract(node_values[grid.activelink_tonode],
                       node_values[grid.activelink_fromnode], out=out)
