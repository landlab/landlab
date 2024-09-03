#!/usr/bin/env python3
"""Functions to set up a finite-volume solution matrix for a landlab grid."""

import numpy as np
from scipy.sparse import csc_matrix

from ._matrix import fill_right_hand_side
from ._matrix import get_matrix_diagonal_elements
from ._matrix import get_matrix_diagonal_elements_with_coef


def get_core_node_at_node(grid):
    """Get node ids as numbered by core nodes.

    Get the core node ID for each node of a grid. If a node is not a core
    node, then use -1.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.

    Returns
    -------
    ndarray of int
        Ids of each of the grid's core nodes.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.utils import get_core_node_at_node

    >>> grid = RasterModelGrid((4, 5))
    >>> get_core_node_at_node(grid)
    array([-1, -1, -1, -1, -1,
           -1,  0,  1,  2, -1,
           -1,  3,  4,  5, -1,
           -1, -1, -1, -1, -1])

    >>> grid.status_at_node[13] = grid.BC_NODE_IS_FIXED_VALUE
    >>> grid.status_at_node[2] = grid.BC_NODE_IS_CLOSED
    >>> get_core_node_at_node(grid)
    array([-1, -1, -1, -1, -1,
           -1,  0,  1,  2, -1,
           -1,  3,  4, -1, -1,
           -1, -1, -1, -1, -1])
    """
    core_node_at_node = -np.ones(grid.number_of_nodes, dtype=int)
    core_node_at_node[grid.core_nodes] = np.arange(grid.number_of_core_nodes, dtype=int)
    return core_node_at_node


def get_matrix_entries(grid, coef_at_link=None):
    """Get entries of a sparse matrix.

    Parameters
    ----------
    grid : RasterModelGrid, HexModelGrid
        A landlab grid.
    coef_at_link : ndarray
        Coefficients at links used to construct the matrix.

    Returns
    -------
    tuple of (data, (row_ind, col_inds))
        Values of matrix elements along with their corresponding row and column
        index.
    """
    core2core = grid.link_with_node_status(
        status_at_tail=grid.BC_NODE_IS_CORE, status_at_head=grid.BC_NODE_IS_CORE
    )
    fv2core = grid.link_with_node_status(
        status_at_tail=grid.BC_NODE_IS_FIXED_VALUE, status_at_head=grid.BC_NODE_IS_CORE
    )
    core2fv = grid.link_with_node_status(
        status_at_tail=grid.BC_NODE_IS_CORE, status_at_head=grid.BC_NODE_IS_FIXED_VALUE
    )

    core_node_at_node = get_core_node_at_node(grid)

    nodes_at_c2fv_link = grid.nodes_at_link[core2fv]
    nodes_at_fv2c_link = grid.nodes_at_link[fv2core]

    core_nodes_at_c2c_link = core_node_at_node[grid.nodes_at_link[core2core]]
    core_nodes_at_c2fv_link = core_node_at_node[nodes_at_c2fv_link]
    core_nodes_at_fv2c_link = core_node_at_node[nodes_at_fv2c_link]

    n_core_nodes = grid.number_of_core_nodes

    values = np.zeros(n_core_nodes + 2 * len(core2core), dtype=float)
    row_inds = np.empty(n_core_nodes + 2 * len(core2core), dtype=int)
    col_inds = np.empty(n_core_nodes + 2 * len(core2core), dtype=int)

    diagonal_values = values[:n_core_nodes]
    diagonal_rows = row_inds[:n_core_nodes]
    diagonal_cols = col_inds[:n_core_nodes]

    upper_values = values[n_core_nodes : n_core_nodes + len(core2core)]
    upper_rows = row_inds[n_core_nodes : n_core_nodes + len(core2core)]
    upper_cols = col_inds[n_core_nodes : n_core_nodes + len(core2core)]

    lower_values = values[n_core_nodes + len(core2core) :]
    lower_rows = row_inds[n_core_nodes + len(core2core) :]
    lower_cols = col_inds[n_core_nodes + len(core2core) :]

    if coef_at_link is None:
        get_matrix_diagonal_elements(
            core_nodes_at_c2c_link,
            core_nodes_at_c2fv_link,
            core_nodes_at_fv2c_link,
            diagonal_values,
        )
        upper_values.fill(1.0)
    else:
        get_matrix_diagonal_elements_with_coef(
            core_nodes_at_c2c_link,
            core_nodes_at_c2fv_link,
            core_nodes_at_fv2c_link,
            coef_at_link[core2core],
            coef_at_link[core2fv],
            coef_at_link[fv2core],
            diagonal_values,
        )
        upper_values[:] = coef_at_link[core2core]

    diagonal_rows[:] = np.arange(n_core_nodes)
    diagonal_cols[:] = diagonal_rows

    upper_rows[:] = core_nodes_at_c2c_link[:, 0]
    upper_cols[:] = core_nodes_at_c2c_link[:, 1]

    lower_values[:] = upper_values
    lower_rows[:] = upper_cols
    lower_cols[:] = upper_rows

    return values, (row_inds, col_inds)


def get_core_node_matrix(grid, value_at_node, coef_at_link=None):
    """A matrix for core nodes and a right-hand side vector.

    Construct and return a matrix for the core nodes, plus a right-hand side vector
    containing values based on the input array `value_at_node`. Optionally,
    `coef_at_link` provides coefficients for each link (default is 1.0).

    Parameters
    ----------
    grid : RasterModelGrid, HexModelGrid
        A landlab grid.
    value_at_node : ndarray
        Values defined at nodes used to construct the right-hand side vector.
    coef_at_link : ndarray, optional
        Coefficents at links used to construct the matrix. If not provided,
        use 1.0.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.utils import get_core_node_matrix

    >>> grid = RasterModelGrid((4, 5))
    >>> grid.status_at_node[13] = grid.BC_NODE_IS_FIXED_VALUE
    >>> grid.status_at_node[2] = grid.BC_NODE_IS_CLOSED

    >>> vals = np.arange(
    ...     grid.number_of_nodes, dtype=np.double
    ... )  # made-up state variable array

    >>> mat, rhs = get_core_node_matrix(grid, vals)
    >>> mat.toarray()
    array([[-4.,  1.,  0.,  1.,  0.],
           [ 1., -3.,  1.,  0.,  1.],
           [ 0.,  1., -4.,  0.,  0.],
           [ 1.,  0.,  0., -4.,  1.],
           [ 0.,  1.,  0.,  1., -4.]])
    >>> rhs
    array([[ -6.],
           [  0.],
           [-25.],
           [-26.],
           [-30.]])

    >>> coefs = np.arange(grid.number_of_links, dtype=np.double)  # coefficient array
    >>> mat, rhs = get_core_node_matrix(grid, vals, coef_at_link=coefs)
    >>> mat.toarray()
    array([[-38.,  10.,   0.,  14.,   0.],
           [ 10., -36.,  11.,   0.,  15.],
           [  0.,  11., -46.,   0.,   0.],
           [ 14.,   0.,   0., -74.,  19.],
           [  0.,  15.,   0.,  19., -78.]])
    >>> rhs
    array([[ -6.],
           [  0.],
           [-25.],
           [-26.],
           [-30.]])
    """
    value_at_node = np.broadcast_to(value_at_node, grid.number_of_nodes)
    if coef_at_link is not None:
        coef_at_link = np.broadcast_to(coef_at_link, grid.number_of_links)

    (values, (row_inds, col_inds)) = get_matrix_entries(grid, coef_at_link=coef_at_link)

    mat = csc_matrix(
        (values, (row_inds, col_inds)),
        shape=(grid.number_of_core_nodes, grid.number_of_core_nodes),
    )

    fv2core = grid.link_with_node_status(
        status_at_tail=grid.BC_NODE_IS_FIXED_VALUE, status_at_head=grid.BC_NODE_IS_CORE
    )
    core2fv = grid.link_with_node_status(
        status_at_tail=grid.BC_NODE_IS_CORE, status_at_head=grid.BC_NODE_IS_FIXED_VALUE
    )

    core_node_at_node = get_core_node_at_node(grid)
    nodes_at_c2fv_link = grid.nodes_at_link[core2fv]
    nodes_at_fv2c_link = grid.nodes_at_link[fv2core]

    rhs = np.zeros(grid.number_of_core_nodes, dtype=float)
    fill_right_hand_side(
        nodes_at_c2fv_link, nodes_at_fv2c_link, core_node_at_node, value_at_node, rhs
    )

    return mat, rhs.reshape((-1, 1))
