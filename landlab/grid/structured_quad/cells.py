import numpy as np

from . import nodes


def number_of_cells(shape):
    """Number of cells is a structured quad grid.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    int :
        Number of cells in the grid.

    Examples
    --------
    >>> number_of_cells((3, 4))
    2
    """
    return np.prod(np.array(shape) - 2)


def shape_of_cells(shape):
    """Shape of cells in a structured quad grid.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    int :
        Shape of cells in grid.

    Examples
    --------
    >>> shape_of_cells((3, 4))
    (1, 2)
    """
    return tuple(np.array(shape) - 2)


def node_id_at_cells(shape):
    """Node ID at each cell.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    int :
        ID of node associated with each cell.

    Examples
    --------
    >>> node_id_at_cells((3, 4))
    array([[5, 6]])
    """
    node_ids = nodes.node_ids(shape)
    return node_ids[1:-1, 1:-1].copy().reshape(shape_of_cells(shape))
