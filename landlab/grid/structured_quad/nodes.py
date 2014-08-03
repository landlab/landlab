import numpy as np


def number_of_nodes(shape):
    """Number of nodes is a structured quad grid.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    int :
        Number of nodes in the grid.

    Examples
    --------
    >>> number_of_nodes((3, 4))
    12
    """
    return np.prod(shape)


def corners(shape):
    """IDs of corner nodes.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    (4, ) ndarray :
        IDs of the corner nodes.

    Examples
    --------
    >>> corners((3, 4))
    array([ 0,  3,  8, 11])
    """
    node_count = number_of_nodes(shape)
    return np.array([0, shape[1] - 1, node_count - shape[1], node_count - 1])