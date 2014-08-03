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