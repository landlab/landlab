import numpy as np


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
    6
    """
    return np.prod(np.array(shape) - 1)