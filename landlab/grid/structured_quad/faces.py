from . import links


def number_of_faces(shape):
    """Number of faces in a structured quad grid.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    int :
        Number of faces in grid.

    Examples
    --------
    >>> from landlab.grid.structured_quad.faces import number_of_faces
    >>> number_of_faces((3, 4))
    7
    """
    if len(shape) != 2:
        raise ValueError('shape must be size 2')

    if min(shape) > 2:
        return ((shape[0] - 1) * (shape[1] - 2) +
                (shape[0] - 2) * (shape[1] - 1))
    else:
        return 0
