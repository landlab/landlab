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
    17
    """
    return links.number_of_links(shape)
