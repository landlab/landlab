import numpy as np

from . import links
from ..base import BAD_INDEX_VALUE


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


def face_at_link(shape, bad_index_value=BAD_INDEX_VALUE):
    """Get face for each link of a grid.

    Examples
    --------
    >>> from landlab.grid.structured_quad.faces import face_at_link
    >>> face_at_link((4, 5), bad_index_value=-1)
    ...     # doctest: +NORMALIZE_WHITESPACE
    array([-1, -1, -1, -1, -1,  0,  1,  2, -1,
            3,  4,  5,  6, -1,  7,  8,  9, -1,
           10, 11, 12, 13, -1, 14, 15, 16, -1,
           -1, -1, -1, -1])
    """
    from .c_faces import _face_at_link

    out = np.full(links.number_of_links(shape), bad_index_value, dtype=int)
    _face_at_link(shape, out)

    return out
