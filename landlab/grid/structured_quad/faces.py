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

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

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


def link_at_face(shape):
    """Get the link at each face.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Examples
    --------
    >>> from landlab.grid.structured_quad.faces import link_at_face
    >>> link_at_face((3, 5)) # doctest: +NORMALIZE_WHITESPACE
    array([  5,  6,  7,
             9, 10, 11, 12,
            14, 15, 16])

    Notes
    -----
    The algorithm is based on the following considerations:

    1. To get the link ID associated with a face ID, we start with the face
       ID and add something to it.
    2. Face 0 crosses the second vertical link. The first NC-1 links are
       horizontal, and run along the bottom of the grid. Then there's a
       vertical link on the lower-left side of the grid. Thus the link that
       crosses face 0 will be ID number NC, where NC is the number of cols.
    3. The offset between face ID and link ID increases by 1 every time we
       finish a row of horizontal faces (NC-2 of them) with vertical links,
       and again when we finish a row of vertical faces (NC-1 of them) with
       horizontal links. Together, a "row" of horizontal and vertical faces
       makes up 2NC-3 faces; this is the variable "fpr" (for "faces per 
       row") below.
    4. The quantity 2 * (face_ids // fpr) increases by 2 for every "full"
       (horizontal plus vertical) row of faces.
    5. The quantity ((face_ids % fpr) >= (NC-2)) increases by 1 every time
       we shift from a horizontal to a vertical row of faces (because
       each "full row" has NC-2 horizontal faces, then NC-1 vertical ones.)
    6. So, to find the offset, we add the "basic" offset, NC, to the "full
       row" offset, 2*(face_ids//fpr), and the "mid-row" offset, 
       ((face_ids % fpr)>=(NC-2)).
    """
    faces = np.arange(number_of_faces(shape))

    faces_per_row = (2 * shape[1]) - 3
    links = (shape[1] + (2 * (faces // faces_per_row)) +
             ((faces % faces_per_row) >= (shape[1] - 2)) + faces)

    return links
