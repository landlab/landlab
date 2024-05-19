cimport cython
from cython.parallel cimport prange


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_nodes_at_face(
    shape,
    cython.integral [:, :] nodes_at_face,
):
    """Get nodes on either side of a face.

    Parameters
    ----------
    shape : tuple of int
        Shape of the grid as `(n_rows, n_cols)`.
    nodes_at_face : ndarray of int, shape `(n_faces, 2)`
        Buffer into which to place node identifiers.
    """
    cdef int face
    cdef int node
    cdef int row
    cdef int col
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]

    # Horizontal faces first
    for row in prange(n_rows - 1, nogil=True, schedule="static"):
        face = row * ((n_cols - 1) + (n_cols - 2))
        for col in range(1, n_cols - 1):
            node = row * n_cols + col
            nodes_at_face[face, 0] = node
            nodes_at_face[face, 1] = node + n_cols
            face = face + 1

    # Vertical faces next
    for row in prange(1, n_rows - 1, nogil=True, schedule="static"):
        face = row * (n_cols - 2) + (row - 1) * (n_cols - 1)
        # face = row * ((n_cols - 2) + (n_cols - 1))
        for col in range(n_cols - 1):
            node = row * n_cols + col
            nodes_at_face[face, 0] = node
            nodes_at_face[face, 1] = node + 1
            face = face + 1
