cimport cython
from cython.parallel cimport prange

ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_nodes_at_face(
    shape,
    id_t [:, :] nodes_at_face,
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
    cdef long vertical_faces_per_row = n_cols - 1
    cdef long horizontal_faces_per_row = n_cols - 2
    cdef long faces_per_row = vertical_faces_per_row + horizontal_faces_per_row

    for row in prange(n_rows - 2, nogil=True, schedule="static"):
        face = row * faces_per_row
        node = row * n_cols + 1

        for col in range(horizontal_faces_per_row):
            nodes_at_face[face, 0] = node
            nodes_at_face[face, 1] = node + n_cols

            node = node + 1
            face = face + 1

        node = (row + 1) * n_cols
        for col in range(vertical_faces_per_row):
            nodes_at_face[face, 0] = node
            nodes_at_face[face, 1] = node + 1

            node = node + 1
            face = face + 1

    face = (n_rows - 2) * faces_per_row
    node = (n_rows - 2) * n_cols + 1
    for col in prange(horizontal_faces_per_row, nogil=True, schedule="static"):
        nodes_at_face[face + col, 0] = node + col
        nodes_at_face[face + col, 1] = node + col + n_cols
