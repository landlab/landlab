import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.int
ctypedef np.int_t DTYPE_t

@cython.boundscheck(False)
def get_node_at_cell(shape, np.ndarray[DTYPE_t, ndim=1] node_at_cell):
    """Get node contained in a cell.

    Parameters
    ----------
    shape : tuple of int
        Shape of the grid as `(n_rows, n_cols)`.
    node_at_cell : ndarray of int
        Buffer into which to place node identifiers.
    """
    cdef int cell
    cdef int cell_rows = shape[0] - 2
    cdef int cell_cols = shape[1] - 2
    cdef int node_cols = shape[1]
    cdef int row_offset
    cdef int row
    cdef int col

    cell = 0
    row_offset = shape[1] + 1
    for row in range(cell_rows):
        for col in range(cell_cols):
            node_at_cell[cell] = row_offset + col
            cell += 1
        row_offset += node_cols


@cython.boundscheck(False)
def get_nodes_at_face(shape, np.ndarray[DTYPE_t, ndim=2] nodes_at_face):
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
    face = 0
    for row in range(n_rows - 1):
        for col in range(1, n_cols - 1):
            node = row * n_cols + col
            nodes_at_face[face, 0] = node
            nodes_at_face[face, 1] = node + n_cols
            face += 1
        face += n_cols - 1

    # Vertical faces next
    face = n_cols - 2
    for row in range(1, n_rows - 1):
        for col in range(n_cols - 1):
            node = row * n_cols + col
            nodes_at_face[face, 0] = node
            nodes_at_face[face, 1] = node + 1
            face += 1
        face += n_cols - 2
