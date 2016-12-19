import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.int
ctypedef np.int_t DTYPE_t


@cython.boundscheck(False)
def fill_node_at_cell(shape, np.ndarray[DTYPE_t, ndim=1] node_at_cell):
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
