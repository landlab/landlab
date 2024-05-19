cimport cython
from cython.parallel cimport prange


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_node_at_cell(
    shape,
    cython.integral [:] node_at_cell,
):
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

    for row in prange(cell_rows, nogil=True, schedule="static"):
        cell = row * cell_cols
        row_offset = (row + 1) * node_cols + 1
        for col in range(cell_cols):
            node_at_cell[cell] = row_offset + col
            cell = cell + 1
