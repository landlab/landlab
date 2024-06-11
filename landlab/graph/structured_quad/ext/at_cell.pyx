cimport cython
from cython.parallel cimport prange

ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_node_at_cell(
    shape,
    id_t[:] node_at_cell,
):
    """Get node contained in a cell.

    Parameters
    ----------
    shape : tuple of int
        Shape of the grid as `(n_rows, n_cols)`.
    node_at_cell : ndarray of int
        Buffer into which to place node identifiers.
    """
    cdef long n_rows = shape[0]
    cdef long n_cols = shape[1]
    cdef long cells_per_row = shape[1] - 2
    cdef long first_node
    cdef long cell
    cdef int row
    cdef int col

    for row in prange(1, n_rows - 1, nogil=True, schedule="static"):
        cell = (row - 1) * cells_per_row
        first_node = row * n_cols + 1
        for col in range(cells_per_row):
            node_at_cell[cell] = first_node + col
            cell = cell + 1

    return node_at_cell
