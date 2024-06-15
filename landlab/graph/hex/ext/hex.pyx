cimport cython

from cython.parallel import prange

from libc.stdlib cimport free
from libc.stdlib cimport malloc

ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_xy_of_node_hex_horizontal(
    shape,
    cython.floating [:] x_of_node,
    cython.floating [:] y_of_node,
):
    """Get x and y coordinates for each node."""
    cdef long n_rows = shape[0]
    cdef long longest_row = n_rows // 2
    cdef long row
    cdef long offset
    cdef long *size_of_row
    cdef long *offset_to_row
    cdef double x0

    try:
        size_of_row = <long*>malloc(sizeof(long) * n_rows)
        offset_to_row = <long*>malloc(sizeof(long) * (n_rows + 1))

        size_of_row[0] = shape[1]
        for row in prange(1, longest_row + 1, nogil=True, schedule="static"):
            size_of_row[row] = size_of_row[row - 1] + 1
        for row in prange(longest_row + 1, n_rows, nogil=True, schedule="static"):
            size_of_row[row] = size_of_row[row - 1] - 1

        offset_to_row[0] = 0
        for row in range(1, n_rows + 1):
            offset_to_row[row] = offset_to_row[row - 1] + size_of_row[row - 1]

        for row in prange(longest_row + 1, nogil=True, schedule="static"):
            x0 = longest_row * 0.5 - row * 0.5

            for offset in range(offset_to_row[row], offset_to_row[row + 1]):
                x_of_node[offset] = x0 + offset - offset_to_row[row]
                y_of_node[offset] = row
            # x0 -= .5

        for row in prange(longest_row + 1, n_rows, nogil=True, schedule="static"):
            x0 = (row - longest_row) * 0.5

            for offset in range(offset_to_row[row], offset_to_row[row + 1]):
                x_of_node[offset] = x0 + offset - offset_to_row[row]
                y_of_node[offset] = row
    finally:
        free(offset_to_row)
        free(size_of_row)


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_xy_of_node_hex_vertical(
    shape,
    cython.floating [:] x_of_node,
    cython.floating [:] y_of_node,
):
    cdef long n_nodes = len(x_of_node)
    cdef int n_cols = shape[1]
    cdef int longest_column = n_cols // 2
    cdef int size_of_longest_column = longest_column + shape[0]
    cdef int n_rows = 2 * size_of_longest_column - 1
    cdef long size_of_odd_row = n_cols // 2
    cdef long size_of_even_row = n_cols - size_of_odd_row
    cdef long n_middle_rows = 2 * shape[0] - 1
    cdef long n_top_rows = n_cols // 2
    cdef long row
    cdef long n
    cdef long *size_of_row
    cdef long *offset_to_row
    cdef float x0
    cdef long offset

    try:
        size_of_row = <long*>malloc(sizeof(long) * n_rows)
        offset_to_row = <long*>malloc(sizeof(long) * (n_rows + 1))

        # Bottom and top rows
        for row in prange(n_top_rows, nogil=True, schedule="static"):
            size_of_row[row] = row + 1
            size_of_row[n_rows - row - 1] = row + 1

        # Middle rows
        for row in prange(0, n_middle_rows, 2, nogil=True, schedule="static"):
            size_of_row[row + n_top_rows] = size_of_even_row
            size_of_row[row + 1 + n_top_rows] = size_of_odd_row

        offset_to_row[0] = 0
        for row in range(1, n_rows + 1):
            offset_to_row[row] = offset_to_row[row - 1] + size_of_row[row - 1]

        for row in prange(0, n_rows, 2, nogil=True, schedule="static"):
            x0 = - (size_of_row[row] - 1) // 2
            for offset in range(offset_to_row[row], offset_to_row[row + 1]):
                x_of_node[offset] = x0 + offset - offset_to_row[row]
                y_of_node[offset] = row

        for row in prange(1, n_rows, 2, nogil=True, schedule="static"):
            x0 = - size_of_row[row] // 2 + .5
            for offset in range(offset_to_row[row], offset_to_row[row + 1]):
                x_of_node[offset] = x0 + offset - offset_to_row[row]
                y_of_node[offset] = row

        x0 = x_of_node[offset_to_row[n_top_rows]]
        for n in range(n_nodes):
            x_of_node[n] -= x0
    finally:
        free(offset_to_row)
        free(size_of_row)


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_xy_of_node_rect_vertical(
    shape,
    cython.floating [:] x_of_node,
    cython.floating [:] y_of_node,
):
    """Get x and y coordinates for each node."""
    cdef long n_cols = shape[1]
    cdef long size_of_odd_row = n_cols // 2
    cdef long size_of_even_row = n_cols - size_of_odd_row
    cdef long n_rows = 2 * shape[0]
    cdef long row
    cdef long offset
    cdef long n

    # even rows
    for row in prange(0, n_rows, 2, nogil=True, schedule="static"):
        offset = row // 2 * n_cols
        for n in range(size_of_even_row):
            x_of_node[offset + n] = n * 2
            y_of_node[offset + n] = row / 2.

    # odd rows
    offset = size_of_even_row
    for row in prange(1, n_rows, 2, nogil=True, schedule="static"):
        offset = size_of_even_row + (row - 1) // 2 * n_cols

        for n in range(size_of_odd_row):
            x_of_node[offset + n] = 2 * n + 1
            y_of_node[offset + n] = row / 2.


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_xy_of_node_rect_horizontal(
    shape,
    cython.floating [:] x_of_node,
    cython.floating [:] y_of_node,
):
    """Get x and y coordinates for each node."""
    cdef int stride = 2 * shape[1]
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int row
    cdef int col
    cdef int node

    for row in range(0, n_rows, 2):
        node = stride * row // 2
        for col in range(n_cols):
            x_of_node[node] = col
            y_of_node[node] = row

            node = node + 1

    for row in range(1, n_rows, 2):
        node = stride * (row - 1) // 2 + n_cols
        for col in range(n_cols):
            x_of_node[node] = col + .5
            y_of_node[node] = row

            node = node + 1
