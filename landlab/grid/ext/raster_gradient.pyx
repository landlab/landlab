cimport cython
from cython.parallel cimport prange

ctypedef fused float_or_int:
    cython.integral
    cython.floating


@cython.boundscheck(False)
@cython.wraparound(False)
def calc_diff_at_link(
    shape,
    const float_or_int[:] value_at_node,
    cython.floating[:] out,
):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int row, col
    cdef int link, node
    cdef int n_links = (2 * n_cols - 1) * n_rows - n_cols
    cdef int links_per_row = 2 * n_cols - 1

    for row in prange(0, n_rows - 1, nogil=True, schedule="static"):
        # Horizontal links
        link = row * links_per_row
        node = row * n_cols
        for col in range(n_cols - 1):
            out[link + col] = value_at_node[node + 1] - value_at_node[node]
            node = node + 1

        # Vertical links
        node = row * n_cols
        for col in range(n_cols - 1, links_per_row):
            out[link + col] = value_at_node[node + n_cols] - value_at_node[node]
            node = node + 1

    # The last row of horizontal links
    node = (n_rows - 1) * n_cols
    for link in range((n_rows - 1) * links_per_row, n_links):
        out[link] = value_at_node[node + 1] - value_at_node[node]
        node = node + 1


@cython.boundscheck(False)
@cython.wraparound(False)
def calc_grad_at_link(
    shape,
    xy_spacing,
    const float_or_int[:] value_at_node,
    cython.floating[:] out,
):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef double dx = xy_spacing[0]
    cdef double dy = xy_spacing[1]

    cdef int row, col
    cdef int link, node
    cdef int n_links = (2 * n_cols - 1) * n_rows - n_cols
    cdef int links_per_row = 2 * n_cols - 1
    cdef double inv_dx = 1.0 / dx
    cdef double inv_dy = 1.0 / dy

    for row in prange(0, n_rows - 1, nogil=True, schedule="static"):
        # Horizontal links
        link = row * links_per_row
        node = row * n_cols
        for col in range(n_cols - 1):
            out[link + col] = (value_at_node[node + 1] - value_at_node[node]) * inv_dx
            node = node + 1

        # Vertical links
        node = row * n_cols
        for col in range(n_cols - 1, links_per_row):
            out[link + col] = inv_dy * (
                value_at_node[node + n_cols] - value_at_node[node]
            )
            node = node + 1

    # The last row of horizontal links
    node = (n_rows - 1) * n_cols
    for link in range((n_rows - 1) * links_per_row, n_links):
        out[link] = (value_at_node[node + 1] - value_at_node[node]) * inv_dx
        node = node + 1
