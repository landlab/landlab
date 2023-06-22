cimport cython

from cython.parallel cimport prange

ctypedef fused float_or_int:
    cython.integral
    cython.floating

@cython.boundscheck(False)
@cython.wraparound(False)
def map_max_of_link_nodes_to_link(
    float_or_int[:] out,
    const float_or_int[:] value_at_node,
    shape,
):
    cdef int node, link
    cdef int row, col
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int links_per_row = 2 * n_cols - 1

    for row in prange(n_rows, nogil=True, schedule="static"):
        link = row * links_per_row
        node = row * n_cols
        for col in range(n_cols - 1):
            out[link] = max(value_at_node[node], value_at_node[node + 1])
            link = link + 1
            node = node + 1

    for row in prange(n_rows - 1, nogil=True, schedule="static"):
        link = row * links_per_row + n_cols - 1
        node = row * n_cols
        for col in range(n_cols):
            out[link] = max(value_at_node[node], value_at_node[node + n_cols])
            link = link + 1
            node = node + 1
