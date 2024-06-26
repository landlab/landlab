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
    cdef long n_rows = shape[0]
    cdef long nodes_per_row = shape[1]
    cdef long horizontal_links_per_row = nodes_per_row - 1
    cdef long vertical_links_per_row = nodes_per_row
    cdef long links_per_row = horizontal_links_per_row + vertical_links_per_row
    cdef long link
    cdef long node
    cdef long row
    cdef long first_link
    cdef long first_node

    for row in prange(n_rows - 1, nogil=True, schedule="static"):
        first_link = row * links_per_row
        first_node = row * nodes_per_row

        node = first_node
        for link in range(first_link, first_link + horizontal_links_per_row):
            out[link] = max(value_at_node[node], value_at_node[node + 1])

            node = node + 1

        first_link = first_link + horizontal_links_per_row
        node = first_node
        for link in range(first_link, first_link + vertical_links_per_row):
            out[link] = max(
                value_at_node[node], value_at_node[node + nodes_per_row]
            )

            node = node + 1

    with nogil:
        first_link = (n_rows - 1) * links_per_row
        first_node = (n_rows - 1) * nodes_per_row

        node = first_node
        for link in range(first_link, first_link + horizontal_links_per_row):
            out[link] = max(value_at_node[node], value_at_node[node + 1])

            node = node + 1
