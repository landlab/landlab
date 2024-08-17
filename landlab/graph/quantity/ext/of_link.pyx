cimport cython
from cython.parallel cimport prange

ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
@cython.wraparound(False)
def calc_midpoint_of_link(
    const id_t [:, :] nodes_at_link,
    const cython.floating [:] x_of_node,
    const cython.floating [:] y_of_node,
    cython.floating [:, :] xy_of_link,
):
    cdef int link
    cdef int n_links = nodes_at_link.shape[0]
    cdef int link_tail
    cdef int link_head

    for link in prange(n_links, nogil=True, schedule="static"):
        link_tail = nodes_at_link[link, 0]
        link_head = nodes_at_link[link, 1]

        xy_of_link[link][0] = (x_of_node[link_tail] + x_of_node[link_head]) * .5
        xy_of_link[link][1] = (y_of_node[link_tail] + y_of_node[link_head]) * .5
