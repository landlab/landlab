cimport cython
from cython.parallel cimport prange
from libc.stdint cimport int8_t

ctypedef fused float_or_int:
    cython.floating
    cython.integral
    int8_t

ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
@cython.wraparound(False)
def get_rightmost_edge_at_patch(
    const id_t [:, :] links_at_patch,
    const cython.floating [:, :] xy_of_link,
    id_t [:] edge,
):
    cdef int n_patches = links_at_patch.shape[0]
    cdef int n_cols = links_at_patch.shape[1]
    cdef int patch
    cdef int link
    cdef int n
    cdef int max_n
    cdef double max_x

    for patch in prange(n_patches, nogil=True, schedule="static"):
        link = links_at_patch[patch, 0]
        max_x, max_n = xy_of_link[link][0], 0

        for n in range(1, n_cols):
            link = links_at_patch[patch, n]
            if link == -1:
                break
            if xy_of_link[link][0] > max_x:
                max_x, max_n = xy_of_link[link][0], n
        edge[patch] = max_n


cdef id_t find_common_node(
    const id_t * link_a,
    const id_t * link_b,
) noexcept nogil:
    if link_a[0] == link_b[0] or link_a[0] == link_b[1]:
        return link_a[0]
    elif link_a[1] == link_b[0] or link_a[1] == link_b[1]:
        return link_a[1]
    else:
        raise ValueError("links are not connected")


cdef long all_nodes_at_patch(
    const id_t * links_at_patch,
    long max_links,
    const id_t * nodes_at_link,
    long * out,
) noexcept nogil:
    cdef long n_links = max_links
    cdef long link
    cdef long i
    cdef long n_nodes = 0

    while links_at_patch[n_links - 1] == -1:
        n_links -= 1

    for i in range(n_links):
        link = links_at_patch[i]

        out[n_nodes] = nodes_at_link[link * 2]
        out[n_nodes + 1] = nodes_at_link[link * 2 + 1]

        n_nodes += 2

    return n_links


cdef void order_nodes_at_patch(
    const id_t * all_nodes_at_patch,
    id_t * out,
    const long n_vertices,
):
    cdef long i
    cdef long vertex

    out[0] = all_nodes_at_patch[1]
    for vertex in range(n_vertices - 1):
        i = vertex * 2
        while all_nodes_at_patch[i] != out[vertex]:
            i += 1
        if i % 2 == 0:
            out[vertex + 1] = all_nodes_at_patch[i + 1]
        else:
            out[vertex + 1] = all_nodes_at_patch[i - 1]


@cython.boundscheck(False)
@cython.wraparound(False)
def get_nodes_at_patch(
    const id_t [:, :] links_at_patch,
    const id_t [:, :] nodes_at_link,
    id_t [:, :] nodes_at_patch,
):
    cdef int n_patches = links_at_patch.shape[0]
    cdef int max_links_at_patch = links_at_patch.shape[1]
    cdef int patch

    for patch in prange(n_patches, nogil=True, schedule="static"):
        _nodes_at_patch(
            &links_at_patch[patch, 0],
            max_links_at_patch,
            &nodes_at_link[0, 0],
            &nodes_at_patch[patch, 0],
        )


cdef long _nodes_at_patch(
    const id_t * links_at_patch,
    const long max_links,
    const id_t * nodes_at_link,
    id_t * out,
) noexcept nogil:
    cdef long n_links = max_links
    cdef long link, next_link
    cdef long i

    while links_at_patch[n_links - 1] == -1:
        n_links -= 1

    next_link = links_at_patch[0]
    for i in range(0, n_links - 1):
        link, next_link = next_link, links_at_patch[i + 1]

        out[i] = find_common_node(
            &nodes_at_link[link * 2], &nodes_at_link[next_link * 2]
        )

    link, next_link = links_at_patch[n_links - 1], links_at_patch[0]
    out[n_links - 1] = find_common_node(
        &nodes_at_link[link * 2], &nodes_at_link[next_link * 2]
    )

    return n_links
