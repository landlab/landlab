import numpy as np
cimport numpy as np
cimport cython

from libc.stdlib cimport malloc, free, qsort

from ...sort.ext.argsort cimport unique_int


@cython.boundscheck(True)
def get_rightmost_edge_at_patch(
    np.ndarray[long, ndim=2, mode="c"] links_at_patch,
    np.ndarray[double, ndim=2, mode="c"] xy_of_link,
    np.ndarray[long, ndim=1, mode="c"] edge):
    cdef int n_patches = links_at_patch.shape[0]
    cdef int n_cols = links_at_patch.shape[1]
    cdef int patch
    cdef int link
    cdef int n
    cdef int max_n
    cdef double max_x

    for patch in range(n_patches):
        link = links_at_patch[patch, 0]
        max_x, max_n = xy_of_link[link][0], 0

        for n in range(1, n_cols):
            link = links_at_patch[patch, n]
            if link == -1:
                break
            if xy_of_link[link][0] > max_x:
                max_x, max_n = xy_of_link[link][0], n
        edge[patch] = max_n


cdef find_common_node(long * link_a, long * link_b):
    if link_a[0] == link_b[0] or link_a[0] == link_b[1]:
        return link_a[0]
    elif link_a[1] == link_b[0] or link_a[1] == link_b[1]:
        return link_a[1]
    else:
        raise ValueError('links are not connected')


cdef all_nodes_at_patch(long * links_at_patch, long max_links,
                        long * nodes_at_link, long * out):
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


cdef order_nodes_at_patch(long * all_nodes_at_patch, long * out,
                          long n_vertices):
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


@cython.boundscheck(True)
def get_nodes_at_patch(np.ndarray[long, ndim=2, mode="c"] links_at_patch,
                       np.ndarray[long, ndim=2, mode="c"] nodes_at_link,
                       np.ndarray[long, ndim=2, mode="c"] nodes_at_patch):
    cdef int n_patches = links_at_patch.shape[0]
    cdef int max_links_at_patch = links_at_patch.shape[1]
    cdef int patch
    cdef int link
    cdef int i
    cdef int n_links
    cdef int n_unique
    cdef long * all_nodes = <long *>malloc(2 * links_at_patch.shape[1] * sizeof(long))

    try:
        for patch in range(n_patches):
            n_links = _nodes_at_patch(
                &links_at_patch[patch, 0], max_links_at_patch,
                &nodes_at_link[0, 0], &nodes_at_patch[patch, 0])
    finally:
        free(all_nodes)


cdef _nodes_at_patch(long * links_at_patch, long max_links,
                     long * nodes_at_link, long * out):
    cdef long n_links = max_links
    cdef long link, next_link, prev_link
    cdef long i

    while links_at_patch[n_links - 1] == -1:
        n_links -= 1

    next_link = links_at_patch[0]
    for i in range(0, n_links - 1):
        link, next_link = next_link, links_at_patch[i + 1]

        out[i] = find_common_node(&nodes_at_link[link * 2],
                                  &nodes_at_link[next_link * 2])

    link, next_link = links_at_patch[n_links - 1], links_at_patch[0]
    out[n_links - 1] = find_common_node(&nodes_at_link[link * 2],
                                        &nodes_at_link[next_link * 2])

    return n_links
