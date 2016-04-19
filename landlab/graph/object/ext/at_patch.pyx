import numpy as np
cimport numpy as np
cimport cython

from libc.stdlib cimport malloc, free, qsort


# from ...ext.argsort cimport unique_int
from ...sort.ext.argsort cimport unique_int


cdef find_common_node(long * link_a, long * link_b):
    if link_a[0] == link_b[0] or link_a[0] == link_b[1]:
        return link_a[0]
    elif link_a[1] == link_b[0] or link_a[1] == link_b[1]:
        return link_a[1]
    else:
        return ValueError('links are not connected')


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

        # out[n_nodes] = nodes_at_link[link][0]
        # out[n_nodes + 1] = nodes_at_link[link][1]

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

    # if out[0] != out[n_vertices - 1]:
    #     raise ValueError('not a closed patch')


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

            # n_links = all_nodes_at_patch(
            #     &links_at_patch[patch, 0], max_links_at_patch,
            #     &nodes_at_link[0, 0], all_nodes)

            # try:
            #   order_nodes_at_patch(all_nodes, &nodes_at_patch[patch, 0], n_links)
            # except ValueError:
            #   print nodes_at_patch[patch]
            #   raise
    finally:
        free(all_nodes)


cdef _nodes_at_patch(long * links_at_patch, long max_links,
                     long * nodes_at_link, long * out):
    cdef long n_links = max_links
    cdef long link, next_link
    cdef long i

    while links_at_patch[n_links - 1] == -1:
        n_links -= 1

    prev_link, link = links_at_patch[n_links - 1], links_at_patch[0]
    out[0] = find_common_node(&nodes_at_link[prev_link * 2],
                              &nodes_at_link[link * 2])

    for i in range(1, n_links):
        prev_link, link = link, links_at_patch[i]
        # link = links_at_patch[i]
        # next_link = links_at_patch[i + 1]

        out[i] = find_common_node(&nodes_at_link[prev_link * 2],
                                  &nodes_at_link[link * 2])

    # link = links_at_patch[n_links - 1]
    # next_link = links_at_patch[0]
    # out[n_links - 1] = find_common_node(&nodes_at_link[link * 2],
    #                                     &nodes_at_link[next_link * 2])

    return n_links


# @cython.boundscheck(False)
def __get_nodes_at_patch(np.ndarray[long, ndim=2, mode="c"] links_at_patch,
                       np.ndarray[long, ndim=2, mode="c"] nodes_at_link,
                       np.ndarray[long, ndim=2, mode="c"] nodes_at_patch):
    cdef int n_patches = links_at_patch.shape[0]
    cdef int max_links_at_patch = links_at_patch.shape[1]
    cdef int patch
    cdef int link
    cdef int i
    cdef int n_nodes
    cdef int n_unique
    cdef int * unique_nodes = <int *>malloc(2 * links_at_patch.shape[1] * sizeof(int))
    cdef long * all_nodes = <long *>malloc(2 * links_at_patch.shape[1] * sizeof(long))

    try:
        for patch in range(n_patches):
            n_nodes = 0
            for i in range(max_links_at_patch):
                link = links_at_patch[patch, i]
                if link == -1:
                    break

                all_nodes[n_nodes] = nodes_at_link[link, 0]
                all_nodes[n_nodes + 1] = nodes_at_link[link, 1]

                n_nodes += 2

            n_unique = unique_int(all_nodes, n_nodes, unique_nodes)
            for i in range(n_unique):
                nodes_at_patch[patch, i] = unique_nodes[i]
    finally:
        free(all_nodes)
        free(unique_nodes)
