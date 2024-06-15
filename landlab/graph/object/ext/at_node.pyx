cimport cython
cimport numpy as np
from cython.parallel cimport parallel
from cython.parallel cimport prange
from libc.stdint cimport int8_t
from libc.stdlib cimport free
from libc.stdlib cimport malloc
import numpy as np

ctypedef fused float_or_int:
    cython.floating
    cython.integral
    long long
    int8_t

ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
@cython.wraparound(False)
cdef long _find_links_at_node(
    const int node,
    const id_t[:, :] nodes_at_link,
    id_t[:] links_at_node,
    int8_t[:] link_dirs_at_node,
) noexcept nogil:
    """Find links touching a node and their directions.

    Parameters
    ----------
    node : int
        A node ID.
    nodes_at_link : ndarray of int, shape `(n_links, 2)`
        Nodes at link tail and head.
    links_at_node : ndarray of int, shape `(max_links_per_node, )`
        Buffer to hold link identifiers for links around node.
    link_dirs_at_node : ndarray of int, shape `(max_links_per_node, )`
        Buffer to hold link directions for links around node.

    Returns
    -------
    int
        The number of links found.
    """
    cdef long link = 0
    cdef long n_links_found = 0
    cdef long max_links_at_node = links_at_node.shape[0]
    cdef long n_links = nodes_at_link.shape[0]

    while n_links_found < max_links_at_node and link < n_links:
        if nodes_at_link[link, 0] == node:
            links_at_node[n_links_found] = link
            link_dirs_at_node[n_links_found] = -1
            n_links_found += 1
        elif nodes_at_link[link, 1] == node:
            links_at_node[n_links_found] = link
            link_dirs_at_node[n_links_found] = 1
            n_links_found += 1

        link += 1

    return n_links_found


@cython.boundscheck(False)
@cython.wraparound(False)
def get_links_at_node(
    const id_t[:, :] nodes_at_link,
    id_t[:, :] links_at_node,
    int8_t[:, :] link_dirs_at_node,
):
    """Get links touching each node and their directions.

    Parameters
    ----------
    nodes_at_link : ndarray of int, shape `(n_links, 2)`
        Node identifiers for each link tail and head.
    links_at_node : ndarray of int, shape `(n_nodes, max_nodes_per_link)`
        Buffer to hold link identifiers for each node.
    link_dirs_at_node : ndarray of int, shape `(n_nodes, max_nodes_per_link)`
        Buffer to hold link directions for each node.
    """
    cdef int node
    cdef int n_nodes = links_at_node.shape[0]

    for node in prange(n_nodes, nogil=True, schedule="static"):
        _find_links_at_node(
            node, nodes_at_link, links_at_node[node], link_dirs_at_node[node]
        )


@cython.boundscheck(False)
@cython.wraparound(False)
def reorder_rows(
    float_or_int[:, :] value_at_row,
    const id_t[:, :] sorted_cols,
):
    cdef long n_rows = value_at_row.shape[0]
    cdef long n_cols = value_at_row.shape[1]
    cdef long row
    cdef long col
    cdef np.ndarray[float_or_int, ndim=2] out = np.empty_like(value_at_row)

    for row in prange(n_rows, nogil=True, schedule="static"):
        for col in range(n_cols):
            out[row, col] = value_at_row[row, sorted_cols[row, col]]

    return out


@cython.boundscheck(False)
@cython.wraparound(False)
def reorder_rows_inplace(
    float_or_int[:, :] value_at_row,
    const id_t[:, :] sorted_cols,
):
    cdef int n_rows = value_at_row.shape[0]
    cdef int n_cols = value_at_row.shape[1]
    cdef int row
    cdef int col
    cdef float_or_int *scratch

    with nogil, parallel():
        scratch = <float_or_int *>malloc(n_cols * sizeof(float_or_int))

        for row in prange(n_rows, schedule="static"):
            for col in range(n_cols):
                scratch[col] = value_at_row[row, sorted_cols[row, col]]
            for col in range(n_cols):
                value_at_row[row, col] = scratch[col]

        free(scratch)


@cython.boundscheck(False)
@cython.wraparound(False)
def reorder_links_at_node(
    id_t[:, :] links_at_node,
    const id_t[:, :] sorted_links,
):
    cdef int n_nodes = links_at_node.shape[0]
    cdef int n_links_per_node = links_at_node.shape[1]
    cdef int i
    cdef int node
    cdef int *buffer = <int *>malloc(n_links_per_node * sizeof(int))

    try:
        for node in range(n_nodes):
            for i in range(n_links_per_node):
                buffer[i] = links_at_node[node, sorted_links[node, i]]
            for i in range(n_links_per_node):
                links_at_node[node, i] = buffer[i]
    finally:
        free(buffer)


@cython.boundscheck(False)
@cython.wraparound(False)
def reorder_link_dirs_at_node(
    int8_t[:, :] link_dirs_at_node,
    const id_t[:, :] sorted_links,
):
    cdef int n_nodes = link_dirs_at_node.shape[0]
    cdef int n_links_per_node = link_dirs_at_node.shape[1]
    cdef int i
    cdef int node
    cdef int *buffer = <int *>malloc(n_links_per_node * sizeof(int))

    try:
        for node in range(n_nodes):
            for i in range(n_links_per_node):
                buffer[i] = link_dirs_at_node[node, sorted_links[node, i]]
            for i in range(n_links_per_node):
                link_dirs_at_node[node, i] = buffer[i]
    finally:
        free(buffer)
