cimport cython
from cython.parallel cimport prange

from libc.stdint cimport int8_t
from libc.stdlib cimport free, malloc


ctypedef fused float_or_int:
    cython.floating
    cython.integral
    int8_t

ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
@cython.wraparound(False)
def find_links_at_node(
    const int node,
    const id_t [:, :] nodes_at_link,
    id_t [:] links_at_node,
    int8_t [:] link_dirs_at_node,
):
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
    cdef int n_links_found
    cdef int max_links_at_node = links_at_node.shape[0]
    cdef int n_links = nodes_at_link.shape[0]

    with nogil:
        n_links_found = _find_links_at_node(
            node,
            &nodes_at_link[0, 0],
            n_links,
            &links_at_node[0],
            &link_dirs_at_node[0],
            max_links_at_node,
        )

    return n_links_found


cdef int _find_links_at_node(
    const int node,
    const id_t * nodes_at_link,
    const int n_links,
    id_t * links_at_node,
    int8_t * link_dirs_at_node,
    const int max_links_at_node,
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
    cdef int link = 0
    cdef int n_links_found = 0
    cdef int offset

    while n_links_found < max_links_at_node and link < n_links:
        offset = 2 * link
        if nodes_at_link[offset] == node:
            links_at_node[n_links_found] = link
            link_dirs_at_node[n_links_found] = -1
            n_links_found += 1
        elif nodes_at_link[offset + 1] == node:
            links_at_node[n_links_found] = link
            link_dirs_at_node[n_links_found] = 1
            n_links_found += 1

        link += 1

    return n_links_found


@cython.boundscheck(False)
@cython.wraparound(False)
def get_links_at_node(
    const id_t [:, :] nodes_at_link,
    id_t [:, :] links_at_node,
    int8_t [:, :] link_dirs_at_node,
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
    cdef int max_links_at_node = links_at_node.shape[0]
    cdef int n_links = nodes_at_link.shape[0]

    for node in prange(n_nodes, nogil=True, schedule="static"):
        _find_links_at_node(
            node,
            &nodes_at_link[0, 0],
            n_links,
            &links_at_node[node, 0],
            &link_dirs_at_node[node, 0],
            max_links_at_node,
        )


@cython.boundscheck(False)
@cython.wraparound(False)
def reorder_links_at_node(
    float_or_int [:, :] links_at_node,
    const cython.integral [:, :] sorted_links,
):
    cdef int n_rows = links_at_node.shape[0]
    cdef int n_cols = links_at_node.shape[1]
    cdef int row
    cdef int col
    cdef float_or_int *buffer = <float_or_int *>malloc(
        n_rows * n_cols * sizeof(float_or_int)
    )
    cdef float_or_int *row_buffer

    for row in prange(n_rows, nogil=True, schedule="static"):
        row_buffer = buffer + row * n_cols

        for col in range(n_cols):
            row_buffer[col] = links_at_node[row, sorted_links[row, col]]
        for col in range(n_cols):
            links_at_node[row, col] = row_buffer[col]

    free(buffer)


@cython.boundscheck(False)
@cython.wraparound(False)
def reorder_link_dirs_at_node(
    int8_t [:, :] link_dirs_at_node,
    cython.integral [:, :] sorted_links,
):
    cdef int n_nodes = link_dirs_at_node.shape[0]
    cdef int n_links_per_node = link_dirs_at_node.shape[1]
    cdef int i
    cdef int node
    cdef int *buffer = <int *>malloc(n_links_per_node * sizeof(int))

    with nogil:
        try:
            for node in range(n_nodes):
                for i in range(n_links_per_node):
                    buffer[i] = link_dirs_at_node[node, sorted_links[node, i]]
                for i in range(n_links_per_node):
                    link_dirs_at_node[node, i] = buffer[i]
        finally:
            free(buffer)
