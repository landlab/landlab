cimport cython
from libc.stdint cimport int8_t
from libc.stdint cimport int32_t
from libc.stdint cimport int64_t
from libc.stdlib cimport calloc
from libc.stdlib cimport free

ctypedef fused id_t:
    int32_t
    int64_t


@cython.boundscheck(False)
@cython.wraparound(False)
def _get_links_at_node(
    const id_t[:, ::1] nodes_at_link,
    id_t[:, ::1] links_at_node,
    int8_t[:, ::1] link_dirs_at_node,
):
    """Get links touching each node and their directions.

    Parameters
    ----------
    nodes_at_link : ndarray of int, shape (n_links, 2)
        Node for each link tail and head.
    links_at_node : ndarray of int, shape (n_nodes, max_links_per_node)
        Buffer to hold link identifiers for each node.
    link_dirs_at_node : ndarray of int8, shape (n_nodes, max_links_per_node)
        Buffer to hold link directions for each node.
    """
    cdef Py_ssize_t n_links = nodes_at_link.shape[0]
    cdef Py_ssize_t n_nodes = links_at_node.shape[0]
    cdef Py_ssize_t max_links_per_node = links_at_node.shape[1]
    cdef Py_ssize_t count
    cdef Py_ssize_t link
    cdef id_t tail, head
    cdef Py_ssize_t* n_links_at_node = <Py_ssize_t*>calloc(n_nodes, sizeof(Py_ssize_t))

    if n_links_at_node == NULL:
        raise MemoryError("malloc failed in _get_links_at_node")

    try:
        with nogil:
            for link in range(n_links):
                tail = nodes_at_link[link, 0]
                head = nodes_at_link[link, 1]

                count = n_links_at_node[tail]
                if count < max_links_per_node:
                    links_at_node[tail, count] = link
                    link_dirs_at_node[tail, count] = -1
                    n_links_at_node[tail] = count + 1

                count = n_links_at_node[head]
                if count < max_links_per_node:
                    links_at_node[head, count] = link
                    link_dirs_at_node[head, count] = 1
                    n_links_at_node[head] = count + 1
    finally:
        free(n_links_at_node)
