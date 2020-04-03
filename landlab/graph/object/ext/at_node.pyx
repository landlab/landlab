import numpy as np
cimport numpy as np
cimport cython

from libc.stdlib cimport malloc, free


ctypedef np.int_t DTYPE_t
ctypedef np.int8_t INT8TYPE_t


@cython.boundscheck(False)
def find_links_at_node(DTYPE_t node,
                       np.ndarray[DTYPE_t, ndim=2] nodes_at_link,
                       np.ndarray[DTYPE_t, ndim=1] links_at_node,
                       np.ndarray[INT8TYPE_t, ndim=1] link_dirs_at_node):
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
    cdef int max_links_at_node = links_at_node.shape[0]
    cdef int n_links = nodes_at_link.shape[0]

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
def get_links_at_node(np.ndarray[DTYPE_t, ndim=2] nodes_at_link,
                      np.ndarray[DTYPE_t, ndim=2] links_at_node,
                      np.ndarray[INT8TYPE_t, ndim=2] link_dirs_at_node):
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

    for node in range(n_nodes):
        find_links_at_node(node, nodes_at_link, links_at_node[node],
                           link_dirs_at_node[node])


@cython.boundscheck(False)
def reorder_links_at_node(np.ndarray[DTYPE_t, ndim=2] links_at_node,
                          np.ndarray[DTYPE_t, ndim=2] sorted_links):
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
def reorder_link_dirs_at_node(np.ndarray[INT8TYPE_t, ndim=2] link_dirs_at_node,
                              np.ndarray[DTYPE_t, ndim=2] sorted_links):
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
