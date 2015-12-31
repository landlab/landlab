import numpy as np
cimport numpy as np
cimport cython


DTYPE = np.int
ctypedef np.int_t DTYPE_t


@cython.boundscheck(False)
def _find_links_at_node(node,
                        np.ndarray[DTYPE_t, ndim=2] nodes_at_link,
                        np.ndarray[DTYPE_t, ndim=1] links_at_node,
                        np.ndarray[DTYPE_t, ndim=1] link_dirs_at_node):
  cdef int link
  cdef int n_links_found
  cdef int n_links
  cdef int max_links_at_node

  max_links_at_node = links_at_node.shape[0]
  n_links = nodes_at_link.shape[0]

  link = 0
  n_links_found = 0
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
def _setup_links_at_node(np.ndarray[DTYPE_t, ndim=2] nodes_at_link,
                         np.ndarray[DTYPE_t, ndim=2] links_at_node,
                         np.ndarray[DTYPE_t, ndim=2] link_dirs_at_node):
  cdef int node
  cdef int n_nodes

  n_nodes = links_at_node.shape[0]

  for node in range(n_nodes):
    _find_links_at_node(node, nodes_at_link, links_at_node[node],
                        link_dirs_at_node[node])


@cython.boundscheck(False)
def _setup_node_at_cell(shape,
                        np.ndarray[DTYPE_t, ndim=1] node_at_cell):
  cdef int cell
  cdef int cell_rows
  cdef int cell_cols
  cdef int row_offset

  cell_rows = shape[0] - 2
  cell_cols = shape[1] - 2

  cell = 0
  row_offset = shape[1] + 1
  for row in range(cell_rows):
    for col in range(cell_cols):
      node_at_cell[cell] = row_offset + col
      cell += 1
    row_offset += shape[1]
