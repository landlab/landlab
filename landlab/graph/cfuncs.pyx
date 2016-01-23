import numpy as np
cimport numpy as np
cimport cython

from libc.stdlib cimport malloc, free, qsort
from libc.math cimport atan2


DTYPE = np.int
ctypedef np.int_t DTYPE_t

from libc.stdlib cimport malloc, free, qsort


ctypedef np.int_t INT_t
ctypedef np.float_t FLOAT_t

cdef struct Sorter:
    INT_t index
    FLOAT_t value


cdef extern from "stdlib.h":
    ctypedef void const_void "const void"
    void qsort(void *base, int nmemb, int size,
            int(*compar)(const_void *, const_void *)) nogil


cdef int _compare(const_void *a, const_void *b):
    cdef double v = ((<Sorter*>a)).value - ((<Sorter*>b)).value
    if v < 0:
        return -1
    else:
        return 1


cdef void _argsort(double * data, int n_elements, Sorter * order):
    cdef int i

    for i in range(n_elements):
        order[i].index = i
        order[i].value = data[i]

    qsort(<void*> order, n_elements, sizeof(Sorter), _compare)


cdef argsort(double * data, int n_elements, int * out):
    cdef Sorter *sorted_struct = <Sorter*>malloc(n_elements * sizeof(Sorter))

    try:
        _argsort(data, n_elements, sorted_struct)
        for i in range(n_elements):
            out[i] = sorted_struct[i].index
    finally:
        free(sorted_struct)


cdef argsort_inplace(double * data, int n_elements, int * out):
    cdef Sorter *sorted_struct = <Sorter*>malloc(n_elements * sizeof(Sorter))

    try:
        _argsort(data, n_elements, sorted_struct)
        for i in range(n_elements):
            out[i] = sorted_struct[i].index
            data[i] = sorted_struct[i].value
    finally:
        free(sorted_struct)


@cython.boundscheck(False)
def _find_links_at_node(DTYPE_t node,
                        np.ndarray[DTYPE_t, ndim=2] nodes_at_link,
                        np.ndarray[DTYPE_t, ndim=1] links_at_node,
                        np.ndarray[DTYPE_t, ndim=1] link_dirs_at_node):
  """Find directions of links touching a node.

  Parameters
  ----------
  node : int
      A node ID.
  nodes_at_link : ndarray of int, shape `(n_links, 2)`
      Nodes at link tail and head.
  links_at_node : ndarray of int, shape `(max_links_per_node, )`
      Buffer to hold link IDs for links around node.
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
  cdef int cell_rows = shape[0] - 2
  cdef int cell_cols = shape[1] - 2
  cdef int node_cols = shape[1]
  cdef int row_offset
  cdef int row
  cdef int col


  cell = 0
  row_offset = shape[1] + 1
  for row in range(cell_rows):
    for col in range(cell_cols):
      node_at_cell[cell] = row_offset + col
      cell += 1
    row_offset += node_cols


@cython.boundscheck(False)
def _remap_nodes_at_link(np.ndarray[DTYPE_t, ndim=2] nodes_at_link,
                         np.ndarray[DTYPE_t, ndim=1] new_nodes):
  cdef int link
  cdef int n_links = len(nodes_at_link)

  for link in range(n_links):
    nodes_at_link[link, 0] = new_nodes[nodes_at_link[link, 0]]
    nodes_at_link[link, 1] = new_nodes[nodes_at_link[link, 1]]


@cython.boundscheck(False)
def _remap_links_at_patch(np.ndarray[DTYPE_t, ndim=1] links_at_patch,
                          np.ndarray[DTYPE_t, ndim=1] new_links):
    cdef int link
    cdef int n_links = links_at_patch.shape[0]

    for link in range(n_links):
        links_at_patch[link] = new_links[links_at_patch[link]]


@cython.boundscheck(False)
def _calc_center_of_patch(np.ndarray[DTYPE_t, ndim=1] links_at_patch,
                          np.ndarray[DTYPE_t, ndim=1] offset_to_patch,
                          np.ndarray[np.float_t, ndim=2] xy_at_link,
                          np.ndarray[np.float_t, ndim=2] xy_at_patch):
    cdef int patch
    cdef int link
    cdef int i
    cdef int offset
    cdef int n_links
    cdef int n_patches = len(xy_at_patch)
    cdef float x
    cdef float y

    for patch in range(n_patches):
        offset = offset_to_patch[patch]
        n_links = offset_to_patch[patch + 1] - offset
        x = 0.
        y = 0.
        for i in range(offset, offset + n_links):
            link = links_at_patch[i]
            x += xy_at_link[link, 0]
            y += xy_at_link[link, 1]
        xy_at_patch[patch, 0] = x / n_links
        xy_at_patch[patch, 1] = y / n_links


@cython.boundscheck(False)
def _resort_patches(np.ndarray[DTYPE_t, ndim=1] links_at_patch,
                    np.ndarray[DTYPE_t, ndim=1] offset_to_patch,
                    np.ndarray[DTYPE_t, ndim=1] sorted_patches):
    cdef int i
    cdef int patch
    cdef int offset
    cdef int n_links
    cdef int n_patches = len(sorted_patches)
    cdef int *new_offset = <int *>malloc(len(offset_to_patch) * sizeof(int))
    cdef int *new_patches = <int *>malloc(len(links_at_patch) * sizeof(int))

    try:
        new_offset[0] = 0
        for patch in range(n_patches):
            offset = offset_to_patch[sorted_patches[patch]]
            n_links = offset_to_patch[sorted_patches[patch] + 1] - offset

            new_offset[patch + 1] = new_offset[patch] + n_links
            for i in range(n_links):
                new_patches[new_offset[patch] + i] = links_at_patch[offset + i]

        for i in range(len(links_at_patch)):
            links_at_patch[i] = new_patches[i]
        for i in range(len(offset_to_patch)):
            offset_to_patch[i] = new_offset[i]
    finally:
        free(new_offset)
        free(new_patches)


@cython.boundscheck(False)
def _setup_links_at_patch(np.ndarray[DTYPE_t, ndim=1] links_at_patch,
                          np.ndarray[DTYPE_t, ndim=1] offset_to_patch,
                          np.ndarray[DTYPE_t, ndim=2] out):
    cdef int i
    cdef int link
    cdef int patch
    cdef int offset
    cdef int n_links
    cdef int n_patches = len(offset_to_patch) - 1

    for patch in range(n_patches):
        offset = offset_to_patch[patch]
        n_links = offset_to_patch[patch + 1] - offset

        link = 0
        for i in range(offset, offset + n_links):
          out[patch, link] = links_at_patch[i]
          link += 1


@cython.boundscheck(False)
def _reorder_links_at_node(np.ndarray[DTYPE_t, ndim=2] links_at_node,
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


cdef calc_centroid(double * points, np.int_t n_points, double * out):
    cdef int i
    cdef double xc = 0.
    cdef double yc = 0.

    for i in range(0, 2 * n_points, 2):
        xc += points[i]
        yc += points[i + 1]
    print xc, yc
    xc /= n_points
    yc /= n_points

    out[0] = xc
    out[1] = yc


cdef calc_spoke_angles(double * hub, double * spokes, np.int_t n_spokes,
                       double * angles):
    cdef int i
    cdef double x0 = hub[0]
    cdef double y0 = hub[1]
    cdef double * spoke = spokes
    cdef double two_pi = 2. * np.pi

    for i in range(n_spokes):
        x = spoke[0]
        y = spoke[1]

        angles[i] = atan2(y - y0, x - x0)
        if angles[i] < 0.:
            angles[i] += two_pi
        print angles[i] * 180. / 3.14
        spoke += 2


cdef argsort_by_angle_around_centroid(double * points,
                                      np.int_t n_points,
                                      int * out):
    cdef double *hub = <double *>malloc(2 * sizeof(double))
    cdef double *angles = <double *>malloc(n_points * sizeof(double))

    try:
        calc_centroid(points, n_points, hub)
        print hub[0], hub[1]
        calc_spoke_angles(hub, points, n_points, angles)
        argsort(angles, n_points, out)
    finally:
        free(angles)
        free(hub)


@cython.boundscheck(False)
def reorder_links_at_patch(np.ndarray[DTYPE_t, ndim=1] links_at_patch,
                           np.ndarray[DTYPE_t, ndim=1] offset_to_patch,
                           np.ndarray[np.float_t, ndim=2] xy_of_link):
    cdef int n_links = xy_of_link.shape[0]
    cdef int n_patches = len(offset_to_patch) - 1
    cdef int link
    cdef int patch
    cdef int offset
    cdef int i
    cdef double *angles = <double *>malloc(n_links * sizeof(double))
    cdef double *points = <double *>malloc(2 * n_links * sizeof(double))
    cdef int *ordered = <int *>malloc(n_links * sizeof(int))
    cdef int *link_buffer = <int *>malloc(n_links * sizeof(int))

    try:
      for patch in range(n_patches):
          offset = offset_to_patch[patch]
          n_links = offset_to_patch[patch + 1] - offset

          for i in range(n_links):
              link = links_at_patch[offset + i]
              points[2 * i] = xy_of_link[link][0]
              points[2 * i + 1] = xy_of_link[link][1]

          argsort_by_angle_around_centroid(points, n_links, ordered)

          for i in range(n_links):
              link_buffer[i] = links_at_patch[offset + ordered[i]]

          for i in range(n_links):
              links_at_patch[offset + i] = link_buffer[i]

    finally:
        free(link_buffer)
        free(ordered)
        free(points)
        free(angles)


@cython.boundscheck(False)
def argsort_spoke_angles(np.ndarray[double, ndim=2, mode="c"] points,
                         np.ndarray[int, ndim=1] out):
    """Sort spokes by angle around a hub.

    Parameters
    ----------
    points : ndarray of float, shape `(n_points, 2)`
        Coordinates of points as (*x*, *y*).
    out : ndarray of int, shape `(n_points, )`
        Indices of sorted points.

    Returns
    -------
    ndarray of int, shape `(n_points, )`
        Indices of sorted points.
    """
    argsort_by_angle_around_centroid(&points[0, 0], points.shape[0], &out[0])
    return out
