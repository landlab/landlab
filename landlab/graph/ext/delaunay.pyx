import numpy as np
cimport numpy as np
cimport cython
from libc.stdlib cimport malloc, free


DTYPE = np.int
ctypedef np.int_t DTYPE_t


@cython.boundscheck(False)
def _setup_links_at_patch(np.ndarray[DTYPE_t, ndim=2] nodes_at_patch,
                          np.ndarray[DTYPE_t, ndim=2] tri_neighbors,
                          np.ndarray[DTYPE_t, ndim=2] nodes_at_link,
                          np.ndarray[DTYPE_t, ndim=2] links_at_patch):
  cdef int i
  cdef int link
  cdef int neighbor
  cdef int n_patches = len(nodes_at_patch)
  cdef int *tri_done = <int *>malloc(n_patches * sizeof(int))
  cdef int *links_per_patch = <int *>malloc(n_patches * sizeof(int))

  if not tri_done or not links_per_patch:
    raise MemoryError(
      'unable to allocate {bytes} bytes'.format(bytes=n_patches * sizeof(int)))

  try:
    for tri in range(n_patches):
      tri_done[tri] = 0
      links_per_patch[tri] = 0

    link = 0
    for tri in range(n_patches):
      for i in (0, 1, 2):
        neighbor = tri_neighbors[tri, i]

        if neighbor == -1 or not tri_done[neighbor]:
          nodes_at_link[link, 0] = nodes_at_patch[tri, (i + 1) % 3]
          nodes_at_link[link, 1] = nodes_at_patch[tri, (i + 2) % 3]

          links_at_patch[tri, links_per_patch[tri]] = link
          links_per_patch[tri] += 1

          if neighbor >= 0:
            links_at_patch[neighbor, links_per_patch[neighbor]] = link
            links_per_patch[neighbor] += 1
            tri_done[tri] = True

          link += 1
  finally:
    free(links_per_patch)
    free(tri_done)
