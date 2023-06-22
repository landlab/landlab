import numpy as np

cimport cython
cimport numpy as np
from cython.parallel cimport prange

DTYPE = int
ctypedef np.int_t DTYPE_t


@cython.boundscheck(False)
def neighbors_at_link(
    np.ndarray[DTYPE_t, ndim=1] links,
    shape,
    np.ndarray[DTYPE_t, ndim=2] out
):
  cdef int stride
  cdef int n_links
  cdef int link
  cdef int i
  cdef is_top, is_bottom, is_left, is_right

  stride = 2 * shape[1] - 1
  n_links = (shape[0] - 1) * shape[1] + shape[0] * (shape[1] - 1)

  for i in range(links.shape[0]):
    link = links[i]

    is_top = link > (n_links - stride)
    is_bottom = link < stride
    is_left = link % stride == 0 or (link + shape[1]) % stride == 0
    is_right = (link - (shape[1] - 2)) % stride == 0 or (link + 1) % stride == 0

    if not is_right:
      out[i, 0] = link + 1

    if not is_top:
      out[i, 1] = link + stride

    if not is_left:
      out[i, 2] = link - 1

    if not is_bottom:
      out[i, 3] = link - stride



@cython.boundscheck(False)
@cython.wraparound(False)
def fill_parallel_links_at_link(
    cython.integral[:, :] out,
    shape,
):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int links_per_row = 2 * shape[1] - 1
    cdef int row, col
    cdef int link

    # Horizontal links
    for row in prange(n_rows, nogil=True, schedule="static"):
        link = row * links_per_row
        out[link, 0] = -1
        out[link, 1] = link + 1
        link = link + 1
        for col in range(1, n_cols - 2):
            out[link, 0] = link - 1
            out[link, 1] = link + 1
            link = link + 1
        out[link, 0] = link - 1
        out[link, 1] = -1

    # First row of vertical links
    link = n_cols - 1
    for col in range(n_cols):
        out[link, 0] = -1
        out[link, 1] = link + links_per_row
        link = link + 1

    # Vertical links
    for row in prange(1, n_rows - 1, nogil=True, schedule="static"):
        link = row * links_per_row + n_cols - 1
        for col in range(n_cols):
            out[link, 0] = link - links_per_row
            out[link, 1] = link + links_per_row
            link = link + 1

    # Last row of vertical links
    link = (n_rows - 2) * links_per_row + n_cols - 1
    for col in range(n_cols):
        out[link, 0] = link - links_per_row
        out[link, 1] = -1
        link = link + 1


@cython.boundscheck(False)
@cython.wraparound(False)
def sum_parallel_links(
    cython.numeric[:] out,
    const cython.numeric[:] value_at_link,
    shape,
):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int links_per_row = 2 * shape[1] - 1
    cdef int row, col
    cdef int link

    for row in prange(n_rows, nogil=True, schedule="static"):
        link = row * links_per_row + 1
        for col in range(1, n_cols - 2):
            out[link] = value_at_link[link - 1] + value_at_link[link + 1]
            link = link + 1

    for row in prange(1, n_rows - 2, nogil=True, schedule="static"):
        link = row * links_per_row + n_cols - 1
        for col in range(n_cols):
            out[link] = value_at_link[link - links_per_row] + value_at_link[link + links_per_row]
            link = link + 1
