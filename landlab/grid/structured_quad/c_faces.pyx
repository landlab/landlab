import numpy as np
cimport numpy as np
cimport cython


DTYPE = np.int
ctypedef np.int_t DTYPE_t


@cython.boundscheck(False)
def _face_at_link(shape, np.ndarray[DTYPE_t, ndim=1] out):
  cdef int stride
  cdef int link
  cdef int face
  cdef int row
  cdef int col

  stride = 2 * shape[1] - 1
  n_rows = shape[0] - 1

  face = 0
  for col in range(shape[1], stride - 1):
    out[col] = face
    face += 1

  for row in range(1, n_rows):
    link = row * stride

    for col in range(0, shape[1] - 1):
      out[link + col] = face
      face += 1

    link = row * stride + shape[1] - 1
    for col in range(1, shape[1] - 1):
      out[link + col] = face
      face += 1
