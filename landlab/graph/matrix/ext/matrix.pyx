import numpy as np
cimport numpy as np
cimport cython

from libc.string cimport memcpy, memmove
from libc.stdlib cimport malloc, free

DTYPE = np.int
ctypedef np.int_t DTYPE_t


cdef roll(void * values, size_t n_values, size_t size, long shift):
    cdef size_t offset
    cdef void * src = values
    cdef void * dst
    cdef void * end
    cdef void * buff

    if shift == 0:
        return
    elif shift < 0:
        offset = (n_values + shift) % n_values
    else:
        offset = shift % n_values

    dst = <char *>values + offset * size
    end = <char *>src + (n_values - offset) * size
    buff = malloc(offset * size)

    memcpy(buff, end, offset * size)
    memmove(dst, src, (n_values - offset) * size)
    memcpy(src, buff, offset * size)

    free(buff)


@cython.boundscheck(False)
def roll_id_matrix_rows(np.ndarray[DTYPE_t, ndim=2] matrix,
                        np.ndarray[DTYPE_t, ndim=1] shift):
    cdef int n_rows = matrix.shape[0]
    cdef int n_cols = matrix.shape[1]
    cdef int row
    cdef int n

    for row in range(n_rows):
        n = 0
        while n < n_cols and matrix[row, n] != -1:
            n += 1

        roll(&matrix[row, 0], n, matrix.itemsize, shift[row])
