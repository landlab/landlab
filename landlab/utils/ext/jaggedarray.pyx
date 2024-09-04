import numpy as np

cimport cython
cimport numpy as np
from libc.stdint cimport int64_t
from libc.string cimport memcpy


cdef void _pad_jaggedarray(
    void * data,
    int64_t * offset,
    size_t n_rows,
    size_t size,
    void * buff,
    size_t n_cols,
) noexcept nogil:
    cdef int row
    cdef void * dst = buff
    cdef void * src = data

    for row in range(n_rows):
        dst = <char *>buff + n_cols * row * size
        src = <char *>data + offset[row] * size

        memcpy(dst, src, size * (offset[row + 1] - offset[row]))


@cython.boundscheck(False)
@cython.wraparound(False)
def unravel(
    np.ndarray data,
    np.ndarray[int64_t, ndim=1, mode="c"] offset,
    np.ndarray out,
):
    _pad_jaggedarray(
        data.data, &offset[0], len(offset) - 1, data.itemsize, out.data, out.shape[1]
    )
