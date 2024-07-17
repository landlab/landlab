cimport cython
from cython.parallel cimport prange
from libc.stdlib cimport free
from libc.stdlib cimport malloc
from libc.string cimport memcpy
from libc.string cimport memmove


cdef void roll(
    void * values,
    size_t n_values,
    size_t size,
    long shift,
) noexcept nogil:
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
@cython.wraparound(False)
def roll_id_matrix_rows(
    cython.integral [:, :] matrix,
    cython.integral [:] shift,
):
    cdef int n_rows = matrix.shape[0]
    cdef int n_cols = matrix.shape[1]
    cdef int itemsize = matrix.itemsize
    cdef int row
    cdef int n

    for row in prange(n_rows, nogil=True, schedule="static"):
        n = 0
        while n < n_cols and matrix[row, n] != -1:
            n = n + 1

        roll(&matrix[row, 0], n, itemsize, shift[row])
