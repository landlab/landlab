cimport cython
from libc.stdint cimport uint8_t

from cython.parallel import prange

ctypedef fused data_t:
    cython.integral
    long long
    cython.floating


@cython.boundscheck(False)
@cython.wraparound(False)
def _unravel(
    data_t[::1] data,
    const Py_ssize_t[::1] offset,
    data_t[:, ::1] out,
):
    cdef Py_ssize_t n_rows = offset.shape[0] - 1
    cdef Py_ssize_t col
    cdef Py_ssize_t row
    cdef Py_ssize_t start
    cdef Py_ssize_t stop

    for row in prange(n_rows, nogil=True, schedule="static"):
        start = offset[row]
        stop = offset[row + 1]
        for col in range(stop - start):
            out[row, col] = data[start + col]


@cython.boundscheck(False)
@cython.wraparound(False)
def _padded_row_contains(
    const Py_ssize_t[:, ::1] padded_array,
    const Py_ssize_t[::1] size_of_row,
    const Py_ssize_t bad_val,
    uint8_t[::1] out,
):
    cdef Py_ssize_t n_rows = padded_array.shape[0]
    cdef Py_ssize_t col
    cdef Py_ssize_t row

    for row in prange(n_rows, nogil=True, schedule="static"):
        out[row] = 0
        for col in range(size_of_row[row]):
            if padded_array[row, col] == bad_val:
                out[row] = 1
                break
