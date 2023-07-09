import numpy as np

cimport cython
cimport numpy as cnp
from cython.parallel cimport prange
from libc.stdlib cimport abs


@cython.boundscheck(False)
@cython.wraparound(False)
def subside_loads(
    cython.floating [:, :] w,
    const cython.floating [:, :] r,
    const long n_rows,
    const long n_cols,
    const cython.floating [:] loads,
    const cython.integral [:] row_of_load,
    const cython.integral [:] col_of_load,
    const long n_loads,
    const double alpha,
    const double gamma_mantle,
):
    cdef long load_row, load_col
    cdef long row, col
    cdef long d_row
    cdef long load
    cdef double c
    cdef double inv_c = 1. / (2. * np.pi * gamma_mantle * alpha ** 2.)

    for row in prange(n_rows, nogil=True, schedule="static"):
        for load in range(n_loads):
            c = loads[load] * inv_c
            load_row = row_of_load[load]
            load_col = col_of_load[load]
            d_row = abs(load_row - row)
            for col in range(n_cols):
                w[row, col] = w[row, col] - c * r[d_row, abs(load_col - col)]
