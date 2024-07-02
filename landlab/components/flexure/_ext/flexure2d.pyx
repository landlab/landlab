cimport cython
from cython.parallel cimport prange
from libc.math cimport M_PI
from libc.stdlib cimport labs

ctypedef fused index_t:
    cython.integral
    long long
    unsigned int
    unsigned long
    unsigned long long


@cython.boundscheck(False)
@cython.wraparound(False)
def subside_loads(
    cython.floating [:, :] w,
    const cython.floating [:, :] r,
    const cython.floating [:] loads,
    const index_t [:] row_of_load,
    const index_t [:] col_of_load,
    const double alpha,
    const double gamma_mantle,
):
    cdef long n_rows = w.shape[0]
    cdef long n_cols = w.shape[1]
    cdef long n_loads = loads.shape[0]
    cdef long load_row
    cdef long load_col
    cdef long row
    cdef long col
    cdef long d_row
    cdef long load
    cdef double c
    cdef double inv_c = 1. / (2. * M_PI * gamma_mantle * alpha ** 2.)

    for row in prange(n_rows, nogil=True, schedule="static"):
        for load in range(n_loads):
            c = loads[load] * inv_c
            load_row = row_of_load[load]
            load_col = col_of_load[load]
            d_row = labs(load_row - row)
            for col in range(n_cols):
                w[row, col] = w[row, col] - c * r[d_row, abs(load_col - col)]
