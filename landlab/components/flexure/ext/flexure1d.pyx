import numpy as np
cimport numpy as np
cimport cython

from libc.math cimport fabs, exp, cos, sin
# from libc.stdlib cimport abs


DTYPE = np.double
ctypedef np.double_t DTYPE_t


@cython.boundscheck(False)
def subside_load_1d(np.ndarray[DTYPE_t, ndim=1] x,
                    np.ndarray[DTYPE_t, ndim=2] loads,
                    DTYPE_t alpha,
                    DTYPE_t rigidity,
                    np.ndarray[DTYPE_t, ndim=2] out):
    cdef int n_rows = loads.shape[0]
    cdef int row

    for row in range(n_rows):
        subside_row(&x[0], len(x), &loads[row, 0], alpha, rigidity,
                    &out[row, 0])


@cython.cdivision(True)
cdef double * subside_row(double * x, int n_points, double * loads,
                          double alpha, double rigidity, double * out) nogil:
    cdef int col
    cdef double load

    for col in range(n_points):
        subside_point_load(x, n_points, x[col], loads[col], alpha, rigidity,
                           out)

    return out


@cython.cdivision(True)
cdef double * subside_point_load(double * x, int n_points, double x_at_load,
                                 double load, double alpha, double rigidity,
                                 double * out) nogil:
    cdef float c = load * alpha ** 3. / (8. * rigidity)
    cdef float dx
    cdef int col

    for col in range(n_points):
        dx = fabs(x[col] - x_at_load) / alpha
        out[col] += c * exp(- dx) * (cos(dx) + sin(dx))

    return out
