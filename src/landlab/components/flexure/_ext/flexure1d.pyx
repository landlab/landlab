cimport cython
from cython.parallel cimport prange
from libc.math cimport cos
from libc.math cimport exp
from libc.math cimport fabs
from libc.math cimport sin


@cython.boundscheck(False)
@cython.wraparound(False)
def subside_load_1d(
    cython.floating [:] x,
    cython.floating [:, :] loads,
    double alpha,
    double rigidity,
    cython.floating [:, :] out,
):
    cdef long n_rows = loads.shape[0]
    cdef long row

    for row in prange(n_rows, nogil=True, schedule="static"):
        _subside_row(x, loads[row], alpha, rigidity, out[row])


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _subside_row(
    cython.floating [:] x,
    cython.floating [:] loads,
    double alpha,
    double rigidity,
    cython.floating [:] out,
) noexcept nogil:
    cdef long n_points = len(x)
    cdef long col

    for col in range(n_points):
        _subside_point_load(x, x[col], loads[col], alpha, rigidity, out)


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _subside_point_load(
    cython.floating [:] x,
    double x_at_load,
    double load,
    double alpha,
    double rigidity,
    cython.floating [:] out,
) noexcept nogil:
    cdef long n_points = len(x)
    cdef float c = load * alpha ** 3. / (8. * rigidity)
    cdef float dx
    cdef long col

    for col in range(n_points):
        dx = fabs(x[col] - x_at_load) / alpha
        out[col] += c * exp(- dx) * (cos(dx) + sin(dx))
