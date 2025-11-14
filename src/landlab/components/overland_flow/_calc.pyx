cimport cython
from cython.parallel cimport prange
from libc.math cimport copysign
from libc.math cimport fabs
from libc.math cimport fmin
from libc.stdint cimport int32_t
from libc.stdint cimport int64_t

ctypedef fused id_t:
    int32_t
    int64_t


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def adjust_unstable_discharge(
    cython.floating [::1] q_at_link,
    const cython.floating [::1] h_at_link,
    *,
    const double dx,
    const double dt,
    const id_t [::1] where,
    cython.floating [::1] out,
):
    cdef long n_links = len(where)
    cdef long i
    cdef long link
    cdef double dx_over_dt = dx / dt
    cdef double factor_of_safety = 0.2
    cdef double q_threshold
    cdef double q
    cdef double h

    for i in prange(n_links, nogil=True, schedule="static"):
        link = where[i]
        q = q_at_link[link]
        h = h_at_link[link]

        q_threshold = (h / 4.0) * dx_over_dt * factor_of_safety
        out[link] = copysign(fmin(fabs(q), q_threshold), q)

    return (<object>out).base
