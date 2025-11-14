cimport cython
from cython.parallel cimport prange
from libc.math cimport copysign
from libc.math cimport fabs
from libc.math cimport fmin
from libc.math cimport sqrt
from libc.stdint cimport int32_t
from libc.stdint cimport int64_t

ctypedef fused id_t:
    int32_t
    int64_t


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def adjust_supercritical_discharge(
    cython.floating [::1] q_at_link,
    const cython.floating [::1] h_at_link,
    *,
    const double g,
    const double froude,
    const id_t [::1] where,
):
    cdef long n_links = len(where)
    cdef long i
    cdef long link
    cdef double factor = froude * sqrt(g)
    cdef double max_q
    cdef double q
    cdef double h

    for i in prange(n_links, nogil=True, schedule="static"):
        link = where[i]
        q = q_at_link[link]
        h = h_at_link[link]

        max_q = factor * h * sqrt(h)
        q_at_link[link] = copysign(fmin(fabs(q), max_q), q)
