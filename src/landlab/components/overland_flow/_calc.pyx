cimport cython
from cython.parallel cimport prange
from libc.stdint cimport int32_t
from libc.stdint cimport int64_t

ctypedef fused id_t:
    int32_t
    int64_t


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def adjust_discharge_for_dry_links(
    const cython.floating [::1] h_at_link,
    const cython.floating [::1] q_at_link,
    *,
    const id_t [::1] where,
    cython.floating [::1] out,
):
    cdef Py_ssize_t n_links = where.shape[0]
    cdef Py_ssize_t link
    cdef Py_ssize_t i

    for i in prange(n_links, nogil=True, schedule="static"):
        link = where[i]
        if h_at_link[link] <= 0.0:
            out[link] = 0.0

    return (<object>out).base
