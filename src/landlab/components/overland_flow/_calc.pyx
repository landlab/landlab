cimport cython
from cython.parallel cimport prange
from libc.math cimport fabs
from libc.math cimport pow
from libc.stdint cimport int32_t
from libc.stdint cimport int64_t

ctypedef fused id_t:
    int32_t
    int64_t


@cython.boundscheck(False)
@cython.wraparound(False)
def calc_discharge_at_links(
    const cython.floating [::1] q_at_link,
    const cython.floating [::1] q_mean_at_link,
    const cython.floating [::1] h_at_link,
    const cython.floating [::1] water_slope_at_link,
    const cython.floating [::1] mannings_at_link,
    *,
    const double g,
    const double dt,
    const id_t [::1] where,
    cython.floating [::1] out,
):
    cdef Py_ssize_t n_links = where.shape[0]
    cdef Py_ssize_t link
    cdef Py_ssize_t i
    cdef double h_to_seven_thirds
    cdef double numerator
    cdef double denominator
    cdef double mannings
    cdef double h
    cdef double g_times_dt = g * dt
    cdef double seven_thirds = 7.0 / 3.0

    for i in prange(n_links, nogil=True, schedule="static"):
        link = where[i]

        if h_at_link[link] > 0.0:
            mannings = mannings_at_link[link]
            h = h_at_link[link]
            h_to_seven_thirds = pow(h, seven_thirds)

            numerator = (
                q_mean_at_link[link] - g_times_dt * h * water_slope_at_link[link]
            )

            denominator = (
                1.0
                + g_times_dt * mannings * mannings * fabs(q_at_link[link])
                / h_to_seven_thirds
            )

            out[link] = numerator / denominator
