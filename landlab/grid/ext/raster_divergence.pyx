import numpy as np

cimport cython
cimport numpy as cnp
from cython.parallel cimport prange

ctypedef fused float_or_int:
    cython.integral
    cython.floating


@cython.boundscheck(False)
@cython.wraparound(False)
def calc_flux_div_at_node(
    shape,
    xy_spacing,
    const float_or_int[:] value_at_link,
    cython.floating[:] out,
):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef double dx = xy_spacing[0]
    cdef double dy = xy_spacing[1]
    cdef int links_per_row = 2 * n_cols - 1
    cdef double inv_area_of_cell = 1.0 / (dx * dy)
    cdef int row, col
    cdef int node, link

    for row in prange(1, n_rows - 1, nogil=True, schedule="static"):
        node = row * n_cols
        link = row * links_per_row

        for col in range(1, n_cols - 1):
            out[node + col] = (
                dy * (value_at_link[link + 1] - value_at_link[link])
                + dx * (value_at_link[link + n_cols] - value_at_link[link - n_cols + 1])
            ) * inv_area_of_cell
            link = link + 1
