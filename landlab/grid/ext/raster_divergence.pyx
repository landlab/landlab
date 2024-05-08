cimport cython
from cython.parallel cimport prange

ctypedef fused float_or_int:
    cython.integral
    cython.floating


@cython.boundscheck(False)
@cython.wraparound(False)
def _calc_flux_div_at_node(
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


@cython.boundscheck(False)
@cython.wraparound(False)
def _calc_net_face_flux_at_cell(
    shape,
    xy_spacing,
    const float_or_int[:] unit_flux_at_face,
    cython.floating[:] out,
):
    cdef double dx = xy_spacing[0]
    cdef double dy = xy_spacing[1]
    cdef int n_cols = shape[1] - 2
    cdef int n_rows = shape[0] - 2
    cdef int n_cells = n_rows * n_cols
    cdef int cell
    cdef int col
    cdef int face
    cdef int row

    for cell in prange(n_cells, nogil=True, schedule="static"):
        out[cell] = 0.0

    for row in prange(n_rows, nogil=True, schedule="static"):
        cell = row * n_cols
        face = n_cols + row * (2 * n_cols + 1)

        for col in range(n_cols):
            out[cell + col] = dy * (
                unit_flux_at_face[face] - unit_flux_at_face[face + 1]
            ) + dx * (
                unit_flux_at_face[face - n_cols] - unit_flux_at_face[face + n_cols + 1]
            )
            face = face + 1
