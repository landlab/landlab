cimport cython

from cython.parallel import prange

from libc.stdint cimport uint8_t

ctypedef fused id_t:
    cython.integral
    long long


cdef int _id_array_contains(const id_t *array, long size, long bad_id) noexcept nogil:
    cdef long n
    for n in range(size):
        if array[n] == bad_id:
            return 1
    return 0


@cython.boundscheck(False)
@cython.wraparound(False)
def id_array_contains(
    const id_t [:, :] corners_at_cell,
    const id_t [:] n_corners_at_cell,
    long bad_val,
    uint8_t [:] out,
):
    cdef long n_cells = corners_at_cell.shape[0]
    cdef long cell

    for cell in prange(n_cells, nogil=True, schedule="static"):
        out[cell] = _id_array_contains(
            &corners_at_cell[cell, 0],
            n_corners_at_cell[cell],
            bad_val,
        )
