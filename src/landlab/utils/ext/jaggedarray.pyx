cimport cython

from cython.parallel import prange

# https://cython.readthedocs.io/en/stable/src/userguide/fusedtypes.html
ctypedef fused id_t:
    cython.integral
    long long


ctypedef fused float_or_int:
    cython.integral
    long long
    cython.floating


@cython.boundscheck(False)
@cython.wraparound(False)
def unravel(
    float_or_int [:] data,
    const id_t [:] offset,
    float_or_int [:, :] out,
):
    cdef long n_rows = len(offset) - 1
    cdef long col
    cdef long row
    cdef long offset_to_row
    cdef long values_per_row

    for row in prange(n_rows, nogil=True, schedule="static"):
        offset_to_row = offset[row]
        values_per_row = offset[row + 1] - offset[row]
        for col in range(values_per_row):
            out[row, col] = data[offset_to_row + col]
