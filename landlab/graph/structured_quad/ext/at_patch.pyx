cimport cython
from cython.parallel cimport prange

ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_links_at_patch(
    shape,
    id_t [:, :] links_at_patch,
):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int links_per_row = 2 * n_cols - 1
    cdef int patches_per_row = n_cols - 1
    cdef int row
    cdef int link
    cdef int patch
    cdef int col

    for row in prange(n_rows - 1, nogil=True, schedule="static"):
        link = row * links_per_row + n_cols
        patch = row * patches_per_row
        for col in range(n_cols - 1):
            links_at_patch[patch, 0] = link
            links_at_patch[patch, 1] = link + n_cols - 1
            links_at_patch[patch, 2] = link - 1
            links_at_patch[patch, 3] = link - n_cols

            patch = patch + 1
            link = link + 1
