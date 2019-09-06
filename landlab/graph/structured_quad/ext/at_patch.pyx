import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.int
ctypedef np.int_t DTYPE_t


@cython.boundscheck(False)
def fill_links_at_patch(shape, np.ndarray[DTYPE_t, ndim=2] links_at_patch):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int links_per_row = 2 * n_cols - 1
    cdef int patches_per_row = n_cols - 1
    cdef int row
    cdef int link
    cdef int patch
    cdef int col

    for row in range(n_rows - 1):
        link = row * links_per_row + n_cols
        patch = row * patches_per_row
        for col in range(n_cols - 1):
            links_at_patch[patch, 0] = link
            links_at_patch[patch, 1] = link + n_cols - 1
            links_at_patch[patch, 2] = link - 1
            links_at_patch[patch, 3] = link - n_cols

            patch += 1
            link += 1
