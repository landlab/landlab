import numpy as np
cimport numpy as np
cimport cython


DTYPE = np.int
ctypedef np.int_t DTYPE_t


@cython.boundscheck(False)
def fill_links_at_patch(np.ndarray[DTYPE_t, ndim=1] links_at_patch,
                        np.ndarray[DTYPE_t, ndim=1] offset_to_patch,
                        np.ndarray[DTYPE_t, ndim=2] out):
    cdef int i
    cdef int link
    cdef int patch
    cdef int offset
    cdef int n_links
    cdef int n_patches = len(offset_to_patch) - 1

    for patch in range(n_patches):
        offset = offset_to_patch[patch]
        n_links = offset_to_patch[patch + 1] - offset

        link = 0
        for i in range(offset, offset + n_links):
          out[patch, link] = links_at_patch[i]
          link += 1
