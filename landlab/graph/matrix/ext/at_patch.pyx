cimport cython
from cython.parallel cimport prange


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_links_at_patch(
    cython.integral [:] links_at_patch,
    cython.integral [:] offset_to_patch,
    cython.integral [:, :] out,
):
    cdef int i
    cdef int link
    cdef int patch
    cdef int offset
    cdef int n_links
    cdef int n_patches = len(offset_to_patch) - 1

    for patch in prange(n_patches, nogil=True, schedule="static"):
        offset = offset_to_patch[patch]
        n_links = offset_to_patch[patch + 1] - offset

        link = 0
        for i in range(offset, offset + n_links):
            out[patch, link] = links_at_patch[i]
            link = link + 1
