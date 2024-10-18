cimport cython

import numpy as np
from cython.parallel import prange

cimport numpy as np
from libc.stdlib cimport free
from libc.stdlib cimport malloc


cdef extern from "math.h":
    double atan2(double y, double x) nogil


from .spoke_sort import sort_spokes_at_wheel

# https://cython.readthedocs.io/en/stable/src/userguide/fusedtypes.html
ctypedef fused id_t:
    cython.integral
    long long

ctypedef fused integral_out_t:
    cython.integral
    long long


@cython.boundscheck(False)
@cython.wraparound(False)
def reverse_one_to_one(
    const id_t [:] mapping,
    id_t [:] out,
):
    cdef int n_elements = mapping.shape[0]
    cdef int index
    cdef int id_

    for index in prange(n_elements, nogil=True, schedule="static"):
        id_ = mapping[index]
        if id_ >= 0:
            out[id_] = index


@cython.boundscheck(False)
@cython.wraparound(False)
def reverse_one_to_many(
    const id_t [:, :] mapping,
    id_t [:, :] out,
):
    cdef int n_elements = mapping.shape[0]
    cdef int n_cols = mapping.shape[1]
    cdef int out_rows = out.shape[0]
    cdef int index
    cdef int id_
    cdef int *count = <int *>malloc(out_rows * sizeof(int))

    try:
        for index in range(out_rows):
            count[index] = 0

        for index in range(n_elements):
            for col in range(n_cols):
                id_ = mapping[index, col]
                if id_ >= 0:
                    out[id_, count[id_]] = index
                    count[id_] += 1
    finally:
        free(count)


@cython.boundscheck(False)
@cython.wraparound(False)
def remap_graph_element(
    id_t [:] elements,
    const id_t [:] old_to_new,
):
    """Remap elements in an array in place.

    Parameters
    ----------
    elements : ndarray of int
        Identifiers of elements.
    old_to_new : ndarray of int
        Mapping from the old identifier to the new identifier.
    """
    cdef int n_elements = elements.shape[0]
    cdef int i

    for i in prange(n_elements, nogil=True, schedule="static"):
        elements[i] = old_to_new[elements[i]]


@cython.boundscheck(False)
@cython.wraparound(False)
def remap_graph_element_ignore(
    id_t [:] elements,
    const id_t [:] old_to_new,
    long bad_val,
):
    """Remap elements in an array in place, ignoring bad values.

    Parameters
    ----------
    elements : ndarray of int
        Identifiers of elements.
    old_to_new : ndarray of int
        Mapping from the old identifier to the new identifier.
    bad_val : int
        Ignore values in the input array when remapping.
    """
    cdef long n_elements = elements.shape[0]
    cdef long i

    for i in prange(n_elements, nogil=True, schedule="static"):
        if elements[i] != bad_val:
            elements[i] = old_to_new[elements[i]]


@cython.boundscheck(False)
@cython.wraparound(False)
def reorder_patches(
    id_t [:] links_at_patch,
    id_t [:] offset_to_patch,
    const id_t [:] sorted_patches,
):
    cdef int i
    cdef int patch
    cdef int offset
    cdef int n_links
    cdef int n_patches = len(sorted_patches)
    cdef int *new_offset = <int *>malloc(len(offset_to_patch) * sizeof(int))
    cdef int *new_patches = <int *>malloc(len(links_at_patch) * sizeof(int))

    try:
        new_offset[0] = 0
        for patch in range(n_patches):
            offset = offset_to_patch[sorted_patches[patch]]
            n_links = offset_to_patch[sorted_patches[patch] + 1] - offset

            new_offset[patch + 1] = new_offset[patch] + n_links
            for i in range(n_links):
                new_patches[new_offset[patch] + i] = links_at_patch[offset + i]

        for i in range(len(links_at_patch)):
            links_at_patch[i] = new_patches[i]
        for i in range(len(offset_to_patch)):
            offset_to_patch[i] = new_offset[i]
    finally:
        free(new_offset)
        free(new_patches)


@cython.boundscheck(False)
@cython.wraparound(False)
def calc_center_of_patch(
    const id_t [:] links_at_patch,
    const id_t [:] offset_to_patch,
    const cython.floating [:, :] xy_at_link,
    cython.floating [:, :] xy_at_patch,
):
    cdef int patch
    cdef int link
    cdef int i
    cdef int offset
    cdef int n_links
    cdef int n_patches = len(xy_at_patch)
    cdef double x
    cdef double y

    for patch in range(n_patches):
        offset = offset_to_patch[patch]
        n_links = offset_to_patch[patch + 1] - offset
        x = 0.
        y = 0.
        for i in range(offset, offset + n_links):
            link = links_at_patch[i]
            x += xy_at_link[link, 0]
            y += xy_at_link[link, 1]
        xy_at_patch[patch, 0] = x / n_links
        xy_at_patch[patch, 1] = y / n_links


@cython.boundscheck(False)
@cython.wraparound(False)
def reorder_links_at_patch(
    id_t [:] links_at_patch,
    id_t [:] offset_to_patch,
    cython.floating [:, :] xy_of_link,
):
    cdef int n_patches = len(offset_to_patch) - 1

    xy_of_patch = np.empty((n_patches, 2), dtype=float)
    calc_center_of_patch(
        links_at_patch, offset_to_patch, xy_of_link, xy_of_patch
    )
    sort_spokes_at_wheel(
        links_at_patch, offset_to_patch, xy_of_patch, xy_of_link
    )


cdef void reverse_order(id_t * array, long size) noexcept nogil:
    cdef long i
    cdef long temp

    for i in range(size // 2):
        temp = array[i]
        array[i] = array[(size - 1) - i]
        array[(size - 1) - i] = temp


@cython.boundscheck(False)
@cython.wraparound(False)
def reverse_element_order(
    id_t [:, :] links_at_patch,
    const id_t [:] patches,
):
    cdef long n_patches = patches.shape[0]
    cdef long max_links = links_at_patch.shape[1]
    cdef long patch
    cdef long n
    cdef long i

    for i in prange(n_patches, nogil=True, schedule="static"):
        patch = patches[i]
        n = 1
        while n < max_links:
            if links_at_patch[patch, n] == -1:
                break
            n = n + 1
        reverse_order(&links_at_patch[patch, 1], n - 1)
