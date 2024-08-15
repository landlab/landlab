cimport cython

import numpy as np
from cython.parallel import prange

cimport numpy as np
from libc.stdlib cimport free
from libc.stdlib cimport malloc


cdef extern from "math.h":
    double atan2(double y, double x) nogil


from .spoke_sort import sort_spokes_at_wheel

ctypedef np.int64_t DTYPE_t
ctypedef np.uint8_t uint8

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

    for i in range(n_elements):
        elements[i] = old_to_new[elements[i]]


@cython.boundscheck(False)
@cython.wraparound(False)
def remap_graph_element_ignore(
    id_t [:] elements,
    const id_t [:] old_to_new,
    int bad_val,
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
    cdef int n_elements = elements.shape[0]
    cdef int i

    for i in range(n_elements):
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


cdef void _offset_to_sorted_blocks(
    const id_t [:] array,
    # id_t *array,
    long len,
    long stride,
    integral_out_t [:] offset,
    # long long [:] offset,
    # id_t *offset,
    long n_values,
) noexcept:
    cdef long i
    cdef long value
    cdef long first_non_negative

    first_non_negative = len * stride
    for i in range(0, len * stride, stride):
        if array[i] >= 0:
            first_non_negative = i
            break

    offset[0] = first_non_negative
    for value in range(1, n_values):
        offset[value] = offset[value - 1]
        for i in range(offset[value], len * stride, stride):
            if array[i] >= value:
                offset[value] = i
                break
        else:
            offset[value] = len * stride

    for value in range(n_values):
        offset[value] = offset[value] // stride


@cython.boundscheck(False)
@cython.wraparound(False)
def offset_to_sorted_block(
    id_t [:, :] sorted_ids,
    integral_out_t [:] offset_to_block,
    # np.ndarray[DTYPE_t, ndim=2, mode="c"] sorted_ids not None,
    # np.ndarray[DTYPE_t, ndim=1, mode="c"] offset_to_block not None,
):
    cdef long n_ids = sorted_ids.shape[0]
    cdef long n_blocks = offset_to_block.shape[0]
    # cdef id_t [:] sorted_ids_flat = sorted_ids
    # cdef id_t [:] sorted_ids_flat = sorted_ids.reshape((-1,))
    cdef id_t *ptr = &sorted_ids[0, 0]
    cdef id_t [:] sorted_ids_flat= <id_t[:sorted_ids.size]>ptr

    _offset_to_sorted_blocks(
        sorted_ids_flat,
        n_ids,
        sorted_ids.shape[1],
        offset_to_block,
        n_blocks,
    )


@cython.boundscheck(False)
@cython.wraparound(False)
def pair_isin(
    # np.ndarray[DTYPE_t, ndim=2, mode="c"] src_pairs not None,
    # np.ndarray[DTYPE_t, ndim=2, mode="c"] pairs not None,
    # # np.ndarray[uint8, ndim=1, mode="c"] out not None,
    # np.ndarray[uint8, ndim=1, mode="c", cast=True] out not None,

    const id_t [:, :] src_pairs,
    const id_t [:, :] pairs,
    uint8 [:] out,
):
    # cdef long n
    cdef long pair
    cdef long n_pairs = pairs.shape[0]
    cdef long n_values = src_pairs.shape[0]
    # cdef long *data = <long *>malloc(n_values * sizeof(long))
    cdef long [:] data = np.full(n_values, 1, dtype=long)
    # cdef const id_t [:] src_pairs_flat = src_pairs.reshape((-1,))
    cdef const id_t *ptr = &src_pairs[0, 0]
    cdef const id_t [:] src_pairs_flat = <const id_t[:src_pairs.size]>ptr
    cdef SparseMatrixInt mat

    # for n in range(n_values):
    #     data[n] = 1
    try:
        # mat = sparse_matrix_alloc_with_tuple(&src_pairs[0, 0], data, n_values, 0)
        mat = sparse_matrix_alloc_with_tuple(src_pairs_flat, data, n_values, 0)
        for pair in range(n_pairs):
            out[pair] = sparse_matrix_get_or_transpose(
                mat, pairs[pair, 0], pairs[pair, 1]
            )
    finally:
        pass
        # free(data)


@cython.boundscheck(False)
@cython.wraparound(False)
def map_pairs_to_values(
    # np.ndarray[DTYPE_t, ndim=2, mode="c"] src_pairs not None,
    # np.ndarray[DTYPE_t, ndim=1, mode="c"] data not None,
    # np.ndarray[DTYPE_t, ndim=2, mode="c"] pairs not None,
    # np.ndarray[DTYPE_t, ndim=1, mode="c"] out not None,
    const id_t [:, :] src_pairs,
    const integral_out_t [:] data,
    const id_t [:, :] pairs,
    integral_out_t [:] out,
):
    cdef long pair
    cdef long n_pairs = out.shape[0]
    cdef long n_values = data.shape[0]
    # cdef const id_t [:] src_pairs_flat = src_pairs.reshape((-1,))
    cdef const id_t *ptr = &src_pairs[0, 0]
    cdef const id_t [:] src_pairs_flat = <const id_t[:src_pairs.size]>ptr
    cdef SparseMatrixInt mat

    # mat = sparse_matrix_alloc_with_tuple(&src_pairs[0, 0], &data[0], n_values, -1)
    mat = sparse_matrix_alloc_with_tuple(src_pairs_flat, data, n_values, -1)

    for pair in range(n_pairs):
        out[pair] = sparse_matrix_get_or_transpose(mat, pairs[pair, 0], pairs[pair, 1])


@cython.boundscheck(False)
@cython.wraparound(False)
def map_rolling_pairs_to_values(
    # np.ndarray[DTYPE_t, ndim=2, mode="c"] src_pairs not None,
    # np.ndarray[DTYPE_t, ndim=1, mode="c"] data not None,
    # np.ndarray[DTYPE_t, ndim=2, mode="c"] pairs not None,
    # np.ndarray[DTYPE_t, ndim=1, mode="c"] size_of_row not None,
    # np.ndarray[DTYPE_t, ndim=2, mode="c"] out not None,
    const id_t [:, :] src_pairs,
    const integral_out_t [:] data,
    const id_t [:, :] pairs,
    const id_t [:] size_of_row,
    integral_out_t [:, :] out,
):
    cdef long n_values = data.shape[0]
    cdef long n_pairs = pairs.shape[0]
    cdef long pair
    # cdef const id_t [:] src_pairs_flat = src_pairs.ravel()
    cdef const id_t *ptr = &src_pairs[0, 0]
    cdef const id_t [:] src_pairs_flat = <const id_t[:2 * len(src_pairs)]>ptr
    cdef SparseMatrixInt mat

    # mat = sparse_matrix_alloc_with_tuple(&src_pairs[0, 0], &data[0], n_values, -1)
    mat = sparse_matrix_alloc_with_tuple(src_pairs_flat, data, n_values, -1)

    for pair in range(n_pairs):
        _map_rolling_pairs(mat, pairs[pair, :], out[pair, :], size_of_row[pair])


cdef _map_rolling_pairs(
    SparseMatrixInt mat,
    const id_t [:] pairs,
    integral_out_t [:] out,
    long size,
):
    cdef long n

    if size > 0:
        for n in range(size - 1):
            out[n] = sparse_matrix_get_or_transpose(mat, pairs[n], pairs[n + 1])

        n = size - 1
        out[n] = sparse_matrix_get_or_transpose(mat, pairs[n], pairs[0])


cdef struct SparseMatrixInt:
    long long *values
    long n_values
    long long *offset_to_row
    long long *col
    long col_start
    long col_stride
    long n_rows
    long n_cols
    long no_val


cdef SparseMatrixInt sparse_matrix_alloc_with_tuple(
    # DTYPE_t *rows_and_cols,
    # DTYPE_t *values,
    # long n_values,
    # long no_val,
    const id_t [:] rows_and_cols,
    const integral_out_t [:] values,
    long n_values,
    long no_val,
):
    cdef long n_rows
    cdef long n_cols
    cdef long max_row = 0
    cdef long max_col = 0
    cdef long i
    cdef SparseMatrixInt mat
    # cdef long long *offset
    cdef long long [:] _col
    cdef long long [:] _values
    # cdef long long * offset = <long long *>malloc((n_rows + 1) * sizeof(long long))
    # cdef long long [:] offset_view = offset
    cdef long long [:] offset

    for i in range(0, n_values * 2, 2):
        if rows_and_cols[i] > max_row:
            max_row = rows_and_cols[i]
        if rows_and_cols[i + 1] > max_col:
            max_col = rows_and_cols[i + 1]
    n_rows = max_row + 1
    n_cols = max_col + 1

    offset = np.empty(n_rows + 1, dtype=long)
    _col = np.array(rows_and_cols, dtype=long)
    _values = np.array(values, dtype=long)
    # offset = <long long *>malloc((n_rows + 1) * sizeof(long long))
    # col = <long long *>malloc((2 * n_values) * sizeof(long long))
    # values = <long long *>malloc(n_values * sizeof(long long))

    _offset_to_sorted_blocks(rows_and_cols, n_values, 2, offset, n_rows + 1)

    mat.values = &_values[0]
    mat.n_values = n_values
    mat.offset_to_row = &offset[0]
    # mat.col = rows_and_cols
    mat.col = &_col[0]
    mat.col_start = 1
    mat.col_stride = 2
    mat.n_rows = n_rows
    mat.n_cols = n_cols
    mat.no_val = no_val

    return mat


cdef sparse_matrix_free(SparseMatrixInt mat):
    free(mat.offset_to_row)
    free(mat.col)
    free(mat.values)


cdef long sparse_matrix_get_or_transpose(SparseMatrixInt mat, long row, long col):
    cdef long val
    val = sparse_matrix_get(mat, row, col)
    if val == mat.no_val:
        val = sparse_matrix_get(mat, col, row)
    return val


cdef long sparse_matrix_get(SparseMatrixInt mat, long row, long col):
    cdef long start
    cdef long stop
    cdef long n
    cdef long i

    if row < 0:
        return mat.no_val
    elif row >= mat.n_rows:
        return mat.no_val
    elif col < 0:
        return mat.no_val
    elif col >= mat.n_cols:
        return mat.no_val

    start = mat.offset_to_row[row]
    stop = mat.offset_to_row[row + 1]

    i = mat.col_start + start * mat.col_stride
    for n in range(start, stop):
        if mat.col[i] == col:
            return mat.values[n]
        i += mat.col_stride

    return mat.no_val
