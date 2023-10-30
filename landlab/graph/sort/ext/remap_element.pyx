import numpy as np

cimport cython
cimport numpy as np
from libc.stdlib cimport free, malloc


cdef extern from "math.h":
    double atan2(double y, double x) nogil


from .spoke_sort import sort_spokes_at_wheel

from .argsort cimport argsort_int

DTYPE = int
ctypedef np.int_t DTYPE_t
ctypedef np.uint8_t uint8


@cython.boundscheck(False)
def reverse_one_to_one(np.ndarray[DTYPE_t, ndim=1] mapping,
                       np.ndarray[DTYPE_t, ndim=1] out):
    cdef int n_elements = mapping.shape[0]
    cdef int index
    cdef int id_

    for index in range(n_elements):
        id_ = mapping[index]
        if id_ >= 0:
          out[id_] = index


@cython.boundscheck(False)
def reverse_one_to_many(np.ndarray[DTYPE_t, ndim=2] mapping,
                        np.ndarray[DTYPE_t, ndim=2] out):
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
def remap_graph_element(np.ndarray[DTYPE_t, ndim=1] elements,
                        np.ndarray[DTYPE_t, ndim=1] old_to_new):
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
def remap_graph_element_ignore(np.ndarray[DTYPE_t, ndim=1] elements,
                               np.ndarray[DTYPE_t, ndim=1] old_to_new,
                               DTYPE_t bad_val):
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
def reorder_patches(np.ndarray[DTYPE_t, ndim=1] links_at_patch,
                    np.ndarray[DTYPE_t, ndim=1] offset_to_patch,
                    np.ndarray[DTYPE_t, ndim=1] sorted_patches):
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
def calc_center_of_patch(np.ndarray[DTYPE_t, ndim=1] links_at_patch,
                         np.ndarray[DTYPE_t, ndim=1] offset_to_patch,
                         np.ndarray[np.float_t, ndim=2] xy_at_link,
                         np.ndarray[np.float_t, ndim=2] xy_at_patch):
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
def reorder_links_at_patch(np.ndarray[DTYPE_t, ndim=1] links_at_patch,
                           np.ndarray[DTYPE_t, ndim=1] offset_to_patch,
                           np.ndarray[np.float_t, ndim=2] xy_of_link):
    cdef int n_patches = len(offset_to_patch) - 1

    xy_of_patch = np.empty((n_patches, 2), dtype=float)
    calc_center_of_patch(links_at_patch, offset_to_patch, xy_of_link,
                         xy_of_patch)
    sort_spokes_at_wheel(links_at_patch, offset_to_patch, xy_of_patch,
                         xy_of_link)


cdef reverse_order(long * array, long size):
    cdef long i
    cdef long temp

    for i in range(size // 2):
        temp = array[i]
        array[i] = array[(size - 1) - i]
        array[(size - 1) - i] = temp


@cython.boundscheck(False)
def reverse_element_order(
    np.ndarray[long, ndim=2] links_at_patch,
    np.ndarray[long, ndim=1] patches
):
    cdef long n_patches = patches.shape[0]
    cdef long max_links = links_at_patch.shape[1]
    cdef long patch
    cdef long n
    cdef long i

    for i in range(n_patches):
        patch = patches[i]
        n = 1
        while n < max_links:
            if links_at_patch[patch, n] == -1:
                break
            n += 1
        reverse_order(&links_at_patch[patch, 1], n - 1)


cdef _offset_to_sorted_blocks(
    DTYPE_t *array,
    long len,
    long stride,
    DTYPE_t *offset,
    long n_values,
):
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
    np.ndarray[DTYPE_t, ndim=2, mode="c"] sorted_ids not None,
    np.ndarray[DTYPE_t, ndim=1, mode="c"] offset_to_block not None,
):
    cdef long n_ids = sorted_ids.shape[0]
    cdef long n_blocks = offset_to_block.shape[0]

    _offset_to_sorted_blocks(
        &sorted_ids[0, 0],
        n_ids,
        sorted_ids.shape[1],
        &offset_to_block[0],
        n_blocks,
    )


@cython.boundscheck(False)
@cython.wraparound(False)
def pair_isin(
    np.ndarray[DTYPE_t, ndim=2, mode="c"] src_pairs not None,
    np.ndarray[DTYPE_t, ndim=2, mode="c"] pairs not None,
    # np.ndarray[uint8, ndim=1, mode="c"] out not None,
    np.ndarray[uint8, ndim=1, mode="c", cast=True] out not None,
):
    cdef long n
    cdef long pair
    cdef long n_pairs = pairs.shape[0]
    cdef long n_values = src_pairs.shape[0]
    cdef DTYPE_t *data = <DTYPE_t *>malloc(n_values * sizeof(DTYPE_t))
    cdef SparseMatrixInt mat

    for n in range(n_values):
        data[n] = 1
    try:
        mat = sparse_matrix_alloc_with_tuple(&src_pairs[0, 0], data, n_values, 0)
        for pair in range(n_pairs):
            out[pair] = sparse_matrix_get_or_transpose(
                mat, pairs[pair, 0], pairs[pair, 1]
            )
    finally:
        free(data)


@cython.boundscheck(False)
@cython.wraparound(False)
def map_pairs_to_values(
    np.ndarray[DTYPE_t, ndim=2, mode="c"] src_pairs not None,
    np.ndarray[DTYPE_t, ndim=1, mode="c"] data not None,
    np.ndarray[DTYPE_t, ndim=2, mode="c"] pairs not None,
    np.ndarray[DTYPE_t, ndim=1, mode="c"] out not None,
):
    cdef long pair
    cdef long n_pairs = out.shape[0]
    cdef long n_values = data.shape[0]
    cdef long val
    cdef SparseMatrixInt mat

    mat = sparse_matrix_alloc_with_tuple(&src_pairs[0, 0], &data[0], n_values, -1)

    for pair in range(n_pairs):
        out[pair] = sparse_matrix_get_or_transpose(mat, pairs[pair, 0], pairs[pair, 1])


@cython.boundscheck(False)
@cython.wraparound(False)
def map_rolling_pairs_to_values(
    np.ndarray[DTYPE_t, ndim=2, mode="c"] src_pairs not None,
    np.ndarray[DTYPE_t, ndim=1, mode="c"] data not None,
    np.ndarray[DTYPE_t, ndim=2, mode="c"] pairs not None,
    np.ndarray[DTYPE_t, ndim=1, mode="c"] size_of_row not None,
    np.ndarray[DTYPE_t, ndim=2, mode="c"] out not None,
):
    cdef long n_values = data.shape[0]
    cdef long n_pairs = pairs.shape[0]
    cdef long pair
    cdef SparseMatrixInt mat

    mat = sparse_matrix_alloc_with_tuple(&src_pairs[0, 0], &data[0], n_values, -1)

    for pair in range(n_pairs):
        _map_rolling_pairs(mat, &pairs[pair, 0], &out[pair, 0], size_of_row[pair])


cdef _map_rolling_pairs(SparseMatrixInt mat, DTYPE_t *pairs, DTYPE_t *out, long size):
    cdef long n
    cdef long val

    if size > 0:
        for n in range(size - 1):
            out[n] = sparse_matrix_get_or_transpose(mat, pairs[n], pairs[n + 1])

        n = size - 1
        out[n] = sparse_matrix_get_or_transpose(mat, pairs[n], pairs[0])


cdef struct SparseMatrixInt:
    DTYPE_t *values
    long n_values
    DTYPE_t *offset_to_row
    DTYPE_t *col
    long col_start
    long col_stride
    long n_rows
    long n_cols
    long no_val


cdef SparseMatrixInt sparse_matrix_alloc_with_tuple(
    DTYPE_t *rows_and_cols,
    DTYPE_t *values,
    long n_values,
    long no_val,
):
    cdef long n_rows
    cdef long n_cols
    cdef long max_row = 0
    cdef long max_col = 0
    cdef long i
    cdef SparseMatrixInt mat
    cdef DTYPE_t *offset

    for i in range(0, n_values * 2, 2):
        if rows_and_cols[i] > max_row:
            max_row = rows_and_cols[i]
        if rows_and_cols[i + 1] > max_col:
            max_col = rows_and_cols[i + 1]
    n_rows = max_row + 1
    n_cols = max_col + 1

    offset = <DTYPE_t *>malloc((n_rows + 1) * sizeof(DTYPE_t))

    _offset_to_sorted_blocks(rows_and_cols, n_values, 2, offset, n_rows + 1)

    mat.values = values
    mat.n_values = n_values
    mat.offset_to_row = offset
    mat.col = rows_and_cols
    mat.col_start = 1
    mat.col_stride = 2
    mat.n_rows = n_rows
    mat.n_cols = n_cols
    mat.no_val = no_val

    return mat


cdef sparse_matrix_free(SparseMatrixInt mat):
    free(mat.offset_to_row)


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
        # if i >= mat.n_values * 2:
        #     print("error")

        if mat.col[i] == col:
            # if n >= mat.n_values:
            #     print("ERROR")
            return mat.values[n]
        i += mat.col_stride

    return mat.no_val
