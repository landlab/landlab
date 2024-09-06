cimport cython
from libc.stdint cimport uint8_t
from libc.stdlib cimport free
from libc.stdlib cimport malloc

# https://cython.readthedocs.io/en/stable/src/userguide/fusedtypes.html
ctypedef fused id_t:
    cython.integral
    long long

ctypedef fused integral_out_t:
    cython.integral
    long long


cdef struct SparseMatrixInt:
    long *values
    long n_values
    long *offset_to_row
    long *col
    long col_start
    long col_stride
    long n_rows
    long n_cols
    long no_val


cdef SparseMatrixInt sparse_matrix_alloc_with_tuple(
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
    cdef long * col_ptr
    cdef long * values_ptr
    cdef long * offset_ptr
    cdef long [:] offset_array

    for i in range(0, n_values * 2, 2):
        if rows_and_cols[i] > max_row:
            max_row = rows_and_cols[i]
        if rows_and_cols[i + 1] > max_col:
            max_col = rows_and_cols[i + 1]
    n_rows = max_row + 1
    n_cols = max_col + 1

    offset_ptr = <long *>malloc((n_rows + 1) * sizeof(long))
    offset_array = <long [:n_rows + 1]>offset_ptr

    col_ptr = <long *>malloc((2 * n_values) * sizeof(long))
    values_ptr = <long *>malloc(n_values * sizeof(long))
    for i in range(2 * n_values):
        col_ptr[i] = rows_and_cols[i]
    for i in range(n_values):
        values_ptr[i] = values[i]

    _offset_to_sorted_blocks(rows_and_cols, n_values, 2, offset_array, n_rows + 1)

    mat.values = values_ptr
    mat.n_values = n_values
    mat.offset_to_row = offset_ptr
    mat.col = col_ptr
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


cdef void _offset_to_sorted_blocks(
    const id_t [:] array,
    long len,
    long stride,
    integral_out_t [:] offset,
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
):
    cdef long n_ids = sorted_ids.shape[0]
    cdef long n_blocks = offset_to_block.shape[0]
    cdef id_t *ptr = &sorted_ids[0, 0]
    cdef id_t [:] sorted_ids_flat = <id_t[:sorted_ids.size]>ptr

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
    const id_t [:, :] src_pairs,
    const id_t [:, :] pairs,
    uint8_t [:] out,
):
    cdef long pair
    cdef long n_pairs = pairs.shape[0]
    cdef long n_values = src_pairs.shape[0]
    cdef long * data_ptr = <long *>malloc(n_values * sizeof(long))
    cdef long [:] data_array = <long [:n_values]>data_ptr
    cdef const id_t *ptr = &src_pairs[0, 0]
    cdef const id_t [:] src_pairs_flat = <const id_t[:src_pairs.size]>ptr
    cdef SparseMatrixInt mat

    for n in range(n_values):
        data_array[n] = 1

    try:
        mat = sparse_matrix_alloc_with_tuple(src_pairs_flat, data_array, n_values, 0)
        for pair in range(n_pairs):
            out[pair] = sparse_matrix_get_or_transpose(
                mat, pairs[pair, 0], pairs[pair, 1]
            )
    finally:
        free(data_ptr)


@cython.boundscheck(False)
@cython.wraparound(False)
def map_pairs_to_values(
    const id_t [:, :] src_pairs,
    const integral_out_t [:] data,
    const id_t [:, :] pairs,
    integral_out_t [:] out,
):
    cdef long pair
    cdef long n_pairs = out.shape[0]
    cdef long n_values = data.shape[0]
    cdef const id_t *ptr = &src_pairs[0, 0]
    cdef const id_t [:] src_pairs_flat = <const id_t[:src_pairs.size]>ptr
    cdef SparseMatrixInt mat

    mat = sparse_matrix_alloc_with_tuple(src_pairs_flat, data, n_values, -1)

    for pair in range(n_pairs):
        out[pair] = sparse_matrix_get_or_transpose(mat, pairs[pair, 0], pairs[pair, 1])


@cython.boundscheck(False)
@cython.wraparound(False)
def map_rolling_pairs_to_values(
    const id_t [:, :] src_pairs,
    const integral_out_t [:] data,
    const id_t [:, :] pairs,
    const id_t [:] size_of_row,
    integral_out_t [:, :] out,
):
    cdef long n_values = data.shape[0]
    cdef long n_pairs = pairs.shape[0]
    cdef long pair
    cdef const id_t *ptr = &src_pairs[0, 0]
    cdef const id_t [:] src_pairs_flat = <const id_t[:2 * len(src_pairs)]>ptr
    cdef SparseMatrixInt mat

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
