cimport cython
cimport numpy as np
from cython.parallel cimport parallel
from cython.parallel cimport prange
from libc.stdint cimport int8_t
from libc.stdlib cimport free
from libc.stdlib cimport malloc
import numpy as np

# https://cython.readthedocs.io/en/stable/src/userguide/fusedtypes.html
ctypedef fused float_or_int:
    cython.floating
    cython.integral
    long long
    int8_t

ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
@cython.wraparound(False)
def reorder_rows(
    float_or_int[:, :] value_at_row,
    const id_t[:, :] sorted_cols,
):
    cdef long n_rows = value_at_row.shape[0]
    cdef long n_cols = value_at_row.shape[1]
    cdef long row
    cdef long col
    cdef np.ndarray[float_or_int, ndim=2] out = np.empty_like(value_at_row)

    for row in prange(n_rows, nogil=True, schedule="static"):
        for col in range(n_cols):
            out[row, col] = value_at_row[row, sorted_cols[row, col]]

    return out


@cython.boundscheck(False)
@cython.wraparound(False)
def reorder_rows_inplace(
    float_or_int[:, :] value_at_row,
    const id_t[:, :] sorted_cols,
):
    cdef int n_rows = value_at_row.shape[0]
    cdef int n_cols = value_at_row.shape[1]
    cdef int row
    cdef int col
    cdef float_or_int *scratch

    with nogil, parallel():
        scratch = <float_or_int *>malloc(n_cols * sizeof(float_or_int))

        for row in prange(n_rows, schedule="static"):
            for col in range(n_cols):
                scratch[col] = value_at_row[row, sorted_cols[row, col]]
            for col in range(n_cols):
                value_at_row[row, col] = scratch[col]

        free(scratch)


@cython.boundscheck(False)
@cython.wraparound(False)
def reorder_links_at_node(
    id_t[:, :] links_at_node,
    const id_t[:, :] sorted_links,
):
    cdef int n_nodes = links_at_node.shape[0]
    cdef int n_links_per_node = links_at_node.shape[1]
    cdef int i
    cdef int node
    cdef int *buffer = <int *>malloc(n_links_per_node * sizeof(int))

    try:
        for node in range(n_nodes):
            for i in range(n_links_per_node):
                buffer[i] = links_at_node[node, sorted_links[node, i]]
            for i in range(n_links_per_node):
                links_at_node[node, i] = buffer[i]
    finally:
        free(buffer)


@cython.boundscheck(False)
@cython.wraparound(False)
def reorder_link_dirs_at_node(
    int8_t[:, :] link_dirs_at_node,
    const id_t[:, :] sorted_links,
):
    cdef int n_nodes = link_dirs_at_node.shape[0]
    cdef int n_links_per_node = link_dirs_at_node.shape[1]
    cdef int i
    cdef int node
    cdef int *buffer = <int *>malloc(n_links_per_node * sizeof(int))

    try:
        for node in range(n_nodes):
            for i in range(n_links_per_node):
                buffer[i] = link_dirs_at_node[node, sorted_links[node, i]]
            for i in range(n_links_per_node):
                link_dirs_at_node[node, i] = buffer[i]
    finally:
        free(buffer)
