import numpy as np
cimport numpy as np
cimport cython

DTYPE_INT = np.int
ctypedef np.int_t DTYPE_INT_t

DTYPE_FLOAT = np.double
ctypedef np.double_t DTYPE_FLOAT_t


def fill_matrix(np.ndarray[DTYPE_INT_t, ndim=1] core2core,
                np.ndarray[DTYPE_INT_t, ndim=1] core2fv,
                np.ndarray[DTYPE_INT_t, ndim=1] fv2core,
                np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_tail,
                np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_head,
                np.ndarray[DTYPE_INT_t, ndim=1] matrow,
                np.ndarray[DTYPE_FLOAT_t, ndim=1] value,
                np.ndarray[DTYPE_FLOAT_t, ndim=2] rhs,
                object mat,
                ):
    """Fill matrix `mat` with values reflecting grid adjacency, and right-side
    vector `rhs` with items from the input `value` array."""
    cdef int ln, t, h

    for ln in core2core:
        t = node_at_link_tail[ln]
        h = node_at_link_head[ln]
        mat[matrow[t], matrow[t]] -= 1.0
        mat[matrow[h], matrow[h]] -= 1.0
        mat[matrow[t], matrow[h]] = 1.0
        mat[matrow[h], matrow[t]] = 1.0

    # Handle core-to-fv links
    for ln in core2fv:
        t = node_at_link_tail[ln]
        h = node_at_link_head[ln]
        mat[matrow[t], matrow[t]] -= 1
        rhs[matrow[t]] -= value[h]

    # Handle fv-to-core links
    for ln in fv2core:
        t = node_at_link_tail[ln]
        h = node_at_link_head[ln]
        mat[matrow[h], matrow[h]] -= 1
        rhs[matrow[h]] -= value[t]

    return mat, rhs


def fill_matrix_with_coefficients(np.ndarray[DTYPE_INT_t, ndim=1] core2core,
                                  np.ndarray[DTYPE_INT_t, ndim=1] core2fv,
                                  np.ndarray[DTYPE_INT_t, ndim=1] fv2core,
                                  np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_tail,
                                  np.ndarray[DTYPE_INT_t, ndim=1] node_at_link_head,
                                  np.ndarray[DTYPE_INT_t, ndim=1] matrow,
                                  np.ndarray[DTYPE_FLOAT_t, ndim=1] value,
                                  np.ndarray[DTYPE_FLOAT_t, ndim=2] rhs,
                                  np.ndarray[DTYPE_FLOAT_t, ndim=1] coef,
                                  object mat,
                                  ):
    """Fill matrix `mat` with values reflecting coefficients related to grid
    adjacency and process coefficient values, and right-side vector `rhs` with
    items from the input `value` array."""
    cdef int ln, t, h

    for ln in core2core:
        t = node_at_link_tail[ln]
        h = node_at_link_head[ln]
        mat[matrow[t], matrow[t]] -= 1.0
        mat[matrow[h], matrow[h]] -= 1.0
        mat[matrow[t], matrow[h]] = 1.0
        mat[matrow[h], matrow[t]] = 1.0

    # Handle core-to-fv links
    for ln in core2fv:
        t = node_at_link_tail[ln]
        h = node_at_link_head[ln]
        mat[matrow[t], matrow[t]] -= 1
        rhs[matrow[t]] -= value[h]

    # Handle fv-to-core links
    for ln in fv2core:
        t = node_at_link_tail[ln]
        h = node_at_link_head[ln]
        mat[matrow[h], matrow[h]] -= 1
        rhs[matrow[h]] -= value[t]

    return mat, rhs
