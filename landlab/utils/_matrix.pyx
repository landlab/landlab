import numpy as np
cimport numpy as np
cimport cython

DTYPE_INT = np.int
ctypedef np.int_t DTYPE_INT_t

DTYPE_FLOAT = np.double
ctypedef np.double_t DTYPE_FLOAT_t


@cython.boundscheck(False)
@cython.wraparound(False)
def get_matrix_diagonal_elements_with_coef(
    np.ndarray[DTYPE_INT_t, ndim=2] core_nodes_at_c2c_link,
    np.ndarray[DTYPE_INT_t, ndim=2] core_nodes_at_c2fv_link,
    np.ndarray[DTYPE_INT_t, ndim=2] core_nodes_at_fv2c_link,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] coef_at_c2c_link,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] coef_at_c2fv_link,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] coef_at_fv2c_link,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] data,
):
    cdef int tail, head
    cdef int link
    cdef int n_links
    cdef double coef

    n_links = len(core_nodes_at_c2c_link)
    for link in range(n_links):
        tail = core_nodes_at_c2c_link[link, 0]
        head = core_nodes_at_c2c_link[link, 1]
        coef = coef_at_c2c_link[link]
        data[tail] -= coef
        data[head] -= coef

    n_links = len(core_nodes_at_c2fv_link)
    for link in range(n_links):
        tail = core_nodes_at_c2fv_link[link, 0]
        data[tail] -= coef_at_c2fv_link[link]

    n_links = len(core_nodes_at_fv2c_link)
    for link in range(n_links):
        head = core_nodes_at_fv2c_link[link, 1]
        data[head] -= coef_at_fv2c_link[link]


@cython.boundscheck(False)
@cython.wraparound(False)
def get_matrix_diagonal_elements(
    np.ndarray[DTYPE_INT_t, ndim=2] core_nodes_at_c2c_link,
    np.ndarray[DTYPE_INT_t, ndim=2] core_nodes_at_c2fv_link,
    np.ndarray[DTYPE_INT_t, ndim=2] core_nodes_at_fv2c_link,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] data,
):
    cdef int tail, head
    cdef int link
    cdef int n_links

    n_links = len(core_nodes_at_c2c_link)
    for link in range(n_links):
        tail = core_nodes_at_c2c_link[link, 0]
        head = core_nodes_at_c2c_link[link, 1]
        data[tail] -= 1.0
        data[head] -= 1.0

    n_links = len(core_nodes_at_c2fv_link)
    for link in range(n_links):
        tail = core_nodes_at_c2fv_link[link, 0]
        data[tail] -= 1.0

    n_links = len(core_nodes_at_fv2c_link)
    for link in range(n_links):
        head = core_nodes_at_fv2c_link[link, 1]
        data[head] -= 1.0


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_right_hand_side(
    np.ndarray[DTYPE_INT_t, ndim=2] nodes_at_c2fv_link,
    np.ndarray[DTYPE_INT_t, ndim=2] nodes_at_fv2c_link,
    np.ndarray[DTYPE_INT_t, ndim=1] core_node_at_node,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] value_at_node,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] out,
):
    cdef int tail, head

    for tail, head in nodes_at_c2fv_link:
        out[core_node_at_node[tail]] -= value_at_node[head]

    for tail, head in nodes_at_fv2c_link:
        out[core_node_at_node[head]] -= value_at_node[tail]
