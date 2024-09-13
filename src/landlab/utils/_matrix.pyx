cimport cython

# https://cython.readthedocs.io/en/stable/src/userguide/fusedtypes.html
ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
@cython.wraparound(False)
def get_matrix_diagonal_elements_with_coef(
    const id_t [:, :] core_nodes_at_c2c_link,
    const id_t [:, :] core_nodes_at_c2fv_link,
    const id_t [:, :] core_nodes_at_fv2c_link,
    const cython.floating [:] coef_at_c2c_link,
    const cython.floating [:] coef_at_c2fv_link,
    const cython.floating [:] coef_at_fv2c_link,
    cython.floating [:] data,
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
    const id_t [:, :] core_nodes_at_c2c_link,
    const id_t [:, :] core_nodes_at_c2fv_link,
    const id_t [:, :] core_nodes_at_fv2c_link,
    cython.floating [:] data,
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
    const id_t [:, :] nodes_at_c2fv_link,
    const id_t [:, :] nodes_at_fv2c_link,
    const id_t [:] core_node_at_node,
    const cython.floating [:] value_at_node,
    cython.floating [:] out,
):
    cdef int tail, head

    for tail, head in nodes_at_c2fv_link:
        out[core_node_at_node[tail]] -= value_at_node[head]

    for tail, head in nodes_at_fv2c_link:
        out[core_node_at_node[head]] -= value_at_node[tail]
