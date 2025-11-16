cimport cython
from cython.parallel cimport prange
from libc.stdint cimport int32_t
from libc.stdint cimport int64_t

ctypedef fused id_t:
    int32_t
    int64_t


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def calc_grad_at_link(
    const cython.floating [::1] value_at_node,
    *,
    const cython.floating [::1] length_of_link,
    const id_t [:, ::1] nodes_at_link,
    const id_t [::1] where,
    cython.floating [::1] out,
):
    """Calculate the gradient of node-values at specified links.

    The function computes finite differences between the head and tail
    nodes of each link, divided by the link length. Only links listed in
    ``where`` are updated; all other elements in ``out`` are left
    unchanged.

    Parameters
    ----------
    value_at_node : array_like of float
        Node-centered values.
    length_of_link : array_like of float
        Length of each link.
    nodes_at_link : array_like of int, shape (n_links, 2)
        Head and tail node for each link. The first column is the
        tail node and the second column is the head node.
    where : array_like of int
        Indices of links to update.
    out : ndarray of float
        Output buffer of shape ``(n_links,)``. Results are written in
        place to the elements indexed by ``where``.

    Returns
    -------
    out : ndarray of float
        Array of gradients at links corresponding to those indexed
        by ``where``. This is the same object as the input ``out`` array.
    """
    cdef long n_links = len(where)
    cdef long i
    cdef long link
    cdef long head
    cdef long tail

    for i in prange(n_links, nogil=True, schedule="static"):
        link = where[i]
        tail = nodes_at_link[link, 0]
        head = nodes_at_link[link, 1]

        out[link] = (
            value_at_node[head] - value_at_node[tail]
        ) / length_of_link[link]

    return (<object>out).base


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def zero_out_dry_links(
    const cython.floating [::1] h_at_link,
    *,
    const id_t [::1] where,
    cython.floating [::1] out,
):
    cdef Py_ssize_t n_links = where.shape[0]
    cdef Py_ssize_t link
    cdef Py_ssize_t i

    for i in prange(n_links, nogil=True, schedule="static"):
        link = where[i]
        if h_at_link[link] <= 0.0:
            out[link] = 0.0

    return (<object>out).base
