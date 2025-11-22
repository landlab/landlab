cimport cython
from cython.parallel cimport prange
from libc.stdint cimport int32_t
from libc.stdint cimport int64_t
from libc.stdint cimport uint8_t

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
def calc_weighted_mean_of_parallel_links(
    const cython.floating[::1] value_at_link,
    const id_t[:, ::1] parallel_links_at_link,
    *,
    const uint8_t[::1] status_at_link,
    const id_t[::1] where,
    const double theta,
    cython.floating[::1] out,
):
    """Calculate a weighted mean of each link and its parallel neighbors.

    For each link in ``where``, compute a weighted combination of that link's
    value and the values of its left and right parallel links. The combination
    weight is controlled by ``theta`` such that::

        out[link] = theta * value_at_link[link] + (1 - theta) / 2 * (
            value_at_link[left] + value_at_link[right]
        )

    Links with invalid or inactive neighbors retain their original value.

    Parameters
    ----------
    value_at_link : array_like of float
        Values defined at links. Must have length equal to the number of links.
    parallel_links_at_link : ndarray of int, shape (n_links, 2)
        IDs of the two parallel links (left, right) associated with each link.
        A value of ``-1`` indicates that no parallel link exists.
    status_at_link : ndarray of uint8
        Link status codes. A value of ``4`` indicates an inactive link that
        should not contribute to the weighted mean.
    where : ndarray of int
        Indices of links for which to update ``out``. Each link in ``where``
        is processed exactly once. Must not contain duplicates.
    theta : float
        Weighting coefficient between 0 and 1. A value of 1.0 means the link's
        own value is kept unchanged; a value of 0.0 means a simple mean of the
        two parallel links.
    out : ndarray of float
        Output buffer into which results are written. Must be the same shape
        and dtype as ``value_at_link``.

    Returns
    -------
    out : ndarray of float
        Array of means at links corresponding to those indexed
        by ``where``. This is the same object as the input ``out`` array.
    """
    cdef long n_links = len(where)
    cdef long link
    cdef long left
    cdef long right
    cdef long i
    cdef long n_neighbors
    cdef double total
    cdef double ONE_MINUS_THETA = (1.0 - theta)
    cdef int LINK_IS_INACTIVE = 4
    cdef int LINK_IS_MISSING = -1

    for i in prange(n_links, nogil=True, schedule="static"):
        link = where[i]
        left = parallel_links_at_link[link, 0]
        right = parallel_links_at_link[link, 1]

        total = 0.0
        n_neighbors = 0

        if left != LINK_IS_MISSING and status_at_link[left] != LINK_IS_INACTIVE:
            n_neighbors = n_neighbors + 1
            total = total + value_at_link[left]

        if right != LINK_IS_MISSING and status_at_link[right] != LINK_IS_INACTIVE:
            n_neighbors = n_neighbors + 1
            total = total + value_at_link[right]

        if n_neighbors == 2:
            out[link] = theta * value_at_link[link] + ONE_MINUS_THETA * total * 0.5
        elif n_neighbors == 1:
            out[link] = theta * value_at_link[link] + ONE_MINUS_THETA * total
        else:
            out[link] = value_at_link[link]

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
    """Zero values for links with non-positive depth.

    Sets ``out[link] = 0.0`` for each link where the corresponding
    entry in ``h_at_link`` is less than or equal to zero, using the
    indices provided in ``where``. Links with positive depth are
    left unchanged.

    Parameters
    ----------
    h_at_link : array_like of float, shape (n_links,)
        Water depth at each link.
    where : array_like of int
        Indices of links to operate on. Usually a subset of links,
        such as those active in the current routing step.
    out : ndarray of float, shape (n_links,)
        Output array into which results are written. Must be
        pre-initialized, entries at dry links are overwritten;
        all other entries are left unchanged.

    Returns
    -------
    out : ndarray of float
        The updated ``out`` array (same object as input), with values
        set to zero at dry links.
    """
    cdef Py_ssize_t n_links = where.shape[0]
    cdef Py_ssize_t link
    cdef Py_ssize_t i

    for i in prange(n_links, nogil=True, schedule="static"):
        link = where[i]
        if h_at_link[link] <= 0.0:
            out[link] = 0.0

    return (<object>out).base
