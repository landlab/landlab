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
def calc_bates_flow_height(
    const cython.floating [::1] z_at_node,
    const cython.floating [::1] h_at_node,
    *,
    const id_t [:, ::1] nodes_at_link,
    const id_t [::1] where,
    cython.floating [::1] out,
):
    """Calculate flow depth at links using the Bates formulation.

    The function computes, for each link in ``where``, the maximum of the
    water-surface elevation (bed elevation plus water depth) between the two
    nodes connected by that link, minus the maximum bed elevation at those
    same nodes.

        h_at_link = max(h_tail + z_tail, h_head + z_head) - max(z_tail, z_head)

    Parameters
    ----------
    z_at_node : ndarray of float, shape (n_nodes,)
        Bed elevation at each node.
    h_at_node : ndarray of float, shape (n_nodes,)
        Water depth at each node.
    nodes_at_link : ndarray of int, shape (n_links, 2)
        IDs of the tail and head nodes for each link.
    where : ndarray of int, shape (n_selected,)
        Indices of links to update.
    out : ndarray of float, shape (n_links,)
        Output array to store the computed flow heights. Only entries indexed
        by ``where`` are modified; all other values remain unchanged.

    Returns
    -------
    out : ndarray of float, shape (n_links,)
        The same array passed as ``out``, with updated values at
        selected links.
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

        out[link] = max(
            h_at_node[tail] + z_at_node[tail],
            h_at_node[head] + z_at_node[head],
        ) - max(z_at_node[tail], z_at_node[head])

    return (<object>out).base
