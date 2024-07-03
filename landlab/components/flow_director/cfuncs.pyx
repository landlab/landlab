cimport cython


@cython.boundscheck(False)
def adjust_flow_receivers(
    const cython.integral [:] src_nodes,
    const cython.integral [:] dst_nodes,
    const cython.floating [:] z,
    const cython.floating [:] link_slope,
    const cython.integral [:] active_links,
    cython.integral [:] receiver,
    cython.integral [:] receiver_link,
    cython.floating [:] steepest_slope,
):
    """Adjust flow receivers based on link slopes and steepest gradients.

    Parameters
    ----------
    src_nodes : array_like
        Ordered upstream node ids.
    dst_nodes : array_like
        Node ids of nodes receiving flow.
    z : array_like
        Node elevations.
    link_slope : array_like
        Link gradients.
    active_links : array_like
        Link IDs for active links.
    receiver : array_like
        Flow-receiver link IDs.
    steepest_slope : array_like
        Gradient of steepest descent from nodes.
    """
    cdef unsigned int n_nodes = src_nodes.shape[0]
    cdef int src_id
    cdef int dst_id

    for i in range(n_nodes):
        src_id = src_nodes[i]
        dst_id = dst_nodes[i]

        if (z[src_id] > z[dst_id]) and (link_slope[i] > steepest_slope[src_id]):
            receiver[src_id] = dst_id
            steepest_slope[src_id] = link_slope[i]
            receiver_link[src_id] = active_links[i]
        elif (z[dst_id] > z[src_id]) and (-link_slope[i] > steepest_slope[dst_id]):
            receiver[dst_id] = src_id
            steepest_slope[dst_id] = - link_slope[i]
            receiver_link[dst_id] = active_links[i]
