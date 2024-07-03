cimport cython


cpdef _thresholder(
    const cython.integral [:] stack,
    const cython.integral [:] link_to_rcvr,
    const cython.integral [:] receivers,
    const cython.floating [:] linkLengths,
    cython.floating [:] el,
    const double slope_thres,
):
    """
    Calcualte D8 flow dirs
    stack: the flow upstream node order
    link_to_rcvr: Link to receiver
    receivers: receivers
    linkLengths: length of links
    el: topographic elevation
    """
    cdef int node

    for node in stack:
        dist_to_rcvr = linkLengths[link_to_rcvr[node]]
        el[node] = min(
            el[node], el[receivers[node]] + slope_thres * dist_to_rcvr
        )
