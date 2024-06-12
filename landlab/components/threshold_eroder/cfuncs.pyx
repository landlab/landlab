import numpy as np

cimport numpy as np

DTYPE_INT = int
ctypedef np.int_t DTYPE_INT_t

DTYPE_FLOAT = np.float64
ctypedef np.float64_t DTYPE_FLOAT_t


cpdef _thresholder(
    np.ndarray[DTYPE_INT_t, ndim=1] stack,
    np.ndarray[DTYPE_INT_t, ndim=1] link_to_rcvr,
    np.ndarray[DTYPE_INT_t, ndim=1] receivers,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] linkLengths,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] el,
    DTYPE_FLOAT_t slope_thres,
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
        el[node] = np.minimum(
            el[node], el[receivers[node]] + slope_thres * dist_to_rcvr
        )
