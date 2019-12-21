import numpy as np
cimport numpy as np
cimport cython


DTYPE_FLOAT = np.double
ctypedef np.double_t DTYPE_FLOAT_t

DTYPE_INT = np.int
#ctypedef np.longlong_t DTYPE_INT_t
ctypedef np.int_t DTYPE_INT_t


@cython.boundscheck(False)
def adjust_flow_receivers(np.ndarray[DTYPE_INT_t, ndim=1] src_nodes,
                          np.ndarray[DTYPE_INT_t, ndim=1] dst_nodes,
                          np.ndarray[DTYPE_FLOAT_t, ndim=1] z,
                          np.ndarray[DTYPE_FLOAT_t, ndim=1] link_slope,
                          np.ndarray[DTYPE_INT_t, ndim=1] active_links,
                          np.ndarray[DTYPE_INT_t, ndim=1] receiver,
                          np.ndarray[DTYPE_INT_t, ndim=1] receiver_link,
                          np.ndarray[DTYPE_FLOAT_t, ndim=1] steepest_slope):
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
    cdef unsigned int n_nodes = src_nodes.size
    cdef int src_id
    cdef int dst_id

    for i in range(n_nodes):
        src_id = src_nodes[i]
        dst_id = dst_nodes[i]

        if z[src_id] > z[dst_id] and link_slope[i] > steepest_slope[src_id]:
            receiver[src_id] = dst_id
            steepest_slope[src_id] = link_slope[i]
            receiver_link[src_id] = active_links[i]
        elif z[dst_id] > z[src_id] and -link_slope[i] > steepest_slope[dst_id]:
            receiver[dst_id] = src_id
            steepest_slope[dst_id] = - link_slope[i]
            receiver_link[dst_id] = active_links[i]
