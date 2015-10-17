import numpy as np
cimport numpy as np
cimport cython

#from libc.math cimport fabs


DTYPE_FLOAT = np.double
ctypedef np.double_t DTYPE_FLOAT_t

DTYPE_INT = np.int
ctypedef np.int_t DTYPE_INT_t


cdef extern from "math.h":
    double fabs(double x) nogil
    double pow(double x, double y) nogil


@cython.boundscheck(False)
def erode_avoiding_pits(np.ndarray[DTYPE_INT_t, ndim=1] src_nodes,
                        np.ndarray[DTYPE_INT_t, ndim=1] dst_nodes,
                        np.ndarray[DTYPE_FLOAT_t, ndim=1] node_z,
                        np.ndarray[DTYPE_FLOAT_t, ndim=1] node_dz):
    """Erode node elevations while avoiding creating pits.

    Parameters
    ----------
    src_nodes : array_like
        Ordered upstream node ids.
    dst_nodes : array_like
        Node ids of nodes receiving flow.
    node_z : array_like
        Node elevations.
    node_dz : array_like
        Node erosion rates.
    """
    cdef unsigned int n_nodes = src_nodes.size
    cdef double z_src_after
    cdef double z_dst_after
    cdef unsigned int src_id
    cdef unsigned int dst_id
    cdef unsigned int i

    for i in range(n_nodes):
        src_id = src_nodes[i]
        dst_id = dst_nodes[src_id]

        z_src_after = node_z[src_id] - node_dz[src_id]
        z_dst_after = node_z[dst_id] - node_dz[dst_id]

        if z_src_after < z_dst_after:
            node_dz[src_id] = (node_z[src_id] - z_dst_after) * 0.999999


@cython.boundscheck(False)
def erode_with_alpha(np.ndarray[DTYPE_INT_t, ndim=1] src_nodes,
                     np.ndarray[DTYPE_INT_t, ndim=1] dst_nodes,
                     np.ndarray[DTYPE_FLOAT_t, ndim=1] alpha,
                     np.ndarray[DTYPE_FLOAT_t, ndim=1] z):
    """Erode node elevations based on a scaling factor.

    Parameters
    ----------
    src_nodes : array_like
        Ordered upstream node ids.
    dst_nodes : array_like
        Node ids of nodes receiving flow.
    alpha : array_like
        Erosion factor.
    z : array_like
        Node elevations.
    """
    cdef unsigned int n_nodes = src_nodes.size
    cdef unsigned int src_id
    cdef unsigned int dst_id
    cdef unsigned int i

    for i in range(n_nodes):
        src_id = src_nodes[i]
        dst_id = dst_nodes[src_id]
        if src_id != dst_id:
            z[src_id] = ((z[src_id] + alpha[src_id] * z[dst_id]) /
                         (1.0 + alpha[src_id]))


def erode_with_link_alpha(np.ndarray[DTYPE_INT_t, ndim=1] src_nodes,
                          np.ndarray[DTYPE_INT_t, ndim=1] dst_nodes,
                          np.ndarray[DTYPE_FLOAT_t, ndim=1] alpha,
                          DTYPE_FLOAT_t n,
                          np.ndarray[DTYPE_FLOAT_t, ndim=1] z):
    """Erode node elevations using alpha scaled by link length.

    Parameters
    ----------
    src_nodes : array_like
        Ordered upstream node ids.
    dst_nodes : array_like
        Node ids of nodes receiving flow.
    alpha : array_like
        Erosion factor scaled by link length to the *n - 1*.
    n : float
        Exponent.
    z : array_like
        Node elevations.
    """
    cdef unsigned int n_nodes = src_nodes.size
    cdef unsigned int src_id
    cdef unsigned int dst_id
    cdef unsigned int i
    cdef double z_diff
    cdef double prev_z
    cdef double next_z
    cdef double f

    for i in range(n_nodes):
        src_id = src_nodes[i]
        dst_id = dst_nodes[src_id]

        if src_id != dst_id:
            next_z = z[src_id]

            while True:
                prev_z = next_z;

                z_diff = prev_z - z[dst_id]
                f = alpha[src_id] * pow(z_diff, n - 1.)
                next_z = prev_z - ((prev_z - z[src_id] + f * z_diff) /
                                   (1. + n * f))

                if fabs((next_z - prev_z) / next_z) < 1.48e-08:
                    break

            z[src_id] = next_z
