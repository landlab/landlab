import numpy as np
cimport numpy as np
cimport cython
from scipy.optimize import newton
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


def erode_with_link_alpha_varthresh(np.ndarray[DTYPE_INT_t, ndim=1] src_nodes,
                                    np.ndarray[DTYPE_INT_t, ndim=1] dst_nodes,
                                    np.ndarray[DTYPE_FLOAT_t, ndim=1] threshsxdt,
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
    threshsxdt : array_like
        Incision thresholds at nodes multiplied by the timestep.
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
    cdef double threshxdt
    cdef double z_diff
    cdef double prev_z
    cdef double next_z
    cdef double f
    cdef double excess_thresh

    for i in range(n_nodes):
        src_id = src_nodes[i]
        dst_id = dst_nodes[src_id]
        threshxdt = threshsxdt[i]

        if src_id != dst_id and z[src_id] > z[dst_id]:
            next_z = z[src_id]
            prev_z = 0.
            niter = 0

            while True:
                niter += 1
                if niter > 40:
                    print('===============')
                    print(z[src_id])
                    print(z[dst_id])
                    print(next_z)
                    print(z_diff)
                    print(f)
                    print(alpha[src_id])
                assert niter < 50, 'failure to converge in SP solver'
                z_diff = next_z - z[dst_id]
                f = alpha[src_id] * pow(z_diff, n - 1.)
                excess_thresh = f * z_diff - threshxdt
                if excess_thresh < 0.:
                    excess_thresh = 0.
                next_z = next_z - ((next_z - z[src_id] + excess_thresh) /
                                   (1. + n * f))
                if next_z < z[dst_id]:
                    next_z = z[dst_id] + 1.e-15  # maintain connectivity
                if next_z != 0.:
                    if fabs((next_z - prev_z)/next_z) < 1.48e-08 or n == 1.:
                        break
                else:
                    break

                prev_z = next_z;

            if next_z < z[src_id]:
                z[src_id] = next_z


def erode_with_link_alpha_fixthresh(np.ndarray[DTYPE_INT_t, ndim=1] src_nodes,
                                    np.ndarray[DTYPE_INT_t, ndim=1] dst_nodes,
                                    DTYPE_FLOAT_t threshxdt,
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
    threshxdt : float
        The (spatially uniform) incision threshold multiplied by the timestep.
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
    cdef unsigned int niter
    cdef double z_diff
    cdef double prev_z
    cdef double next_z
    cdef double f
    cdef double excess_thresh

    for i in range(n_nodes):
        src_id = src_nodes[i]
        dst_id = dst_nodes[src_id]

        if src_id != dst_id and z[src_id] > z[dst_id]:
            next_z = z[src_id]
            prev_z = 0.
            niter = 0

            while True:
                niter += 1
                if niter > 40:
                    print('===============')
                    print(z[src_id])
                    print(z[dst_id])
                    print(next_z)
                    print(z_diff)
                    print(f)
                    print(alpha[src_id])
                assert niter < 50, 'failure to converge in SP solver'
                z_diff = next_z - z[dst_id]
                f = alpha[src_id] * pow(z_diff, n - 1.)
                excess_thresh = f * z_diff - threshxdt
                if excess_thresh < 0.:
                    excess_thresh = 0.
                next_z = next_z - ((next_z - z[src_id] + excess_thresh) /
                                   (1. + n * f))
                if next_z < z[dst_id]:
                   next_z = z[dst_id] + 1.e-15  # maintain connectivity
                if next_z != 0.:
                    if fabs((next_z - prev_z)/next_z) < 1.48e-08 or n == 1.:
                        break
                else:
                    break

                prev_z = next_z;

            if next_z < z[src_id]:
                z[src_id] = next_z
                
def smooth_stream_power_eroder_solver(np.ndarray[DTYPE_INT_t, ndim=1] src_nodes,
                                      np.ndarray[DTYPE_INT_t, ndim=1] dst_nodes,
                                      np.ndarray[DTYPE_FLOAT_t, ndim=1] z,
                                      np.ndarray[DTYPE_FLOAT_t, ndim=1] alpha,
                                      np.ndarray[DTYPE_FLOAT_t, ndim=1] gamma,
                                      np.ndarray[DTYPE_FLOAT_t, ndim=1] delta):
    """
    Erode node elevations for SmoothStreamPower eroder.  

    Parameters
    ----------
    src_nodes : array_like
        Ordered upstream node ids.
    dst_nodes : array_like
        Node ids of nodes receiving flow.
    z : array_like
        Node elevations.
    alpha : array_like
        Parameter = K A^m dt / L (nondimensional) at all nodes
    gamma : array_like
        Parameter = omega_c * dt (dimension of L, because omega_c [=] L/T)
    delta : array_like
        Parameter = K A^m / (L * wc) [=] L^{-1} (so d * z [=] [-])
    """
    # define internal variables
    cdef unsigned int n_nodes = src_nodes.size
    cdef unsigned int src_id
    cdef unsigned int dst_id
    cdef unsigned int i
    cdef double epsilon
    
    # loop through nodes from bottom of the stream network. 
    for i in range(len(src_nodes)):
        
        # get the node ids of the source and reciever. 
        src_id = src_nodes[i]
        dst_id = dst_nodes[src_id]
                    
        # if the node doesn't drain to iteself and the source is above the 
        # reciever, continue. 
        if src_id != dst_id and z[src_id] > z[dst_id]:
            
            # calculate epsilon
            epsilon = (alpha[src_id] * z[dst_id]
                       + gamma[src_id] + z[src_id])
            
            # calculate new z using newton's method. 
            z[src_id] = newton(smooth_stream_power_erosion_equation, 
                               z[src_id],
                               fprime=smooth_stream_power_erosion_prime,
                               args=(alpha[src_id],
                                     z[dst_id],
                                     gamma[src_id],
                                     delta[src_id],
                                     epsilon)) 
        # nothing needs to be returned as this method updates z. 
            
def smooth_stream_power_erosion_equation(DTYPE_FLOAT_t x, 
                                         DTYPE_FLOAT_t a, 
                                         DTYPE_FLOAT_t b, 
                                         DTYPE_FLOAT_t c, 
                                         DTYPE_FLOAT_t d, 
                                         DTYPE_FLOAT_t e):
    """Equation for elevation of a node at timestep t+1 for SmoothStreamPower.

    Parameters
    ----------
    x : float
        Value of new elevation
    a : float
        Parameter = K A^m dt / L (nondimensional)
    b : float
        Elevation of downstream node, z_j
    c : float
        Parameter = omega_c * dt (dimension of L, because omega_c [=] L/T)
    d : float
        Parameter = K A^m / (L * wc) [=] L^{-1} (so d * z [=] [-])
    e : float
        z(t) + a z_j + (wc * dt)
    """
    cdef double f
    
    f = x * (1.0 + a) + c * np.exp(-d * (x - b)) - e
    
    return f


def smooth_stream_power_erosion_prime(DTYPE_FLOAT_t x, 
                                      DTYPE_FLOAT_t a, 
                                      DTYPE_FLOAT_t b, 
                                      DTYPE_FLOAT_t c, 
                                      DTYPE_FLOAT_t d,
                                      DTYPE_FLOAT_t e):
    
    """Derivative of the equation for elevation of a node for SmoothStreamPower.
    
    Parameters
    ----------
    x : float
        Value of new elevation
    a : float
        Parameter = K A^m dt / L
    b : float
        Elevation of downstream node, z_j
    c : float
        Parameter = omega_c * dt (dimension of L, because omega_c [=] L/T)
    d : float
        Parameter = K A^m / (L * wc) [=] L^{-1} (so d * z [=] [-])
    e : n/a
        Placeholder; not used but must be input for solver.
    """
    cdef double f

    f = (1.0 + a) - c * d * np.exp(-d * (x - b))
    
    return f
