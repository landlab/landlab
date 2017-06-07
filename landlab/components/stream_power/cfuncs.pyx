import numpy as np
cimport numpy as np
cimport cython

#from libc.math cimport fabs

# as best as KRB can tell, the function _brentq is in c.
# not sure what, if anything I need to do to 
from scipy.optimize._zeros import _brentq as brentq

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
                
def brent_method_erode_variable_threshold(np.ndarray[DTYPE_INT_t, ndim=1] src_nodes,
                                          np.ndarray[DTYPE_INT_t, ndim=1] dst_nodes,
                                          np.ndarray[DTYPE_FLOAT_t, ndim=1] threshsxdt,
                                          np.ndarray[DTYPE_FLOAT_t, ndim=1] alpha,
                                          DTYPE_FLOAT_t n,
                                          np.ndarray[DTYPE_FLOAT_t, ndim=1] z):
    """Erode node elevations using Brent's method for stability.

    The alpha value is given as 
    
    alpha = delta_t*K * (rainfall_intensity*A)**m/(delta_x**n)
    
    It will be multiplied by the value:
        (z_node(t) - z_downstream(t+delta_t))**(n-1)
        
    to become the alpha value defined below in the function used in erode_fun
    used in the root-finding operation. 
    
    Parameters
    ----------
    src_nodes : array_like
        Ordered upstream node ids.
    dst_nodes : array_like
        Node ids of nodes receiving flow.
    threshsxdt : array_like
        Incision thresholds at nodes multiplied by the timestep.
    alpha : array_like
        Erosion factor. 
    n : float
        Exponent.
    z : array_like
        Node elevations.
    """     

    cdef unsigned int n_nodes = src_nodes.size
    cdef unsigned int src_id
    cdef unsigned int dst_id
    cdef unsigned int i
    cdef double z_old
    cdef double z_downstream
    cdef double thresholddt
    cdef double z_diff_old
    cdef double z_diff
    cdef double alpha_param
    cdef double beta_param
    cdef double check_function
    cdef double x
                
    # Solve using Brent's method
    for i in range(n_nodes):
        
        # get IDs for source and reciever nodes
        src_id = src_nodes[i]
        dst_id = dst_nodes[src_id]
        
        # if a node does not flow to itself, continue
        if src_id != dst_id:
             
            # Get values for z at present node and present time, 
            # and z downstream at t + delta t (which should have just been
            # solved for)
            z_old = z[src_id]
            z_downstream = z[dst_id]
            
            # Get the threshold value. 
            thresholddt = threshsxdt[src_id]
       
            # calculate the difference between z_old and z_downstream
            z_diff_old = z_old - z_downstream
            
            # if flow was reversed, z_diff_old may be negative, which can
            # lead to Runtime or Floating PointErrors. This is dealt with 
            #in part by setting alpha to zero. But to prevent these errors, 
            # here we also set z_diff_old, alpha, and beta to zero. 
            # this will result in a value of 
            
            if z_diff_old >=0:

                # using z_diff_old, calculate the alpha paramter of Braun and
                # Willet by calculating alpha times z
            
                alpha_param = alpha[src_id] * pow(z_diff_old, n-1.0)
                                
                beta_param = thresholddt / z_diff_old
                    
                # check if the threshold has been exceeded:
                check_function = erode_fn(1, alpha_param, beta_param, n)
                
                if check_function <= 0: 
                    # if the threshold was not exceeded 
                    # do not change the elevation
                    # this means that the maximum possible slope value 
                    # does not produce stream power needed to exceed the erosion
                    # threshold
                    pass
                else:
                    # if the threshold was exceeded, then there will be a zero
                    # between x = 0 and x= 1
                    
                    # solve using brentq, which requires a zero to exist 
                    # in between the two end values
                    
                    # if n is 1, solution has an analytical solution. 
                    if n != 1.0:
                        
                        x = brentq(erode_fn, 
                                           0.0, 1.0, 
                                           1e-12,4.4408920985006262e-16, 
                                           100, 
                                           (alpha_param, beta_param, n),
                                           False,
                                           True)
                        
#                        x = brentq(erode_fn, 
#                                   0.0,                                  
#                                   1.0,
#                                   args=(alpha_param, beta_param, n),
#                                   maxiter=200)
                    else:
                        x = (1.0 + beta_param)/(1.0 + alpha_param)
                    # just in case, 
                    if x>0:
                        z[src_id] = z_downstream + x * (z_old - z_downstream)
                    else:
                        z[src_id] = z_downstream + 1.0e-15
                        
def brent_method_erode_fixed_threshold(np.ndarray[DTYPE_INT_t, ndim=1] src_nodes,
                                       np.ndarray[DTYPE_INT_t, ndim=1] dst_nodes,
                                       DTYPE_FLOAT_t threshsxdt,
                                       np.ndarray[DTYPE_FLOAT_t, ndim=1] alpha,
                                       DTYPE_FLOAT_t n,
                                       np.ndarray[DTYPE_FLOAT_t, ndim=1] z):
    
    """Erode node elevations.
    
    The alpha value is given as 
    
    alpha = delta_t*K * (rainfall_intensity*A)**m/(delta_x**n)
    
    It will be multiplied by the value:
        (z_node(t) - z_downstream(t+delta_t))**(n-1)
        
    to become the alpha value defined below in the function used in erode_fun
    used in the root-finding operation. 
    
    Parameters
    ----------
    src_nodes : array_like
        Ordered upstream node ids.
    dst_nodes : array_like
        Node ids of nodes receiving flow.
    threshsxdt : float
        Incision thresholds at nodes multiplied by the timestep.
    alpha : array_like
        Erosion factor.
    n : float
        Exponent.
    z : array_like
        Node elevations.
    """     

    cdef unsigned int n_nodes = src_nodes.size
    cdef unsigned int src_id
    cdef unsigned int dst_id
    cdef unsigned int i
    cdef double z_old
    cdef double z_downstream
    cdef double z_diff_old
    cdef double z_diff
    cdef double alpha_param
    cdef double beta_param
    cdef double check_function
    cdef double x
                
    # Solve using Brent's method
    for i in range(n_nodes):
        
        # get IDs for source and reciever nodes
        src_id = src_nodes[i]
        dst_id = dst_nodes[src_id]
        
        # if a node does not flow to itself, continue
        if src_id != dst_id:
             
            # Get values for z at present node and present time, 
            # and z downstream at t + delta t (which should have just been
            # solved for)
            z_old = z[src_id]
            z_downstream = z[dst_id]
            
            # Get the threshold value. 
            # its constant so we already have it. 
       
            # calculate the difference between z_old and z_downstream
            z_diff_old = z_old - z_downstream
            
            # if flow was reversed, z_diff_old may be negative, which can
            # lead to Runtime or Floating PointErrors. This is dealt with 
            #in part by setting alpha to zero. But to prevent these errors, 
            # here we also set z_diff_old, alpha, and beta to zero. 
            # this will result in a value of 
            
            if z_diff_old >=0:

                # using z_diff_old, calculate the alpha paramter of Braun and
                # Willet by calculating alpha times z
            
                alpha_param = alpha[src_id] * pow(z_diff_old, n-1.0)
                                
                beta_param = threshsxdt / z_diff_old
                    
                # check if the threshold has been exceeded:
                check_function = erode_fn(1, alpha_param, beta_param, n)
                
                if check_function <= 0: 
                    # if the threshold was not exceeded 
                    # do not change the elevation
                    # this means that the maximum possible slope value 
                    # does not produce stream power needed to exceed the erosion
                    # threshold
                    pass
                else:
                    # if the threshold was exceeded, then there will be a zero
                    # between x = 0 and x= 1
                    
                    # solve using brentq, which requires a zero to exist 
                    # in between the two end values
                    
                    # if n is 1, solution has an analytical solution. 
                    if n != 1.0:
                        
                        x = brentq(erode_fn, 
                                           0.0, 1.0, 
                                           1e-12,4.4408920985006262e-16, 
                                           100, 
                                           (alpha_param, beta_param, n),
                                           False,
                                           True)
                        
#                        x = brentq(erode_fn, 
#                                   0.0,                                  
#                                   1.0,
#                                   args=(alpha_param, beta_param, n),
#                                   maxiter=200)
                    else:
                        x = (1.0 + beta_param)/(1.0 + alpha_param)
                    # just in case, 
                    if x>0:
                        z[src_id] = z_downstream + x * (z_old - z_downstream)
                    else:
                        z[src_id] = z_downstream + 1.0e-15
                        
                        
def erode_fn(DTYPE_FLOAT_t x, 
             DTYPE_FLOAT_t alpha, 
             DTYPE_FLOAT_t beta, 
             DTYPE_FLOAT_t n):
    """Evaluates the solution to the water-depth equation.

    Called by scipy.brentq() to find solution for $x$ using Brent's method.

    Parameters
    ----------
    x : float
        normalized elevation, see below.
    alpha : float
        alpha parameter, see below.
    beta : float
        beta parameter, see below.
    n : float
        n exponent
    check : boolean, default is True
        flag to determine if a ValueError should be thrown if a check 
        identifies that the threshold value is high enough such that no erosion
        will occur.
        
    
    This equation represents the implicit solution for normalized topographic 
    elevation $x$ at the  next time step. This solution is inspired by the 
    Appendix of Braun and Willet (2012) but was generalized to include an a 
    threshold value such that if the threshold is not exceeded, no erosion will
    occur. 
    
    Consider stream power erosion under the equation:

        E = K * (rainfall_intensity*A)**m * S**n - threshold_sp,
    
    on a grid with link delta_x and for a timestep of delta_t. 
    
    When iterating from downstream to upstream in the drainage stack, at a 
    given node at time = t+delta_t, the value of the node at time = t, and the
    value of the downstream node at time t+delta_t is known. 
    
    Define 
    x = (z_node(t+delta_t) - z_downstream(t+delta_t))/(z_node(t) - z_downstream(t+delta_t)).
    
    A discretized version of the stream power equation above yeilds the equation
    
    f = x - 1 + alpha*(x**n) - beta
    
    where 
    
    alpha = delta_t*K * (rainfall_intensity*A)**m/(delta_x**n) * (z_node(t) - z_downstream(t+delta_t))**(n-1)
    
    and 
    
    beta = threshold_sp*delta_t/(z_node(t) - z_downstream(t+delta_t))
    
    Finding the root of f provides the implicit solution for the stream power
    equation. 
    
    If f(x=1) = 0, then no erosion occurs as potential erosion is cancelled by
    the erosion threshold and the topography at the given node does not change. 
    
    If f(x=0) = 0, then the topography at the given node becomes that of the 
    the downstream mode. 
    
    If the threshold term, beta is zero, this equation collapses to the form 
    given by Braun and Willet (2012).
    
    When the threshold term is greater than zero, it is possible that the no
    erosion will occur. In this case, an evaluation of f(x=1) will yeild a 
    negative number. 
    
    """
    cdef double f
        
    f = x - 1.0 + (alpha * (x ** n)) - beta

    return f