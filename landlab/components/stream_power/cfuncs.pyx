import numpy as np
cimport numpy as np
cimport cython
from scipy.optimize import newton
#from libc.math cimport fabs

# suspect that the function _brentq is in c and thus this is the most effective
# method for using the brentq method in cython.
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
    # define internally used variables.
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

    # Loop through nodes.
    for i in range(n_nodes):

        # get IDs for source and reciever nodes
        src_id = src_nodes[i]
        dst_id = dst_nodes[src_id]

        # if a node does not flow to itself, and the source node is above the
        # destination node
        if src_id != dst_id and z[src_id] > z[dst_id]:

            # Get values for z at present node and present time,
            # and z downstream at t + delta t (which should have been
            # previously solved for)
            z_old = z[src_id]
            z_downstream = z[dst_id]

            # Get the threshold value. In this function, it is spatially variable
            thresholddt = threshsxdt[src_id]

            # calculate the difference between z_old and z_downstream
            z_diff_old = z_old - z_downstream

            # using z_diff_old, calculate the alpha paramter of Braun and
            # Willet by calculating alpha times z

            alpha_param = alpha[src_id] * pow(z_diff_old, n-1.0)

            # Calculate the beta parameter that accounts for the possible
            # presence of a threshold.
            beta_param = thresholddt / z_diff_old

            # check if the threshold has been exceeded by passing a value of
            # x = 1 to the erode_fn. If this returns a value of less than
            # zero, this means that the the maximum possible slope value  does
            # not produce stream power needed to exceed the erosion threshold
            check_function = erode_fn(1, alpha_param, beta_param, n)

            # if the threshold was not exceeded do not change the elevation,
            # otherwise calculate the erosion rate
            if check_function > 0:
                # if the threshold was exceeded, then there will be a zero
                # between x = 0 and x= 1

                # solve using brentq, which requires a zero to exist
                # in between the two end values

                # if n is 1, finding x has an analytical solution. Otherwise,
                # use the the numerical solution given by root finding
                if n != 1.0:

                    # The threshold values passed here are the defaults if one
                    # were to import brentq from scipy.optimize
                    x = brentq(erode_fn,
                                       0.0, 1.0,
                                       1e-12,4.4408920985006262e-16,
                                       100,
                                       (alpha_param, beta_param, n),
                                       False,
                                       True)

                else:
                    # Analytical solution
                    x = (1.0 + beta_param)/(1.0 + alpha_param)

                # If x is provided as a value greater than zero, calculate
                # z at t=t+delta_t useing the values of x, z_downstream and
                # z_old as given by the definition of x (see erode_fn for
                # details). If x is equal to zero, set it as just slightly
                # higher than x_downstream.
                if x>0:
                    z[src_id] = z_downstream + x * (z_old - z_downstream)
                else:
                    z[src_id] = z_downstream + 1.0e-15

                # Nothing is returned from this function as it serves to update
                # the array z.


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
    # define internally used variables.
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

    # Loop through nodes.
    for i in range(n_nodes):

        # get IDs for source and receiver nodes
        src_id = src_nodes[i]
        dst_id = dst_nodes[src_id]

        # if a node does not flow to itself, and the source node is above the
        # destination node
        if src_id != dst_id and z[src_id] > z[dst_id]:

            # Get values for z at present node and present time,
            # and z downstream at t + delta t (which should have been
            # previously solved for)
            z_old = z[src_id]
            z_downstream = z[dst_id]

            # Get the threshold value. In this function, it is constant, so we
            # already have it.
            # threshsxdt = threshsxdt

            # calculate the difference between z_old and z_downstream
            z_diff_old = z_old - z_downstream

            # using z_diff_old, calculate the alpha paramter of Braun and
            # Willet by calculating alpha times z

            alpha_param = alpha[src_id] * pow(z_diff_old, n-1.0)

            # Calculate the beta parameter that accounts for the possible
            # presence of a threshold.
            beta_param = threshsxdt / z_diff_old
            # check if the threshold has been exceeded by passing a value of
            # x = 1 to the erode_fn. If this returns a value of less than
            # zero, this means that the the maximum possible slope value  does
            # not produce stream power needed to exceed the erosion threshold
            check_function = erode_fn(1, alpha_param, beta_param, n)

            # if the threshold was not exceeded do not change the elevation,
            # otherwise calculate the erosion rate
            if check_function > 0:
                # if the threshold was exceeded, then there will be a zero
                # between x = 0 and x= 1

                # solve using brentq, which requires a zero to exist
                # in between the two end values

                # if n is 1, finding x has an analytical solution. Otherwise,
                # use the the numerical solution given by root finding
                if n != 1.0:

                    # The threshold values passed here are the defaults if one
                    # were to import brentq from scipy.optimize
                    x = brentq(erode_fn,
                                       0.0, 1.0,
                                       1e-12,4.4408920985006262e-16,
                                       100,
                                       (alpha_param, beta_param, n),
                                       False,
                                       True)

                else:
                    # Analytical solution
                    x = (1.0 + beta_param)/(1.0 + alpha_param)

                # If x is provided as a value greater than zero, calculate
                # z at t=t+delta_t useing the values of x, z_downstream and
                # z_old as given by the definition of x (see erode_fn for
                # details). If x is equal to zero, set it as just slightly
                # higher than x_downstream.
                if x>0:
                    z[src_id] = z_downstream + x * (z_old - z_downstream)
                else:
                    z[src_id] = z_downstream + 1.0e-15

                # Nothing is returned from this function as it serves to update
                # the array z.


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

    If the threshold term, beta, is zero, this equation collapses to the form
    given by Braun and Willet (2012).

    When the threshold term is greater than zero, it is possible that no
    erosion will occur. In this case, an evaluation of f(x=1) will yeild a
    negative number.

    If for the values of alpha and beta provided, f(x=1)<0, then this function
    has no root and use in a solver such as Brent's method in which a solution
    interval is required will fail.

    It is recommended that before using this method in a solver, that the
    function be evaluated with x=1 to determine if it any erosion occured.

    """
    cdef double f

    f = x - 1.0 + (alpha * (x ** n)) - beta

    return f


def smooth_stream_power_eroder_solver(np.ndarray[DTYPE_INT_t, ndim=1] src_nodes,
                                      np.ndarray[DTYPE_INT_t, ndim=1] dst_nodes,
                                      np.ndarray[DTYPE_FLOAT_t, ndim=1] z,
                                      np.ndarray[DTYPE_FLOAT_t, ndim=1] alpha,
                                      np.ndarray[DTYPE_FLOAT_t, ndim=1] gamma,
                                      np.ndarray[DTYPE_FLOAT_t, ndim=1] delta):
    """Erode node elevations using Newtons Method for smoothed Stream Power. "

    This method takes three parameters, alpha, gamma, and delta. 

    alpha = K A^m dt / L
    
    delta = K A^m / (L * wc)
    
    gamma = omega_c * dt 

    This method will use the new_elev and new_elev_prime equations. 

    Parameters
    ----------
    src_nodes : array_like
        Ordered upstream node ids.
    dst_nodes : array_like
        Node ids of nodes receiving flow.
    alpha : array_like
        Erosion equation parameter. 
    gamma : array_like
        Erosion equation parameter. 
    delta : array_like
        Erosion equation parameter. 
    z : array_like
        Node elevations.
    """
    cdef unsigned int n_nodes = src_nodes.size
    cdef unsigned int src_id
    cdef unsigned int dst_id
    cdef unsigned int i

    cdef double epilon



    for i in range(len(src_nodes)):
        src_id = src_nodes[i]
        dst_id = dst_nodes[src_id]

        if src_id != dst_id and z[src_id] > z[dst_id]:

            # calculate epsilon
            epsilon = (alpha[src_id] * z[dst_id]
                       + gamma[src_id] + z[src_id])

            # calculate new z
            z[src_id] = newton(new_elev, z[src_id],
                             fprime=new_elev_prime,
                             args=(alpha[src_id],
                                   z[dst_id],
                                   gamma[src_id],
                                   delta[src_id],
                                   epsilon))


def new_elev(DTYPE_FLOAT_t x,
             DTYPE_FLOAT_t a,
             DTYPE_FLOAT_t b,
             DTYPE_FLOAT_t c,
             DTYPE_FLOAT_t d,
             DTYPE_FLOAT_t e):
    """Equation for elevation of a node at timestep t+1.

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


def new_elev_prime(DTYPE_FLOAT_t x,
                   DTYPE_FLOAT_t a,
                   DTYPE_FLOAT_t b,
                   DTYPE_FLOAT_t c,
                   DTYPE_FLOAT_t d,
                   DTYPE_FLOAT_t e):
    """Equation for elevation of a node at timestep t+1.

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
        Placeholder; not used
    """
    cdef double f

    f = (1.0 + a) - c * d * np.exp(-d * (x - b))

    return f


cpdef DTYPE_FLOAT_t sed_flux_fn_gen_genhump(DTYPE_FLOAT_t rel_sed_flux_in,
                                            DTYPE_FLOAT_t kappa,
                                            DTYPE_FLOAT_t nu,
                                            DTYPE_FLOAT_t c,
                                            DTYPE_FLOAT_t phi,
                                            DTYPE_FLOAT_t norm):
    """
    Returns f(qs,qc) assuming a generalized humped function per
    Hobley et al., 2011.

    Note that this function permits values outside those that are physically
    meaningful for the relative sediment flux, i.e., qs/qc < 0. or > 1.

    Parameters
    ----------
    rel_sed_flux_in : float
        The relative sediment flux; carried flux divided by capacity.
    kappa, nu, c, phi : float
        Parameters controlling the shape of the generalized nonlinear sediment
        flux curve.
    norm : float
        Parameter to normalize the curve such that its maximum is 1.
    
    Returns
    -------
    fqs : float
        The erosion efficiency, f(q_s).
    """
    return (
        norm*kappa*(rel_sed_flux_in**nu + c)*np.exp(-phi*rel_sed_flux_in))


cpdef DTYPE_FLOAT_t sed_flux_fn_gen_lindecl(DTYPE_FLOAT_t rel_sed_flux_in,
                                            DTYPE_FLOAT_t kappa,
                                            DTYPE_FLOAT_t nu,
                                            DTYPE_FLOAT_t c,
                                            DTYPE_FLOAT_t phi,
                                            DTYPE_FLOAT_t norm):
    """
    Returns f(qs,qc) assuming a linear decline model (see e.g. Gasparini
    et al., 2006).
    kappa, nu, c, phi, & norm are all dummy variables here.

    Note that this function permits values outside those that are physically
    meaningful for the relative sediment flux, i.e., qs/qc < 0. or > 1.
    
    Parameters
    ----------
    rel_sed_flux_in : float
        The relative sediment flux; carried flux divided by capacity.
    kappa, nu, c, phi, norm : float
        Unused placeholder parameters.
    
    Returns
    -------
    fqs : float
        The erosion efficiency, f(q_s).
    """
    return 1.-rel_sed_flux_in


cpdef DTYPE_FLOAT_t sed_flux_fn_gen_almostparabolic(
                                                DTYPE_FLOAT_t rel_sed_flux_in,
                                                DTYPE_FLOAT_t kappa,
                                                DTYPE_FLOAT_t nu,
                                                DTYPE_FLOAT_t c,
                                                DTYPE_FLOAT_t phi,
                                                DTYPE_FLOAT_t norm):
    """
    Returns f(qs,qc) assuming the almost parabolic humped model (see
    Gasparini et al., 2006).
    kappa, nu, c, phi, & norm are all dummy variables here.

    Note that this function permits values outside those that are physically
    meaningful for the relative sediment flux, i.e., qs/qc < 0. or > 1.

    Parameters
    ----------
    rel_sed_flux_in : float
        The relative sediment flux; carried flux divided by capacity.
    kappa, nu, c, phi, norm : float
        Unused placeholder parameters.
    
    Returns
    -------
    fqs : float
        The erosion efficiency, f(q_s).
    """
    return np.where(rel_sed_flux_in > 0.1,
                    1. - 4.*(rel_sed_flux_in-0.5)**2.,
                    2.6*rel_sed_flux_in+0.1)


cpdef DTYPE_FLOAT_t sed_flux_fn_gen_const(DTYPE_FLOAT_t rel_sed_flux_in,
                                          DTYPE_FLOAT_t kappa,
                                          DTYPE_FLOAT_t nu,
                                          DTYPE_FLOAT_t c,
                                          DTYPE_FLOAT_t phi,
                                          DTYPE_FLOAT_t norm):
    """
    Returns 1, and thus no sed flux effects.
    kappa, nu, c, phi, & norm are all dummy variables here.

    Note that this function permits values outside those that are physically
    meaningful for the relative sediment flux, i.e., qs/qc < 0. or > 1.

    Parameters
    ----------
    rel_sed_flux_in : float
        The relative sediment flux; carried flux divided by capacity.
    kappa, nu, c, phi, norm : float
        Unused placeholder parameters.
    
    Returns
    -------
    fqs : float
        The erosion efficiency, f(q_s).
    """
    return 1.


cpdef void get_sed_flux_function_pseudoimplicit_bysedout(
        DTYPE_FLOAT_t sed_in_bydt,
        DTYPE_FLOAT_t trans_cap_vol_out_bydt,
        DTYPE_FLOAT_t prefactor_for_volume_bydt,
        DTYPE_FLOAT_t cell_area,
        sed_flux_fn_gen,
        DTYPE_FLOAT_t kappa, DTYPE_FLOAT_t nu, DTYPE_FLOAT_t c,
        DTYPE_FLOAT_t phi, DTYPE_FLOAT_t norm,
        DTYPE_INT_t pseudoimplicit_repeats,
        np.ndarray[DTYPE_FLOAT_t, ndim=1] out_array):
    """
    This function uses a pseudoimplicit method to calculate the sediment
    flux function for a node, and also returns dz/dt and the rate of
    sediment output from the node. This version stabilises the sediment
    out of the node, rather than the sed flux function itself.

    Note that this method now operates in PER TIME units; this was not
    formerly the case.

    Parameters
    ----------
    sed_in_bydt : float
        Total rate of incoming sediment, sum(Q_s_in)/dt
    trans_cap_vol_out_bydt : float
        Volumetric transport capacity as a rate (i.e., m**3/s) on outgoing
        link
    prefactor_for_volume_bydt : float
        Equal to K*A**m*S**n * cell_area.
    cell_area : float
        The area of the cell.
    sed_flux_fn_gen : function
        Function to calculate the sed flux function. Takes inputs
        rel_sed_flux_in, kappa, nu, c, phi, norm, where last 5 are dummy
        unless type is generalized_humped.
    kappa, nu, c, phi, norm : float
        Params for the sed flux function. Values if generalized_humped,
        zero otherwise.
    pseudoimplicit_repeats : int
        Maximum number of loops to perform with the pseudoimplicit
        iterator, seeking a stable solution. Convergence is typically
        rapid.
        out_array : array of floats
            Array to be filled, containing:
            dzbydt: Rate of change of substrate elevation,
            vol_pass_rate: Q_s/dt on the outgoing link,
            rel_sed_flux: f(Q_s/Q_c),
            error_in_sed_flux_fn: Measure of how well converged rel_sed_flux is
    """
    cdef unsigned int i
    cdef double rel_sed_flux_in
    cdef double last_rel_sed_flux
    cdef double rel_sed_flux
    cdef double sed_flux_fn
    cdef double sed_vol_added_bydt
    cdef double new_sed_vol_added_bydt
    cdef double prop_added
    cdef double error_in_sed_vol_added
    cdef double excess_trans_capacity

    excess_trans_capacity = trans_cap_vol_out_bydt - sed_in_bydt
    if excess_trans_capacity < 1.e-10:   #  can be -ve, note
        out_array[0] = 0.
        out_array[1] = trans_cap_vol_out_bydt
        out_array[2] = 1.  # arbitrary; probably more stable in later use
        out_array[3] = trans_cap_vol_out_bydt
    else:
        rel_sed_flux_in = sed_in_bydt / trans_cap_vol_out_bydt
        if rel_sed_flux_in > 1.:
            rel_sed_flux_in = 1.
        last_sed_flux = rel_sed_flux_in
        sed_flux_fn = sed_flux_fn_gen(
            rel_sed_flux_in, kappa, nu, c, phi, norm)
        sed_vol_added_bydt = 0.  # prefactor_for_volume_bydt * sed_flux_fn

        for i in range(pseudoimplicit_repeats):
            prop_added = sed_vol_added_bydt / trans_cap_vol_out_bydt
            rel_sed_flux = prop_added + last_rel_sed_flux
            if rel_sed_flux < 0.:
                rel_sed_flux = 0.
            if rel_sed_flux > 1.:
                rel_sed_flux = 1.
            rel_sed_flux = 0.5 * (rel_sed_flux + rel_sed_flux_in)

            sed_flux_fn = sed_flux_fn_gen(
                rel_sed_flux, kappa, nu, c, phi, norm)
            new_sed_vol_added_bydt = prefactor_for_volume_bydt * sed_flux_fn
            if new_sed_vol_added_bydt > excess_trans_capacity:
                new_sed_vol_added_bydt = excess_trans_capacity
            error_in_sed_vol_added = abs(
                new_sed_vol_added_bydt - sed_vol_added_bydt
            )  # absolute, as ratios crash at low erosion rates
            if error_in_sed_vol_added < 1.e-3:
                break
            last_rel_sed_flux = rel_sed_flux
            sed_vol_added_bydt = new_sed_vol_added_bydt
            
        # note that the method will silently terminate even if we still have
        # bad convergence. Note this is very rare.

        out_array[0] = new_sed_vol_added_bydt / cell_area
        out_array[1] = sed_in_bydt + new_sed_vol_added_bydt  # sed passed
        out_array[2] = (
            0.5 * (out_array[1] + sed_in_bydt) / trans_cap_vol_out_bydt
        )
        out_array[3] = error_in_sed_vol_added


cpdef void iterate_sde_downstream(
                np.ndarray[DTYPE_INT_t, ndim=1] s_in,
                np.ndarray[DTYPE_FLOAT_t, ndim=1] cell_areas,
                np.ndarray[DTYPE_FLOAT_t, ndim=1] hillslope_sediment,
                np.ndarray[DTYPE_FLOAT_t, ndim=1] hillslope_sediment_flux,
                np.ndarray[DTYPE_FLOAT_t, ndim=1] river_volume_flux_into_node,
                np.ndarray[DTYPE_FLOAT_t, ndim=1] transport_capacities,
                np.ndarray[DTYPE_FLOAT_t, ndim=1] erosion_prefactor_withS,
                np.ndarray[DTYPE_FLOAT_t, ndim=1] rel_sed_flux,
                np.ndarray[dtype=np.int8_t, ndim=1] is_it_TL,
                np.ndarray[DTYPE_FLOAT_t, ndim=1] vol_drop_rate,
                np.ndarray[DTYPE_INT_t, ndim=1] flow_receiver,
                DTYPE_INT_t pseudoimplicit_repeats,
                np.ndarray[DTYPE_FLOAT_t, ndim=1] dzbydt,
                sed_flux_fn_gen,
                DTYPE_FLOAT_t kappa,
                DTYPE_FLOAT_t nu,
                DTYPE_FLOAT_t c,
                DTYPE_FLOAT_t phi,
                DTYPE_FLOAT_t norm):
    """
    Iterates down a drainage network, redistributing sediment and solving
    the sediment flux dependent incision equations.

    Parameters
    ----------
    s_in : array
        The upstream node order
    cell_areas : array
        The areas of all cells in the grid
    hillslope_sediment : array
        Depth of sediment in the channel at each node.
    hillslope_sediment_flux : array
        The existing volume of sediment on the channel bed at a node,
        expressed as volume per unit time of the timestep. This turns
        the accumulated sediment on the bed into a virtual sediment supply
        from upstream, such that at the end of the step the same depth of
        sediment would be present at the node if no other transport
        occurred.
    river_volume_flux_into_node : array
        Total ""true" river flux coming into node from upstream.
    transport_capacities : array
        The bedload transport capacity at each node, expressed as a flux.
    erosion_prefactor_withS : array
        Equal to K * A**m * S**n at nodes
    rel_sed_flux : array
        The sediment flux as a function of the transport capacity.
    is_it_TL : boolean array
        Describes whether the sediment transported at the node is at
        capacity or not.
    vol_drop_rate : array
        Flux of sediment ending up on the bed during transport.
    flow_receiver : array
        The downstream node ID.
    pseudoimplicit_repeats : int
        Maximum number of loops to perform with the pseudoimplicit
        iterator, seeking a stable solution. Convergence is typically
        rapid.
    dzbydt : array
        The rate of change of *bedrock* surface elevation.
    sed_flux_fn_gen : function
        Function to calculate the sed flux function. Takes inputs
        rel_sed_flux_in, kappa, nu, c, phi, norm, where last 5 are dummy
        unless type is generalized_humped.
    kappa, nu, c, phi, norm : float
        Params for the sed flux function. Values if generalized_humped,
        zero otherwise.
    """
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] out_array = np.empty(4, dtype=float)
    cdef unsigned int i
    cdef double cell_area
    cdef double flood_depth_flux
    cdef double sed_flux_into_this_node_bydt
    cdef double node_capacity
    cdef double vol_prefactor_bydt
    cdef double vol_pass_rate
    cdef double depth_sed_in
    
    for i in s_in[::-1]:  # work downstream
        cell_area = cell_areas[i]
        sed_flux_into_this_node_bydt = (
            hillslope_sediment_flux[i] +
            river_volume_flux_into_node[i])
        node_capacity = transport_capacities[i]
        # ^we work in volume discharge, not volume per se here

        if sed_flux_into_this_node_bydt < node_capacity:
            # ^note incision is forbidden at capacity
            vol_prefactor_bydt = erosion_prefactor_withS[i]*cell_area
            get_sed_flux_function_pseudoimplicit_bysedout(
                    sed_flux_into_this_node_bydt,
                    node_capacity,
                    vol_prefactor_bydt, cell_area,
                    sed_flux_fn_gen,
                    kappa, nu, c, phi, norm,
                    pseudoimplicit_repeats, out_array)
            dzbydt[i] = -out_array[0]
            # ^minus returns us to the correct sign convention
            vol_pass_rate = out_array[1]
            rel_sed_flux[i] = out_array[2]
            # error_in_sed_flux = out_array[3]
        else:
            is_it_TL[i] = 1
            rel_sed_flux[i] = 1.
            dzbydt[i] = 0.
            vol_pass_rate = node_capacity
            vol_drop_rate[i] = sed_flux_into_this_node_bydt - vol_pass_rate
        
        assert vol_drop_rate[i] >= 0.
        river_volume_flux_into_node[flow_receiver[i]] += vol_pass_rate
