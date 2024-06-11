import numpy as np

cimport numpy as np

from scipy.optimize import newton

# suspect that the function _brentq is in c and thus this is the most effective
# method for using the brentq method in cython.
from scipy.optimize._zeros import _brentq as brentq

# from libc.math cimport fabs


DTYPE_FLOAT = np.double
ctypedef np.double_t DTYPE_FLOAT_t

DTYPE_INT = int
ctypedef np.int_t DTYPE_INT_t


cdef extern from "math.h":
    double fabs(double x) nogil
    double pow(double x, double y) nogil


def brent_method_erode_variable_threshold(
    np.ndarray[DTYPE_INT_t, ndim=1] src_nodes,
    np.ndarray[DTYPE_INT_t, ndim=1] dst_nodes,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] threshsxdt,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] alpha,
    DTYPE_FLOAT_t n,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] z,
):
    """Erode node elevations using Brent's method for stability.

    The alpha value is given as

    alpha = delta_t*K * (A)**m/(delta_x**n)

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
    cdef unsigned int n_nodes = src_nodes.shape[0]
    cdef unsigned int src_id
    cdef unsigned int dst_id
    cdef unsigned int i
    cdef double z_old
    cdef double z_downstream
    cdef double thresholddt
    cdef double z_diff_old
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
                    x = brentq(
                        erode_fn,
                        0.0,
                        1.0,
                        1e-12,
                        4.4408920985006262e-16,
                        100,
                        (alpha_param, beta_param, n),
                        False,
                        True,
                    )

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


def brent_method_erode_fixed_threshold(
    np.ndarray[DTYPE_INT_t, ndim=1] src_nodes,
    np.ndarray[DTYPE_INT_t, ndim=1] dst_nodes,
    DTYPE_FLOAT_t threshsxdt,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] alpha,
    DTYPE_FLOAT_t n,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] z,
):
    """Erode node elevations.

    The alpha value is given as::

        alpha = delta_t*K * (A)**m/(delta_x**n)

    It will be multiplied by the value::

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
    cdef unsigned int n_nodes = src_nodes.shape[0]
    cdef unsigned int src_id
    cdef unsigned int dst_id
    cdef unsigned int i
    cdef double z_old
    cdef double z_downstream
    cdef double z_diff_old
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
                    x = brentq(
                        erode_fn,
                        0.0,
                        1.0,
                        1e-12,
                        4.4408920985006262e-16,
                        100,
                        (alpha_param, beta_param, n),
                        False,
                        True,
                    )

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


def erode_fn(
    DTYPE_FLOAT_t x,
    DTYPE_FLOAT_t alpha,
    DTYPE_FLOAT_t beta,
    DTYPE_FLOAT_t n,
):
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

    Consider stream power erosion under the equation::

        E = K * (A)**m * S**n - threshold_sp,

    on a grid with link delta_x and for a timestep of delta_t.

    When iterating from downstream to upstream in the drainage stack, at a
    given node at time = t+delta_t, the value of the node at time = t, and the
    value of the downstream node at time t+delta_t is known.

    Define::

        x = (
            z_node(t + delta_t) - z_downstream(t + delta_t)
        ) / (z_node(t) - z_downstream(t + delta_t))

    A discretized version of the stream power equation above yeilds the equation::

        f = x - 1 + alpha*(x**n) - beta

    where::

        alpha = delta_t * K * A ** m / (delta_x ** n) * (
            z_node(t) - z_downstream(t + delta_t)
        ) ** (n - 1)

    and::

        beta = threshold_sp * delta_t / (z_node(t) - z_downstream(t + delta_t))

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

    f = x - 1.0 + (alpha * pow(x, n)) - beta

    return f


def smooth_stream_power_eroder_solver(
    np.ndarray[DTYPE_INT_t, ndim=1] src_nodes,
    np.ndarray[DTYPE_INT_t, ndim=1] dst_nodes,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] z,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] alpha,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] gamma,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] delta,
):
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
    cdef unsigned int src_id
    cdef unsigned int dst_id
    cdef unsigned int i

    for i in range(len(src_nodes)):
        src_id = src_nodes[i]
        dst_id = dst_nodes[src_id]

        if src_id != dst_id and z[src_id] > z[dst_id]:

            # calculate epsilon
            epsilon = (alpha[src_id] * z[dst_id] + gamma[src_id] + z[src_id])

            # calculate new z
            z[src_id] = newton(
                new_elev,
                z[src_id],
                fprime=new_elev_prime,
                args=(
                    alpha[src_id],
                    z[dst_id],
                    gamma[src_id],
                    delta[src_id],
                    epsilon,
                )
            )


def new_elev(
    DTYPE_FLOAT_t x,
    DTYPE_FLOAT_t a,
    DTYPE_FLOAT_t b,
    DTYPE_FLOAT_t c,
    DTYPE_FLOAT_t d,
    DTYPE_FLOAT_t e,
):
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


def new_elev_prime(
    DTYPE_FLOAT_t x,
    DTYPE_FLOAT_t a,
    DTYPE_FLOAT_t b,
    DTYPE_FLOAT_t c,
    DTYPE_FLOAT_t d,
    DTYPE_FLOAT_t e,
):
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
