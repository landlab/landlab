cimport cython
from cython.parallel cimport prange
from libc.math cimport cbrt
from libc.math cimport copysign
from libc.math cimport fabs
from libc.math cimport fmin
from libc.math cimport sqrt
from libc.stdint cimport int32_t
from libc.stdint cimport int64_t
from libc.stdint cimport uint8_t

from landlab.grid.linkstatus import LinkStatus

ctypedef fused id_t:
    int32_t
    int64_t


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def adjust_unstable_discharge(
    cython.floating [::1] q_at_link,
    const cython.floating [::1] h_at_link,
    *,
    const double dx,
    const double dt,
    const id_t [::1] where,
    cython.floating [::1] out,
):
    """Clamp discharge to maintain stability under a shallow-water CFL constraint.

    Limit the magnitude of discharge such that it does not exceed a scaled
    Courant–type stability threshold. The maximum allowable discharge is
    determined as::

        q_max = factor * (h / 4) * (dx / dt)

    where ``dx`` is the grid spacing, ``dt`` is the timestep, ``h`` is water depth at
    the link, and ``factor`` is a safety coefficient (currently set to ``0.8``). The
    value of ``q_at_link`` is clamped to ``[-q_max, +q_max]`` while preserving its sign.

    Parameters
    ----------
    q_at_link : ndarray of float, shape (n_links,)
        Water discharge at each link.
    h_at_link : ndarray of float, shape (n_links,)
        Water depth at each link.
    dx : float
        Grid spacing (must be > 0).
    dt : double
        Timestep (must be > 0).
    where : ndarray of int
        Indices of links for which discharge should be adjusted.
    out : ndarray of float, shape (n_links,)
        Output array into which adjusted discharge values are written. May be the
        same memory as ``q_at_link`` for in-place updates.

    Returns
    -------
    ndarray of float
        The adjusted discharge array.
    """
    cdef Py_ssize_t n_links = where.shape[0]
    cdef Py_ssize_t i
    cdef id_t link
    cdef double dx_over_dt = dx / dt
    cdef double factor_of_safety = 0.8
    cdef double q_threshold
    cdef double q
    cdef double h

    for i in prange(n_links, nogil=True, schedule="static"):
        link = where[i]
        q = q_at_link[link]
        h = h_at_link[link]

        q_threshold = 0.25 * h * dx_over_dt * factor_of_safety
        out[link] = copysign(fmin(fabs(q), q_threshold), q)

    return (<object>out).base


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


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def adjust_supercritical_discharge(
    cython.floating [::1] q_at_link,
    const cython.floating [::1] h_at_link,
    *,
    const double g,
    const double froude,
    const id_t [::1] where,
    cython.floating [::1] out,
):
    """Reduce discharge at selected links to enforce supercritical-flow limits.

    Clamp discharge values to satisfy the stability condition based on
    the specified Froude number. For each link, discharge is limited to

        abs(q) <= froude * h * sqrt(g * h)

    which corresponds to maintaining a Froude number that does not exceed the
    specified threshold.

    Parameters
    ----------
    q_at_link : cython.floating[::1]
        Discharge at each link.
    h_at_link : cython.floating[::1]
        Water depth at each link.
    g : double
        Gravitational acceleration (must be positive).
    froude : double
        Maximum allowed Froude number (must be non-negative). Values of
        discharge that would imply a higher Froude number are reduced.
    where : cython integral array
        Array of link indices at which discharge should be adjusted.
    out : cython.floating[::1]
        Output array to receive adjusted values. `out` may be the same
        array as `q_at_link`, in which case changes are made in-place.

    Returns
    -------
    ndarray of float
        The modified discharge values.
    """
    cdef long n_links = len(where)
    cdef long i
    cdef long link
    cdef double factor = froude * sqrt(g)
    cdef double q
    cdef double h

    for i in prange(n_links, nogil=True, schedule="static"):
        link = where[i]
        q = q_at_link[link]
        h = h_at_link[link]

        out[link] = copysign(fmin(fabs(q), factor * h * sqrt(h)), q)

    return (<object>out).base


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def calc_discharge_at_links(
    const cython.floating [::1] q_at_link,
    const cython.floating [::1] q_mean_at_link,
    const cython.floating [::1] h_at_link,
    const cython.floating [::1] water_slope_at_link,
    const cython.floating [::1] mannings_at_link,
    *,
    const double g,
    const double dt,
    const id_t [::1] where,
    cython.floating [::1] out,
):
    """Calculate discharge values at links.

    For each link index in ``where``, the discharge is computed using::

        q_new = (
            q_mean - g * dt * h * slope
        ) / (
            1 + g * dt * mannings**2 * abs(q) / h**(7 / 3)
        )

    Note that this formula is only applied when ``h > 0.0``. For negative or zero
    flow depths, the ``out`` array is unmodified. This is Equation 23 [1]_.

    Parameters
    ----------
    q_at_link : cython.floating[::1]
        Current discharge values at links.
    q_mean_at_link : cython.floating[::1]
        Averaged discharge values.
    h_at_link : cython.floating[::1]
        Water depth at each link.
    water_slope_at_link : cython.floating[::1]
        Local water surface slope at each link.
    mannings_at_link : cython.floating[::1]
        Manning roughness coefficient.
    g : double
        Gravitational acceleration.
    dt : double
        Model timestep.
    where : cython integral array
        Indices of links to update.
    out : cython.floating[::1]
        Output array where results are written. Can be the same array as
        ``q_at_link`` for in-place modification.

    Returns
    -------
    ndarray
        Array of updated discharge values.

    References
    ----------
    .. [1] De Almeida, Gustavo AM, et al. "Improving the stability
       of a simple formulation of the shallow water equations for 2‐D
       flood modeling." Water Resources Research 48.5 (2012).
       https://doi.org/10.1029/2011WR011570
    """
    cdef Py_ssize_t n_links = where.shape[0]
    cdef Py_ssize_t link
    cdef Py_ssize_t i
    cdef double h_to_seven_thirds
    cdef double numerator
    cdef double denominator
    cdef double mannings
    cdef double h
    cdef double g_times_dt = g * dt
    cdef double h_min = 1e-6

    for i in prange(n_links, nogil=True, schedule="static"):
        link = where[i]
        h = h_at_link[link]

        if h > h_min:
            mannings = mannings_at_link[link]
            h_to_seven_thirds = h * h * cbrt(h)

            numerator = (
                q_mean_at_link[link] - g_times_dt * h * water_slope_at_link[link]
            )

            denominator = (
                1.0
                + g_times_dt * mannings * mannings * fabs(q_at_link[link])
                / h_to_seven_thirds
            )
            out[link] = numerator / denominator
        elif h > 0.0:
            out[link] = 0.0

    return (<object>out).base


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def calc_grad_at_link(
    const cython.floating [::1] value_at_node,
    *,
    const cython.floating [::1] length_of_link,
    const id_t [:, ::1] nodes_at_link,
    const id_t [::1] where,
    cython.floating [::1] out,
):
    """Calculate the gradient of node-values at specified links.

    The function computes finite differences between the head and tail
    nodes of each link, divided by the link length. Only links listed in
    ``where`` are updated; all other elements in ``out`` are left
    unchanged.

    Parameters
    ----------
    value_at_node : array_like of float
        Node-centered values.
    length_of_link : array_like of float
        Length of each link.
    nodes_at_link : array_like of int, shape (n_links, 2)
        Head and tail node for each link. The first column is the
        tail node and the second column is the head node.
    where : array_like of int
        Indices of links to update.
    out : ndarray of float
        Output buffer of shape ``(n_links,)``. Results are written in
        place to the elements indexed by ``where``.

    Returns
    -------
    out : ndarray of float
        Array of gradients at links corresponding to those indexed
        by ``where``. This is the same object as the input ``out`` array.
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

        out[link] = (
            value_at_node[head] - value_at_node[tail]
        ) / length_of_link[link]

    return (<object>out).base


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def calc_weighted_mean_of_parallel_links(
    const cython.floating[::1] value_at_link,
    const id_t[:, ::1] parallel_links_at_link,
    *,
    const uint8_t[::1] status_at_link,
    const id_t[::1] where,
    const double theta,
    cython.floating[::1] out,
):
    """Calculate a weighted mean of each link and its parallel neighbors.

    For each link in ``where``, compute a weighted combination of that link's
    value and the values of its left and right parallel links. The combination
    weight is controlled by ``theta`` such that::

        out[link] = theta * value_at_link[link] + (1 - theta) / 2 * (
            value_at_link[left] + value_at_link[right]
        )

    Links with invalid or inactive neighbors retain their original value.

    Parameters
    ----------
    value_at_link : array_like of float
        Values defined at links. Must have length equal to the number of links.
    parallel_links_at_link : ndarray of int, shape (n_links, 2)
        IDs of the two parallel links (left, right) associated with each link.
        A value of ``-1`` indicates that no parallel link exists.
    status_at_link : ndarray of uint8
        Link status codes. A value of ``LinkStatus.INACTIVE`` indicates an
        inactive link that should not contribute to the weighted mean.
    where : ndarray of int
        Indices of links for which to update ``out``. Each link in ``where``
        is processed exactly once. Must not contain duplicates.
    theta : float
        Weighting coefficient between 0 and 1. A value of 1.0 means the link's
        own value is kept unchanged; a value of 0.0 means a simple mean of the
        two parallel links.
    out : ndarray of float
        Output buffer into which results are written. Must be the same shape
        and dtype as ``value_at_link``.

    Returns
    -------
    out : ndarray of float
        Array of means at links corresponding to those indexed
        by ``where``. This is the same object as the input ``out`` array.
    """
    cdef long n_links = len(where)
    cdef long link
    cdef long left
    cdef long right
    cdef long i
    cdef long n_neighbors
    cdef double total
    cdef double ONE_MINUS_THETA = (1.0 - theta)
    cdef int LINK_IS_INACTIVE = LinkStatus.INACTIVE
    cdef int LINK_IS_MISSING = -1

    for i in prange(n_links, nogil=True, schedule="static"):
        link = where[i]
        left = parallel_links_at_link[link, 0]
        right = parallel_links_at_link[link, 1]

        total = 0.0
        n_neighbors = 0

        if left != LINK_IS_MISSING and status_at_link[left] != LINK_IS_INACTIVE:
            n_neighbors = n_neighbors + 1
            total = total + value_at_link[left]

        if right != LINK_IS_MISSING and status_at_link[right] != LINK_IS_INACTIVE:
            n_neighbors = n_neighbors + 1
            total = total + value_at_link[right]

        if n_neighbors == 2:
            out[link] = theta * value_at_link[link] + ONE_MINUS_THETA * total * 0.5
        elif n_neighbors == 1:
            out[link] = theta * value_at_link[link] + ONE_MINUS_THETA * total
        else:
            out[link] = value_at_link[link]

    return (<object>out).base


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def zero_out_dry_links(
    const cython.floating [::1] h_at_link,
    *,
    const id_t [::1] where,
    cython.floating [::1] out,
):
    """Zero values for links with non-positive depth.

    Sets ``out[link] = 0.0`` for each link where the corresponding
    entry in ``h_at_link`` is less than or equal to zero, using the
    indices provided in ``where``. Links with positive depth are
    left unchanged.

    Parameters
    ----------
    h_at_link : array_like of float, shape (n_links,)
        Water depth at each link.
    where : array_like of int
        Indices of links to operate on. Usually a subset of links,
        such as those active in the current routing step.
    out : ndarray of float, shape (n_links,)
        Output array into which results are written. Must be
        pre-initialized, entries at dry links are overwritten;
        all other entries are left unchanged.

    Returns
    -------
    out : ndarray of float
        The updated ``out`` array (same object as input), with values
        set to zero at dry links.
    """
    cdef Py_ssize_t n_links = where.shape[0]
    cdef Py_ssize_t link
    cdef Py_ssize_t i

    for i in prange(n_links, nogil=True, schedule="static"):
        link = where[i]
        if h_at_link[link] <= 0.0:
            out[link] = 0.0

    return (<object>out).base
