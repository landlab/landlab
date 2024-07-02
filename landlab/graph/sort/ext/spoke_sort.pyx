cimport cython

from cython.parallel import prange

from libc.math cimport M_PI
from libc.math cimport atan2
from libc.stdlib cimport free
from libc.stdlib cimport malloc

from .argsort cimport argsort_flt

ctypedef fused id_t:
    cython.integral
    long long


cdef void _calc_spoke_angles(
    cython.floating * hub,
    cython.floating * spokes,
    int n_spokes,
    cython.floating * angles,
) noexcept nogil:
    cdef int i
    cdef double x0 = hub[0]
    cdef double y0 = hub[1]
    cdef cython.floating * spoke = spokes
    cdef double two_pi = 2. * M_PI
    cdef double x
    cdef double y

    for i in range(n_spokes):
        x = spoke[0]
        y = spoke[1]

        angles[i] = atan2(y - y0, x - x0)
        if angles[i] < 0.:
            angles[i] += two_pi
        spoke += 2


cdef void _argsort_spokes_around_hub(
    id_t * spokes,
    int n_spokes,
    cython.floating * xy_of_spoke,
    cython.floating * xy_of_hub,
    id_t * ordered,
) noexcept nogil:
    cdef int point
    cdef int spoke
    cdef cython.floating * points = <cython.floating *>malloc(
        2 * n_spokes * sizeof(cython.floating)
    )
    cdef cython.floating * angles = <cython.floating *>malloc(
        n_spokes * sizeof(cython.floating)
    )

    try:
        point = 0
        for spoke in range(n_spokes):
            points[point] = xy_of_spoke[2 * spokes[spoke]]
            points[point + 1] = xy_of_spoke[2 * spokes[spoke] + 1]
            point += 2

        _calc_spoke_angles(xy_of_hub, points, n_spokes, angles)
        argsort_flt(angles, n_spokes, ordered)
    finally:
        free(angles)
        free(points)


cdef void _sort_spokes_around_hub(
    id_t * spokes,
    int n_spokes,
    cython.floating * xy_of_spoke,
    cython.floating * xy_of_hub,
) noexcept nogil:
    cdef int spoke
    cdef id_t * ordered = <id_t *>malloc(n_spokes * sizeof(id_t))
    cdef id_t * temp = <id_t *>malloc(n_spokes * sizeof(id_t))

    try:
        _argsort_spokes_around_hub(spokes, n_spokes, xy_of_spoke, xy_of_hub, ordered)

        for spoke in range(n_spokes):
            temp[spoke] = spokes[ordered[spoke]]

        for spoke in range(n_spokes):
            spokes[spoke] = temp[spoke]
    finally:
        free(temp)
        free(ordered)


@cython.boundscheck(False)
@cython.wraparound(False)
def calc_spoke_angles(
    cython.floating [:] hub,
    cython.floating [:] spokes,
    cython.floating [:] angles,
):
    cdef int n_spokes = spokes.shape[0] // 2
    _calc_spoke_angles(&hub[0], &spokes[0], n_spokes, &angles[0])


@cython.boundscheck(False)
@cython.wraparound(False)
def sort_spokes_around_hub(
    id_t [:] spokes,
    cython.floating [:, :] xy_of_spoke,
    cython.floating [:] xy_of_hub,
):
    cdef int n_spokes = spokes.shape[0]

    _sort_spokes_around_hub(&spokes[0], n_spokes, &xy_of_spoke[0, 0], &xy_of_hub[0])


@cython.boundscheck(False)
@cython.wraparound(False)
def argsort_points_around_hub(
    cython.floating [:, :] points,
    cython.floating [:] hub,
    id_t [:] out,
):
    """Sort spokes by angle around a hub.

    Parameters
    ----------
    points : ndarray of float, shape `(n_points, 2)`
        Coordinates of points as (*x*, *y*).
    out : ndarray of int, shape `(n_points, )`
        Indices of sorted points.

    Returns
    -------
    ndarray of int, shape `(n_points, )`
        Indices of sorted points.
    """
    cdef int n_points = points.shape[0]
    cdef cython.floating *angles = <cython.floating *>malloc(
        n_points * sizeof(cython.floating)
    )

    try:
        _calc_spoke_angles(&hub[0], &points[0, 0], n_points, angles)
        argsort_flt(angles, n_points, &out[0])
    finally:
        free(angles)

    return out


@cython.boundscheck(False)
@cython.wraparound(False)
def argsort_spokes_at_wheel(
    id_t [:] spokes_at_wheel,
    id_t [:] offset_to_wheel,
    cython.floating [:, :] xy_of_hub,
    cython.floating [:, :] xy_of_spoke,
    id_t [:] ordered,
):
    cdef long n_wheels = len(offset_to_wheel) - 1
    cdef long i
    cdef id_t n_spokes
    cdef id_t * wheel
    cdef id_t * order

    for i in prange(n_wheels, nogil=True, schedule="static"):
        n_spokes = offset_to_wheel[i + 1] - offset_to_wheel[i]
        wheel = &spokes_at_wheel[offset_to_wheel[i]]
        order = &ordered[offset_to_wheel[i]]

        _argsort_spokes_around_hub(
            wheel, n_spokes, &xy_of_spoke[0, 0], &xy_of_hub[i, 0], order
        )


@cython.boundscheck(False)
@cython.wraparound(False)
def sort_spokes_at_wheel(
    id_t [:] spokes_at_wheel,
    id_t [:] offset_to_wheel,
    cython.floating [:, :] xy_of_hub,
    cython.floating [:, :] xy_of_spoke,
):
    """Sort spokes about multiple hubs.

    Parameters
    ----------
    spokes_at_wheel : ndarray of int
        Spokes for each wheel.
    offset_to_wheel : ndarray of int
        Offset into *spokes_at_wheel* for each wheel.
    xy_of_hub : ndarray of float, shape `(n_hubs, 2)`
        Coordinates of each hub as `(x, y)`.
    xy_of_spoke : ndarray of float, shape `(n_spokes, 2)`
        Coordinates of the end of each spoke as `(x, y)`.
    """
    cdef int n_wheels = len(offset_to_wheel) - 1
    cdef int i
    cdef int n_spokes
    cdef int spoke
    cdef id_t * wheel

    for i in prange(n_wheels, nogil=True, schedule="static"):
        wheel = &spokes_at_wheel[offset_to_wheel[i]]
        n_spokes = offset_to_wheel[i + 1] - offset_to_wheel[i]
        for spoke in range(n_spokes):
            if wheel[spoke] == -1:
                n_spokes = spoke
                break

        _sort_spokes_around_hub(wheel, n_spokes, &xy_of_spoke[0, 0], &xy_of_hub[i, 0])
