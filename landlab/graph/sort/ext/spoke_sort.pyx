import numpy as np
cimport numpy as np
cimport cython

from libc.stdlib cimport malloc, free, qsort
from libc.math cimport atan2

from .argsort cimport argsort


cdef _calc_spoke_angles(double * hub, double * spokes, np.int_t n_spokes,
                        double * angles):
    cdef int i
    cdef double x0 = hub[0]
    cdef double y0 = hub[1]
    cdef double * spoke = spokes
    cdef double two_pi = 2. * np.pi
    cdef double x
    cdef double y

    for i in range(n_spokes):
        x = spoke[0]
        y = spoke[1]

        angles[i] = atan2(y - y0, x - x0)
        if angles[i] < 0.:
            angles[i] += two_pi
        spoke += 2


cdef _argsort_spokes_around_hub(long * spokes, int n_spokes,
                                double * xy_of_spoke, double * xy_of_hub,
                                int * ordered):
    cdef int point
    cdef int spoke
    cdef double * points = <double *>malloc(2 * n_spokes * sizeof(double))
    cdef double * angles = <double *>malloc(n_spokes * sizeof(double))
    # cdef int * ordered = <int *>malloc(n_spokes * sizeof(int))
    # cdef int * temp = <int *>malloc(n_spokes * sizeof(int))
    
    try:
        point = 0
        for spoke in range(n_spokes):
            points[point] = xy_of_spoke[2 * spokes[spoke]]
            points[point + 1] = xy_of_spoke[2 * spokes[spoke] + 1]
            point += 2

        _calc_spoke_angles(xy_of_hub, points, n_spokes, angles)
        argsort(angles, n_spokes, ordered)

        # for spoke in range(n_spokes):
        #     temp[spoke] = spokes[ordered[spoke]]

        # for spoke in range(n_spokes):
        #     spokes[spoke] = temp[spoke]
    finally:
        free(angles)
        # free(temp)
        # free(ordered)
        free(points)


cdef _sort_spokes_around_hub(long * spokes, int n_spokes, double * xy_of_spoke,
                             double * xy_of_hub):
    cdef int point
    cdef int spoke
    # cdef double * points = <double *>malloc(2 * n_spokes * sizeof(double))
    # cdef double * angles = <double *>malloc(n_spokes * sizeof(double))
    cdef int * ordered = <int *>malloc(n_spokes * sizeof(int))
    cdef int * temp = <int *>malloc(n_spokes * sizeof(int))
    
    try:
        _argsort_spokes_around_hub(spokes, n_spokes, xy_of_spoke, xy_of_hub,
                                   ordered)

        # point = 0
        # for spoke in range(n_spokes):
        #     points[point] = xy_of_spoke[2 * spokes[spoke]]
        #     points[point + 1] = xy_of_spoke[2 * spokes[spoke] + 1]
        #     point += 2

        # _calc_spoke_angles(xy_of_hub, points, n_spokes, angles)
        # argsort(angles, n_spokes, ordered)

        for spoke in range(n_spokes):
            temp[spoke] = spokes[ordered[spoke]]

        for spoke in range(n_spokes):
            spokes[spoke] = temp[spoke]
    finally:
        # free(angles)
        free(temp)
        free(ordered)
        # free(points)


@cython.boundscheck(False)
def calc_spoke_angles(np.ndarray[double, ndim=1, mode="c"] hub,
                      np.ndarray[double, ndim=1, mode="c"] spokes,
                      np.ndarray[double, ndim=1, mode="c"] angles):
    cdef int n_spokes = spokes.shape[0] // 2
    _calc_spoke_angles(&hub[0], &spokes[0], n_spokes, &angles[0])


@cython.boundscheck(False)
def sort_spokes_around_hub(np.ndarray[long, ndim=1, mode="c"] spokes,
                           np.ndarray[double, ndim=2, mode="c"] xy_of_spoke,
                           np.ndarray[double, ndim=1, mode="c"] xy_of_hub):
    cdef int n_spokes = spokes.shape[0]

    _sort_spokes_around_hub(&spokes[0], n_spokes, &xy_of_spoke[0, 0],
                            &xy_of_hub[0])


@cython.boundscheck(False)
def argsort_points_around_hub(np.ndarray[double, ndim=2, mode="c"] points,
                              np.ndarray[double, ndim=1, mode="c"] hub,
                              np.ndarray[int, ndim=1] out):
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
    cdef double *angles = <double *>malloc(n_points * sizeof(double))

    try:
        _calc_spoke_angles(&hub[0], &points[0, 0], n_points, angles)
        argsort(angles, n_points, &out[0])
    finally:
      free(angles)

    return out


@cython.boundscheck(False)
def argsort_spokes_at_wheel(np.ndarray[long, ndim=1, mode="c"] spokes_at_wheel,
                            np.ndarray[long, ndim=1, mode="c"] offset_to_wheel,
                            np.ndarray[double, ndim=2, mode="c"] xy_of_hub,
                            np.ndarray[double, ndim=2, mode="c"] xy_of_spoke,
                            np.ndarray[int, ndim=1, mode="c"] ordered):
    cdef int n_wheels = len(offset_to_wheel) - 1
    cdef int i
    cdef int n_spokes
    cdef long * wheel
    cdef int * order

    wheel = &spokes_at_wheel[0]
    order = &ordered[0]
    for i in range(n_wheels):
        n_spokes = offset_to_wheel[i + 1] - offset_to_wheel[i]

        _argsort_spokes_around_hub(wheel, n_spokes, &xy_of_spoke[0, 0],
                                   &xy_of_hub[i, 0], order)

        order += n_spokes
        wheel += n_spokes


@cython.boundscheck(False)
def sort_spokes_at_wheel(np.ndarray[long, ndim=1, mode="c"] spokes_at_wheel,
                         np.ndarray[long, ndim=1, mode="c"] offset_to_wheel,
                         np.ndarray[double, ndim=2, mode="c"] xy_of_hub,
                         np.ndarray[double, ndim=2, mode="c"] xy_of_spoke):
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
    cdef long * wheel

    # wheel = &spokes_at_wheel[0]
    for i in range(n_wheels):
        # wheel = &spokes_at_wheel[0] + offset_to_wheel[i]
        wheel = &spokes_at_wheel[offset_to_wheel[i]]
        n_spokes = offset_to_wheel[i + 1] - offset_to_wheel[i]
        # n_spokes = spokes_per_wheel[i]
        for spoke in range(n_spokes):
            if wheel[spoke] == -1:
                n_spokes = spoke
                break

        _sort_spokes_around_hub(wheel, n_spokes, &xy_of_spoke[0, 0],
                                &xy_of_hub[i, 0])
        # wheel += n_spokes
