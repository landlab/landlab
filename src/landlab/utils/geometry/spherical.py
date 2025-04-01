"""
spherical_geometry: routines for various spherical geometry calculations,
such as cartesian-spherical coordinate conversion, arc length,
triangle area, etc.

Greg Tucker, University of Colorado Boulder, first version October 2023
"""

import numpy as np


def arc_length(p1, p2, r):
    """
    Calculate and return great-circle arc length in radians between two points p1 & p2.

    Parameters
    ----------
    p1, p2 : array of float (3, )
        (x, y, z) coordinates of each point
    r : float
        Radius
    """
    return np.arccos(np.dot(p1, p2) / r**2)


def cartesian_to_spherical(x, y, z):
    """
    Return spherical coordinates corresponding to given (x,y,z) coordinates.

    Parameters
    ----------
    x, y, z : ndarray of float
        Cartesian coordinates

    Examples
    --------
    >>> import numpy as np
    >>> x = np.array([0, 1, 1, 2, 0, -1, -1, 0])
    >>> y = np.array([1, 0, 1, 0, 1, 0, -1, 0])
    >>> z = np.array([1, 1, 0, 0, -1, 1, 0, 1])
    >>> r, p, t = cartesian_to_spherical(x, y, z)
    >>> np.round(r, 2)
    array([1.41, 1.41, 1.41, 2.  , 1.41, 1.41, 1.41, 1.  ])
    >>> np.round(p / np.pi, 2)
    array([0.5 , 0.  , 0.25, 0.  , 0.5 , 1.  , 1.25, 0.  ])
    >>> np.round(t / np.pi, 2)
    array([0.25, 0.25, 0.5 , 0.5 , 0.75, 0.25, 0.5 , 0.  ])
    """
    r = np.sqrt(x * x + y * y + z * z)
    nonzero = np.logical_or(x != 0.0, y != 0.0)
    phi = np.zeros(len(x))
    phi[nonzero] = np.arccos(x[nonzero] / np.sqrt(x[nonzero] ** 2 + y[nonzero] ** 2))
    phi[y < 0.0] = 2 * np.pi - phi[y < 0.0]
    theta = np.arccos(z / r)
    return r, phi, theta


def rotate_around_y_axis(x, z, angle):
    """
    Calculate positions of points with given (x,z) coordinates after rotation around
    y axis by given angle.

    Parameters
    ----------
    x, z : float or ndarray of float
        (x, z) coordinates of one or points
    angle : float
        angle of rotation, radians

    Examples
    --------
    >>> import numpy as np
    >>> x = np.array([1.0, 1.0, 2.0, -1.0, -1.0])
    >>> z = np.array([1.0, 0.0, 0.0, 1.0, 0.0])
    >>> ang = np.array([-np.pi / 4, -np.pi / 2, -np.pi / 2, np.pi / 4, np.pi / 2])
    >>> xr, zr = rotate_around_y_axis(x, z, ang)
    >>> np.round(np.abs(xr), 2)
    array([0., 0., 0., 0., 0.])
    >>> np.round(zr, 2)
    array([1.41, 1.  , 2.  , 1.41, 1.  ])
    """
    xry = x * np.cos(angle) + z * np.sin(angle)
    zry = -x * np.sin(angle) + z * np.cos(angle)
    return xry, zry


def rotate_around_z_axis(x, y, angle):
    """
    Calculate positions of points with given (x,y) coordinates after rotation around
    z axis by given angle.

    Parameters
    ----------
    x, y : float or ndarray of float
        (x, y) coordinates of one or points
    angle : float
        angle of rotation, radians

    Examples
    --------
    >>> import numpy as np
    >>> x = np.array([0.0, 1.0, 2.0, -1.0])
    >>> y = np.array([1.0, 0.0, 0.0, -1.0])
    >>> ang = np.array([-np.pi / 2, -np.pi / 2, np.pi / 2, 3 * np.pi / 4])
    >>> xr, yr = rotate_around_z_axis(x, y, ang)
    >>> np.round(np.abs(xr), 2)
    array([1.  , 0. , 0.  , 1.41])
    >>> np.round(yr, 2)
    array([ 0., -1.,  2., -0.])
    """
    xrz = x * np.cos(angle) - y * np.sin(angle)
    yrz = x * np.sin(angle) + y * np.cos(angle)
    return xrz, yrz


def rotate_zy(x, y, z, phi, theta):
    """
    Rotate points around z axis then y axis.

    Parameters
    ----------
    x, y, z : ndarray of float
        (x,y,z) coordinates of points to be rotated
    phi : float
        Angle for rotation about z axis, radians
    theta : float
        Angle for rotation about y axis, radians

    Examples
    --------
    >>> import numpy as np
    >>> x = np.array([1.0, -1.0])
    >>> y = np.array([1.0, 0.0])
    >>> z = np.array([0.0, 1.0])
    >>> phi = np.array([-np.pi / 4, np.pi])
    >>> theta = np.array([-np.pi / 2, -np.pi / 4])
    >>> xr, yr, zr = rotate_zy(x, y, z, phi, theta)
    >>> np.round(xr, 2)
    array([0., 0.])
    >>> np.abs(np.round(yr, 2))
    array([0., 0.])
    >>> np.round(zr, 2)
    array([1.41, 1.41])
    """
    xrz, yrz = rotate_around_z_axis(x, y, phi)
    xrzy, zry = rotate_around_y_axis(xrz, z, theta)
    return xrzy, yrz, zry


def radial_length_of_sphertri_sides(p0, p1, p2, r=1.0):
    """
    Calculate and return the radial length of the 3 sides of a spherical triangle
    with given cartesian coords.

    Parameters
    ----------
    p0, p1, p2 : ndarray of float (3, )
        (x,y,z) coordinates of spherical triangle's three points
    r : float (optional)
        Sphere radius (default 1.0)

    Examples
    --------

    A triangular patch on an icosahedron

    >>> import numpy as np
    >>> t = (1.0 + 5.0**0.5) / 2.0
    >>> p0 = np.array([-1.0, t, 0.0])
    >>> p1 = np.array([-t, 0.0, 1.0])
    >>> p2 = np.array([0.0, 1.0, t])
    >>> a, b, c = radial_length_of_sphertri_sides(p0, p1, p2, np.sqrt(np.sum(p0**2)))
    >>> (int(10 * a), int(10 * b), int(10 * c))
    (11, 11, 11)
    """
    if r != 1.0:
        p0 = p0.copy() / r
        p1 = p1.copy() / r
        p2 = p2.copy() / r
    a = np.arccos(np.dot(p1, p2))
    b = np.arccos(np.dot(p0, p2))
    c = np.arccos(np.dot(p0, p1))
    return a, b, c


def spher_angle_from_sides(s, a, b):
    """
    Calculate and return the angle between two 2 sides of a spherical
    triangle with semiperimeter (half sum of radial side lengths) s, and
    radial lengths of the adjacent sides a and b.

    Uses the half-sine rule of spherical trigonometry.

    Parameters
    ----------
    s : float
        Semiperimeter
    a, b : float
        Radial lengths of two sides, radians

    Examples
    --------

    Sides of a spherical triangle in an icosahedron:

    >>> a = 1.1071487177940904
    >>> int(10 * spher_angle_from_sides(1.5 * a, a, a) / np.pi)
    4
    """
    return 2.0 * np.arcsin(
        np.sqrt((np.sin(s - a) * np.sin(s - b)) / (np.sin(a) * np.sin(b)))
    )


def angles_of_sphertri(a, b, c):
    """
    Calculate and return angles of the spherical triangle with radial side lengths a, b, and c.

    Parameters
    ----------
    a, b, c : float
        Lengths of 3 sides of triangle, radians

    Examples
    --------
    >>> import numpy as np
    >>> a = b = c = 1.1071487177940904
    >>> A, B, C = angles_of_sphertri(a, b, c)
    >>> (int(10 * A / np.pi), int(10 * B / np.pi), int(10 * C / np.pi))
    (4, 4, 4)
    """
    s = 0.5 * (a + b + c)
    A = spher_angle_from_sides(s, b, c)
    B = spher_angle_from_sides(s, a, c)
    C = spher_angle_from_sides(s, a, b)
    return A, B, C


def area_of_sphertri(p0, p1, p2, R):
    """
    Calculate and return area of a spherical triangle with 3D coordinates p0, p1, p2,
    and radius R.

    Uses Girard's Theorem for the area of a spherical triangle (equal to the
    "spherical excess" times R^2).

    Parameters
    ----------
    p0, p1, p2 : array of float (3, )
        (x, y, z) coordinates for each of three points
    R : float
        Radius

    Examples
    --------

    Triangle on an icosahedron should be 1/20th of the total area of 4 pi R^2

    >>> import numpy as np
    >>> t = (1.0 + 5.0**0.5) / 2.0
    >>> R0 = np.sqrt(1 + t**2)
    >>> p0 = np.array([-1.0, t, 0.0]) / R0
    >>> p1 = np.array([-t, 0.0, 1.0]) / R0
    >>> p2 = np.array([0.0, 1.0, t]) / R0
    >>> 20 * area_of_sphertri(p0, p1, p2, 1.0) / np.pi
    4.0
    """
    a, b, c = radial_length_of_sphertri_sides(p0, p1, p2, R)
    A, B, C = angles_of_sphertri(a, b, c)
    return (A + B + C - np.pi) * R * R
