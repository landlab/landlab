#!/usr/bin/env python

import numpy as np
import scipy.special

_POISSON = 0.25

_N_PROCS = 4


def get_flexure_parameter(h, E, n_dim, gamma_mantle=33000.0):
    """
    Calculate the flexure parameter based on some physical constants. *h* is
    the Effective elastic thickness of Earth's crust (m), *E* is Young's
    Modulus, and *n_dim* is the number of spatial dimensions for which the
    flexure parameter is used. The number of dimension must be either 1, or
    2.

    Examples
    --------
    >>> from __future__ import print_function
    >>> from landlab.components.flexure import get_flexure_parameter

    >>> eet = 65000.
    >>> youngs = 7e10
    >>> alpha = get_flexure_parameter(eet, youngs, 1)
    >>> print('%.3f' % round(alpha, 3))
    119965.926

    >>> alpha = get_flexure_parameter(eet, youngs, 2)
    >>> print('%.2f' % alpha)
    84828.72
    """
    D = E * pow(h, 3) / 12.0 / (1.0 - pow(_POISSON, 2))

    if n_dim not in (1, 2):
        raise ValueError("n_dim must be either 1 or 2")

    if n_dim == 2:
        alpha = pow(D / gamma_mantle, 0.25)
    else:
        alpha = pow(4.0 * D / gamma_mantle, 0.25)

    return alpha


def _calculate_distances(locs, coords):
    r = pow(coords[0][:, np.newaxis] - locs[0], 2)
    r += pow(coords[1][:, np.newaxis] - locs[1], 2)
    return np.sqrt(r, out=r)


def _calculate_deflections(load, locs, coords, alpha, out=None, gamma_mantle=33000.0):
    c = -load / (2.0 * np.pi * gamma_mantle * pow(alpha, 2.0))
    r = _calculate_distances(locs, coords) / alpha

    scipy.special.kei(r, out=r)
    np.multiply(r, c[np.newaxis, :], out=r)
    return np.sum(r, axis=1, out=out)


def subside_point_load(load, loc, coords, params=None, out=None):
    """Calculate deflection at points due a point load.

    Calculate deflections on a grid, defined by the points in the *coords*
    tuple, due to a point load of magnitude *load* applied at *loc*.

    *x* and *y* are the x and y coordinates of each node of the solution
    grid (in meters). The scalars *eet* and *youngs* define the crustal
    properties.

    Parameters
    ----------
    load : float
        Magnitude of the point load.
    loc : float or tuple
        Location of the load as either a scalar or as (*x*, *y*)
    coords : ndarray
        Array of points to calculate deflections at
    params : dict-like
        Physical parameters used for deflection calculation. Valid keys are
        - *eet*: Effective elastic thickness
        - *youngs*: Young's modulus
    out : ndarray, optional
        Array to put deflections into.

    Returns
    -------
    out : ndarray
        Array of deflections.

    Examples
    --------

    >>> from landlab.components.flexure import subside_point_load

    >>> params = dict(eet=65000., youngs=7e10)
    >>> load = 1e9

    Define a unifrom rectilinear grid.

    >>> x = np.arange(0, 10000, 100.)
    >>> y = np.arange(0, 5000, 100.)
    >>> (x, y) = np.meshgrid(x, y)
    >>> x.shape = (x.size, )
    >>> y.shape = (y.size, )

    Calculate deflections due to a load applied at position (5000., 2500.).

    >>> import six
    >>> x = np.arange(0, 10000, 1000.)
    >>> y = np.arange(0, 5000, 1000.)
    >>> (x, y) = np.meshgrid(x, y)
    >>> x.shape = (x.size, )
    >>> y.shape = (y.size, )
    >>> dz = subside_point_load(load, (5000., 2500.), (x, y), params=params)
    >>> print('%.5g' % round(dz.sum(), 9))
    2.6267e-05
    >>> six.print_(round(dz.min(), 9))
    5.24e-07
    >>> six.print_(round(dz.max(), 9))
    5.26e-07

    >>> dz = subside_point_load((1e9, 1e9), ((5000., 5000.), (2500., 2500.)),
    ...                         (x, y), params=params)
    >>> six.print_(round(dz.min(), 9) / 2.)
    5.235e-07
    >>> six.print_(round(dz.max(), 9) / 2.)
    5.265e-07
    """
    params = params or dict(eet=6500.0, youngs=7.0e10)
    eet, youngs = params["eet"], params["youngs"]
    gamma_mantle = params.get("gamma_mantle", 33000.0)

    load = np.asarray(load).reshape((-1,))
    loc = np.asarray(loc).reshape((-1, len(load)))
    coords = np.asarray(coords)
    if coords.ndim == 1:
        coords = np.expand_dims(coords, axis=0)

    n_dim = len(loc)
    if n_dim not in (1, 2):
        raise ValueError("number of dimension must be 1 or 2")
    if len(coords) != n_dim:
        raise ValueError("number of dimensions in coordinates doesn't match loc")

    if out is None:
        out = np.empty(coords[0].size, dtype=np.float)

    alpha = get_flexure_parameter(eet, youngs, n_dim, gamma_mantle=gamma_mantle)

    if n_dim == 2:
        _calculate_deflections(
            load, loc, coords, alpha, out=out, gamma_mantle=gamma_mantle
        )
    else:
        x, x0 = np.meshgrid(loc[0], coords[0])
        c = load / (2.0 * alpha * gamma_mantle)
        r = abs(x - x0) / alpha
        out[:] = (c * np.exp(-r) * (np.cos(r) + np.sin(r))).sum(axis=1)

    return out
