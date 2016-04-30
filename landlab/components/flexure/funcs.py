#!/usr/bin/env python

import numpy as np
import scipy.special
from multiprocessing import Pool

_POISSON = .25

_N_PROCS = 4


def get_flexure_parameter(h, E, n_dim, gamma_mantle=33000.):
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
    D = E * pow(h, 3) / 12. / (1. - pow(_POISSON, 2))

    assert(n_dim == 1 or n_dim == 2)

    if n_dim == 2:
        alpha = pow(D / gamma_mantle, .25)
    else:
        alpha = pow(4. * D / gamma_mantle, .25)

    return alpha


def _calculate_distances(locs, coords):
    if isinstance(locs[0], (float, int)):
        return np.sqrt(pow(coords[0] - locs[0], 2) +
                       pow(coords[1] - locs[1], 2))
    else:
        r = pow(coords[0][:, np.newaxis] - locs[0], 2)
        r += pow(coords[1][:, np.newaxis] - locs[1], 2)
        return np.sqrt(r, out=r)


def _calculate_deflections(load, locs, coords, alpha, out=None,
                           gamma_mantle=33000.):
    c = - load / (2. * np.pi * gamma_mantle * pow(alpha, 2.))
    r = _calculate_distances(locs, coords) / alpha

    if isinstance(c, (float, int)):
        return np.multiply(scipy.special.kei(r), c, out=out)
    else:
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
    params = params or dict(eet=6500., youngs=7.e10)
    eet, youngs = params['eet'], params['youngs']
    gamma_mantle = params.get('gamma_mantle', 33000.)

    assert(len(loc) in [1, 2])
    assert(len(coords) == len(loc))
    assert(len(coords[0].shape) == 1)

    if not isinstance(load, (int, float, np.ndarray)):
        load = np.array(load)

    if out is None:
        out = np.empty(coords[0].size, dtype=np.float)

    alpha = get_flexure_parameter(eet, youngs, len(loc),
                                  gamma_mantle=gamma_mantle)

    if len(loc) == 2:
        _calculate_deflections(load, loc, coords, alpha, out=out,
                               gamma_mantle=gamma_mantle)
    else:
        c = load / (2. * alpha * gamma_mantle)
        r = abs(coords[0] - loc[0]) / alpha
        out[:] = c * np.exp(-r) * (np.cos(r) + np.sin(r))

    return out


def subside_point_loads(loads, locs, coords, params=None, deflection=None,
                        n_procs=1):
    """Calculate deflection at points due multiple point loads.

    Calculate lithospheric deflections due to *loads* at coordinates
    specified by the *locs* tuple. *coords* is a tuple that gives the
    coordinates of each point where deflections are calculated; *locs* is
    positions of the applied loads. Since this function calculates the 1D
    or 2D flexure equation, *coords* and *locs* must have either one or two
    elements.

    Parameters
    ----------
    load : array_like
        Magnitude of the point loads.
    loc : tuple of (loc_x, loc_y)
        Load locations.
    coords : ndarray
        Array of points to calculate deflections at
    params : dict-like
        Physical parameters used for deflection calculation. Valid keys are
        - *eet*: Effective elastic thickness
        - *youngs*: Young's modulus
        - *gamma_mantle*: Specific weight of the mantle
    out : ndarray, optional
        Array to put deflections into.

    Returns
    -------
    out : ndarray
        Array of deflections.
    """
    params = params or dict(eet=6500., youngs=7.e10)
    eet, youngs = params['eet'], params['youngs']
    gamma_mantle = params.get('gamma_mantle', 33000.)

    if deflection is None:
        deflection = np.empty(coords[0].size, dtype=np.float)

    assert(len(coords) in [1, 2])
    assert(len(locs) == len(coords))
    assert(loads.size == locs[0].size)

    if n_procs > 1:
        _subside_in_parallel(deflection, loads, locs, coords, eet, youngs,
                             gamma_mantle, n_procs=n_procs)
    else:
        for index in loads.nonzero()[0]:
            loc = [dim.flat[index] for dim in locs]
            deflection += subside_point_load(loads.flat[index], loc,
                                             coords, eet, youngs,
                                             gamma_mantle)
    return deflection


def _subside_point_load_helper(args):
    return subside_point_load(*args)


def _subside_in_parallel(dz, loads, locs, coords, eet, youngs, gamma_mantle,
                         n_procs=4):
    args = []
    for index in loads.nonzero()[0]:
        loc = (locs[0].flat[index], locs[1].flat[index])
        args.append((loads.flat[index], loc, coords, eet, youngs,
                     gamma_mantle))

    pool = Pool(processes=n_procs)

    results = pool.map(_subside_point_load_helper, args)
    for result in results:
        try:
            dz += result
        except ValueError:
            result.shape = dz.shape
            dz += result


if __name__ == '__main__':
    import doctest
    doctest.testmod()
