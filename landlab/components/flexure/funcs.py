#!/usr/bin/env python

import numpy as np
import scipy.special
from multiprocessing import Pool

_RHO_MANTLE = 3300.
_GRAVITY = 9.81
_POISSON = .25

_N_PROCS = 4


def get_flexure_parameter(h, E, n_dim):
    """
    Calculate the flexure parameter based on some physical constants. *h* is
    the Effective elastic thickness of Earth's crust (m), *E* is Young's
    Modulus, and *n_dim* is the number of spatial dimensions for which the
    flexure parameter is used. The number of dimension must be either 1, or
    2.

    === Example ===

    >>> from landlab.components.flexure.funcs import get_flexure_parameter

    >>> eet = 65000.
    >>> youngs = 7e10
    >>> alpha = get_flexure_parameter(eet, youngs, 1)
    >>> print round(alpha,3)
    120542.629

    >>> alpha = get_flexure_parameter(eet, youngs, 2)
    >>> print round(alpha,3)
    85236.51
    """
    D = E * pow(h, 3) / 12. / (1. - pow(_POISSON, 2))
    rho_m = _RHO_MANTLE

    assert(n_dim == 1 or n_dim == 2)

    if n_dim == 2:
        alpha = pow(D / (rho_m * _GRAVITY), .25)
    else:
        alpha = pow(4. * D / (rho_m * _GRAVITY), .25)

    return alpha


def subside_point_load(x, y, load, eet, youngs, i_load, j_load):
    """
    Calculate deflections on a grid due to a point load of magnitude *load*
    applied at the index *i_load*, *j_load*.

    *x* and *y* are the x and y coordinates of each node of the solution
    grid (in meters). The scalars *eet* and *youngs* define the crustal
    properties.

    === Example ===

    >>> from landlab.components.flexure.funcs import subside_point_load

    >>> EET = 65000.
    >>> Youngs = 7e10
    >>> load = 1e9
    >>> x = np.arange(0, 10000, 100.)
    >>> y = np.arange(0, 5000, 100.)
    >>> (x,y) = np.meshgrid(x, y)
    >>> dz = subside_point_load(x, y, load, EET, Youngs, 25, 50)
    >>> print round(dz.sum(), 9)
    0.002652179
    >>> print round(dz.min(), 9)
    5.29e-07
    >>> print round(dz.max(), 9)
    5.31e-07
    """
    alpha = get_flexure_parameter(eet, youngs, x.ndim)

    if x.ndim == 2:
        x_0 = x[i_load][j_load]
        y_0 = y[i_load][j_load]
        c = load / (2. * np.pi * _RHO_MANTLE * _GRAVITY * pow(alpha, 2.))
        r = np.sqrt(pow(x - x_0, 2) + pow(y - y_0, 2)) / alpha
        dz = - c * scipy.special.kei(r)
    elif x.ndim == 1:
        x_0 = x[i_load]
        c = load / (2. * alpha * _RHO_MANTLE * _GRAVITY)
        r = abs(x - x_0)/alpha
        dz = c * np.exp(-r) * (np.cos(r) + np.sin(r))

    return dz


def subside_grid(dz, x, y, load, eet, youngs, parallel=True):
    if dz.ndim == 2:
        if parallel:
            _subside_in_parallel(dz, x, y, load, eet, youngs)
        else:
            load_locs = scipy.where(load > 0)
            for (i, j) in zip(*load_locs):
                dz += subside_point_load(x, y, load[i, j], eet, youngs, i, j)
    else:
        for i in xrange(load.shape[0]):
            if abs(load[i]) > 0:
                dz += subside_point_load(x, None, load[i], eet, youngs,
                                         i, None)
    return dz


def _subside_point_load_helper(args):
    return subside_point_load(*args)


def _subside_in_parallel(dz, x, y, load, eet, youngs, n_procs=4):
    load_locs = scipy.where(load > 0)

    args = []
    for (i, j) in zip(*load_locs):
        args.append((x, y, load[i, j], eet, youngs, i, j))
    pool = Pool(processes=n_procs)

    results = pool.map(_subside_point_load_helper, args)
    for result in results:
        dz += result


def calc_def(func, args):
    return func(*args)


def worker(input, output):
    for func, args in iter(input.get, 'STOP'):
        result = calc_def(func, args)
        output.put(result)


if __name__ == '__main__':
    main()
