#! /usr/bin/env python
import numpy as np


def quad(x,a,b,c):
   y = a*x*x + b*x + c
   return y;

def find_shoreline2(x, z, sea_level=0., kind='cubic'):
    """Find the shoreline of a profile.

    Parameters
    ----------
    x : array of float
        X-positions of profile.
    z : array of float
        Elevations along the profile.
    sea_level : float, optional
   
        Elevation of sea level.
    kind : str, optional
        Interpolation method used to find shoreline. Values are the same
        as those used for `scipy.interpolate.interp1d`. Default is `'cubic'`.

    Returns
    -------
    float
        X-position of the shoreline.

    Examples
    --------
    >>> from landlab.components.submarine_diffusion.utils import find_shoreline
    >>> import numpy as np

    Create a linearly-dipping profile.

    >>> x = np.arange(10.)
    >>> z = - x + 5.

    Find the shoreline.

    >>> find_shoreline(x, z, kind='linear')
    5.0
    >>> find_shoreline(x, z, sea_level=.25, kind='linear')
    4.75

    If sea level is higher/lower than the max/min elevation, return
    the first/last *x* value.

    >>> find_shoreline(x, z, sea_level=100., kind='linear')
    0.0
    >>> find_shoreline(x, z, sea_level=-100., kind='linear')
    9.0
    """
    from scipy.interpolate import interp1d
    from scipy.optimize import bisect

    (below_water, ) = np.where(z < sea_level)

    if len(below_water) == 0:
        x_of_shoreline = x[-1]
    elif len(below_water) == len(x):
        x_of_shoreline = x[0]
    else:
        index_at_shore = below_water[0]
        l1 = index_at_shore - 1
        l2 = index_at_shore - 2
        l3 = index_at_shore - 3
        a = ((z[l2]-z[l3])*x[l1] + (z[l1]-z[l2])*x[l3] + (z[l3]-z[l1])*x[l2])/ (
                (x[l3]-x[l2])*x[l1]*x[l1] + (-x[l3]*x[l3] + x[l2]*x[l2])*x[l1] + (x[l2]*x[l3]*x[l3]-x[l2]*x[l2]*x[l3]))
        b = ((-z[l2]+z[l3])*x[l1]*x[l1] + (-z[l1]*x[l3]+z[l1]*x[l2])*x[l1] + (
                                               (-z[l1]+z[l2])*x[l3]*x[l3] + z[l1]*x[l2]*x[l3] - z[l3]*x[l2]*x[l2]))/ (
                (x[l3]-x[l2])*x[l1]*x[l1] + (-x[l3]*x[l3] + x[l2]*x[l2])*x[l1] + (x[l2]*x[l3]*x[l3]-x[l2]*x[l2]*x[l3]))
        c = ((z[l2]*x[l3]-z[l3]*x[l2])*x[l1]*x[l1] + ((z[l1]-z[l2])*x[l3]*x[l3] - (
                                                      z[l1]*x[l2]*x[l3])+ z[l3]*x[l2]*x[l2])*x[l1])/ (
                (x[l3]-x[l2])*x[l1]*x[l1] + (-x[l3]*x[l3] + x[l2]*x[l2])*x[l1] + (x[l2]*x[l3]*x[l3]-x[l2]*x[l2]*x[l3]))
        yy = z
        yy[index_at_shore: len(yy)]   = quad(x[index_at_shore: len(yy)],a,b,c)
        func = interp1d(x, yy - sea_level, kind='quadratic')
        x_of_shoreline1 = bisect(func, x[index_at_shore - 1], x[index_at_shore])
        print (x_of_shoreline1)
        l1 = index_at_shore 
        l2 = index_at_shore + 1
        l3 = index_at_shore + 2
        a = ((z[l2]-z[l3])*x[l1] + (z[l1]-z[l2])*x[l3] + (z[l3]-z[l1])*x[l2])/ (
                (x[l3]-x[l2])*x[l1]*x[l1] + (-x[l3]*x[l3] + x[l2]*x[l2])*x[l1] + (x[l2]*x[l3]*x[l3]-x[l2]*x[l2]*x[l3]))
        b = ((-z[l2]+z[l3])*x[l1]*x[l1] + (-z[l1]*x[l3]+z[l1]*x[l2])*x[l1] + (
                                               (-z[l1]+z[l2])*x[l3]*x[l3] + z[l1]*x[l2]*x[l3] - z[l3]*x[l2]*x[l2]))/ (
                (x[l3]-x[l2])*x[l1]*x[l1] + (-x[l3]*x[l3] + x[l2]*x[l2])*x[l1] + (x[l2]*x[l3]*x[l3]-x[l2]*x[l2]*x[l3]))
        c = ((z[l2]*x[l3]-z[l3]*x[l2])*x[l1]*x[l1] + ((z[l1]-z[l2])*x[l3]*x[l3] - (
                                                      z[l1]*x[l2]*x[l3])+ z[l3]*x[l2]*x[l2])*x[l1])/ (
                (x[l3]-x[l2])*x[l1]*x[l1] + (-x[l3]*x[l3] + x[l2]*x[l2])*x[l1] + (x[l2]*x[l3]*x[l3]-x[l2]*x[l2]*x[l3]))
        yy = z
        yy[1:index_at_shore - 1]   = quad(x[1:index_at_shore - 1],a,b,c)
        x_of_shoreline2 = bisect(func, x[index_at_shore - 1], x[index_at_shore])
        print x_of_shoreline2
        shore = (x_of_shoreline1 + x_of_shoreline2)/2
    return shore