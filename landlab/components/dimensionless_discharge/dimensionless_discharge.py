#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate dimensionless dischange of stream sections based on Tang (2019)
"""

import numpy as np
from numpy.core.records import array


_WATER_DENSITY = 997.9 # kg/m^3
_GRAVITY = 9.8  # m/s^2


class DimensionlessDischange(Component):
    r"""Component that calculates dimensionless dischange of stream segments.

    The DimensionlessDischarge component calculates the dimensionless discharge 
    value for every stream segment in a watershed. It also computes a threshold
    which is used predict whether or not a segment would have a debris flow. If 
    the dimensionless discharge value is greater than the threshold value, we 
    predict that a debris flow would occure in this area.

    Parameters
    ----------
    soil_density : float list
        Density of soil in watershed (kg/m^3)
    flux : float list
        flux in each stream segment (m/s)
    d50 : float list
        soil grain size (m)
    number_segments : int
        number of stream segments that will be used
    C : float
        
    N : float
       
    stream_segment_slopes : float list
        slope of each stream segment

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> rg = RasterModelGrid((4, 5), 10.0)
    >>> kw = KinwaveOverlandFlowModel(rg)
    >>> kw.vel_coef
    100.0
    >>> rg.at_node['surface_water__depth']
    array([ 0.,  0.,  0.,  0.,  0.,
            0.,  0.,  0.,  0.,  0.,
            0.,  0.,  0.,  0.,  0.,
            0.,  0.,  0.,  0.,  0.])

    

    References
    ----------
   
    """

    _name = "DimensionlessDischangeModel"

    _unit_agnostic = False

    _info = {
    "soil_density": {
        "dtype": list[float],
        "intent": "in",
        "optional": False,
        "units": "kg/m^2",
        "doc": "density of soil in stream segment",
    },
     "water__velocity": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/s",
            "mapping": "link",
            "doc": "flow velocity component in the direction of the link",
        },
    }

    

    def __init__():
        pass