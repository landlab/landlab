#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate dimensionless dischange of stream sections based on Tang (2019)
"""

import math
from landlab import Component


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
        "dimensionless_discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "none",
            "mapping": "node",
            "doc": "Dimensionless discharge value for a stream segment.",
        },
        "dimensionless_discharge_above_threshold": {
            "dtype": bool,
            "intent": "out",
            "optional": False,
            "units": "none",
            "mapping": "node",
            "doc": "Dimensionless discharge value for a stream segment.",
        }
    }
    
    def __init__(self, grid, soil_density=1330, d50=[], C=12.0, N=0.85, stream_slopes=[], flux=[[]]):
        """Initialize the DimensionlessDischange.

        Parameters
        ----------
        soil_density : float, required (defaults to empty list [])
            density of soil (kg/m^3)
        d50 : list[float], required (defaults to 1 hour)
            soil partical size (m)
        C : float, optional (defaults to 12.0)
            
        N : float, defaults to 0.85
            
        stream_slope: list[float], required (defaults to empty list [])
            Slope of each segment in the stream
        stream_flux: list[list[float]], required (default to empty 2D list)
            Flux value calculated for each stream segment
        """

        super().__init__(grid)

        # Store parameters and do unit conversion
        self._current_time = 0
        self._iteration = 0
        self._soil_density = soil_density
        self._d50 = d50
        self._C = C
        self._N = N
        self._stream_slopes = stream_slopes
        self._flux = flux
        self._dimensionless_discharge = self.grid.at_node["dimensionless_discharge"]
        self._dimensionless_discharge_above_threshold = self.grid.at_node["dimensionless_discharge_above_threshold"]

        #set threshold values for each segment
        self._dimensionless_discharge_threshold_value = [[ 0 for i in range(len(self._flux[1]))] for in range(len(self._flux))]
        for i in range(len(self._dimensionless_discharge_threshold_value)):
            for j in range(len(i)):
                self._dimensionless_discharge_threshold_value[i][j] = C/(math.tan(self._stream_slopes[i])**N)

    def run_one_step(self, dt):
        for i in range(self._iteration):
            self._dimensionless_discharge[self.iteration][i] = self._flux[0][i]/math.sqrt(((self._soil_density
        -_WATER_DENSITY)/_WATER_DENSITY)*_GRAVITY*(self._d50[i]**3))
            if self._dimensionless_discharge[self.iteration][i] >= self._dimensionless_discharge_threshold_value[i]:
                self._dimensionless_discharge_above_threshold[self.iteration][i] = 1

        self.iteration+=1
        self._current_time += dt

        