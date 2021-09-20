#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate dimensionless dischange of stream sections based on Tang (2019)
"""

import math
from landlab import Component


_WATER_DENSITY = 997.9 # kg/m^3
_GRAVITY = 9.8  # m/s^2


class DimensionlessDischarge(Component):
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
    number_segments : int
        number of stream segments that will be used
    C : float
        
    N : float

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
        },
        "flux": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/s",
            "mapping": "node",
            "doc": "flux in each stream segment",
        },
        "d50": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "soil grain size",
        },
        "stream_slopes": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "None",
            "mapping": "node",
            "doc": "slope of each stream segment",
        }
    }
    
    def __init__(self, grid, soil_density=1330, C=12.0, N=0.85):
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
        self._soil_density = soil_density
        self._C = C
        self._N = N
        self._stream_slopes = self.grid.at_node["stream_slopes"]

        #set threshold values for each segment
        dimensionless_discharge = self.grid.add_zeros('node', 'dimensionless_discharge')
        dimensionless_discharge_above_threshold = self.grid.add_zeros('node', 'dimensionless_discharge_above_threshold')
        dimensionless_discharge_threshold_value =self.grid.add_zeros('node', 'dimensionless_discharge_threshold_value')
        for i in range(len(self.grid.at_node["dimensionless_discharge_threshold_value"])):
            self.grid.at_node["dimensionless_discharge_threshold_value"][i] = C/(math.tan(self._stream_slopes[i])**N)

    def run_one_step(self, dt):
        for i in range(len(self.grid.at_node["dimensionless_discharge"])):
            self.grid.at_node["dimensionless_discharge"][i] = self.grid.at_node["flux"][i]/math.sqrt(((self._soil_density
        -_WATER_DENSITY)/_WATER_DENSITY)*_GRAVITY*(self.grid.at_node["d50"][i]**3))
            if self.grid.at_node["dimensionless_discharge"][i] >= self.grid.at_node["dimensionless_discharge_threshold_value"][i]:
                self.grid.at_node["dimensionless_discharge_above_threshold"][i] = 1
            else: 
                self.grid.at_node["dimensionless_discharge_above_threshold"][i] = 0

        self._current_time += dt

        