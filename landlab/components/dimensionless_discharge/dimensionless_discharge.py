#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate dimensionless dischange of stream sections based on Tang (2019)
"""

import math
from landlab import Component
from landlab.utils.return_array import return_array_at_node
import numpy as np


class DimensionlessDischarge(Component):
    r"""Component that calculates dimensionless dischange of stream 
    segments.

   The dimensionless discharge model calculates the unitless discharge 
   value of streams and can be used to help determine locations of 
   debris flows. It uses an equation from Tang et al. (2019) to 
   calculate the dimensionless discharge as well as the threshold for 
   whether a debris flow will occur for a specified location. 

    Parameters
    ----------
    soil_density : float list
        Density of soil in watershed (kg/m^3)
    number_segments : int
        number of stream segments that will be used
    water_density = float
        density of water in waterhed (kg/m^3)
    C : float
        Empirically derived constant.
    N : float
        Empirically derived constant

    Examples
    --------
    >>> from landlab.components import DimensionlessDischarge
    >>> from landlab import RasterModelGrid
    >>> import random
    >>> watershed_grid = RasterModelGrid((3, 3))
    >>> flux = watershed_grid.add_ones('node', 'flux')
    >>> d50 = watershed_grid.add_ones('node', 'd50')
    >>> watershed_grid.at_node['dem_values'] = np.array([[1.1, 2, 3, 4, 2, 3, 4, 5, 3]])
    >>> dd = DimensionlessDischarge(watershed_grid)
    >>> dd.run_one_step(.5)
    >>> print(watershed_grid.at_node['dimensionless_discharge'])
    [ 0.55372743  0.55372743  0.55372743  0.55372743  0.55372743  
            0.55372743  0.55372743  0.55372743  0.55372743]

    

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
        "dem_values": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "slope of each stream segment",
        }
    }
    
    def __init__(self, grid, soil_density = 1330,  \
        water_density = 997.9, C = 12.0, N = 0.85):
        """Initialize the DimensionlessDischange.

        Parameters
        ----------
        soil_density : float, required (defaults to empty list [])
            density of soil (kg/m^3)
        C : float, optional (defaults to 12.0)
            
        N : float, defaults to 0.85
            
        stream_slope: list[float], required (defaults to empty list [])
            Slope of each segment in the stream
        stream_flux: list[list[float]], required (default to empty 2D list)
            Flux value calculated for each stream segment
        """

        super().__init__(grid)

        # Store parameters
        self._current_time = 0
        self._soil_density = soil_density
        self._C = [C]*grid.number_of_nodes
        self._N = [N]*grid.number_of_nodes
        # change DEM values into slope of a stream segment
        self._stream_slopes = grid.calc_slope_at_node(elevs="dem_values")
        self.water_density = water_density
        self.gravity = 9.8 

        #set threshold values for each segment
        dimensionless_discharge = self.grid.add_zeros('node', 
            'dimensionless_discharge')
        dimensionless_discharge_above_threshold = self.grid.add_zeros('node', 
            'dimensionless_discharge_above_threshold')
        self.grid.at_node["dimensionless_discharge_above_threshold"] = \
            np.array([[False]*self.grid.number_of_nodes])
        dimensionless_discharge_threshold_value = self.grid.add_zeros('node', 
            'dimensionless_discharge_threshold_value')
        self.grid.at_node["dimensionless_discharge_threshold_value"] = \
            self._C / (self._stream_slopes ** self._N)

    def run_one_step(self, dt):
        
        self.grid.at_node["dimensionless_discharge"] = \
            self.grid.at_node["flux"] / np.sqrt(((self._soil_density -
                self.water_density) / self.water_density) *
                 self.gravity * (self.grid.at_node["d50"] ** 3))
        
        self.grid.at_node["dimensionless_discharge_above_threshold"] = [False if \
            self.grid.at_node["dimensionless_discharge"][i] >= \
            self.grid.at_node["dimensionless_discharge_threshold_value"][i] \
            else True for i in range(self.grid.number_of_nodes)]
       

        self._current_time += dt

        