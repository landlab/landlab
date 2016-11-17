#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 08:32:48 2016

@author: RCGlade
"""

from landlab import Component, CLOSED_BOUNDARY
import numpy as np



class ExponentialWeathering(Component):
    
    """
    This component implements exponential weathering of bedrock on hillslopes. 
    Uses exponential soil production function in the style of Heimsath 1997.
        
    Parameters
    ----------
    grid: ModelGrid
        Landlab ModelGrid object
    wstar: float
	characteristic weathering depth
    wnot: float
	maximum weathering rate for bare bedrock 
        
 

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import ExponentialWeathering
    >>> mg = RasterModelGrid((5, 5))
    >>> soilz = mg.add_zeros('node', 'soil__depth')
    >>> soilrate = mg.add_ones('node', 'weathering__rate')
    >>> expw = ExponentialWeathering(mg)
    >>> expw.exponentialweather()
    >>> np.allclose(mg.at_node['weathering__rate'], 1.)
    True
    """

    _name = 'ExponentialWeathering'

    _input_var_names = (
        'soil__depth',
    )
    
    _output_var_names = (
        'weathering__rate',
    )
        
    _var_units = {
        'soil__depth' : 'm',
        'weathering__rate' : 'm/yr',
    }
    
    _var_mapping = {
        'soil__depth' : 'node',
        'weathering__rate' : 'node',
        
    }
        
    _var_doc = {
        'soil__depth':
                'depth of soil/weather bedrock',
        'weathering__rate':
                'rate of soil production at nodes',

    }

    def __init__(self, grid, wnot=1, wstar=1, **kwds):
        
                
        #Store grid and parameters
        self._grid = grid
        self.wstar = wstar
        self.wnot = wnot
        
        #create fields
    
        
        #soil depth
        if 'soil__depth' in grid.at_node:
            self.depth = grid.at_node['soil__depth']
        else:
            self.depth = grid.add_zeros('node','soil__depth')

            
        #weathering rate
        if 'weathering__rate' in grid.at_node:
            self.weather = grid.at_node['weathering__rate']
        else:
            self.weather = grid.add_zeros('node','weathering__rate')

        self._active_nodes = self.grid.status_at_node != CLOSED_BOUNDARY
        
    
    def exponentialweather(self, current_time=0.0, **kwds):
        """Calculate soil flux for a time period 'dt'.
        """
        
        #weather
        self.weather[self._active_nodes] = (self.wnot*np.exp(-self.depth[self._active_nodes]/self.wstar))

        
        
