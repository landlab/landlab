#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 18:30:44 2018

@author: barnhark
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create a block of rock with different properties. 

@author: barnhark
"""
import numpy as np
from landlab.layers import RockBlock


class LayeredRockBlock(RockBlock):

    """Create LayeredRockBlock

   
    
    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.layers import RockBlock
    >>> mg = RasterModelGrid(3, 3)
    >>> depths = []
    >>> ids = []
    >>> attrs = {'K_sp': {1: 0.001
    ...                    2: 0.0001},
    ...          'D': {1: 0.01, 
    ...                2: 0.001}}
    >>> rb = RockBlock(grid, depths, ids, attrs)
    
    
    """

    _name = 'LayeredRockBlock'

    _cite_as = """ """

    

    def __init__(self, grid, z0s, ids, attrs, x0=0, y0=0,  function=lambda x, y: 0):
        """Initialize the flexure component.

        Parameters
        ----------
        grid : RasterModelGrid
            A grid.
        x
        y
        zs
        ids
        
        
        """
      
        self._grid = grid
        
        if np.asarray(z0s).size != np.asarray(ids).size:
            msg = 'size of zs and ids must be the same'
            raise ValueError(msg)
            
        if np.any(np.diff(z0s) < 0):
            msg = 'bad order'
            raise ValueError(msg)
            
        z_surf = function(self._grid.x_of_node - y0, self._grid.y_of_node - y0)
        
        layer_thicknesses = []
        layer_ids = []
        
        num_layers = np.asarray(z0s).size
        
        last_layer_elev = np.zeros(self._grid.number_of_nodes)
        
        # create layers (here listed from the top to the bottom.)
        for i in range(num_layers):
            
            layer_depth = z_surf + z0s[i]
            layer_depth[layer_depth<0] = 0
            
            layer_thickness = layer_depth.copy() - last_layer_elev.copy()
            
            last_layer_elev = layer_depth.copy()
            
            layer_thicknesses.append(layer_thickness)
            layer_ids.append(ids[i] * np.ones(z_surf.size))

        super(LayeredRockBlock, self).__init__(grid, layer_thicknesses, layer_ids, attrs)
       