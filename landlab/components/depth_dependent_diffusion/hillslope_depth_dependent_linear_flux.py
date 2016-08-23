# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 08:32:48 2016

@author: RCGlade
"""

from landlab import Component
import numpy as np
from landlab import INACTIVE_LINK, CLOSED_BOUNDARY

class DepthDependentDiffuser(Component):
    
    """
    This component implements a depth and slope dependent linear diffusion rule
    in the style of Johnstone and Hilley (2014). 

        
    Parameters
    ----------
    grid: ModelGrid
        Landlab ModelGrid object
    k: float
        Hillslope efficiency, m/yr
    hstar: float
        characteristic transport soil depth, m
    


    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import ExponentialWeathering, DepthDependentDiffuser
    >>> mg = RasterModelGrid((5, 5))
    >>> soilTh = mg.add_zeros('node', 'soil__depth')
    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> BRz = mg.add_zeros('node', 'bedrock__elevation')
    >>> expweath = ExponentialWeathering(mg)
    >>> DDdiff = DepthDependentDiffuser(mg)
    >>> expweath.exponentialweather()
    >>> np.allclose(mg.at_node['weathering__rate'], 1.)
    True
    >>> DDdiff.soilflux(2.)
    >>> np.allclose(mg.at_node['topographic__elevation'], 0.)
    True
    >>> np.allclose(mg.at_node['bedrock__elevation'], -2.)
    True
    >>> np.allclose(mg.at_node['soil__depth'], 2.)
    True

    Now with a slope:
    >>> mg = RasterModelGrid((5, 5))
    >>> soilTh = mg.add_zeros('node', 'soil__depth')
    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> BRz = mg.add_zeros('node', 'bedrock__elevation')
    >>> z += mg.node_x.copy()
    >>> BRz += mg.node_x/2.
    >>> soilTh[:] = z - BRz
    >>> expweath = ExponentialWeathering(mg)
    >>> DDdiff = DepthDependentDiffuser(mg)
    >>> expweath.exponentialweather()
    >>> mynodes = mg.nodes[2, :]
    >>> np.allclose(
    ...     mg.at_node['weathering__rate'][mynodes],
    ...     np.array([ 1., 0.60653066, 0.36787944, 0.22313016, 0.13533528]))
    True
    >>> DDdiff.soilflux(2.)
    >>> np.allclose(
    ...     mg.at_node['topographic__elevation'][mynodes],
    ...     np.array([0., 1.47730244, 2.28949856, 3.17558975, 4.]))
    True
    >>> np.allclose(mg.at_node['bedrock__elevation'][mynodes],
    ...     np.array([-2., -0.71306132, 0.26424112, 1.05373968, 1.72932943]))
    True
    >>> np.allclose(mg.at_node['soil__depth'], z - BRz)
    True

    """

    _name = 'DepthDependentDiffuser'

    _input_var_names = (
        'topographic__elevation',
        'soil__depth',
        'weathering__rate',
    )
    
    _output_var_names = (
        'soil__flux',
        'topographic__slope',
        'topographic__elevation',
        'bedrock__elevation',
        'soil__depth',   
    )
        
    _var_units = {
        'topographic__elevation' : 'm',
        'topographic__slope' : 'm/m',
        'soil__depth' : 'm',
        'soil__flux' : 'm^2/yr',
        'weathering__rate' : 'm/yr',
        'bedrock__elevation' : 'm',
    }
    
    _var_mapping = {
        'topographic__elevation' : 'node',
        'topographic__slope' : 'link',
        'soil__depth' : 'node',
        'soil__flux' : 'link',
        'weathering__rate' : 'node',
        'bedrock__elevation' :'node',
    }
        
    _var_doc = {
        'topographic__elevation':
                'elevation of the ground surface',
        'topographic__slope':
                'gradient of the ground surface',
        'soil__depth':
                'depth of soil/weather bedrock',
        'soil__flux':
                'flux of soil in direction of link', 
        'weathering__rate':
                'rate of soil production at nodes',
        'bedrock__elevation':
                'elevation of the bedrock surface',
    }

    def __init__(self,grid,k=1,hstar=1,**kwds):
        
       
        
        #Store grid and parameters
        self._grid=grid
        self.k=k
        self.hstar=hstar
        
        
        
        #create fields
        #elevation
        if 'topographic__elevation' in self.grid.at_node:
            self.elev=self.grid.at_node['topographic__elevation']
        else:
            self.elev=self.grid.add_zeros('node','topographic__elevation')
        
        #slope
        if 'topographic__slope' in self.grid.at_link:
            self.slope=self.grid.at_link['topographic__slope']
        else:
            self.slope=self.grid.add_zeros('link','topographic__slope')
        
        #soil depth
        if 'soil__depth' in self.grid.at_node:
            self.depth=self.grid.at_node['soil__depth']
        else:
            self.depth=self.grid.add_zeros('node','soil__depth')
        
        #soil flux
        if 'soil__flux' in self.grid.at_link:
            self.flux=self.grid.at_link['soil__flux']
        else:
            self.flux=self.grid.add_zeros('link','soil__flux')
            
        #weathering rate
        if 'weathering__rate' in self.grid.at_node:
            self.weather=self.grid.at_node['weathering__rate']
        else:
            self.weather=self.grid.add_zeros('node','weathering__rate')
            
        #bedrock elevation
        if 'bedrock__elevation' in self.grid.at_node:
            self.bedrock=self.grid.at_node['bedrock__elevation']
        else:
            self.bedrock=self.grid.add_zeros('node','bedrock__elevation')

        self._active_nodes = self.grid.status_at_node != CLOSED_BOUNDARY
        
        
    
    def soilflux(self, dt):
        """Calculate soil flux for a time period 'dt'.

        Parameters:
        ----------

        dt: float (time)
            The imposed timestep.

        """


        #update soil depth
        self.grid.at_node['soil__depth'][:] = self.grid.at_node['topographic__elevation'] - self.grid.at_node['bedrock__elevation']
        
        #Calculate soil depth at links.
        H_link=self.grid.map_value_at_max_node_to_link('topographic__elevation','soil__depth')
        
        #Calculate gradients
        slope=self.grid.calc_grad_at_link(self.elev)
        slope[self.grid.status_at_link == INACTIVE_LINK] = 0.
        
        #Calculate flux
        self.flux[:]=-self.k*slope*self.hstar*(1.-np.exp(-H_link/self.hstar))
        
        #Calculate flux divergence
        dqdx=self.grid.calc_flux_div_at_node(self.flux)
        dqdx[self.grid.status_at_node == CLOSED_BOUNDARY] = 0.
        
        #Calculate change in soil depth
        dhdt=self.weather-dqdx
        
        #Calculate soil depth at nodes
        self.depth[self._active_nodes] += dhdt[self._active_nodes]*dt
        
        #prevent negative soil thickness
        self.depth[self.depth < 0.0]=0.0

        #Calculate bedrock elevation
        self.bedrock[self._active_nodes]-=self.weather[self._active_nodes]*dt

        #Update topography
        self.elev[self._active_nodes]=self.depth[self._active_nodes]+self.bedrock[self._active_nodes]


