# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 11:15:21 2016

@author: RCGlade
"""

#Cubic hillslope flux component

from landlab import Component
import numpy as np
from landlab import INACTIVE_LINK, CLOSED_BOUNDARY

class CubicNonLinearDiffuser(Component):
    
    """
    hillslope evolution using a cubic formulation of hillslope flux derived from
    a taylor expansion of the nonlinear flux rule, following Ganti et al., 2012. 
    
    Parameters
    ----------
    grid: ModelGrid
            Landlab ModelGrid object
    soil_transport_coefficient: float
            Hillslope efficiency, m/yr
    slope_crit: float
            Critical slope 
        
    Examples
    --------
    >>> import numpy as np
    >>> import decimal
    >>> from landlab import RasterModelGrid 
    >>> from landlab.plot.imshow import imshow_node_grid
    >>> mg = RasterModelGrid((3, 3))
    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> initial_slope=1.0
    >>> leftmost_elev=1000.
    >>> z[:] = leftmost_elev
    >>> z[:] += (initial_slope * np.amax(mg.x_of_node)) - (initial_slope * mg.x_of_node)
    >>> mg.set_closed_boundaries_at_grid_edges(False, True, False, True)
    >>> cubicflux = CubicNonLinearDiffuser(mg, soil_transport_coefficient=0.5, slope_crit=0.1)
    >>> cubicflux.run_one_step(1.)
    >>> np.allclose(
    ...     mg.at_node['topographic__elevation'],
    ...     np.array([ 1002.,  1001.,  1000.,  1002.,  1001.,  1000.,  1002.,  1001.,
    ...     1000.]))
    True
    """

    _name = 'CubicNonLinearDiffuser'

    _input_var_names = set((
        'topographic__elevation',
    ))

    _output_var_names = set((
        'soil__flux',
        'topographic__slope',
        'topographic__elevation',
    ))

    _var_units = {
        'topographic__elevation' : 'm',
        'topographic__slope' : 'm/m',
        'soil__flux' : 'm^2/yr',
    }

    _var_mapping = {
        'topographic__elevation' : 'node',
        'topographic__slope' : 'link',
        'soil__flux' : 'link',
    }

    _var_doc = {
        'topographic__elevation':
                'elevation of the ground surface',
        'topographic__slope':
                'gradient of the ground surface',
        'soil__flux':
                'flux of soil in direction of link', 
    }

    def __init__(self, grid, soil_transport_coefficient=1., slope_crit=1.,
                 **kwds):
        
        """Initialize CubicNonLinearDiffuser.
        """

        # Store grid and parameters
        self._grid = grid
        self.k = soil_transport_coefficient
        self.slope_crit = slope_crit

        # For efficiency, store the quantity Kd / Sc^2
        self.k_over_slope_crit_sq = (soil_transport_coefficient
                                     / (slope_crit * slope_crit))

        # Create fields:
        #
        # elevation
        if 'topographic__elevation' in self.grid.at_node:
            self.elev = self.grid.at_node['topographic__elevation']
        else:
            self.elev = self.grid.add_zeros('node', 'topographic__elevation')

        # slope gradient
        if 'topographic__slope' in self.grid.at_link:
            self.slope = self.grid.at_link['topographic__slope']
        else:
            self.slope = self.grid.add_zeros('link', 'topographic__slope')

        # soil flux
        if 'soil__flux' in self.grid.at_link:
            self.flux = self.grid.at_link['soil__flux']
        else:
            self.flux = self.grid.add_zeros('link', 'soil__flux')


    def soilflux(self, dt):
        """Calculate soil flux for a time period 'dt'.
        """

        # Calculate gradients
        self.slope[:] = self.grid.calc_grad_at_link(self.elev)
        self.slope[self.grid.status_at_link == INACTIVE_LINK] = 0.

        # Calculate flux
        self.flux[:] = -(self.k * self.slope
                         + (self.k_over_slope_crit_sq 
                            * np.power(self.slope, 3)))

        # Calculate flux divergence
        dqdx = self.grid.calc_flux_div_at_node(self.flux)
        dqdx[self.grid.status_at_node == CLOSED_BOUNDARY] = 0.

        # Update topography
        self.elev -= dqdx


    def run_one_step(self, dt, **kwds):
        """
        Advance cubic soil flux component by one time step of size dt.

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        """
        self.soilflux(dt, **kwds)    
