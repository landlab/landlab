#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 10:13:38 2017

@author: margauxmouchene
"""

# Hillslope diffusion from Carretier et al 2016, after Davy and Lague 2012
# L>=dx to keep seposition flux D dx smaller than incoming sediment flux qs
# Works on regular grid only (dx in definition of transport length L)


from landlab import Component
import numpy as np
from landlab import INACTIVE_LINK, CLOSED_BOUNDARY
from landlab.components import FlowDirectorSteepest
# from landlab.components.flow_director.flow_direction_DN import grid_flow_directions


class TransportLengthHillslopeDiffuser(Component):

    """
    description...

    Construction::
        TransportLengthHillslopeDiffuser(grid, ...)


    Parameters
    ----------
    grid: ModelGrid
            Landlab ModelGrid object
    ...
    Returns
    ----------
    ...

    Examples
    --------
    >>> import numpy as np
    >>>
    """

# TO DO ###########################################
    _name = 'TransportLengthHillslopeDiffuser'

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
        ###############################################


    def __init__(self, grid, erodibility, slope_crit=1.,
                 **kwds):

        """Initialize Diffuser.
        """

        # Store grid and parameters
        self._grid = grid
        self.k = erodibility
        self.slope_crit = slope_crit

        # Create fields:

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

        # deposition (defined at nodes)
        if 'deposition' in self.grid.at_node:
            self.depo = self.grid.at_node['deposition']
        else:
            self.depo = self.grid.add_zeros('node', 'deposition')

        # transferred sediments (not deposited) (defined at nodes)
        if 'transfer' in self.grid.at_node:
            self.trans = self.grid.at_node['transfer']
        else:
            self.trans = self.grid.add_zeros('node', 'transfer')

        # transport length (defined at nodes)
        if 'transport_length' in self.grid.at_node:
            self.L = self.grid.at_node['transport_length']
        else:
            self.L = self.grid.add_zeros('node', 'transport_length')

        # flux in (defined at nodes)
        if 'flux_in' in self.grid.at_node:
            self.flux_in = self.grid.at_node['flux_in']
        else:
            self.flux_in = self.grid.add_zeros('node', 'flux_in')

        # flux out (defined at nodes)
        if 'flux_out' in self.grid.at_node:
            self.flux_out = self.grid.at_node['flux_out']
        else:
            self.flux_out = self.grid.add_zeros('node', 'flux_out')

        # erosion (defined at nodes)
        if 'erosion' in self.grid.at_node:
            self.erosion = self.grid.at_node['erosion']
        else:
            self.erosion = self.grid.add_zeros('node', 'erosion')

#        self.flux_in = np.zeros(self.grid.number_of_nodes)
#        self.flux_out = np.zeros(self.grid.number_of_nodes)
#        self.erosion = np.zeros(self.grid.number_of_nodes)
#        self.depo = np.zeros(self.grid.number_of_nodes)
#        self.trans = np.zeros(self.grid.number_of_nodes)

    def soilflux(self, dt):
        """Calculate soil flux for a time period 'dt'.
        """

        # Downstream steepest slope at node:
        self.steepest = self.grid.at_node['topographic__steepest_slope'] # / 100
        # On each node, node ID of downstream receiver node
        # (on node (i), ID of node that receives flow from node (i)):
        self.receiver = self.grid.at_node['flow__receiver_node']
#
#        from matplotlib.pyplot import figure
#        from landlab.plot import imshow_grid
#        figure(7)
#        im = imshow_grid(self.grid, 'topographic__steepest_slope', plot_name='Steepest Slope')
#        figure(8)
#        im = imshow_grid(self.grid, 'flow__receiver_node', plot_name='Receiver node')


        for i in self.grid.core_nodes:
            # Calculate influx on node i = outflux of nodes whose receiver is i
            self.flux_in[self.receiver[i]] += self.flux_out[i]
            
            # Calculate transport length
            if self.steepest[i] >= self.slope_crit:
                self.L[i] = 1000000000.   # 'Infinite' length, for stability
            else:
                self.L[i] = (self.grid.dx)/(1-(np.power(((self.steepest[i])/self.slope_crit), 2)))

        cores = self.grid.core_nodes
        # Calculate deposition on node
        self.depo[cores] = self.flux_in[cores] / self.L[cores]
        # Calculate transfer over node
        self.trans[cores] = self.flux_in[cores] - self.depo[cores]

            # Calculate erosion on node (positive value)
        for i in self.grid.core_nodes:
            if self.steepest[i] > self.slope_crit:
                self.steepest[i] = self.slope_crit
            else:
                pass
            self.erosion[:] = self.k * self.steepest

            # Calculate outflux
        self.flux_out[:] = self.erosion + self.trans     # Flux out of node

#        from matplotlib.pyplot import figure
#        from landlab.plot import imshow_grid
#        figure()
#        im = imshow_grid(self.grid, self.flux_in, plot_name='Flux in')
#        figure()
#        im = imshow_grid(self.grid, self.L, plot_name='L')
#        figure()
#        im = imshow_grid(self.grid, self.depo, plot_name='Depo')
#        figure()
#        im = imshow_grid(self.grid, self.trans, plot_name='Trans')
#        figure()
#        im = imshow_grid(self.grid, self.erosion, plot_name='Erosion')
#        figure()
#        im = imshow_grid(self.grid, self.flux_out, plot_name='Flux out')

        # Update topography
        cores = self.grid.core_nodes
        self.elev[cores] += (-self.erosion[cores] + self.depo[cores]) * dt 
        # divide by area? / (pow(self.grid.dx, 2)))

        # reset erosion, depo, trans to 0    
        self.erosion[:] = 0.
        self.depo[:] = 0.
        self.trans[:] = 0.
        self.flux_in[:] = 0.


    #self.flux_in[self.grid.boundary_nodes[:]] = 0.
    # flux in boundary nodes is flushed out of grid


    def run_one_step(self, dt, **kwds):
        """
        Advance transport length model soil flux component 
        by one time step of size dt.

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        """
        self.soilflux(dt, **kwds)    
