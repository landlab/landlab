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

        # Flow direction
        self.fdir = FlowDirectorSteepest(self.elev)

    def soilflux(self, dt):
        """Calculate soil flux for a time period 'dt'.
        """

        # Calculate gradients
        self.slope[:] = self.grid.calc_grad_at_link(self.elev)
        self.slope[self.grid.status_at_link == INACTIVE_LINK] = 0.

        # Run flow directior (steepest slope)
        self.fdir.run_one_step()
        # Downstream steepest slope at node:
        self.steepest = self.grid.at_node['topographic__steepest_slope']
        # On each node, node ID of downstream receiver node
        # (on node (i), ID of node that receives flow from node (i)):
        self.receiver = self.grid.at_node['flow__receiver_node']

        for i in (self.grid.core_nodes):
            # Calculate influx on node i = outflux of nodes whose receiver is i
            self.flux_in[self.receiver[i]] += self.flux_out[i]

            # Calculate deposition on node and transfer over node
            if self.flux_in[i] > 0.:
                # Calculate transport length
                self.L = (self.grid.dx)/(1-(self.steepest/self.slope_crit)**2)
                # Calculate deposition
                self.depo[i] = self.flux_in[i]/self.L[i]
                # Calculate transfer
                self.trans[i] = self.flux_in[i]-self.depo[i]
            else:
                self.depo[i] = 0.
                self.trans[i] = 0.

            # Calculate erosion on node
        for i in range(self.grid.core_nodes):
            if self.slope[i] > self.slope_crit:
                self.slope[i] = self.slope_crit
            else:
                pass
            self.erosion[:] = -self.k * self.steepest

            # Calculate outflux
        self.flux_out[:] = self.erosion + self.trans     # Flux out of node

        # Update topography
        self.elev += ((-self.erosion + self.depo) * dt)/(grid.dx*grid.dy)

        # reset erosion, flux_in, flux_out, depo, trans to 0
        self.flux_in[:] = 0.
        self.flux_out[:] = 0.
        self.depo[:] = 0.
        self.trans[:] = 0.


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