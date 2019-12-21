#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 10:13:38 2017

@author: margauxmouchene
"""


import numpy as np

from landlab import Component


class TransportLengthHillslopeDiffuser(Component):

    """
    Hillslope diffusion component in the style of Carretier et al. (2016,
        ESurf), and Davy and Lague (2009)

    dz/dt = - E + D (+ Uplift)
    D = qs / L
    E = k * S
    L = dx / (1 - (S / Sc)^2)

    Works on regular raster-type grid (RasterModelGrid, dx=dy).
    To be coupled with FlowDirectorSteepest for the calculation of steepest
     slope at each timestep.

    Component written by Margaux Mouchene, 2017

    Parameters
    ----------
    grid : ModelGrid
        Landlab ModelGrid object
    erodibility: float
        Erodibility coefficient [L/T]
    slope_crit: float (default=1.)
        Critical slope [L/L]

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowDirectorSteepest
    >>> from landlab.components import TransportLengthHillslopeDiffuser

    Define grid and initial topography:
        - 3x5 grid
        - east and west boundaries are open, north and south are closed
        - Initial topography is plane at base level on the boundaries and
            1m of elevation elsewhere (core)

    >>> mg = RasterModelGrid((5, 5))
    >>> mg.set_closed_boundaries_at_grid_edges(False, True, False, True)
    >>> z = np.array([0., 0., 0., 0., 0.,
    ...               0., 1., 1., 1., 0.,
    ...               0., 1., 1., 1., 0.,
    ...               0., 1., 1., 1., 0.,
    ...               0., 0., 0., 0., 0.])
    >>> _ = mg.add_field('node', 'topographic__elevation', z)

    Instantiate Flow director (steepest slope type) and TL hillslope diffuser

    >>> fdir = FlowDirectorSteepest(mg)
    >>> tl_diff = TransportLengthHillslopeDiffuser(mg, \
                                                   erodibility=0.001,\
                                                   slope_crit=0.6)

    Run the components for ten short timepsteps

    >>> for t in range(10):
    ...     fdir.run_one_step()
    ...     tl_diff.run_one_step(1.)

    Check final topography

    >>> np.allclose(
    ...     mg.at_node['topographic__elevation'],
    ...     np.array([ 0.,  0.        ,  0.        ,  0.        ,  0.,
    ...                0.,  0.96175283,  0.99982519,  0.96175283,  0.,
    ...                0.,  0.96175283,  0.99982519,  0.96175283,  0.,
    ...                0.,  0.96175283,  0.99982519,  0.96175283,  0.,
    ...                0.,  0.        ,  0.        ,  0.        ,  0.]))
    True
    """

    _name = "TransportLengthHillslopeDiffuser"

    _input_var_names = set(
        ("topographic__elevation", "flow__receiver_node", "topographic__steepest_slope")
    )

    _output_var_names = set(
        (
            "topographic__elevation",
            "sediment__deposition_rate",
            "sediment__transfer_rate",
            "sediment__deposition_coeff",
            "sediment__flux_in",
            "sediment__flux_out",
            "sediment__erosion_rate",
        )
    )

    _var_units = {
        "topographic__elevation": "m",
        "flow__receiver_node": "-",
        "topographic__steepest_slope": "m/m",
        "sediment__deposition_rate": "m/yr",
        "sediment__transfer_rate": "m/yr",
        "sediment__deposition_coeff": "-",
        "sediment__flux_in": "m/yr",
        "sediment__flux_out": "m/yr",
        "sediment__erosion_rate": "m/yr",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "flow__receiver_node": "node",
        "topographic__steepest_slope": "node",
        "sediment__deposition_rate": "node",
        "sediment__transfer_rate": "node",
        "sediment__deposition_coeff": "node",
        "sediment__flux_in": "node",
        "sediment__flux_out": "node",
        "sediment__erosion_rate": "node",
    }

    _var_doc = {
        "topographic__elevation": "Elevation of the ground surface",
        "flow__receiver_node": "Node array of receivers (node that receives flow from "
        "current node)",
        "topographic__steepest_slope": "Steepest gradient of the ground surface at each node",
        "sediment__deposition_rate": "Deposition rate on node",
        "sediment__transfer_rate": "Rate of transferred sediment across a node (incoming "
        "sediment - deposited sediment on node)",
        "sediment__deposition_coeff": "Fraction of incoming sediment that is deposited on the node",
        "sediment__flux_in": "Incoming sediment rate on node (=qs/dx)",
        "sediment__flux_out": "Outgoing sediment rate on node = sediment eroded on"
        " node + sediment transported across node from upstream",
        "sediment__erosion_rate": "Erosion rate on node",
    }

    def __init__(self, grid, erodibility, slope_crit=1.0):

        """Initialize Diffuser.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        erodibility: float
            Erodibility coefficient [L/T]
        slope_crit: float (default=1.)
            Critical slope [L/L]
        """
        if grid.at_node["flow__receiver_node"].size != grid.size("node"):
            msg = (
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that TransportLengthHillslopeDiffuser is compatible "
                "with route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )
            raise NotImplementedError(msg)

        # Store grid and parameters
        self._grid = grid
        self.k = erodibility
        self.slope_crit = slope_crit

        # Create fields:
        # Elevation
        if "topographic__elevation" in self.grid.at_node:
            self.elev = self.grid.at_node["topographic__elevation"]
        else:
            self.elev = self.grid.add_zeros("node", "topographic__elevation")

        # Deposition
        if "sediment__deposition_rate" in self.grid.at_node:
            self.depo = self.grid.at_node["sediment__deposition_rate"]
        else:
            self.depo = self.grid.add_zeros("node", "sediment__deposition_rate")

        # Transferred sediments (crossing over node)
        if "sediment__transfer_rate" in self.grid.at_node:
            self.trans = self.grid.at_node["sediment__transfer_rate"]
        else:
            self.trans = self.grid.add_zeros("node", "sediment__transfer_rate")

        # Transport coefficient
        if "sediment__deposition_coeff" in self.grid.at_node:
            self.d_coeff = self.grid.at_node["sediment__deposition_coeff"]
        else:
            self.d_coeff = self.grid.add_zeros("node", "sediment__deposition_coeff")

        # Flux in
        if "sediment__flux_in" in self.grid.at_node:
            self.flux_in = self.grid.at_node["sediment__flux_in"]
        else:
            self.flux_in = self.grid.add_zeros("node", "sediment__flux_in")

        # Flux out
        if "sediment__flux_out" in self.grid.at_node:
            self.flux_out = self.grid.at_node["sediment__flux_out"]
        else:
            self.flux_out = self.grid.add_zeros("node", "sediment__flux_out")

        # Erosion
        if "sediment__erosion_rate" in self.grid.at_node:
            self.erosion = self.grid.at_node["sediment__erosion_rate"]
        else:
            self.erosion = self.grid.add_zeros("node", "sediment__erosion_rate")

    def tldiffusion(self, dt):
        """Calculate hillslope diffusion for a time period 'dt'.

        Parameters:
        grid : ModelGrid
            Landlab ModelGrid object
        dt: float (time)
            The imposed timestep.
        """

        # Reset erosion, depo, trans and flux_in to 0
        self.erosion[:] = 0.0
        self.depo[:] = 0.0
        self.trans[:] = 0.0
        self.flux_in[:] = 0.0

        # Downstream steepest slope at node:
        self.steepest = self.grid.at_node["topographic__steepest_slope"]
        # On each node, node ID of downstream receiver node
        # (on node (i), ID of node that receives flow from node (i)):
        self.receiver = self.grid.at_node["flow__receiver_node"]

        dx = self.grid.dx
        cores = self.grid.core_nodes

        # Calculate influx rate on node i  = outflux of nodes
        # whose receiver is i
        for i in self.grid.core_nodes:
            self.flux_in[self.receiver[i]] += self.flux_out[i]

            # Calculate transport coefficient
            # When S ~ Scrit, d_coeff is set to "infinity", for stability and
            # so that there is no deposition
            if self.steepest[i] >= self.slope_crit:
                self.d_coeff[i] = 1000000000.0
            else:
                self.d_coeff[i] = 1 / (
                    1 - (np.power(((self.steepest[i]) / self.slope_crit), 2))
                )

        # Calculate deposition rate on node
        self.depo[cores] = self.flux_in[cores] / self.d_coeff[cores]

        # Calculate erosion rate on node (positive value)
        # If S > Scrit, erosion is simply set for the slope to return to Scrit
        # Otherwise, erosion is slope times erodibility coefficent
        for i in self.grid.core_nodes:
            if self.steepest[i] > self.slope_crit:
                self.erosion[i] = dx * (self.steepest[i] - self.slope_crit) / (100 * dt)
            else:
                self.erosion[i] = self.k * self.steepest[i]

            # Update elevation
            self.elev[i] += (-self.erosion[i] + self.depo[i]) * dt

        # Calculate transfer rate over node
        self.trans[cores] = self.flux_in[cores] - self.depo[cores]

        # Calculate outflux rate
        self.flux_out[:] = self.erosion + self.trans

    def run_one_step(self, dt):
        """
        Advance transport length-model hillslope diffusion component
        by one time step of size dt and tests for timestep stability.

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        """
        self.tldiffusion(dt)

        # Test code stability for timestep dt
        # Raise unstability error if local slope is reversed by erosion
        # and deposition during a timestep dt
        elev_dif = self.elev - self.elev[self.receiver]
        s = elev_dif[np.where(self.grid.at_node["flow__sink_flag"] == 0)]
        if np.any(s < -1) is True:
            raise ValueError(
                "The component is unstable" " for such a large timestep " "on this grid"
            )
        else:
            pass
