# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 08:32:48 2016

@author: RCGlade
"""

import numpy as np

from landlab import INACTIVE_LINK, Component


class DepthDependentDiffuser(Component):

    """
    This component implements a depth and slope dependent linear diffusion rule
    in the style of Johnstone and Hilley (2014).

    This component will ignore soil thickness located at non-core nodes.

    Parameters
    ----------
    grid: ModelGrid
        Landlab ModelGrid object
    linear_diffusivity: float
        Hillslope diffusivity, m**2/yr
        Equivalent to the soil creep efficiency
        times the soil transport decay depth.
    soil_transport_decay_depth: float
        characteristic transport soil depth, m

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import ExponentialWeatherer
    >>> from landlab.components import DepthDependentDiffuser
    >>> mg = RasterModelGrid((5, 5))
    >>> soilTh = mg.add_zeros('node', 'soil__depth')
    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> BRz = mg.add_zeros('node', 'bedrock__elevation')
    >>> expweath = ExponentialWeatherer(mg)
    >>> DDdiff = DepthDependentDiffuser(mg)
    >>> expweath.calc_soil_prod_rate()
    >>> np.allclose(mg.at_node['soil_production__rate'][mg.core_nodes], 1.)
    True
    >>> DDdiff.soilflux(2.)
    >>> np.allclose(mg.at_node['topographic__elevation'][mg.core_nodes], 0.)
    True
    >>> np.allclose(mg.at_node['bedrock__elevation'][mg.core_nodes], -2.)
    True
    >>> np.allclose(mg.at_node['soil__depth'][mg.core_nodes], 2.)
    True

    Now with a slope:

    >>> mg = RasterModelGrid((3, 5))
    >>> soilTh = mg.add_zeros('node', 'soil__depth')
    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> BRz = mg.add_zeros('node', 'bedrock__elevation')
    >>> z += mg.node_x.copy()
    >>> BRz += mg.node_x/2.
    >>> soilTh[:] = z - BRz
    >>> expweath = ExponentialWeatherer(mg)
    >>> DDdiff = DepthDependentDiffuser(mg)
    >>> expweath.calc_soil_prod_rate()
    >>> np.allclose(
    ...     mg.at_node['soil_production__rate'][mg.core_nodes],
    ...     np.array([ 0.60653066,  0.36787944,  0.22313016]))
    True
    >>> DDdiff.soilflux(2.)
    >>> np.allclose(
    ...     mg.at_node['topographic__elevation'][mg.core_nodes],
    ...     np.array([ 1.47730244,  2.28949856,  3.17558975]))
    True
    >>> np.allclose(mg.at_node['bedrock__elevation'][mg.core_nodes],
    ...     np.array([-0.71306132,  0.26424112,  1.05373968]))
    True
    >>> np.allclose(mg.at_node['soil__depth'], z - BRz)
    True
    """

    _name = "DepthDependentDiffuser"

    _input_var_names = (
        "topographic__elevation",
        "soil__depth",
        "soil_production__rate",
    )

    _output_var_names = (
        "soil__flux",
        "topographic__slope",
        "topographic__elevation",
        "bedrock__elevation",
        "soil__depth",
    )

    _var_units = {
        "topographic__elevation": "m",
        "topographic__slope": "m/m",
        "soil__depth": "m",
        "soil__flux": "m^2/yr",
        "soil_production__rate": "m/yr",
        "bedrock__elevation": "m",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "topographic__slope": "link",
        "soil__depth": "node",
        "soil__flux": "link",
        "soil_production__rate": "node",
        "bedrock__elevation": "node",
    }

    _var_doc = {
        "topographic__elevation": "elevation of the ground surface",
        "topographic__slope": "gradient of the ground surface",
        "soil__depth": "depth of soil/weather bedrock",
        "soil__flux": "flux of soil in direction of link",
        "soil_production__rate": "rate of soil production at nodes",
        "bedrock__elevation": "elevation of the bedrock surface",
    }

    def __init__(self, grid, linear_diffusivity=1.0, soil_transport_decay_depth=1.0):
        """Initialize the DepthDependentDiffuser."""

        # Store grid and parameters
        self._grid = grid
        self.K = linear_diffusivity
        self.soil_transport_decay_depth = soil_transport_decay_depth

        # create fields
        # elevation
        if "topographic__elevation" in self.grid.at_node:
            self.elev = self.grid.at_node["topographic__elevation"]
        else:
            self.elev = self.grid.add_zeros("node", "topographic__elevation")

        # slope
        if "topographic__slope" in self.grid.at_link:
            self.slope = self.grid.at_link["topographic__slope"]
        else:
            self.slope = self.grid.add_zeros("link", "topographic__slope")

        # soil depth
        if "soil__depth" in self.grid.at_node:
            self.depth = self.grid.at_node["soil__depth"]
        else:
            self.depth = self.grid.add_zeros("node", "soil__depth")

        # soil flux
        if "soil__flux" in self.grid.at_link:
            self.flux = self.grid.at_link["soil__flux"]
        else:
            self.flux = self.grid.add_zeros("link", "soil__flux")

        # weathering rate
        if "soil_production__rate" in self.grid.at_node:
            self.soil_prod_rate = self.grid.at_node["soil_production__rate"]
        else:
            self.soil_prod_rate = self.grid.add_zeros("node", "soil_production__rate")

        # bedrock elevation
        if "bedrock__elevation" in self.grid.at_node:
            self.bedrock = self.grid.at_node["bedrock__elevation"]
        else:
            self.bedrock = self.grid.add_zeros("node", "bedrock__elevation")

    def soilflux(self, dt):
        """Calculate soil flux for a time period 'dt'.

        Parameters
        ----------

        dt: float (time)
            The imposed timestep.
        """

        # update soil thickness
        self.grid.at_node["soil__depth"][:] = (
            self.grid.at_node["topographic__elevation"]
            - self.grid.at_node["bedrock__elevation"]
        )

        # Calculate soil depth at links.
        H_link = self.grid.map_value_at_max_node_to_link(
            "topographic__elevation", "soil__depth"
        )

        # Calculate gradients
        slope = self.grid.calc_grad_at_link(self.elev)
        slope[self.grid.status_at_link == INACTIVE_LINK] = 0.0

        # Calculate flux
        self.flux[:] = (
            -self.K * slope * (1.0 - np.exp(-H_link / self.soil_transport_decay_depth))
        )

        # Calculate flux divergence
        dqdx = self.grid.calc_flux_div_at_node(self.flux)

        # Calculate change in soil depth
        dhdt = self.soil_prod_rate - dqdx

        # Calculate soil depth at nodes
        self.depth[self.grid.core_nodes] += dhdt[self.grid.core_nodes] * dt

        # prevent negative soil thickness
        self.depth[self.depth < 0.0] = 0.0

        # Calculate bedrock elevation
        self.bedrock[self.grid.core_nodes] -= (
            self.soil_prod_rate[self.grid.core_nodes] * dt
        )

        # Update topography
        self.elev[self.grid.core_nodes] = (
            self.depth[self.grid.core_nodes] + self.bedrock[self.grid.core_nodes]
        )

    def run_one_step(self, dt):
        """

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        """

        self.soilflux(dt)
