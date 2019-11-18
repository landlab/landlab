# -*- coding: utf-8 -*-
"""
GroundwaterDupuitPercolator Component

@author: G Tucker
"""

import numpy as np

from landlab import Component
from landlab.grid.mappers import map_mean_of_link_nodes_to_link
from landlab.utils import return_array_at_link, return_array_at_node

ACTIVE_LINK = 0


class GroundwaterDupuitPercolator(Component):
    """
    Simulate groundwater flow in a shallow unconfined aquifer.

    The GroundwaterDupuitPercolator uses the Dupuit approximation that the
    hydraulic gradient is equal to the slope of the water table.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> mg = RasterModelGrid((3, 3))
    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> gdp = GroundwaterDupuitPercolator(mg)

    Notes
    -----
    Groundwater discharge per unit length, q, is calculated as:

        q = - K H dw/dx,

    where K is hydraulic conductivity, H is aquifer thickness, w is water table
    height, and x is horizontal distance.

    An explicit forward-in-time finite-volume method is used to implement a
    numerical solution. Flow discharge between neighboring nodes is calculated
    using the average depth at the nodes.

    Note that the current version does NOT handle surface seepage.
    """

    _name = "GroundwaterDupuitPercolator"

    _input_var_names = set(("topographic__elevation", "aquifer_base__elevation"))

    _output_var_names = set(
        (
            "aquifer__thickness",
            "water_table__elevation",
            "hydraulic__gradient",
            "groundwater__specific_discharge",
            "groundwater__velocity",
        )
    )

    _var_units = {
        "topographic__elevation": "m",
        "aquifer_base__elevation": "m",
        "aquifer__thickness": "m",
        "water_table__elevation": "m",
        "hydraulic__gradient": "m/m",
        "groundwater__specific_discharge": "m2/s",
        "groundwater__velocity": "m/s",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "aquifer_base__elevation": "node",
        "aquifer__thickness": "node",
        "water_table__elevation": "node",
        "hydraulic__gradient": "link",
        "groundwater__specific_discharge": "link",
        "groundwater__velocity": "link",
    }

    _var_doc = {
        "topographic__elevation": "elevation of land surface",
        "aquifer_base__elevation": "elevation of impervious layer",
        "aquifer__thickness": "thickness of saturated zone",
        "water_table__elevation": "elevation of water table",
        "hydraulic__gradient": "gradient of water table in link direction",
        "groundwater__specific_discharge": "discharge per width in link dir",
        "groundwater__velocity": "velocity of groundwater in link direction",
    }

    def __init__(self, grid, hydraulic_conductivity=0.01, recharge_rate=1.0e-8):
        """Initialize the GroundwaterDupuitPercolator.

        Parameters
        ----------
        grid: ModelGrid
            Landlab ModelGrid object
        hydraulic_conductivity: float, field name (link), or array of float
            saturated hydraulic conductivity, m/s
            Default = 0.01 m/s
        recharge_rate: float, field name (node), or array of float
            Rate of recharge, m/s
            Default = 1.0e-8 m/s
        """
        # Store grid
        self._grid = grid

        # Shorthand
        self._cores = grid.core_nodes
        self._inactive_links = np.where(grid.status_at_link != ACTIVE_LINK)[0]

        # Convert parameters to fields if needed, and store a reference
        self._K = return_array_at_link(grid, hydraulic_conductivity)
        self._recharge = return_array_at_node(grid, recharge_rate)

        # Create fields:

        if "topographic__elevation" in self._grid.at_node:
            self._elev = self._grid.at_node["topographic__elevation"]
        else:
            self._elev = self._grid.add_ones("node", "topographic__elevation")

        if "aquifer_base__elevation" in self._grid.at_node:
            self._base = self._grid.at_node["aquifer_base__elevation"]
        else:
            self._base = self._grid.add_zeros("node", "aquifer_base__elevation")

        if "water_table__elevation" in self._grid.at_node:
            self._wtable = self._grid.at_node["water_table__elevation"]
        else:
            self._wtable = self._grid.add_zeros("node", "water_table__elevation")

        if "aquifer__thickness" in self._grid.at_node:
            self._thickness = self._grid.at_node["aquifer__thickness"]
        else:
            self._thickness = self._grid.add_zeros("node", "aquifer__thickness")
            self._thickness[:] = self._wtable - self._base

        if "hydraulic__gradient" in self._grid.at_link:
            self._hydr_grad = self._grid.at_link["hydraulic__gradient"]
        else:
            self._hydr_grad = self._grid.add_zeros("link", "hydraulic__gradient")

        if "groundwater__specific_discharge" in self._grid.at_link:
            self._q = self._grid.at_link["groundwater__specific_discharge"]
        else:
            self._q = self._grid.add_zeros("link", "groundwater__specific_discharge")

        if "groundwater__velocity" in self._grid.at_link:
            self._vel = self._grid.at_link["groundwater__velocity"]
        else:
            self._vel = self._grid.add_zeros("link", "groundwater__velocity")

    def run_one_step(self, dt, **kwds):
        """
        Advance component by one time step of size dt.

        Parameters
        ----------
        dt: float (time in seconds)
            The imposed timestep.
        """

        # Calculate hydraulic gradient
        self._hydr_grad[:] = self._grid.calc_grad_at_link(self._wtable)
        self._hydr_grad[self._inactive_links] = 0.0

        # Calculate groundwater velocity
        self._vel[:] = -self._K * self._hydr_grad
        self._vel[self._grid.status_at_link != 0] = 0.0

        # Aquifer thickness at links
        hlink = map_mean_of_link_nodes_to_link(self._grid, "aquifer__thickness")

        # Calculate specific discharge
        self._q[:] = hlink * self._vel

        # Mass balance
        dqdx = self._grid.calc_flux_div_at_node(self._q)
        dhdt = self._recharge - dqdx

        # Update
        self._thickness[self._grid.core_nodes] += dhdt[self._cores] * dt

        # Recalculate water surface height
        self._wtable[self._grid.core_nodes] = (
            self._base[self._cores] + self._thickness[self._cores]
        )
