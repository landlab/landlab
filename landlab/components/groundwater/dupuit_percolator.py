# -*- coding: utf-8 -*-
"""
GroundwaterDupuitPercolator Component

@author: G Tucker
"""

import numpy as np
from landlab import Component
from landlab.utils import return_array_at_node, return_array_at_link
from landlab.grid.mappers import map_mean_of_link_nodes_to_link


ACTIVE_LINK = 0


class GroundwaterDupuitPercolator(Component):
    """
    Simulate groundwater flow in a shallow unconfined aquifer.

    The GroundwaterDupuitPercolator uses the Dupuit approximation that the
    hydraulic gradient is equal to the slope of the water table.

    Parameters
    ----------
    grid: ModelGrid
            Landlab ModelGrid object
    hydraulic_conductivity: float, field name, or array of float
            saturated hydraulic conductivity, m/s
            Default = 0.01 m/s
    recharge_rate: float, field name, or array of float
            Rate of recharge, m/s
            Default = 1.0e-8 m/s

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

    _input_var_names = set(("topographic__elevation",
                            "aquifer_base__elevation"))

    _output_var_names = set(
        ("aquifer__thickness", "water_table__elevation", "hydraulic__gradient",
         "groundwater__specific_discharge", "groundwater__velocity")
    )

    _var_units = {
        "topographic__elevation": "m",
        "aquifer_base__elevation": "m",
        "aquifer__thickness": "m",
        "water_table__elevation": "m",
        "hydraulic__gradient": "m/m",
        "groundwater__specific_discharge": "m2/s",
        "groundwater__velocity": "m/s"
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

    def __init__(self, grid, hydraulic_conductivity=0.01,
                 recharge_rate=1.0e-8):
        """Initialize the GroundwaterDupuitPercolator.

        Parameters
        ----------
        grid: ModelGrid
            Landlab ModelGrid object
        hydraulic_conductivity: float, field name, or array of float
                saturated hydraulic conductivity, m/s
                Default = 0.01 m/s
        recharge_rate: float, field name, or array of float
                Rate of recharge, m/s
                Default = 1.0e-8 m/s
        """
        # Store grid
        self._grid = grid

        # Shorthand
        self.cores = grid.core_nodes
        self.inactive_links = np.where(grid.status_at_link != ACTIVE_LINK)[0]

        # Convert parameters to fields if needed, and store a reference
        self.K = return_array_at_link(grid, hydraulic_conductivity)
        self.recharge = return_array_at_node(grid, recharge_rate)

        # Create fields:

        if "topographic__elevation" in self.grid.at_node:
            self.elev = self.grid.at_node["topographic__elevation"]
        else:
            self.elev = self.grid.add_ones("node", "topographic__elevation")

        if "aquifer_base__elevation" in self.grid.at_node:
            self.base = self.grid.at_node["aquifer_base__elevation"]
        else:
            self.base = self.grid.add_zeros("node", "aquifer_base__elevation")

        if "water_table__elevation" in self.grid.at_node:
            self.wtable = self.grid.at_node["water_table__elevation"]
        else:
            self.wtable = self.grid.add_zeros("node", "water_table__elevation")

        if "aquifer__thickness" in self.grid.at_node:
            self.thickness = self.grid.at_node["aquifer__thickness"]
        else:
            self.thickness = self.grid.add_zeros("node", "aquifer__thickness")
            self.thickness[:] = self.wtable - self.base

        if "hydraulic__gradient" in self.grid.at_link:
            self.hydr_grad = self.grid.at_link["hydraulic__gradient"]
        else:
            self.hydr_grad = self.grid.add_zeros("link", "hydraulic__gradient")

        if "groundwater__specific_discharge" in self.grid.at_link:
            self.q = self.grid.at_link["groundwater__specific_discharge"]
        else:
            self.q = self.grid.add_zeros("link",
                                         "groundwater__specific_discharge")

        if "groundwater__velocity" in self.grid.at_link:
            self.vel = self.grid.at_link["groundwater__velocity"]
        else:
            self.vel = self.grid.add_zeros("link", "groundwater__velocity")

    def run_one_step(self, dt, **kwds):
        """
        Advance component by one time step of size dt.

        Parameters
        ----------
        dt: float (time in seconds)
            The imposed timestep.
        """

        # Calculate hydraulic gradient
        self.hydr_grad[:] = self._grid.calc_grad_at_link(self.wtable)
        self.hydr_grad[self.inactive_links] = 0.0

        # Calculate groundwater velocity
        self.vel[:] = -self.K * self.hydr_grad
        self.vel[self._grid.status_at_link != 0] = 0.0

        # Aquifer thickness at links
        hlink = map_mean_of_link_nodes_to_link(self._grid,
                                               'aquifer__thickness')

        # Calculate specific discharge
        self.q[:] = hlink * self.vel

        # Mass balance
        dqdx = self._grid.calc_flux_div_at_node(self.q)
        dhdt = self.recharge - dqdx

        # Update
        self.thickness[self._grid.core_nodes] += dhdt[self.cores] * dt

        # Recalculate water surface height
        self.wtable[self._grid.core_nodes] = (self.base[self.cores]
                                              + self.thickness[self.cores])
