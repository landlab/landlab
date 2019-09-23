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
    >>> abe = mg.add_zeros('node', 'aquifer_base__elevation')
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

    _info = {
        "aquifer__thickness": {
            "dtype":None,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "thickness of saturated zone",
        },
        "aquifer_base__elevation": {
            "dtype":None,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "elevation of impervious layer",
        },
        "groundwater__specific_discharge": {
            "dtype":None,
            "intent": "out",
            "optional": False,
            "units": "m2/s",
            "mapping": "link",
            "doc": "discharge per width in link dir",
        },
        "groundwater__velocity": {
            "dtype":None,
            "intent": "out",
            "optional": False,
            "units": "m/s",
            "mapping": "link",
            "doc": "velocity of groundwater in link direction",
        },
        "hydraulic__gradient": {
            "dtype":None,
            "intent": "out",
            "optional": False,
            "units": "m/m",
            "mapping": "link",
            "doc": "gradient of water table in link direction",
        },
        "topographic__elevation": {
            "dtype":None,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "elevation of land surface",
        },
        "water_table__elevation": {
            "dtype":None,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "elevation of water table",
        },
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
        super(GroundwaterDupuitPercolator, self).__init__(grid)

        # Shorthand
        self._cores = grid.core_nodes
        self._inactive_links = np.where(grid.status_at_link != ACTIVE_LINK)[0]

        # Convert parameters to fields if needed, and store a reference
        self._K = return_array_at_link(grid, hydraulic_conductivity)
        self._recharge = return_array_at_node(grid, recharge_rate)

        # Create fields:
        self._elev = self._grid.at_node["topographic__elevation"]
        self._base = self._grid.at_node["aquifer_base__elevation"]

        self.initialize_output_fields()

        self._wtable = self._grid.at_node["water_table__elevation"]
        self._thickness = self._grid.at_node["aquifer__thickness"]
        self._thickness[:] = self._wtable - self._base
        self._hydr_grad = self._grid.at_link["hydraulic__gradient"]
        self._q = self._grid.at_link["groundwater__specific_discharge"]
        self._vel = self._grid.at_link["groundwater__velocity"]

    def run_one_step(self, dt):
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
