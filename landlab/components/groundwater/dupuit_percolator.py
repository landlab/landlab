# -*- coding: utf-8 -*-
"""GroundwaterDupuitPercolator Component.

@author: G Tucker, D Litwin
"""

import numpy as np

from landlab import Component
from landlab.grid.base import INACTIVE_LINK
from landlab.grid.mappers import (
    map_max_of_node_links_to_node,
    map_mean_of_link_nodes_to_link,
    map_value_at_max_node_to_link,
)
from landlab.utils import return_array_at_link, return_array_at_node


# regularization functions used to deal with numerical demons of seepage
def _regularize_G(u, reg_factor):
    """Smooths transition of step function with an exponential.

    0<=u<=1.
    """
    return np.exp(-(1 - u) / reg_factor)


def _regularize_R(u):
    """ramp function on u."""
    return u * np.greater_equal(u, 0)


def get_link_hydraulic_conductivity(grid, K):
    """Returns array of hydraulic conductivity on links, allowing for aquifers
    with laterally anisotropic hydraulic conductivity.

    Parameters
    ----------
    K: (2x2) array of floats (m/s)
        The hydraulic conductivity tensor:
        [[Kxx, Kxy],[Kyx,Kyy]]
    """

    u = grid.unit_vector_at_link
    K_link = np.zeros(len(u))
    for i in range(len(u)):
        K_link[i] = np.dot(np.dot(u[i, :], K), u[i, :])
    return K_link


class GroundwaterDupuitPercolator(Component):
    """Simulate groundwater flow in a shallow unconfined aquifer.

    The GroundwaterDupuitPercolator solves the Boussinesq equation for
    flow in an unconfined aquifer over an impermeable aquifer base that may
    have uneven slope. This method uses the Dupuit approximation that the
    hydraulic gradient is equal to the slope of the water table.

    Examples
    --------
    Import the grid class and component

    >>> from landlab import RasterModelGrid
    >>> from landlab.components import GroundwaterDupuitPercolator

    Initialize the grid and component

    >>> grid = RasterModelGrid((10, 10), xy_spacing=10.0)
    >>> elev = grid.add_zeros('node', 'topographic__elevation')
    >>> abe = grid.add_zeros('node', 'aquifer_base__elevation')
    >>> elev[:] = 5.0
    >>> gdp = GroundwaterDupuitPercolator(grid)

    Run component forward. Note all time units in the model are in seconds.

    >>> dt = 1E4
    >>> for i in range(100):
    ...     gdp.run_one_step(dt)

    In an example that generates surface water leakage, the surface water flux
    out of the domain can be calculated only after a FlowAccumulator is run.
    Here is a more advanced example with a sloping aquifer that returns surface water flow.

    >>> from landlab.components import FlowAccumulator

    Set boundary conditions and initialize grid

    >>> grid = RasterModelGrid((5, 41), xy_spacing=10.0)
    >>> grid.set_closed_boundaries_at_grid_edges(True, True, False, True)

    Make a sloping, 3 m thick aquifer, initially fully saturated

    >>> elev = grid.add_zeros('node', 'topographic__elevation')
    >>> elev[:] = grid.x_of_node/100+3
    >>> base = grid.add_zeros('node', 'aquifer_base__elevation')
    >>> base[:] = grid.x_of_node/100
    >>> wt = grid.add_zeros('node', 'water_table__elevation')
    >>> wt[:] = grid.x_of_node/100 + 3

    Initialize components

    >>> gdp = GroundwaterDupuitPercolator(grid, recharge_rate=1E-7)
    >>> fa = FlowAccumulator(grid, runoff_rate='surface_water__specific_discharge')

    Advance timestep

    >>> dt = 1E3
    >>> for i in range(1000):
    ...     gdp.run_one_step(dt)

    Calculate surface water flux out of domain

    >>> fa.run_one_step()
    >>> np.testing.assert_almost_equal(gdp.calc_sw_flux_out(),0.0005077)


    Notes
    -----
    Groundwater discharge per unit length, q, is calculated as:

        q = - K H ( dH/dx cos(alpha) + sin(alpha) ),

    where K is hydraulic conductivity, H is aquifer thickness, alpha is
    the slope angle of the aquifer base, and x is horizontal distance.

    Surface water discharge per unit area, qs, is calculated as:

        qs = G( H/(Z-Zb) ) * R( f - dq/dx)

    where G is a smoothed step function, R is the ramp function, Z is the
    topographic elevation, Zb is the aquifer base elevation, and f is
    the recharge rate.

    An explicit forward-in-time finite-volume method is used to implement a
    numerical solution. Flow discharge between neighboring nodes is calculated
    using the saturated thickness at the up-gradient node.
    """

    _name = "GroundwaterDupuitPercolator"

    _info = {
        "aquifer__thickness": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "thickness of saturated zone",
        },
        "aquifer_base__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "elevation of impervious layer",
        },
        "aquifer_base__gradient": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/m",
            "mapping": "link",
            "doc": "gradient of the aquifer base in the link direction",
        },
        "groundwater__specific_discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m2/s",
            "mapping": "link",
            "doc": "discharge per width in link dir",
        },
        "groundwater__velocity": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/s",
            "mapping": "link",
            "doc": "velocity of groundwater in link direction",
        },
        "hydraulic__gradient": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/m",
            "mapping": "link",
            "doc": "gradient of water table in link direction",
        },
        "surface_water__specific_discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/s",
            "mapping": "node",
            "doc": "rate of seepage to surface",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "water_table__elevation": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "elevation of water table",
        },
        "water_table__velocity": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/s",
            "mapping": "node",
            "doc": "rate of change of water table elevation",
        },
    }

    def __init__(
        self,
        grid,
        hydraulic_conductivity=0.001,
        porosity=0.2,
        recharge_rate=1.0e-8,
        regularization_f=1e-2,
        courant_coefficient=0.01,
    ):
        """Initialize the GroundwaterDupuitPercolator.

        Parameters
        ----------
        grid: ModelGrid
            Landlab ModelGrid object
        hydraulic_conductivity: float, field name, or array of float
            saturated hydraulic conductivity, m/s
            Default = 0.001 m/s
        porosity: float, field name or array of float
            the porosity of the aquifer [-]
            Default = 0.2
        recharge_rate: float, field name, or array of float
            Rate of recharge, m/s
            Default = 1.0e-8 m/s
        regularization_f: float
            factor controlling the smoothness of the transition between
            surface and subsurface flow
        courant_coefficient: float (-)
            The muliplying factor on the condition that the timestep is
            smaller than the minimum link length over groundwater flow
            velocity. This parameter is only used with
            ``run_with_adaptive_time_step_solver`` and must be greater than
            zero.
        """
        super(GroundwaterDupuitPercolator, self).__init__(grid)

        # Shorthand
        self._cores = grid.core_nodes

        # Create fields:

        self._elev = self._grid.at_node["topographic__elevation"]
        self._base = self._grid.at_node["aquifer_base__elevation"]

        self.initialize_output_fields()

        self._wtable = self._grid.at_node["water_table__elevation"]
        self._wtable[grid.closed_boundary_nodes] = 0

        self._thickness = self._grid.at_node["aquifer__thickness"]
        self._thickness[:] = self._wtable - self._base
        self._thickness[grid.closed_boundary_nodes] = 0

        self._hydr_grad = self._grid.at_link["hydraulic__gradient"]
        self._base_grad = self._grid.at_link["aquifer_base__gradient"]

        self._q = self._grid.at_link["groundwater__specific_discharge"]
        self._qs = self._grid.at_node["surface_water__specific_discharge"]

        self._vel = self._grid.at_link["groundwater__velocity"]

        self._dhdt = self._grid.at_node["water_table__velocity"]

        # Convert parameters to fields if needed, and store a reference
        self._K = return_array_at_link(grid, hydraulic_conductivity)
        self._recharge = return_array_at_node(grid, recharge_rate)
        self._n = return_array_at_node(grid, porosity)
        self._n_link = map_mean_of_link_nodes_to_link(self._grid, self._n)
        self._r = regularization_f
        self._S = abs(grid.calc_grad_at_link(self._elev))
        self._S_node = map_max_of_node_links_to_node(grid, self._S)

        # save courant_coefficient (and test)
        self.courant_coefficient = courant_coefficient

    @property
    def courant_coefficient(self):
        """Courant coefficient for adaptive time step.

        Parameters
        ----------
        courant_coefficient: float (-)
            The muliplying factor on the condition that the timestep is
            smaller than the minimum link length over groundwater flow
            velocity. This parameter is only used with
            ``run_with_adaptive_time_step_solver`` and must be greater than
            zero.
        """

    @courant_coefficient.setter
    def courant_coefficient(self, new_val):
        if new_val <= 0:
            raise ValueError("courant_coefficient must be > 0.")
        self._courant_coefficient = new_val

    @property
    def K(self):
        """hydraulic conductivity at link (m/s)"""
        return self._K

    @K.setter
    def K(self, new_val):
        """set hydraulic conductivity at link (m/s)"""
        self._K = new_val

    @property
    def recharge(self):
        """recharge rate (m/s)"""
        return self._recharge

    @recharge.setter
    def recharge(self, new_val):
        """set recharge rate (m/s)"""
        self._recharge = new_val

    @property
    def n(self):
        """porosity of the aquifer (-)"""
        return self._n

    def calc_recharge_flux_in(self):
        """Calculate flux into the domain from recharge.

        Includes recharge that may immediately become saturation excess
        overland flow. (m3/s)
        """
        return np.sum(self._grid.area_of_cell * self._recharge[self._cores])

    def calc_gw_flux_out(self):
        """Groundwater flux through open boundaries may be positive (out of the
        domain) or negative (into the domain).

        This function determines the correct sign for specific discharge
        based upon this convention, and sums the flux across the
        boundary faces. (m3/s)
        """
        # get links at open boundary nodes
        open_nodes = self._grid.open_boundary_nodes
        links_at_open = self._grid.links_at_node[open_nodes]
        link_dirs_at_open = self._grid.active_link_dirs_at_node[open_nodes]

        # find active links at open boundary nodes
        active_links_at_open = links_at_open[link_dirs_at_open != 0]
        active_link_dirs_at_open = link_dirs_at_open[link_dirs_at_open != 0]

        # get groundwater specific discharge at these locations
        q_at_open_links = self._grid.at_link["groundwater__specific_discharge"][
            active_links_at_open
        ]

        # get cell widths at these locations
        faces = self._grid.face_at_link[active_links_at_open]
        face_widths = self._grid.length_of_face[faces]

        # get volume flux out at these locations
        gw_volume_flux_rate_out = (
            q_at_open_links * active_link_dirs_at_open * face_widths
        )

        return np.sum(gw_volume_flux_rate_out)

    def calc_sw_flux_out(self):
        """Surface water flux out of the domain through seepage and saturation
        excess.

        Note that model does not allow for reinfiltration.  (m3/s)
        """
        return np.sum(
            self._grid.at_node["surface_water__discharge"][
                self._grid.open_boundary_nodes
            ]
        )

    def calc_gw_flux_at_node(self):
        """Calculate the sum of the groundwater flux leaving a node.

        (m2/s)
        """
        gw = (
            self._grid.at_link["groundwater__specific_discharge"][
                self._grid.links_at_node
            ]
            * self._grid.link_dirs_at_node
        )
        gw_out = -np.sum(gw * (gw < 0), axis=1)
        return gw_out

        # Old definition of gw flux at node.
        # return map_max_of_node_links_to_node(self._grid,self._grid.dx* abs(self._grid.at_link['groundwater__specific_discharge']))

    def calc_shear_stress_at_node(self, n_manning=0.05):
        """Calculate the shear stress Tau based upon the equations: (N/m2)

            Tau = rho g S d
            d = n q / ( dx S^1/2)^3/5

        where rho is the density of water, g is the gravitational constant,
        S is the slope, d is the water depth calculated with manning's equation,
        n is Manning's n, q is surface water discharge, and dx is the grid cell
        width.

        Parameters
        ----------
        n_manning: float or array of float (-)
            Manning's n at nodes, giving surface roughness.
        """
        rho = 1000  # kg/m3
        g = 9.81  # m/s2
        return (
            rho
            * g
            * self._S_node
            * (
                (n_manning * self._grid.at_node["surface_water__discharge"] / 3600)
                / (self._grid.dx * np.sqrt(self._S_node))
            )
            ** (3 / 5)
        )

    def calc_total_storage(self):
        """calculate the current water storage in the aquifer (m3)"""
        return np.sum(
            self._n[self._cores]
            * self._grid.area_of_cell
            * self._grid.at_node["aquifer__thickness"][self._cores]
        )

    def run_one_step(self, dt):
        """Advance component by one time step of size dt.

        Parameters
        ----------
        dt: float (time in seconds)
            The imposed timestep.
        """

        # Calculate base gradient
        self._base_grad[self._grid.active_links] = self._grid.calc_grad_at_link(
            self._base
        )[self._grid.active_links]

        # Calculate hydraulic gradient
        self._hydr_grad[self._grid.active_links] = self._grid.calc_grad_at_link(
            self._thickness
        )[self._grid.active_links]

        # Calculate groundwater velocity
        self._vel[:] = -self._K * (
            self._hydr_grad * np.cos(np.arctan(abs(self._base_grad)))
            + np.sin(np.arctan(self._base_grad))
        )
        self._vel[self._grid.status_at_link == INACTIVE_LINK] = 0.0

        # Aquifer thickness at links (upwind)
        hlink = map_value_at_max_node_to_link(
            self._grid, "water_table__elevation", "aquifer__thickness"
        )

        # Calculate specific discharge
        self._q[:] = hlink * self._vel

        # Groundwater flux divergence
        dqdx = self._grid.calc_flux_div_at_node(self._q)

        # Determine the relative aquifer thickness, 1 if permeable thickness is 0.
        soil_present = self._elev - self._base > 0.0
        rel_thickness = np.ones_like(self._elev)
        rel_thickness[soil_present] = np.minimum(
            1,
            self._thickness[soil_present]
            / (self._elev[soil_present] - self._base[soil_present]),
        )

        # Calculate surface discharge at nodes
        self._qs[:] = _regularize_G(rel_thickness, self._r) * _regularize_R(
            self._recharge - dqdx
        )

        # Mass balance
        self._dhdt[:] = (1 / self._n) * (self._recharge - dqdx - self._qs)

        # Update
        self._thickness[self._cores] += self._dhdt[self._cores] * dt
        self._thickness[self._thickness < 0] = 0.0

        # Recalculate water surface height
        self._wtable[self._cores] = (
            self._base[self._cores] + self._thickness[self._cores]
        )

    def run_with_adaptive_time_step_solver(self, dt):
        """Advance component by one time step of size dt, subdividing the
        timestep into substeps as necessary to meet a Courant condition
        (defined in the init). Note this method only returns the fluxes at the
        last subtimestep.

        Parameters
        ----------
        dt: float (time in seconds)
            The imposed timestep.
        """

        remaining_time = dt
        self._num_substeps = 0

        while remaining_time > 0.0:
            # Calculate base gradient
            self._base_grad[self._grid.active_links] = self._grid.calc_grad_at_link(
                self._base
            )[self._grid.active_links]

            # Calculate hydraulic gradient
            self._hydr_grad[self._grid.active_links] = self._grid.calc_grad_at_link(
                self._thickness
            )[self._grid.active_links]

            # Calculate groundwater velocity
            self._vel[:] = -self._K * (
                self._hydr_grad * np.cos(np.arctan(abs(self._base_grad)))
                + np.sin(np.arctan(self._base_grad))
            )
            self._vel[self._grid.status_at_link == INACTIVE_LINK] = 0.0

            # Aquifer thickness at links (upwind)
            hlink = map_value_at_max_node_to_link(
                self._grid, "water_table__elevation", "aquifer__thickness"
            )

            # Calculate specific discharge
            self._q[:] = hlink * self._vel

            # Groundwater flux divergence
            dqdx = self._grid.calc_flux_div_at_node(self._q)

            # Determine the relative aquifer thickness, 1 if permeable thickness is 0.
            soil_present = self._elev - self._base > 0.0
            rel_thickness = np.ones_like(self._elev)
            rel_thickness[soil_present] = np.minimum(
                1,
                self._thickness[soil_present]
                / (self._elev[soil_present] - self._base[soil_present]),
            )

            # Calculate surface discharge at nodes
            self._qs[:] = _regularize_G(rel_thickness, self._r) * _regularize_R(
                self._recharge - dqdx
            )

            # Mass balance
            self._dhdt[:] = (1 / self._n) * (self._recharge - dqdx - self._qs)

            # calculate criteria for timestep
            max_vel = max(abs(self._vel / self._n_link))
            grid_dist = min(self._grid.length_of_link)
            substep_dt = np.nanmin(
                [self._courant_coefficient * grid_dist / max_vel, remaining_time]
            )

            # Update
            self._thickness[self._cores] += self._dhdt[self._cores] * substep_dt
            self._thickness[self._thickness < 0] = 0.0

            # Recalculate water surface height
            self._wtable[self._cores] = (
                self._base[self._cores] + self._thickness[self._cores]
            )

            # calculate the time remaining and advance count of substeps
            remaining_time -= substep_dt
            self._num_substeps += 1
