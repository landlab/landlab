# -*- coding: utf-8 -*-
"""
GroundwaterDupuitPercolator Component

@author: G Tucker, D Litwin, K Barnhart
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
    """ Smooths transition of step function with an exponential. 0<=u<=1. """
    return np.exp(-(1 - u) / reg_factor)


def _regularize_R(u):
    """ ramp function on u """
    return u * np.greater_equal(u, 0)


def get_link_hydraulic_conductivity(grid, K):
    """
    Returns array of hydraulic conductivity on links, allowing
    for aquifers with laterally anisotropic hydraulic conductivity

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
    r"""
    Simulate groundwater flow in a shallow unconfined aquifer.

    The GroundwaterDupuitPercolator solves the Boussinesq equation for
    flow in an unconfined aquifer over an impermeable aquifer base and
    calculates groundwater return flow to the surface. This method uses the
    Dupuit approximation that the hydraulic gradient is zero in the direction
    normal to the aquifer base. For this reason, aquifer thickness
    and regolith thickness are calculated in the direction normal to the
    aquifer base. The hydraulic gradient is also calculated relative
    to the aquifer base. For more details, see component documentation
    :ref:`here <dupuit_theory>`.

    Parameters
    ----------
    grid: ModelGrid
            Landlab ModelGrid object
    hydraulic_conductivity: float, field name, or array of float
            saturated hydraulic conductivity, m/s.
            Default = 0.001 m/s
    porosity: float, field name or array of float
            the porosity of the aquifer [-].
            Default = 0.2
    recharge_rate: float, field name, or array of float
            Rate of recharge, m/s.
            Default = 1.0e-8 m/s
    regularization_f: float
            factor controlling the smoothness of the transition between
            surface and subsurface flow.
            Default = 1e-2

    Examples
    --------
    Import the grid class and component

    >>> from landlab import RasterModelGrid
    >>> from landlab.components import GroundwaterDupuitPercolator

    Initialize the grid and component

    >>> grid = RasterModelGrid((10, 10), xy_spacing=10.0)
    >>> elev = grid.add_zeros('node', 'topographic__elevation')
    >>> elev[:] = 5.0
    >>> gdp = GroundwaterDupuitPercolator(grid)

    Run component forward. Note all time units in the model are in seconds.

    >>> dt = 1E4
    >>> for i in range(100):
    ...     gdp.run_one_step(dt)

    When the model generates groundwater return flow, the surface water flux
    out of the domain can be calculated only after a FlowAccumulator is run.
    Below is a more complex model that demonstrates this case.

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
    Below is a summary of the theory and numerical implementation of
    the ``GroundwaterDupuitPercolator``. A complete description can be found
    :ref:`here <dupuit_theory>`.

    Groundwater discharge per unit length, :math:`q`, is calculated as:

    .. math::
        q = -K_{sat} h \big( \nabla z \big) \cos^2 (\alpha)

    where :math:`K_{sat}` is the saturated hydraulic conductivity, :math:`h` is
    the aquifer thickness, and :math:`\alpha` is the slope angle of the aquifer base.

    Surface water discharge per unit area, :math:`q_s`, is calculated as:

    .. math::
        q_s = \mathcal{G}_r \bigg( \frac{h}{d} \bigg) \mathcal{R} \big(-\cos(\alpha) \nabla \cdot q + f \big)

    where :math:`\mathcal{G}_r` is a smoothed step function, :math:`\mathcal{R}` is the ramp function,
    :math:`d` is the regolith thickness, and :math:`f` is the recharge rate.

    The evolution of aquifer thickness is then given by:

    .. math::
        n \frac{\partial h}{\partial t} = f - q_s - \nabla \cdot q

    where :math:`n` is the drainable porosity.

    An explicit forward in time finite volume method is used to implement a
    numerical solution. Flow discharge between neighboring nodes is calculated
    using the saturated thickness at the up-gradient node.

    """

    _name = "GroundwaterDupuitPercolator"

    _input_var_names = set(("topographic__elevation", "aquifer_base__elevation"))

    _output_var_names = set(
        (
            "aquifer__thickness",
            "water_table__elevation",
            "aquifer_base__gradient",
            "hydraulic__gradient",
            "groundwater__specific_discharge",
            "groundwater__velocity",
            "surface_water__specific_discharge",
            "average_surface_water__specific_discharge",
            "",
        )
    )

    _var_units = {
        "topographic__elevation": "m",
        "aquifer_base__elevation": "m",
        "aquifer__thickness": "m",
        "aquifer_base__gradient": "m/m",
        "water_table__elevation": "m",
        "hydraulic__gradient": "m/m",
        "groundwater__specific_discharge": "m2/s",
        "groundwater__velocity": "m/s",
        "surface_water__specific_discharge": "m/s",
        "average_surface_water__specific_discharge": "m/s",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "aquifer_base__elevation": "node",
        "aquifer__thickness": "node",
        "water_table__elevation": "node",
        "hydraulic__gradient": "link",
        "aquifer_base__gradient": "link",
        "groundwater__specific_discharge": "link",
        "groundwater__velocity": "link",
        "surface_water__specific_discharge": "node",
        "average_surface_water__specific_discharge": "node",
    }

    _var_doc = {
        "topographic__elevation": "elevation of land surface",
        "aquifer_base__elevation": "elevation of impermeable base",
        "aquifer__thickness": "vertical thickness of saturated zone",
        "water_table__elevation": "elevation of water table",
        "hydraulic__gradient": "gradient of water table in link direction",
        "aquifer_base__gradient": "gradient of the aquifer base in the link direction",
        "groundwater__specific_discharge": "discharge per width in link dir",
        "groundwater__velocity": "groundwater darcy flux in link direction",
        "surface_water__specific_discharge": "rate of groundwater return flow plus recharge on saturated area",
        "average_surface_water__specific_discharge": "average surface water specific discharge over variable timesteps",
    }

    def __init__(
        self,
        grid,
        hydraulic_conductivity=0.001,
        porosity=0.2,
        recharge_rate=1.0e-8,
        regularization_f=1e-2,
    ):
        """

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
        """
        # Store grid
        self._grid = grid

        # Shorthand
        self._cores = grid.core_nodes

        # Create fields:
        if "topographic__elevation" in self.grid.at_node:
            self._elev = self.grid.at_node["topographic__elevation"]
        else:
            self._elev = self.grid.add_ones("node", "topographic__elevation")

        if "aquifer_base__elevation" in self.grid.at_node:
            self._base = self.grid.at_node["aquifer_base__elevation"]
        else:
            self._base = self.grid.add_zeros("node", "aquifer_base__elevation")

        if "aquifer_base__gradient" in self.grid.at_link:
            self._base_grad = self.grid.at_link["aquifer_base__gradient"]
        else:
            self._base_grad = self.grid.add_zeros("link", "aquifer_base__gradient")
            self._base_grad[:] = self._grid.calc_grad_at_link(self._base)

        if "water_table__elevation" in self.grid.at_node:
            self._wtable = self.grid.at_node["water_table__elevation"]
        else:
            self._wtable = self.grid.add_zeros("node", "water_table__elevation")
        self._wtable[grid.closed_boundary_nodes] = 0

        if "aquifer__thickness" in self.grid.at_node:
            self._thickness = self.grid.at_node["aquifer__thickness"]
        else:
            self._thickness = self.grid.add_zeros("node", "aquifer__thickness")
            self._thickness[:] = self._wtable - self._base
        self._thickness[grid.closed_boundary_nodes] = 0

        if "hydraulic__gradient" in self.grid.at_link:
            self._hydr_grad = self.grid.at_link["hydraulic__gradient"]
        else:
            self._hydr_grad = self.grid.add_zeros("link", "hydraulic__gradient")

        if "groundwater__specific_discharge" in self.grid.at_link:
            self._q = self.grid.at_link["groundwater__specific_discharge"]
        else:
            self._q = self.grid.add_zeros("link", "groundwater__specific_discharge")

        if "groundwater__velocity" in self.grid.at_link:
            self._vel = self.grid.at_link["groundwater__velocity"]
        else:
            self._vel = self.grid.add_zeros("link", "groundwater__velocity")

        if "surface_water__specific_discharge" in self.grid.at_node:
            self._qs = self.grid.at_node["surface_water__specific_discharge"]
        else:
            self._qs = self.grid.add_zeros("node", "surface_water__specific_discharge")

        if "average_surface_water__specific_discharge" in self.grid.at_node:
            self._qsavg = self.grid.at_node["average_surface_water__specific_discharge"]
        else:
            self._qsavg = self.grid.add_zeros(
                "node", "average_surface_water__specific_discharge"
            )

        # Convert parameters to fields if needed, and store a reference
        self._K = return_array_at_link(grid, hydraulic_conductivity)
        self._recharge = return_array_at_node(grid, recharge_rate)
        self._n = return_array_at_node(grid, porosity)
        self._n_link = map_mean_of_link_nodes_to_link(self._grid, self._n)
        self._r = regularization_f

    @property
    def K(self):
        """hydraulic conductivity at link (m/s)"""
        return self._K

    @K.setter
    def K(self, new_val):
        """set hydraulic conductivity at link (m/s)"""
        self._K = return_array_at_link(self._grid, new_val)

    @property
    def recharge(self):
        """recharge rate (m/s)"""
        return self._recharge

    @recharge.setter
    def recharge(self, new_val):
        """set recharge rate (m/s)"""
        self._recharge = return_array_at_node(self._grid, new_val)

    @property
    def n(self):
        """porosity of the aquifer (-)"""
        return self._n

    @n.setter
    def n(self, new_val):
        """set aquifer porosity"""
        self._n = return_array_at_node(self._grid, new_val)
        self._n_link = map_mean_of_link_nodes_to_link(self._grid, self._n)

    @property
    def number_of_substeps(self):
        """
        The numer of substeps used by the run_with_adaptive_time_step_solver
        method in the latest method call.
        """
        if self._num_substeps:
            return self._num_substeps
        else:
            print("The method run_with_adaptive_time_step_solver has not been used")

        return self._num_substeps

    def calc_recharge_flux_in(self):
        """
        Calculate flux into the domain from recharge. Includes recharge that
        may immediately become saturation excess overland flow. (m3/s)
        """
        return np.sum(self._grid.area_of_cell * self._recharge[self._cores])

    def calc_gw_flux_out(self):
        """
        Groundwater flux through open boundaries may be positive (out of
        the domain) or negative (into the domain). This function determines the
        correct sign for specific discharge based upon this convention,
        and sums the flux across the boundary faces. (m3/s)
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
        face_widths = self._grid.width_of_face[faces]

        # get volume flux out at these locations
        gw_volume_flux_rate_out = (
            q_at_open_links * active_link_dirs_at_open * face_widths
        )

        return np.sum(gw_volume_flux_rate_out)

    def calc_sw_flux_out(self):
        """
        Surface water flux out of the domain through seepage and saturation excess.
        Note that model does not allow for reinfiltration.  (m3/s)
        """
        return np.sum(
            self._grid.at_node["surface_water__discharge"][
                self._grid.open_boundary_nodes
            ]
        )

    def calc_gw_flux_at_node(self):
        """
        Calculate the sum of the groundwater flux leaving a node. (m2/s)
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
        r"""
        Calculate the shear stress :math:`\tau` based upon the equations: (N/m2)

        .. math::
            \tau = \rho g S d

        .. math::
            d = \bigg( \frac{n Q}{S^{1/2} dx} \bigg)^{3/2}

        where :math:`\rho` is the density of water, :math:`g` is the gravitational constant,
        :math:`S` is the topographic slope, :math:`d` is the water depth calculated with Manning's equation,
        :math:`n` is Manning's n, :math:`q` is surface water discharge, and :math:`dx` is the grid cell
        width.

        Parameters
        ----------
        n_manning: float or array of float (-)
            Manning's n at nodes, giving surface roughness.

        """
        self._S = abs(self._grid.calc_grad_at_link(self._elev))
        self._S_node = map_max_of_node_links_to_node(self._grid, self._S)
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
        """
        calculate the current water storage in the aquifer (m3)
        """
        return np.sum(
            self._n[self._cores]
            * self._grid.area_of_cell
            * self._grid.at_node["aquifer__thickness"][self._cores]
        )

    def run_one_step(self, dt, **kwds):
        """
        Advance component by one time step of size dt.

        Parameters
        ----------
        dt: float (time in seconds)
            The imposed timestep.
        """

        if (self._wtable > self._elev).any():
            self._wtable[self._wtable > self._elev] = self._elev[
                self._wtable > self._elev
            ]
            self._thickness[self._cores] = (
                self._wtable[self._cores] - self._base[self._cores]
            )

        # Calculate base gradient
        self._base_grad[self._grid.active_links] = self._grid.calc_grad_at_link(
            self._base
        )[self._grid.active_links]
        cosa = np.cos(np.arctan(self._base_grad))

        # Calculate hydraulic gradient
        self._hydr_grad[self._grid.active_links] = (
            self._grid.calc_grad_at_link(self._wtable)[self._grid.active_links]
            * cosa[self._grid.active_links]
        )

        # Calculate groundwater velocity
        self._vel[:] = -self._K * self._hydr_grad
        self._vel[self._grid.status_at_link == INACTIVE_LINK] = 0.0

        # Aquifer thickness at links (upwind)
        hlink = (
            map_value_at_max_node_to_link(
                self._grid, "water_table__elevation", "aquifer__thickness"
            )
            * cosa
        )

        # Calculate specific discharge
        self._q[:] = hlink * self._vel

        # Groundwater flux divergence
        dqdx = self._grid.calc_flux_div_at_node(self._q)

        # Determine the relative aquifer thickness, 1 if permeable thickness is 0.
        soil_present = (self._elev - self._base) > 0.0
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
        dhdt = (1 / self._n) * (self._recharge - self._qs - dqdx)

        # Update
        self._thickness[self._cores] += dhdt[self._cores] * dt
        self._thickness[self._thickness < 0] = 0.0

        # Recalculate water surface height
        self._wtable[self._cores] = (self._base + self._thickness)[self._cores]

    def run_with_adaptive_time_step_solver(self, dt, courant_coefficient=0.1, **kwds):
        """
        Advance component by one time step of size dt, subdividing the timestep
        into substeps as necessary to meet a Courant condition.
        Note this method only returns the fluxes at the last subtimestep.

        Parameters
        ----------
        dt: float (time in seconds)
            The imposed timestep.
        courant_coefficient: float (-)
            The muliplying factor on the condition that the timestep is
            smaller than the minimum link length over linear groundwater flow velocity
        """

        # check water table above surface
        if (self._wtable > self._elev).any():
            self._wtable[self._wtable > self._elev] = self._elev[
                self._wtable > self._elev
            ]
            self._thickness[self._cores] = (
                self._wtable[self._cores] - self._base[self._cores]
            )

        # Calculate base gradient
        self._base_grad[self._grid.active_links] = self._grid.calc_grad_at_link(
            self._base
        )[self._grid.active_links]
        cosa = np.cos(np.arctan(self._base_grad))

        # Initialize reg_thickness, rel_thickness
        reg_thickness = self._elev - self._base
        soil_present = reg_thickness > 0.0
        rel_thickness = np.ones_like(self._elev)

        # Initialize for average surface discharge
        qs_cumulative = np.zeros_like(self._elev)

        # Initialize variable timestep
        remaining_time = dt
        self._num_substeps = 0

        while remaining_time > 0.0:

            # Calculate hydraulic gradient
            self._hydr_grad[self._grid.active_links] = (
                self._grid.calc_grad_at_link(self._wtable)[self._grid.active_links]
                * cosa[self._grid.active_links]
            )

            # Calculate groundwater velocity
            self._vel[:] = -self._K * self._hydr_grad
            self._vel[self._grid.status_at_link == INACTIVE_LINK] = 0.0

            # Aquifer thickness at links (upwind)
            hlink = (
                map_value_at_max_node_to_link(
                    self._grid, "water_table__elevation", "aquifer__thickness"
                )
                * cosa
            )

            # Calculate specific discharge
            self._q[:] = hlink * self._vel

            # Groundwater flux divergence
            dqdx = self._grid.calc_flux_div_at_node(self._q)

            # calculate relative thickness
            rel_thickness[soil_present] = np.minimum(
                1, self._thickness[soil_present] / (reg_thickness[soil_present])
            )

            # Calculate surface discharge at nodes
            self._qs[:] = _regularize_G(rel_thickness, self._r) * _regularize_R(
                self._recharge - dqdx
            )

            # Mass balance
            dhdt = (1 / self._n) * (self._recharge - self._qs - dqdx)

            # calculate criteria for timestep
            max_vel = max(abs(self._vel / self._n_link))
            grid_dist = min(self._grid.length_of_link)
            substep_dt = np.nanmin(
                [courant_coefficient * grid_dist / max_vel, remaining_time]
            )

            # Update
            self._thickness[self._cores] += dhdt[self._cores] * substep_dt
            self._thickness[self._thickness < 0] = 0.0

            # Recalculate water surface height
            self._wtable[self._cores] = (self._base + self._thickness)[self._cores]

            # add cumulative sw discharge in substeps
            qs_cumulative += self._qs * substep_dt

            # calculate the time remaining and advance count of substeps
            remaining_time -= substep_dt
            self._num_substeps += 1

        self._qsavg[:] = qs_cumulative / dt
