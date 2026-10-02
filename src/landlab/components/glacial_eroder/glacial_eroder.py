#!/usr/bin/env python
"""
Component that simulates the growth, flow, and erosion of glaciers.
"""

import numpy as np

from landlab import Component


class GlacialEroder(Component):
    """Simulates the growth, flow, and erosion of glaciers.

    This component can either simulate erosion due to a user defined
    "fixed ice thickess" at each grid node, or compute an internally
    consistent change in ice thickness based on model dynamics. In
    both cases, the topography is eroded by the glacier using the
    "Glacial Erosion Rule". The glacial erosion rule relates the
    erosion rate to the basal sliding velocity of the glacier through:

    ..math::

        E = k * |u_b|^l

    where E is the erosion rate, k is the glacial erosion coefficient,
    u_b is the basal sliding velocity, and l is the sliding exponent.
    The glacial erosion rule is a simplified representation of the complex
    processes that occur at the base of a glacier, and ignores subglacial
    hydrology and bedrock properties.

    The basal sliding velocity is determined using an empirical relationship
    based on the model basal shear stress and a characteristic basal sliding
    velocity at a characteristic basal shear stress (from Kessler et al., 2006).
    The expression for the basal sliding velocity is given by:

    ..math::

        u_b = u_c * exp(1 - tau_b / tau_c)

    where u_b is the basal sliding velocity, u_c is a characteristic basal
    sliding velocity as characteristic basal stress tau_c, and tau_b is
    the modeled basal shear stress. The basal shear stress is gravitationally driven,
    and is computed as:

    ..math::

        tau_b = rho_ice * g * H * sin(theta)

    where rho_ice is the density of ice, g is the acceleration due to gravity,
    H is the ice thickness, and theta is the surface slope angle.

    When dynamically computing ice thickness, the model updates the ice thickness
    at each time step based on the balance between accumulation, ablation, and net
    Accumulation is determined based on an equilibrium line altitude (ELA). When the
    ice surface is below the ELA, no accumulation occurs. When the ice surface is above
    the ELA, accumulation occurs at a rate proportional to the height above the ELA.
    Beyond the ELA, ablation also occurs via a geothermal heat flux, and due to shear
    heating at the base of the glacier due to basal sliding. The migration of the glacier
    across the topography is determined by computing a depth averaged velocity assuming the
    Shallow Ice Approximation and a Glen-Nye power-law rheology. The velocity is calculated
    using:

    ..math::

        u_ice = (2 * A / (n_exp + 2)) * rho_ice * g *
                (grad(H + elevation))^n_exp * H^(n_exp + 1)

    where u_ice is the depth-averaged ice velocity, A is the ice flow coefficient,
    n_exp is the Glen-Nye flow law exponent, rho_ice is the density of ice, g is the
    acceleration due to gravity, H is the ice thickness, and grad(H + elevation) is
    the surface slope of the ice. When simulating dynamic evolution of the ice
    thickness, the net velocity flux at each node is given by the sum of u_ice and
    the basal sliding velocity u_b.
    """

    _name = "GlacialEroder"

    _unit_agnostic = True

    _info = {
        "dynamic_ice__thickness": {
            "dtype": float,
            "intent": "inout",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Dynamic ice thickness",
        },
        "fixed_ice__thickness": {
            "dtype": float,
            "intent": "inout",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Fixed ice thickness",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
    }

    def __init__(
        self,
        grid,
        rho_ice=917.0,
        latent_heat=3.34e5,
        Q_geo=0.065,
        Q_accum=1.0,
        max_Q_accum=1.0,
        Q_ablat=1.0,
        ice_flow_coefficient=2.4e-24,
        gravity=9.81,
        n_exp=3,
        glacial_erosion_coefficient=2.7e-7,
        sliding_exponent=2.0,
        ela_height=5000,
        ela_thickness_band=500,
        characteristic_u=20.0,
        characteristic_tau=1e5,
        courant=1.0,
    ):
        """Initialize the glacial eroder component.

        Parameters
        ----------
        grid : RasterModelGrid
            The Landlab grid.
        rho_ice : float, optional
            Density of ice (kg/m^3). Default is 917.0.
        latent_heat : float, optional
            Latent heat of fusion for ice (J/kg). Default is 3.34e5.
        Q_geo : float, optional
            Geothermal heat flux (W/m^2). Default is 0.065.
        Q_accum : float, optional
            Accumulation heat flux (W/m^2). Default is 1.0.
        max_Q_accum : float, optional
            Maximum accumulation heat flux (W/m^2). Default is 1.0.
        Q_ablat : float, optional
            Ablation heat flux (W/m^2). Default is 1.0.
        ice_flow_coefficient : float, optional
            Ice flow coefficient (Pa^-3 s^-1). Default is 2.4e-24.
        gravity : float, optional
            Gravitational acceleration (m/s^2). Default is 9.81.
        n_exp : float, optional
            Flow law exponent (dimensionless). Default is 3.
        glacial_erosion_coefficient : float, optional
            Glacial erosion coefficient (m^(1-2n) s^(2n-1)). Default is 2.7e-7.
        sliding_exponent : float, optional
            Sliding exponent (dimensionless). Default is 2.0.
        ela_height : float, ndarray, or str
            Equilibrium line altitude (m). If an array or a field name, this must
            correspond to the ELA at each node. If a float, the same value is applied
            to all nodes.
        ela_thickness_band : float, optional
            Thickness band of the equilibrium line altitude (m). Default is 500.
        characteristic_u : float, optional
            Characteristic ice velocity (m/s). Default is 1.0.
        characteristic_tau : float, optional
            Characteristic basal shear stress (Pa). Default is 1.0.
        courant : float, optional
            Courant number for numerical stability (dimensionless). Default is 0.5.
        """
        super().__init__(grid)

        use_fixed_ice = "fixed_ice__thickness" in self._grid.at_node
        use_dynamic_ice = "dynamic_ice__thickness" in self._grid.at_node

        if not use_fixed_ice and not use_dynamic_ice:
            raise ValueError(
                "Must have either fixed or dynamic ice thickness. "
                "For fixed ice, include 'fixed_ice__thickness' at "
                "nodes, and for dynamic ice, include "
                "'dynamic_ice__thickness' at nodes."
            )

        if use_fixed_ice and use_dynamic_ice:
            raise ValueError(
                "Cannot have both fixed and dynamic ice thickness "
                "at the same time. For fixed ice thickness, include "
                "only 'fixed_ice__thickness' at nodes. For dynamic "
                "ice thickness, include only 'dynamic_ice__thickness' "
                "at nodes."
            )

        self._evolve_ice_thickness = use_dynamic_ice
        self._ice_string = (
            "dynamic_ice__thickness" if use_dynamic_ice else "fixed_ice__thickness"
        )

        self._eps = 1e-10

        self.rho_ice = rho_ice
        self.latent_heat = latent_heat
        self._Q_geo = Q_geo
        self._Q_accum = Q_accum
        self._Q_ablat = Q_ablat
        self._max_Q_accum = max_Q_accum
        self.ice_flow_coefficient = ice_flow_coefficient
        self.gravity = gravity
        self._n_exp = n_exp
        self.glacial_erosion_coefficient = glacial_erosion_coefficient
        self._sliding_exponent = sliding_exponent
        self._ela_height = self._validate_equilibrium_line_altitude(grid, ela_height)
        self._ela_thickness_band = ela_thickness_band
        self._characteristic_u = characteristic_u
        self._characteristic_tau = characteristic_tau
        self.courant = courant

    @staticmethod
    def _validate_equilibrium_line_altitude(grid, ela_height):
        if isinstance(ela_height, str):
            if ela_height in grid.at_node:
                h_ela = grid.at_node[ela_height]
            else:
                raise ValueError(
                    f"ela_height {ela_height!r}, it must be defined " "at either nodes."
                )
        elif np.ndim(ela_height) == 0:
            h_ela = float(ela_height)
        else:
            h_ela = np.asarray(ela_height)
            if h_ela.size != grid.number_of_nodes:
                raise ValueError("ela_height must be defined at nodes.")

        return h_ela

    @property
    def rho_ice(self):
        """Density of ice (kg/m^3)."""
        return self._rho_ice

    @rho_ice.setter
    def rho_ice(self, value):
        if value > 0.0:
            self._rho_ice = float(value)
        else:
            raise ValueError("ice density must be greater than zero.")

    @property
    def latent_heat(self):
        """Latent heat of fusion for ice (J/kg)."""
        return self._latent_heat

    @latent_heat.setter
    def latent_heat(self, value):
        if value > 0.0:
            self._latent_heat = float(value)
        else:
            raise ValueError("latent heat must be greater than zero.")

    @property
    def Q_geo(self):
        """Geothermal heat flux (W/m^2)."""
        return self._Q_geo

    @property
    def Q_accum(self):
        """Accumulation heat flux (W/m^2)."""
        return self._Q_accum

    @property
    def Q_ablat(self):
        """Ablation heat flux (W/m^2)."""
        return self._Q_ablat

    @property
    def max_Q_accum(self):
        """Maximum accumulation heat flux (W/m^2)."""
        return self._max_Q_accum

    @property
    def ice_flow_coefficient(self):
        """Ice flow coefficient (Pa^-3 s^-1)."""
        return self._ice_flow_coefficient

    @ice_flow_coefficient.setter
    def ice_flow_coefficient(self, value):
        if value > 0.0:
            self._ice_flow_coefficient = float(value)
        else:
            raise ValueError("ice flow coefficient must be greater than zero.")

    @property
    def gravity(self):
        """Gravitational acceleration (m/s^2)."""
        return self._gravity

    @gravity.setter
    def gravity(self, value):
        if value > 0.0:
            self._gravity = float(value)
        else:
            raise ValueError("gravity must be greater than zero.")

    @property
    def n_exp(self):
        """Flow law exponent (dimensionless)."""
        return self._n_exp

    @property
    def glacial_erosion_coefficient(self):
        """Glacial erosion coefficient (m^(1-2n) s^(2n-1))."""
        return self._glacial_erosion_coefficient

    @glacial_erosion_coefficient.setter
    def glacial_erosion_coefficient(self, value):
        if value > 0.0:
            self._glacial_erosion_coefficient = float(value)
        else:
            raise ValueError("glacial erosion coefficient must be greater than zero.")

    @property
    def sliding_exponent(self):
        """Sliding exponent (dimensionless)."""
        return self._sliding_exponent

    @property
    def ela_height(self):
        """Equilibrium line altitude (m)."""
        return self._ela_height

    @property
    def ela_thickness_band(self):
        """Equilibrium line altitude thickness band (m)."""
        return self._ela_thickness_band

    @property
    def characteristic_u(self):
        """Characteristic sliding velocity (m/s)."""
        return self._characteristic_u

    @property
    def characteristic_tau(self):
        """Characteristic basal stress (Pa)."""
        return self._characteristic_tau

    @property
    def courant(self):
        """Courant number for CFL condition (dimensionless)."""
        return self._courant

    @courant.setter
    def courant(self, value):
        if value > 0.0 and value <= 1.0:
            self._courant = float(value)
        else:
            raise ValueError(
                "courant number must be greater than zero and less "
                "than or equal to one."
            )

    def run_one_step(self, dt):
        """Run the glacial eroder for one timestep, dt.

        Parameters
        ----------
        dt : float
            The imposed timestep.
        """
        remaining = dt
        while remaining > 0.0:
            deformation_velocity = self._calc_ice_deformation_velocity()
            sliding_velocity = self._calc_ice_sliding_velocity()
            substep = min(
                remaining,
                self._calc_cfl_timestep(sliding_velocity, deformation_velocity),
            )

            if self._evolve_ice_thickness:
                ice_thickness_change = self._calc_ice_thickness_change(
                    deformation_velocity, sliding_velocity
                )

            erosion_rate = self._calc_erosion_rate(sliding_velocity)

            if self._evolve_ice_thickness:
                self._grid.at_node[self._ice_string] += ice_thickness_change * substep

            self._grid.at_node["topographic__elevation"] -= erosion_rate * substep

            self._grid.at_node[self._ice_string] = np.maximum(
                self._grid.at_node[self._ice_string], 0.0
            )

            remaining -= substep

    def _calc_cfl_timestep(self, sliding_velocity, deformation_velocity):
        """Calculate the largest stable explicit timestep for advection of ice.

        Returns
        -------
        float
            The CFL-limited timestep.
        """
        max_velocity = np.max(abs(deformation_velocity + sliding_velocity))
        max_link = np.amax(self._grid.length_of_link)

        if max_velocity <= 0.0:
            return np.inf

        return self._courant * max_link / max_velocity

    def _calc_ice_deformation_velocity(self):
        """Calculate the ice internal deformation velocity assuming Glen's flow law
           and the Shallow Ice Approximation.

        Returns
        -------
        ndarray at links
            The ice deformation velocity at each link.
        """
        ice_thickness = self._grid.at_node[self._ice_string]
        elevation = self._grid.at_node["topographic__elevation"]
        ice_surface = elevation + ice_thickness

        prefactor = (
            2
            * self._ice_flow_coefficient
            / (self._n_exp + 2)
            * self._rho_ice
            * self._gravity**self._n_exp
        )
        ice_surface_gradient = self._grid.calc_grad_at_link(ice_surface)
        ice_thickness_at_links = self._grid.map_value_at_max_node_to_link(
            ice_surface, ice_thickness
        )

        ice_thickness_at_links = np.maximum(ice_thickness_at_links, self._eps)

        deformation_velocity = (
            -prefactor
            * ice_surface_gradient**self._n_exp
            * ice_thickness_at_links ** (self._n_exp + 1)
        )

        deformation_velocity[
            self._grid.status_at_link != self._grid.BC_LINK_IS_ACTIVE
        ] = 0.0

        return deformation_velocity

    def _calc_ice_sliding_velocity(self):
        """Calculate the ice sliding velocity.

        Returns
        -------
        ndarray at links
            The ice sliding velocity at each link.
        """
        basal_driving_stress = self._calc_basal_driving_stress()

        sliding_velocity = np.zeros_like(basal_driving_stress)

        # Check to make sure we are not dividing by zero, which can occur when
        # there is no ice, or when the gradient is zero.
        basal_driving_stress_exists = basal_driving_stress > 0
        sliding_velocity[basal_driving_stress_exists] = (
            -self._characteristic_u
            * np.exp(
                1
                - self._characteristic_tau
                / basal_driving_stress[basal_driving_stress_exists]
            )
        )

        sliding_velocity[self._grid.status_at_link != self._grid.BC_LINK_IS_ACTIVE] = (
            0.0
        )

        return sliding_velocity

    def _calc_basal_driving_stress(self):
        """Calculate the basal driving stress on the ice.

        Returns
        -------
        ndarray at links
            The basal driving stress at each link.
        """
        elevation = self._grid.at_node["topographic__elevation"]
        ice_thickness_at_links = self._grid.map_value_at_max_node_to_link(
            elevation,
            self._grid.at_node[self._ice_string],
        )

        ice_thickness_at_links = np.maximum(ice_thickness_at_links, self._eps)

        basal_driving_stress = (
            self._rho_ice
            * self._gravity
            * ice_thickness_at_links
            * abs(self._grid.calc_grad_at_link(elevation))
        )
        return basal_driving_stress

    def _calc_ice_accumulation_term(self):
        """Calculate the ice accumulation term.

        Returns
        -------
        ndarray at nodes
            The ice accumulation term at each node.
        """

        elevation = self._grid.at_node["topographic__elevation"]
        ice_thickness = self._grid.at_node[self._ice_string]
        ice_surface = elevation + ice_thickness

        ice_accum = np.zeros_like(ice_surface)

        mask = (ice_surface >= self._ela_height) & (
            ice_surface <= self._ela_height + self._ela_thickness_band
        )

        ice_accum[ice_surface < self._ela_height] = -self._Q_ablat
        ice_accum[ice_surface > self._ela_height + self._ela_thickness_band] = (
            self._Q_accum
        )
        ice_accum[mask] = (
            np.clip(
                (ice_surface[mask] - self._ela_height) / self._ela_thickness_band,
                0.0,
                1.0,
            )
            * self._Q_accum
        )

        return np.minimum(ice_accum, self._max_Q_accum)

    def _calc_ice_thickness_change(self, deformation_velocity, sliding_velocity):
        """Calculate the change in ice thickness.
        Parameters
        ----------
        deformation_velocity : ndarray at links
            The ice deformational velocity at each link.
        sliding_velocity : ndarray at links
            The ice sliding velocity at each link.

        Returns
        -------
        ndarray at nodes
            The change in ice thickness at each node.
        """
        ice_accumulation = self._calc_ice_accumulation_term()

        net_ice_velocity = deformation_velocity + sliding_velocity
        ice_flux_at_node = self._grid.calc_flux_div_at_node(net_ice_velocity)

        frictional_heat_flux = self._calc_basal_driving_stress() * sliding_velocity

        frictional_heat_flux_at_node = abs(
            self._grid.map_upwind_node_link_mean_to_node(frictional_heat_flux)
        )

        ice_change = (
            -(self._Q_geo + frictional_heat_flux_at_node)
            / (self._rho_ice * self._latent_heat)
            + ice_accumulation
            - ice_flux_at_node
        )
        return ice_change

    def _calc_erosion_rate(self, sliding_velocity):
        """Calculate the erosion rate.

        Parameters
        ----------
        sliding_velocity : ndarray at links
            The ice sliding velocity at each link.

        Returns
        -------
        ndarray at nodes
            The erosion rate at each node.
        """

        # Map the mean ice sliding velocity determined at links to the downwind nodes.
        # This will be used to determine the glacial erosion rate, so taking the max value
        # will result in the largest amount of erosion.
        ice_velocity_at_node = abs(
            self._grid.map_upwind_node_link_mean_to_node(sliding_velocity)
        )
        erosion_rate = self._glacial_erosion_coefficient * (
            ice_velocity_at_node**self._sliding_exponent
        )

        return erosion_rate
