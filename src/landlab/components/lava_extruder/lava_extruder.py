#!/usr/bin/env python
"""
Component that simulates the effusion of lava from a vent and the
subsequent routing and freezing of the lava over the topography.
"""

import numpy as np
import scipy.constants
from requireit import require_between
from requireit import require_greater_than
from requireit import require_nonnegative
from requireit import require_positive

import landlab.components.lava_extruder.lava_helper_functions as lhf
from landlab import Component


class LavaExtruder(Component):
    """Effusion, routing, and freezing of lava over topography.

    This component takes the coordinates of a vent (or vents)
    and simulates effusion of lava from the vent(s) over the topography.
    The component advects the lava over the topography with a Bingham
    visco-plastic rheology, and simulates the loss of thermal energy
    as the lava cools and eventually freezes.

    A Bingham rheology relates the shear stress to the strain rate and
    a yield stress through the expression:

    .. math::

        strain_rate = (shear_stress - yield_stress) / lava_viscosity

    For shear stresses below the yield stress, the lava does not flow because
    the strain rate is 0. The shear stress acting on the lava is gravitational,
    and is equal to:

    .. math::

        shear_stress = lava_density * gravity * lava_thickness * tan(theta)

    where :math:`\theta` is the slope of the topography from the horizontal.
    The strain rate can be related to the velocity of the lava flow through:

    .. math::

        strain_rate = d(velocity) / d(z)

    where :math:`z` is the vertical coordinate measured from the base of the lava
    layer. z = 0 at the base of the lava layer and H at the top of the lava layer.
    Plugging in the strain rate and shear stress expressions into the constitutive
    equation, and rearranging the terms we end up with:

    .. math::

        d(velocity) = (lava_density * gravity * lava_thickness * tan(theta) -
                       yield_stress) / lava_viscosity * d(z)

    which can be integrated to obtain the velocity profile within the lava layer. The
    boundary condition which allows us to determine the constant of integration is
    that at z = 0, the velocity is zero (no slip condition). After integrating and
    averaging the velocity over the thickness of the lava, we end up with the velocity
    expression used to advect the lava in this model:

    .. math::

        velocity = lava_thickness / lava_viscosity *
                   (lava_density * gravity * lava_thickness * tan(theta) / 3 -
                    yield_stress / 2)

    To simulate cooling and freezing of the lava, the component accounts for three
    sources of energy loss from the lava flow:

    1. Radiative heat loss from the surface of the lava when no "lid" is present.
    2. Conductive heat loss to the underlying topography.
    3. Conductive heat loss through a "lid" that forms on top of the lava as it cools.

    Lava is erupted with a sensible heat content that is equal to the eruption
    temperature minus the freezing temperature. While the lava is above the freezing
    temperature, energy is strictly removed from the sensible heat content of the lava.
    When the sensible heat content is zero, energy is removed from the lava through
    freezing, which is proportional to the latent heat of freezing. All cooling and
    freezing calculations assume a bulk energy/thermal balance throughout the entire
    lava layer at each point on the grid. As long as lava remains underneath a lid,
    the lid will continue to insulate the lava from radiative heat and conductive
    heat loss. Once lava is fully frozen, the lid thickness is added to the topography.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import LavaExtruder
    >>> dt = 30.0
    >>> time_to_run = 600.0

    >>> mg = RasterModelGrid((3, 10), xy_spacing=(100.0, 100.0))

    >>> mg.add_zeros("topographic__elevation", at="node")
    >>> mg.add_zeros("lava__thickness", at="node")
    >>> mg.at_node["topographic__elevation"][:] += mg.x_of_node

    >>> mg.set_fixed_value_boundaries_at_grid_edges(True, True, True, True)

    >>> lf = LavaExtruder(
    ...     mg,
    ...     vent_coordinates=[[800, 100]],
    ...     eruptive_volume_flux=10.0,
    ...     viscosity=100.0,
    ...     specific_heat_capacity=300,
    ... )

    >>> elapsed_time = 0.0
    >>> while elapsed_time < time_to_run:
    ...     if elapsed_time + dt > time_to_run:
    ...         dt = time_to_run - elapsed_time
    ...     lf.run_one_step(dt)
    ...     elapsed_time += dt
    ...

    >>> print(np.round(mg.at_node["lava__thickness"][12:20], 3))
    array([0.,    0.,    0.092, 0.127, 0.127, 0.127, 0.127, 0.])

    """

    _name = "LavaExtruder"

    _unit_agnostic = False

    _info = {
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "lava__thickness": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Thickness of the lava layer",
        },
        "lava__temperature": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "K",
            "mapping": "node",
            "doc": "Temperature of the lava layer",
        },
    }

    def __init__(
        self,
        grid,
        *,
        vent_coordinates=None,
        lava_density=2700.0,
        emissivity=0.9,
        specific_heat_capacity=1000.0,
        thermal_conductivity=2.0,
        viscosity=10000.0,
        latent_heat_of_freezing=4.0e5,
        eruptive_volume_flux=0.1,
        gravity=scipy.constants.g,
        yield_stress=1000.0,
        eruption_temperature=1473.0,
        freezing_temperature=800.0,
        ambient_temperature=300.0,
        max_lava_thickness=10.0,
        basal_thickness=20.0,
        max_lid_thickness=5.0,
        min_lid_thickness=1.0e-2,
        courant=1.0,
        min_dt=0.0,
    ):
        """Initialize the LavaExtruder component.

        Parameters
        ----------
        grid : RasterModelGrid
            A grid.
        lava_density : float, optional
            Density of the lava (kg / m^3).
        emissivity : float, optional
            Emissivity of the lava surface (dimensionless).
        specific_heat_capacity : float, optional
            Specific heat capacity of the lava (J / (kg*K)).
        thermal_conductivity : float, optional
            Thermal conductivity of the lava (W / (m*K))
        viscosity : float, optional
            Viscosity of the lava (Pa*s).
        latent_heat_of_freezing : float, optional
            Latent heat of freezing (J / kg).
        eruptive_volume_flux : float, optional
            Volume flux of lava (m^3 / s).
        vent_coordinates : array_like, optional
            X- and Y-coordinates of vent nodes.
        gravity : float, optional
            Acceleration due to gravity (m / s^2).
        yield_stress : float, optional
            Yield stress of the lava (Pa).
        eruption_temperature : float, optional
            Eruption temperature of the lava (K).
        freezing_temperature : float, optional
            Freezing temperature of the lava (K).
        max_lava_thickness : float, optional
            Maximum thickness of the lava flow (m).
        basal_thickness : float, optional
            Thickness of the basal layer (m).
        max_lid_thickness : float, optional
            Maximum thickness of the lid layer (m).
        courant : float, optional
            Courant number for the time-stepping scheme.
        min_dt : float, optional
            Minimum allowed time step (s).
        """

        super().__init__(grid)
        super().initialize_output_fields()

        self._lava_density = require_positive(lava_density, name="lava_density")
        self._emissivity = require_between(
            emissivity, a_min=0.0, a_max=1.0, name="emissivity"
        )
        self._viscosity = require_positive(viscosity, name="viscosity")
        self._specific_heat_capacity = require_positive(
            specific_heat_capacity, name="specific_heat_capacity"
        )
        self._thermal_conductivity = require_nonnegative(
            thermal_conductivity, name="thermal_conductivity"
        )
        self._latent_heat_of_freezing = require_positive(
            latent_heat_of_freezing, name="latent_heat_of_freezing"
        )
        self._eruptive_volume_flux = require_nonnegative(
            eruptive_volume_flux, name="eruptive_volume_flux"
        )
        self._gravity = require_positive(gravity, name="gravity")
        self._yield_stress = require_nonnegative(yield_stress, name="yield_stress")
        self._ambient_temperature = require_nonnegative(
            ambient_temperature, name="ambient_temperature"
        )
        self._freezing_temperature = require_greater_than(
            freezing_temperature, self._ambient_temperature, name="freezing_temperature"
        )
        self._eruption_temperature = require_greater_than(
            eruption_temperature,
            self._freezing_temperature,
            name="eruption_temperature",
        )
        self._max_lava_thickness = require_positive(
            max_lava_thickness, name="max_lava_thickness"
        )
        self._basal_thickness = require_positive(
            basal_thickness, name="basal_thickness"
        )
        self._min_lid_thickness = require_positive(
            min_lid_thickness, name="min_lid_thickness"
        )
        self._max_lid_thickness = require_greater_than(
            max_lid_thickness, self._min_lid_thickness, name="max_lid_thickness"
        )
        self._courant = require_between(
            courant, 0.0, 1.0, inclusive_min=False, inclusive_max=True, name="courant"
        )
        self._min_dt = require_nonnegative(min_dt, name="min_dt")

        # Set the vent nodes based on the provided coordinates
        self._vent_nodes = []
        if vent_coordinates is not None:
            for vent_coordinate in vent_coordinates:
                x, y = vent_coordinate
                if not (
                    self._grid.x_of_node.min() <= x <= self._grid.x_of_node.max()
                    and self._grid.y_of_node.min() <= y <= self._grid.y_of_node.max()
                ):
                    raise ValueError("Vent coordinates must be within the domain.")
                self._vent_nodes.append(self._grid.find_nearest_node(vent_coordinate))

        # Minimum epsilon for determining when lava thickness is effectively zero.
        self._dry_thickness_threshold = 1e-14

        self.freezing_rate = np.zeros(self._grid.number_of_nodes)
        self.lid_thickness = np.zeros(self._grid.number_of_nodes)
        self.lava_thermal_energy = np.zeros(self._grid.number_of_nodes)
        self.vent_flux = np.zeros(self._grid.number_of_nodes)

        nodes_with_lava_initially = self._grid.at_node["lava__thickness"] > 0
        if np.any(nodes_with_lava_initially):
            if np.any(
                self._grid.at_node["lava__temperature"][nodes_with_lava_initially]
                < self._freezing_temperature
            ):
                raise ValueError(
                    "You are setting a non-zero lava thickness without setting "
                    "an initial temperature that is at or above the lava freezing "
                    "temperature. Please prescribe an initial temperature that is "
                    "equal to or greater than the freezing temperature everywhere "
                    "that lava is present."
                )
            self.lava_thermal_energy[nodes_with_lava_initially] = (
                self._grid.at_node["lava__thickness"][nodes_with_lava_initially]
                * (
                    self._grid.at_node["lava__temperature"][nodes_with_lava_initially]
                    - self._freezing_temperature
                )
                * self._specific_heat_capacity
                * self._lava_density
                * self._grid.cell_area_at_node[nodes_with_lava_initially]
            )

    @property
    def lava_density(self):
        """Density of the lava (kg/m^3)."""
        return self._lava_density

    @property
    def emissivity(self):
        """Emissivity of the lava surface (dimensionless)."""
        return self._emissivity

    @property
    def viscosity(self):
        """Viscosity of the lava (Pa*s)."""
        return self._viscosity

    @property
    def specific_heat_capacity(self):
        """Specific heat capacity of the lava (J/(kg*K))."""
        return self._specific_heat_capacity

    @property
    def thermal_conductivity(self):
        """Thermal conductivity of the lava (W/(m*K))."""
        return self._thermal_conductivity

    @property
    def latent_heat_of_freezing(self):
        """Latent heat of freezing (J/kg)."""
        return self._latent_heat_of_freezing

    @property
    def eruptive_volume_flux(self):
        """Volume flux of lava (m^3/s)."""
        return self._eruptive_volume_flux

    @property
    def gravity(self):
        """Acceleration due to gravity (m/s^2)."""
        return self._gravity

    @property
    def yield_stress(self):
        """Yield stress of the lava (Pa)."""
        return self._yield_stress

    @property
    def ambient_temperature(self):
        """Surface temperature (K)."""
        return self._ambient_temperature

    @property
    def freezing_temperature(self):
        """Freezing temperature of the lava (K)."""
        return self._freezing_temperature

    @property
    def eruption_temperature(self):
        """Eruption temperature of the lava (K)."""
        return self._eruption_temperature

    @property
    def max_lava_thickness(self):
        """Maximum thickness of the lava flow (m)."""
        return self._max_lava_thickness

    @property
    def basal_thickness(self):
        """Thickness of the basal layer (m)."""
        return self._basal_thickness

    @property
    def min_lid_thickness(self):
        """Minimum thickness of the lid layer (m)."""
        return self._min_lid_thickness

    @property
    def max_lid_thickness(self):
        """Maximum thickness of the lid layer (m)."""
        return self._max_lid_thickness

    @property
    def courant(self):
        """Courant number for the simulation."""
        return self._courant

    @property
    def min_dt(self):
        """Minimum allowable time step for the simulation."""
        return self._min_dt

    def _calc_lava_velocity(self):
        """Calculate lava velocity."""
        lava_thickness = self._grid.at_node["lava__thickness"]
        lava_surface = self._grid.at_node["topographic__elevation"] + lava_thickness

        lava_surface_grad = self._grid.calc_grad_at_link(lava_surface)
        lava_at_link = self._grid.map_value_at_max_node_to_link(
            lava_surface, lava_thickness
        )

        lava_speed = lhf.calc_lava_speed(
            lava_at_link,
            lava_surface_grad,
            lava_density=self._lava_density,
            gravity=self._gravity,
            viscosity=self._viscosity,
            yield_stress=self._yield_stress,
        )

        lava_speed[self._grid.status_at_link != self._grid.BC_LINK_IS_ACTIVE] = 0.0
        lava_flux = -np.sign(lava_surface_grad) * lava_at_link * lava_speed

        return lava_speed, lava_flux

    def _calc_CFL_advection_timestep(self, lava_speed, dt, sub_dt):
        """Calculate the CFL-limited advection timestep.

        Parameters
        ----------
        lava_speed : ndarray
            The magnitude of the lava velocity at each link.
        dt : float
            The current time step.
        sub_dt : float
            The accumulated sub-time step.
        """

        # Apply the CFL condition to each link with a non-zero velocity, and take the
        # smallest timestep from these.
        link_dt = np.full(self._grid.number_of_links, np.inf)
        current_dt = np.min(
            self._courant
            * np.divide(
                self._grid.length_of_link,
                lava_speed,
                where=lava_speed > 0.0,
                out=link_dt,
            )
        )

        return min(current_dt, dt - sub_dt)

    def _calc_CFL_freezing_timestep(self, conductive_cooling_rate, dt, sub_dt):
        """Calculate the CFL-limited freezing timestep.

        Parameters
        ----------
        conductive_cooling_rate : ndarray
            The conductive cooling rate at each node.
        dt : float
            The current time step.
        sub_dt : float
            The accumulated sub-time step.
        """
        lava_thickness = self._grid.at_node["lava__thickness"]
        cell_areas = self._grid.area_of_cell[self._grid.cell_at_node]

        freezing_nodes = (
            (lava_thickness > self._dry_thickness_threshold)
            & (self.lava_thermal_energy <= 0.0)
            & (conductive_cooling_rate < 0.0)
        )
        freezing_nodes[self._vent_nodes] = False

        current_dt = dt

        if np.any(freezing_nodes):
            available_energy_for_freezing = (
                lava_thickness[freezing_nodes]
                * self._lava_density
                * self._latent_heat_of_freezing
                * cell_areas[freezing_nodes]
            )
            freezing_dt = (
                self._courant
                * available_energy_for_freezing
                / np.abs(conductive_cooling_rate[freezing_nodes])
            )
            current_dt = np.clip(np.min(freezing_dt), self._min_dt, current_dt)

        return min(current_dt, dt - sub_dt)

    def _calc_CFL_radiative_timestep(self, radiative_cooling_rate, dt, sub_dt):
        """Calculate the CFL-limited radiative cooling timestep.

        Parameters
        ----------
        radiative_cooling_rate : ndarray
            The radiative cooling rate at each node.
        dt : float
            The current time step.
        sub_dt : float
            The accumulated sub-time step.
        """
        radiate_thermal_energy = (
            (self.lava_thermal_energy > 0.0)
            & (radiative_cooling_rate < 0.0)
            & (self.lid_thickness == self._min_lid_thickness)
        )

        radiate_thermal_energy[self._vent_nodes] = False

        current_dt = dt

        if np.any(radiate_thermal_energy):
            radiative_dt = np.min(
                self.lava_thermal_energy[radiate_thermal_energy]
                / np.abs(radiative_cooling_rate[radiate_thermal_energy])
            )
            current_dt = np.clip(self._courant * radiative_dt, self._min_dt, current_dt)

        return min(current_dt, dt - sub_dt)

    def _calc_lava_temperature(
        self, lava_thermal_energy, lava_thickness, cell_areas, mask, lava_temperature
    ):

        lava_temperature[mask] = np.minimum(
            lhf.thermal_energy_to_temperature(
                lava_thermal_energy[mask],
                lava_thickness[mask],
                cell_areas[mask],
                lava_density=self._lava_density,
                specific_heat_capacity=self._specific_heat_capacity,
            )
            + self._freezing_temperature,
            self._eruption_temperature,
        )
        lava_temperature[self._vent_nodes] = self._eruption_temperature
        lava_temperature[:] = np.maximum(lava_temperature, self._ambient_temperature)

        return lava_temperature

    def _limit_lava_flux(self, lava_flux, dt):
        """Scale outgoing link fluxes so no node loses more lava than it holds.

        Parameters
        ----------
        lava_flux : ndarray
            The flux of lava at each link.
        dt : float
            The time step over which the flux is applied.

        Returns
        -------
        ndarray
            The limited flux at each link.
        """
        lava_thickness = self._grid.at_node["lava__thickness"]

        # Get the links corresponding to the faces of the grid.
        links = self._grid.link_at_face
        flux = lava_flux[links]

        # Find links with outgoing flux
        donor_links = np.where(
            flux > 0.0,
            self._grid.node_at_link_tail[links],
            self._grid.node_at_link_head[links],
        )

        # Calculate the net outflow from each node.
        net_outflow = np.bincount(
            donor_links,
            weights=np.abs(flux) * self._grid.length_of_face,
            minlength=self._grid.number_of_nodes,
        )

        # Volume rate a node can supply: vent input plus the lava it holds.
        available_lava = (
            self.vent_flux + lava_thickness / dt
        ) * self._grid.cell_area_at_node

        # Scale the outgoing flux so that no node loses more lava than it has available.
        scale = np.ones(self._grid.number_of_nodes)
        np.divide(available_lava, net_outflow, out=scale, where=net_outflow > 0.0)
        scale = np.minimum(scale, 1.0)

        limited_lava_flux = lava_flux.copy()
        limited_lava_flux[links] = flux * scale[donor_links]
        return limited_lava_flux

    def _calc_radiative_cooling_rate(self):
        """Calculate the radiative cooling rate of the lava."""

        lava_thickness = self._grid.at_node["lava__thickness"]
        lava_temperature = self._grid.at_node["lava__temperature"]
        cell_areas = self._grid.area_of_cell[self._grid.cell_at_node]

        has_lava_no_lid = (lava_thickness > self._dry_thickness_threshold) & (
            self.lid_thickness == self._min_lid_thickness
        )

        radiative_cooling_rate = np.zeros_like(lava_temperature)
        radiative_cooling_rate[has_lava_no_lid] = (
            -self._emissivity
            * scipy.constants.Stefan_Boltzmann
            * (lava_temperature[has_lava_no_lid] ** 4 - self._ambient_temperature**4)
            * cell_areas[has_lava_no_lid]
        )
        radiative_cooling_rate[self._vent_nodes] = 0.0

        return radiative_cooling_rate

    def _calc_conductive_cooling_rate(self):
        """Calculate the conductive cooling rate of the lava."""

        lava_thickness = self._grid.at_node["lava__thickness"]
        lava_temperature = self._grid.at_node["lava__temperature"]

        has_lava_and_lid = (lava_thickness > self._dry_thickness_threshold) & (
            self.lid_thickness > self._min_lid_thickness
        )
        cell_areas = self._grid.area_of_cell[self._grid.cell_at_node]

        lava_temperature[:] = self._calc_lava_temperature(
            self.lava_thermal_energy,
            lava_thickness,
            cell_areas,
            has_lava_and_lid,
            lava_temperature,
        )

        basal_energy_loss_rate = -np.maximum(
            0.0,
            lhf.conductive_heat_flux(
                lava_temperature,
                cell_areas,
                self._basal_thickness,
                thermal_conductivity=self._thermal_conductivity,
                ambient_temperature=self._ambient_temperature,
            ),
        )

        self.lid_thickness[:] = np.maximum(self.lid_thickness, self._min_lid_thickness)

        atmos_energy_loss_rate = -np.maximum(
            0.0,
            lhf.conductive_heat_flux(
                lava_temperature,
                cell_areas,
                self.lid_thickness,
                thermal_conductivity=self._thermal_conductivity,
                ambient_temperature=self._ambient_temperature,
            ),
        )

        conductive_cooling_rate = basal_energy_loss_rate + atmos_energy_loss_rate
        conductive_cooling_rate[self._vent_nodes] = 0.0
        return conductive_cooling_rate

    def _calc_advective_energy_change(
        self, limited_dt, lava_flux, elevation, lava_thickness, cell_areas
    ):
        lava_temperature = self._grid.at_node["lava__temperature"]
        self.lava_thermal_energy[self._vent_nodes] += lhf.vent_thermal_energy(
            lava_temperature[self._vent_nodes],
            lava_thickness[self._vent_nodes],
            cell_areas[self._vent_nodes],
            lava_density=self._lava_density,
            specific_heat_capacity=self._specific_heat_capacity,
        )

        lava_surface = elevation + lava_thickness
        temperature_at_link = self._grid.map_value_at_max_node_to_link(
            lava_surface, (lava_temperature - self._freezing_temperature)
        )
        thermal_flux_at_link = (
            self._lava_density
            * self._specific_heat_capacity
            * lava_flux
            * temperature_at_link
        )
        thermal_flux_at_link[
            self._grid.status_at_link != self._grid.BC_LINK_IS_ACTIVE
        ] = 0.0

        thermal_flux_div_at_node = self._grid.calc_flux_div_at_node(
            thermal_flux_at_link
        )

        dUdt_adv = -thermal_flux_div_at_node * cell_areas
        total_advective_energy_change = dUdt_adv * limited_dt

        # Make sure that advection does not make lava_thermal_energy negative
        self.lava_thermal_energy[:] += total_advective_energy_change
        self.lava_thermal_energy[:] = np.maximum(self.lava_thermal_energy, 0.0)
        return self.lava_thermal_energy

    def _radiate_thermal_energy(self, radiative_cooling_change):
        """Update lava thermal energy due to radiative cooling.

        Parameters
        ----------
        radiative_cooling_change : ndarray
            The change in thermal energy due to radiative cooling at each node.

        Returns
        -------
        lava_thermal_energy : ndarray
            The updated lava thermal energy at each node.
        """
        radiate_thermal_enegy = radiative_cooling_change <= 0.0

        thermal_energy_change = np.zeros(self._grid.number_of_nodes)

        thermal_energy_change[radiate_thermal_enegy] = -np.minimum(
            abs(radiative_cooling_change[radiate_thermal_enegy]),
            self.lava_thermal_energy[radiate_thermal_enegy],
        )

        self.lava_thermal_energy[:] += thermal_energy_change
        return self.lava_thermal_energy

    def _conduct_thermal_energy(self, conductive_energy_change):
        conduct_thermal_energy = conductive_energy_change < 0.0

        thermal_energy_change = np.zeros(self._grid.number_of_nodes)

        thermal_energy_change[conduct_thermal_energy] = -np.minimum(
            abs(conductive_energy_change[conduct_thermal_energy]),
            self.lava_thermal_energy[conduct_thermal_energy],
        )
        self.lava_thermal_energy[:] += thermal_energy_change

        conductive_energy_change[conduct_thermal_energy] -= thermal_energy_change[
            conduct_thermal_energy
        ]

        return conductive_energy_change

    def _calc_freezing_rate(
        self, lava_thickness, cell_areas, conductive_energy_change, limited_dt
    ):
        self.freezing_rate[:] = 0.0

        # Determine how much mass energy there is so that we know how much energy
        # can be removed from the system due to freezing at each point.
        total_mass_energy = (
            lava_thickness
            * cell_areas
            * self._lava_density
            * self._latent_heat_of_freezing
        )

        remove_mass_energy = conductive_energy_change < 0.0

        conductive_energy_change[remove_mass_energy] = np.maximum(
            conductive_energy_change[remove_mass_energy],
            -total_mass_energy[remove_mass_energy],
        )

        self.freezing_rate[remove_mass_energy] = -conductive_energy_change[
            remove_mass_energy
        ] / (
            self._lava_density
            * self._latent_heat_of_freezing
            * cell_areas[remove_mass_energy]
            * limited_dt
        )

        return self.freezing_rate

    def _update_lava_thickness(self, lava_thickness, lava_flux, limited_dt):
        lava_flux_div_at_node = self._grid.calc_flux_div_at_node(lava_flux)

        dHdt = self.vent_flux - lava_flux_div_at_node - self.freezing_rate

        lava_thickness[self._grid.core_nodes] += (
            dHdt[self._grid.core_nodes] * limited_dt
        )

        no_lava = lava_thickness <= self._dry_thickness_threshold
        lava_thickness[no_lava] = 0.0

        return lava_thickness

    def _partition_lava_to_lid_and_topography(
        self, elevation, lava_thickness, limited_dt
    ):
        frozen_lava_thickness = self.freezing_rate * limited_dt

        add_to_lid = self.lid_thickness < self._max_lid_thickness

        thickness_to_add_to_lid = np.minimum(
            frozen_lava_thickness, self._max_lid_thickness - self.lid_thickness
        )

        self.lid_thickness[add_to_lid] += thickness_to_add_to_lid[add_to_lid]

        frozen_lava_thickness[add_to_lid] -= thickness_to_add_to_lid[add_to_lid]

        # Add the rest of the frozen lava to the topography
        elevation[self._grid.core_nodes] += frozen_lava_thickness[self._grid.core_nodes]

        # Check to see if any lava is thicker than the maximum allowed thickness.
        too_much_lava = lava_thickness > self._max_lava_thickness

        # Add the difference to the topography
        elevation[too_much_lava] += (
            lava_thickness[too_much_lava] - self._max_lava_thickness
        )
        lava_thickness[too_much_lava] = self._max_lava_thickness

        no_lava = lava_thickness <= self._dry_thickness_threshold
        elevation[no_lava] += self.lid_thickness[no_lava] - self._min_lid_thickness
        self.lid_thickness[no_lava] = self._min_lid_thickness

    def _calc_changes(
        self, limited_dt, lava_flux, conductive_cooling_rate, radiative_cooling_rate
    ):
        """Calculate the changes in lava thermal energy and thickness
           due to advection and cooling.

        Parameters
        ----------
        limited_dt : float
            The limited time step for the current iteration.
        lava_flux : ndarray
            The flux of lava at each link.
        conductive_cooling_rate : ndarray
            The conductive cooling rate at each node.
        radiative_cooling_rate : ndarray
            The radiative cooling rate at each node.
        """

        elevation = self._grid.at_node["topographic__elevation"]
        lava_thickness = self._grid.at_node["lava__thickness"]
        lava_temperature = self._grid.at_node["lava__temperature"]
        cell_areas = self._grid.area_of_cell[self._grid.cell_at_node]

        lava_flux = self._limit_lava_flux(lava_flux, limited_dt)

        self.lava_thermal_energy = self._calc_advective_energy_change(
            limited_dt, lava_flux, elevation, lava_thickness, cell_areas
        )

        radiative_cooling_change = radiative_cooling_rate * limited_dt

        self.lava_thermal_energy = self._radiate_thermal_energy(
            radiative_cooling_change
        )

        conductive_energy_change = conductive_cooling_rate * limited_dt

        conductive_energy_change = self._conduct_thermal_energy(
            conductive_energy_change
        )

        self.freezing_rate = self._calc_freezing_rate(
            lava_thickness, cell_areas, conductive_energy_change, limited_dt
        )

        lava_thickness = self._update_lava_thickness(
            lava_thickness, lava_flux, limited_dt
        )

        self._partition_lava_to_lid_and_topography(
            elevation, lava_thickness, limited_dt
        )

        has_lava = lava_thickness > self._dry_thickness_threshold
        lava_temperature[:] = self._calc_lava_temperature(
            self.lava_thermal_energy,
            lava_thickness,
            cell_areas,
            has_lava,
            lava_temperature,
        )

        # If there is no lava and there is a lid, accrete the lid to the topography.
        no_lava = lava_thickness <= self._dry_thickness_threshold
        lava_temperature[no_lava] = self._ambient_temperature

    def _calc_vent_flux(self):
        """Calculate the lava flux from the vent."""

        self.vent_flux[:] = 0.0
        self.vent_flux[self._vent_nodes] = lhf.vent_flux(
            self._eruptive_volume_flux,
            self._grid.area_of_cell[self._grid.cell_at_node[self._vent_nodes]],
        )

        lava_thickness = self._grid.at_node["lava__thickness"]

        if np.amin(lava_thickness) < self._dry_thickness_threshold:
            lava_thickness[:] = np.maximum(lava_thickness, 0.0)

    def run_one_step(self, dt):
        """Route the lava over the topography, cooling and freezing
        it as necessary for one time step.

        Parameters
        ----------
        dt : float
            Time step over which to route the lava.
        """

        sub_dt = 0.0
        while sub_dt < dt:

            self._calc_vent_flux()

            lava_speed, lava_flux = self._calc_lava_velocity()
            conductive_rate = self._calc_conductive_cooling_rate()
            radiative_rate = self._calc_radiative_cooling_rate()
            advection_dt = self._calc_CFL_advection_timestep(lava_speed, dt, sub_dt)
            freezing_dt = self._calc_CFL_freezing_timestep(conductive_rate, dt, sub_dt)
            radiative_dt = self._calc_CFL_radiative_timestep(radiative_rate, dt, sub_dt)
            current_dt = min(advection_dt, freezing_dt, radiative_dt)
            self._calc_changes(current_dt, lava_flux, conductive_rate, radiative_rate)

            sub_dt += current_dt
