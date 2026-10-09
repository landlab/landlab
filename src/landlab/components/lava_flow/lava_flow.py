#!/usr/bin/env python
"""
Component that simulates the effusion of lava from a vent and the
subsequent routing and freezing of the lava over the topography.
"""

import numpy as np

from landlab import Component


class LavaFlow(Component):
    """Effusion, routing, and freezing of lava over topography.

    This component takes the coordinates of a vent (or vents)
    and simulates effusion of lava from the vent(s) over the topography.
    The component advects the lava over the topography with a Bingham
    visco-plastic rheology, and simulates the loss of thermal energy
    as the lava cools and eventually freezes.

    A Bingham rheology relates the shear stress to the strain rate and
    a yield stress through the expression:

    ..math::

        strain_rate = 1/viscosity * (shear_stress - yield_stress)

    For shear stresses below the yield stress, the lava does not flow because
    the strain rate is 0. The shear_stress acting on the lava is gravitational,
    and is equal to:

    ..math::

        shear_stress = rho_lava * gravity * lava__thickness * sin(theta)

    where :math:`rho_lava` is the density of the lava, :math:`gravity` is the
    gravitational acceleration, :math:`lava__thickness` is the thickness of the
    lava layer, and :math:`theta` is the slope of the topography from the horizontal.
    The strain rate can be related to the velocity of the lava by:

    ..math::

        strain_rate = d(velocity)/dz

    where :math:`z` is the vertical coordinate measured from the base of the lava
    layer. z = 0 at the base of the lava layer and H at the top of the lava layer.
    Plugging in the strain rate and shear stress expressions into the constitutive
    equation, and rearranging the terms we end up with:

    ..math::

        d(velocity) = 1/viscosity * (rho_lava * gravity * lava__thickness *
                      sin(theta) - yield_stress) * dz

    which can be integrated to obtain the velocity profile of the lava layer. The
    boundary condition which allows us to determine the integral constant is that
    at z = 0, the velocity is zero. After integrating and averaging the velocity
    over the thickness of the lava, we end up with the velocity expression used to
    advect the lava in this model:

    ..math::

        velocity = lava_thickness/viscosity * (rho_lava * gravity * lava__thickness *
                   sin(theta) / 3 - yield_stress / 2)

    To simulate cooling and freezing of the lava, the component accounts for three
    sources of energy lossfrom the lava flow:
    1. Radiative heat loss from the surface of the lava when no "lid" is present.
    2. Conductive heat loss to the underlying topography.
    3. Conductive heat loss through a "lid" that forms on top of the lava as it cools.

    Lava is erupted with a sensible heat content which is equal to the eruption
    temperature minus the freezing temperture. While the lava is above the freezing
    temperature, energy is strictly removed from the sensible heat content of the lava.
    When the sensible heat content is zero, energy is removed from the lava through
    freezing, which is proportional to the latent heat of freezing. All cooling and
    freezing calculations assume a bulk energy/thermal balance throughout the entire
    lava layer at each point on the grid. As long as lava remains underneat a lid,
    the lid will continue to insulate the lava from radiative heat and conductive
    heat loss. Once lava is full frozen, the lid thickness is added to the topography.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import LavaFlow
    >>> dt = 30.0
    >>> time_to_run = 600.0

    >>> mg = RasterModelGrid((3, 10), xy_spacing=(100.0, 100.0))

    >>> mg.add_zeros("topographic__elevation", at="node")
    >>> mg.add_zeros("lava__thickness", at="node")
    >>> mg.at_node["topographic__elevation"][:] += mg.x_of_node

    >>> mg.set_fixed_value_boundaries_at_grid_edges(True, True, True, True)

    >>> lf = LavaFlow(
    ...     mg, vent_coordinates=[800, 100], volume_flux=10.0, viscosity=100.0, Cp=300
    ... )

    >>> elapsed_time = 0.0
    >>> while elapsed_time < time_to_run:
    ...     if elapsed_time + dt > time_to_run:
    ...         dt = time_to_run - elapsed_time
    ...     lf.run_one_step(dt)
    ...     elapsed_time += dt
    ...

    >>> print(np.round(mg.at_node["lava__thickness"][12:20], 3))
    array([0.,    0.,    0.086, 0.127, 0.127, 0.127, 0.127, 0.])

    >>> print(np.round(mg.at_node["lava__temperature"][12:20], 1))
    array([300.,   300.,   800.,   840.9,  980.1, 1151.5, 1473.,   300.])
    """

    _name = "LavaFlow"

    _unit_agnostic = True

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
            "doc": "Thickness of the lava layer (m)",
        },
        "lava__temperature": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "K",
            "mapping": "node",
            "doc": "Temperature of the lava layer (K)",
        },
        "lid__thickness": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Thickness of the solidified lid (m)",
        },
        "lava__thermal_energy": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "J",
            "mapping": "node",
            "doc": "Thermal energy of the lava layer (J)",
        },
        "freezing__rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/s",
            "mapping": "node",
            "doc": "Rate of freezing of the lava layer (m/s)",
        },
        "vent__flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m^3/s",
            "mapping": "node",
            "doc": "Volume flux of lava at the vent nodes (m^3/s)",
        },
    }

    def __init__(
        self,
        grid,
        vent_coordinates=None,
        rho_lava=2700.0,
        Cp=1000.0,
        k_lava=2.0,
        viscosity=10000.0,
        L_freeze=4.0e5,
        volume_flux=0.1,
        gravity=9.80665,
        yield_stress=1000.0,
        T_erupt=1473.0,
        T_freeze=800.0,
        T_surf=300.0,
        max_lava_thickness=10.0,
        basal_thickness=20.0,
        max_lid_thickness=5.0,
        min_lid_thickness=1.0e-2,
        courant=1.0,
        min_dt=0.0,
    ):
        """Initialize the flexure component.

        Parameters
        ----------
        grid : RasterModelGrid
            A grid.
        rho_lava : float, optional
            Density of the lava (kg / m^3).
        Cp : float, optional
            Specific heat capacity of the lava (J / (kg*K)).
        k_lava : float, optional
            Thermal conductivity of the lava (W / (m*K))
        viscosity : float, optional
            Viscosity of the lava (Pa*s).
        L_freeze : float, optional
            Latent heat of freezing (J / kg).
        volume_flux : float, optional
            Volume flux of lava (m^3 / s).
        vent_coordinates : array_like, optional
            Indices of vent nodes.
        gravity : float, optional
            Acceleration due to gravity (m / s^2).
        yield_stress : float, optional
            Yield stress of the lava (Pa).
        T_erupt : float, optional
            Eruption temperature of the lava (K).
        T_freeze : float, optional
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

        self.rho_lava = rho_lava
        self.viscosity = viscosity
        self.Cp = Cp
        self.k_lava = k_lava
        self.L_freeze = L_freeze
        self.volume_flux = volume_flux
        self.gravity = gravity
        self.yield_stress = yield_stress
        self.T_surf = T_surf
        self.T_freeze = T_freeze
        self.T_erupt = T_erupt
        self.max_lava_thickness = max_lava_thickness
        self.basal_thickness = basal_thickness
        self.min_lid_thickness = min_lid_thickness
        self.max_lid_thickness = max_lid_thickness
        self.courant = courant
        self.min_dt = min_dt

        # Stores the coefficients used for the lava velocity
        self._coef1 = self.rho_lava * self.gravity / (3.0 * self.viscosity)
        self._coef2 = self.yield_stress / (2.0 * self.viscosity)

        # Store the Stefann-Boltzmann constant for radiative heat loss
        self._sb_constant = 5.670374419e-8

        # Minimum epsilon for determining when lava thickness is effectively zero.
        self._eps = 1e-8

        # Set the vent nodes based on the provided coordinates
        if vent_coordinates is not None:
            self._vent_nodes = self._grid.find_nearest_node(vent_coordinates)
        else:
            self._vent_nodes = []

    @property
    def rho_lava(self):
        """Density of the lava (kg/m^3)."""
        return self._rho_lava

    @rho_lava.setter
    def rho_lava(self, value):
        if value > 0.0:
            self._rho_lava = float(value)
        else:
            raise ValueError("lava density must be greater than zero.")

    @property
    def viscosity(self):
        """Viscosity of the lava (Pa*s)."""
        return self._viscosity

    @viscosity.setter
    def viscosity(self, value):
        if value > 0.0:
            self._viscosity = float(value)
        else:
            raise ValueError("lava viscosity must be greater than zero.")

    @property
    def Cp(self):
        """Specific heat capacity of the lava (J/(kg*K))."""
        return self._Cp

    @Cp.setter
    def Cp(self, value):
        if value > 0.0:
            self._Cp = float(value)
        else:
            raise ValueError("lava specific heat capacity must be greater than zero.")

    @property
    def k_lava(self):
        """Thermal conductivity of the lava (W/(m*K))."""
        return self._k_lava

    @k_lava.setter
    def k_lava(self, value):
        if value > 0.0:
            self._k_lava = float(value)
        else:
            raise ValueError("lava thermal conductivity must be greater than zero.")

    @property
    def L_freeze(self):
        """Latent heat of freezing (J/kg)."""
        return self._L_freeze

    @L_freeze.setter
    def L_freeze(self, value):
        if value > 0.0:
            self._L_freeze = float(value)
        else:
            raise ValueError("latent heat of freezing must be greater than zero.")

    @property
    def volume_flux(self):
        """Volume flux of lava (m^3/s)."""
        return self._volume_flux

    @volume_flux.setter
    def volume_flux(self, value):
        if value > 0.0:
            self._volume_flux = float(value)
        else:
            raise ValueError("lava volume flux must be greater than zero.")

    @property
    def gravity(self):
        """Acceleration due to gravity (m/s^2)."""
        return self._gravity

    @gravity.setter
    def gravity(self, value):
        if value > 0.0:
            self._gravity = float(value)
        else:
            raise ValueError("gravity must be greater than zero.")

    @property
    def yield_stress(self):
        """Yield stress of the lava (Pa)."""
        return self._yield_stress

    @yield_stress.setter
    def yield_stress(self, value):
        if value >= 0.0:
            self._yield_stress = float(value)
        else:
            raise ValueError("yield stress must be greater than or equal to zero.")

    @property
    def T_surf(self):
        """Surface temperature (K)."""
        return self._T_surf

    @T_surf.setter
    def T_surf(self, value):
        if value > 0.0:
            self._T_surf = float(value)
        else:
            raise ValueError("surface temperature must be greater than zero.")

    @property
    def T_freeze(self):
        """Freezing temperature of the lava (K)."""
        return self._T_freeze

    @T_freeze.setter
    def T_freeze(self, value):
        if value > self.T_surf:
            self._T_freeze = float(value)
        else:
            raise ValueError(
                "freezing temperature must be greater than the surface temperature."
            )

    @property
    def T_erupt(self):
        """Eruption temperature of the lava (K)."""
        return self._T_erupt

    @T_erupt.setter
    def T_erupt(self, value):
        if value > self.T_surf and value >= self.T_freeze:
            self._T_erupt = float(value)
        else:
            raise ValueError(
                "eruption temperature must be greater than the "
                "surface temperature and greater than or equal to "
                "the freezing temperature."
            )

    @property
    def max_lava_thickness(self):
        """Maximum thickness of the lava flow (m)."""
        return self._max_lava_thickness

    @max_lava_thickness.setter
    def max_lava_thickness(self, value):
        if value > 0.0:
            self._max_lava_thickness = float(value)
        else:
            raise ValueError("maximum lava thickness must be greater than zero.")

    @property
    def basal_thickness(self):
        """Thickness of the basal layer (m)."""
        return self._basal_thickness

    @basal_thickness.setter
    def basal_thickness(self, value):
        if value > 0.0:
            self._basal_thickness = float(value)
        else:
            raise ValueError("basal thickness must be greater than zero.")

    @property
    def min_lid_thickness(self):
        """Minimum thickness of the lid layer (m)."""
        return self._min_lid_thickness

    @min_lid_thickness.setter
    def min_lid_thickness(self, value):
        if value > 0.0:
            self._min_lid_thickness = float(value)
        else:
            raise ValueError("minimum lid thickness must be greater than zero.")

    @property
    def max_lid_thickness(self):
        """Maximum thickness of the lid layer (m)."""
        return self._max_lid_thickness

    @max_lid_thickness.setter
    def max_lid_thickness(self, value):
        if value >= self.min_lid_thickness:
            self._max_lid_thickness = float(value)
        else:
            raise ValueError(
                "maximum lid thickness must be greater than or "
                "equal to the minimum lid thickness."
            )

    @property
    def courant(self):
        """Courant number for the simulation."""
        return self._courant

    @courant.setter
    def courant(self, value):
        if value > 0.0:
            self._courant = float(value)
        else:
            raise ValueError("Courant number must be greater than zero.")

    @property
    def min_dt(self):
        """Minimum allowable time step for the simulation."""
        return self._min_dt

    @min_dt.setter
    def min_dt(self, value):
        if value >= 0.0:
            self._min_dt = float(value)
        else:
            raise ValueError("Minimum time step must be greater than or equal to zero.")

    def _calc_lava_velocity(self):
        """Calculate lava velocity."""
        lava_thickness = self._grid.at_node["lava__thickness"]
        lava_surface = self._grid.at_node["topographic__elevation"] + lava_thickness

        lava_surface_grad = self._grid.calc_grad_at_link(lava_surface)
        lava_at_link = self._grid.map_value_at_max_node_to_link(
            lava_surface, lava_thickness
        )

        velmag = np.maximum(
            self._coef1 * lava_at_link**2 * np.abs(lava_surface_grad)
            - self._coef2 * lava_at_link,
            0.0,
        )
        velmag[self._grid.status_at_link != self._grid.BC_LINK_IS_ACTIVE] = 0.0
        lava_flux = -np.sign(lava_surface_grad) * lava_at_link * velmag

        return velmag, lava_flux

    def _calc_CFL_advection_timestep(self, velmag, dt, sub_dt):
        """Calculate the CFL-limited advection timestep.

        Parameters
        ----------
        velmag : ndarray
            The magnitude of the lava velocity at each link.
        dt : float
            The current time step.
        sub_dt : float
            The accumulated sub-time step.
        """

        velmax = np.amax(velmag)
        current_dt = dt
        max_link = np.amax(self._grid.length_of_link)

        if velmax > 0.0:
            current_dt = min(self._courant * (max_link / velmax), current_dt)

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
        lava_thermal_energy = self._grid.at_node["lava__thermal_energy"]
        cell_areas = self._grid.area_of_cell[self._grid.cell_at_node]

        freezing_nodes = (
            (lava_thickness > 0.0)
            & (lava_thermal_energy == 0.0)
            & (conductive_cooling_rate < 0.0)
        )
        freezing_nodes[self._vent_nodes] = False

        current_dt = dt

        if np.any(freezing_nodes):
            available_energy_for_freezing = (
                lava_thickness[freezing_nodes]
                * self._rho_lava
                * self._L_freeze
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
        lava_thermal_energy = self._grid.at_node["lava__thermal_energy"]
        has_thermal_energy = (lava_thermal_energy > 0.0) & (
            radiative_cooling_rate < 0.0
        )
        has_thermal_energy[self._vent_nodes] = False

        current_dt = dt

        if np.any(has_thermal_energy):
            radiative_dt = np.min(
                lava_thermal_energy[has_thermal_energy]
                / np.abs(radiative_cooling_rate[has_thermal_energy])
            )
            current_dt = np.clip(self._courant * radiative_dt, self._min_dt, current_dt)

        return min(current_dt, dt - sub_dt)

    def _compute_lava_temperature(
        self, lava_thermal_energy, lava_thickness, cell_areas, mask, lava_temperature
    ):

        lava_temperature[mask] = np.minimum(
            self._T_freeze
            + lava_thermal_energy[mask]
            / (self._rho_lava * self._Cp * cell_areas[mask] * lava_thickness[mask]),
            self._T_erupt,
        )
        lava_temperature[self._vent_nodes] = self._T_erupt
        lava_temperature[:] = np.maximum(lava_temperature, self._T_surf)

        return lava_temperature

    def _calc_radiative_cooling_rate(self):
        """Calculate the radiative cooling rate of the lava."""

        lava_temperature = self._grid.at_node["lava__temperature"]
        lava_thickness = self._grid.at_node["lava__thickness"]
        lid_thickness = self._grid.at_node["lid__thickness"]
        cell_areas = self._grid.area_of_cell[self._grid.cell_at_node]

        has_lava_no_lid = (lava_thickness > 0.0) & (
            lid_thickness == self._min_lid_thickness
        )

        radiative_cooling_rate = np.zeros_like(lava_temperature)
        radiative_cooling_rate[has_lava_no_lid] = (
            -0.9
            * self._sb_constant
            * (lava_temperature[has_lava_no_lid] ** 4 - self._T_surf**4)
            * cell_areas[has_lava_no_lid]
        )
        radiative_cooling_rate[self._vent_nodes] = 0.0

        return radiative_cooling_rate

    def _calc_conductive_cooling_rate(self):
        """Calculate the conductive cooling rate of the lava."""

        lava_temperature = self._grid.at_node["lava__temperature"]
        lava_thermal_energy = self._grid.at_node["lava__thermal_energy"]
        lava_thickness = self._grid.at_node["lava__thickness"]
        lid_thickness = self._grid.at_node["lid__thickness"]

        has_lava_and_lid = (lava_thickness > 0.0) & (
            lid_thickness > self._min_lid_thickness
        )
        cell_areas = self._grid.area_of_cell[self._grid.cell_at_node]

        lava_temperature[:] = self._compute_lava_temperature(
            lava_thermal_energy,
            lava_thickness,
            cell_areas,
            has_lava_and_lid,
            lava_temperature,
        )

        basal_energy_loss_rate = np.minimum(
            0.0,
            -self._k_lava
            * (lava_temperature - self._T_surf)
            / self._basal_thickness
            * cell_areas,
        )

        lid_thickness[:] = np.maximum(lid_thickness, self._min_lid_thickness)
        atmos_energy_loss_rate = np.minimum(
            0.0,
            -self._k_lava
            * (lava_temperature - self._T_surf)
            / lid_thickness
            * cell_areas,
        )

        conductive_cooling_rate = basal_energy_loss_rate + atmos_energy_loss_rate
        conductive_cooling_rate[self._vent_nodes] = 0.0
        return conductive_cooling_rate

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
        lava_temperature = self._grid.at_node["lava__temperature"]
        lava_thermal_energy = self._grid.at_node["lava__thermal_energy"]
        lava_thickness = self._grid.at_node["lava__thickness"]
        freezing_rate = self._grid.at_node["freezing__rate"]
        lid_thickness = self._grid.at_node["lid__thickness"]
        vent_flux = self._grid.at_node["vent__flux"]
        cell_areas = self._grid.area_of_cell[self._grid.cell_at_node]

        lava_thermal_energy[self._vent_nodes] += (
            self._volume_flux
            * self._rho_lava
            * self._Cp
            * (self._T_erupt - self._T_freeze)
            * limited_dt
        )
        has_lava = lava_thickness > 0.0

        lava_surface = elevation + lava_thickness
        temperature_at_link = self._grid.map_value_at_max_node_to_link(
            lava_surface, (lava_temperature - self._T_freeze)
        )
        thermal_flux_at_link = (
            self._rho_lava * self._Cp * lava_flux * temperature_at_link
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
        lava_thermal_energy[:] += total_advective_energy_change
        lava_thermal_energy[:] = np.maximum(lava_thermal_energy, 0.0)

        # Determine how much mass energy there is so that we know how much energy
        # can be removed from the system due to freezing at each point.
        total_mass_energy = (
            lava_thickness * cell_areas * self._rho_lava * self._L_freeze
        )

        thermal_energy_change = np.zeros_like(total_mass_energy)

        radiative_cooling_change = radiative_cooling_rate * limited_dt
        conductive_energy_change = conductive_cooling_rate * limited_dt

        radiate_thermal_enegy = radiative_cooling_change <= 0.0

        thermal_energy_change[radiate_thermal_enegy] = -np.minimum(
            abs(radiative_cooling_change[radiate_thermal_enegy]),
            lava_thermal_energy[radiate_thermal_enegy],
        )

        lava_thermal_energy[:] += thermal_energy_change

        thermal_energy_change[:] = 0.0
        conduct_thermal_energy = conductive_energy_change < 0.0

        thermal_energy_change[conduct_thermal_energy] = -np.minimum(
            abs(conductive_energy_change[conduct_thermal_energy]),
            lava_thermal_energy[conduct_thermal_energy],
        )
        lava_thermal_energy[:] += thermal_energy_change

        conductive_energy_change[conduct_thermal_energy] -= thermal_energy_change[
            conduct_thermal_energy
        ]

        lava_temperature[:] = self._compute_lava_temperature(
            lava_thermal_energy, lava_thickness, cell_areas, has_lava, lava_temperature
        )

        freezing_rate[:] = 0.0

        remove_mass_energy = conductive_energy_change < 0.0
        freezing_rate[remove_mass_energy] = -conductive_energy_change[
            remove_mass_energy
        ] / (
            self._rho_lava
            * self._L_freeze
            * cell_areas[remove_mass_energy]
            * limited_dt
        )

        max_net_outflux = vent_flux + (lava_thickness / limited_dt) - freezing_rate
        lava_flux_div_at_node = self._grid.calc_flux_div_at_node(lava_flux)
        lava_flux_div_at_node_limited = lava_flux_div_at_node.copy()
        lava_flux_div_at_node_limited[self._grid.core_nodes] = np.minimum(
            lava_flux_div_at_node_limited[self._grid.core_nodes],
            max_net_outflux[self._grid.core_nodes],
        )

        dHdt = vent_flux - lava_flux_div_at_node_limited - freezing_rate

        # print(freezing_rate)

        lava_thickness[self._grid.core_nodes] += (
            dHdt[self._grid.core_nodes] * limited_dt
        )

        frozen_lava_thickness = freezing_rate * limited_dt

        add_to_lid = lid_thickness < self._max_lid_thickness

        thickness_to_add_to_lid = np.minimum(
            frozen_lava_thickness, self._max_lid_thickness - lid_thickness
        )

        lid_thickness[add_to_lid] += thickness_to_add_to_lid[add_to_lid]

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

        # If there is no lava and there is a lid, accrete the lid to the topography.
        no_lava = lava_thickness < self._eps
        lava_temperature[no_lava] = self._T_surf

        elevation[no_lava] += lid_thickness[no_lava] - self._min_lid_thickness
        lid_thickness[no_lava] = self._min_lid_thickness

    def _calc_vent_flux(self):
        """Calculate the lava flux from the vent."""

        vent_flux = self._grid.at_node["vent__flux"]
        vent_flux[:] = 0.0
        vent_flux[self._vent_nodes] = (
            self._volume_flux
            / self._grid.area_of_cell[self._grid.cell_at_node[self._vent_nodes]]
        )

        lava_thickness = self._grid.at_node["lava__thickness"]

        if np.amin(lava_thickness) < self._eps:
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

            velmag, lava_flux = self._calc_lava_velocity()
            conductive_rate = self._calc_conductive_cooling_rate()
            radiative_rate = self._calc_radiative_cooling_rate()
            advection_dt = self._calc_CFL_advection_timestep(velmag, dt, sub_dt)
            freezing_dt = self._calc_CFL_freezing_timestep(conductive_rate, dt, sub_dt)
            radiative_dt = self._calc_CFL_radiative_timestep(radiative_rate, dt, sub_dt)
            current_dt = min(advection_dt, freezing_dt, radiative_dt)
            self._calc_changes(current_dt, lava_flux, conductive_rate, radiative_rate)

            sub_dt += current_dt
