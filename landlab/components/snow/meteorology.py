"""Landlab component that calculates the energy fluxes

This component calculates several energy fluxes (e.g., solar radiation,
long wave radiation, sensible heat, latent heat) to estimate the net total energy
flux that can be used as the input for the snow energy balance component.
The code is implemented based on the TopoFlow meteorology component
(by Scott D. Packham) with some adjustments.
https://github.com/peckhams/topoflow36/blob/master/topoflow/components/met_base.py

@author: Tian Gan Oct, 2023

"""

import numpy as np

from landlab import Component


class Meteorology(Component):

    """Simulate snowmelt process using snow energy balance method.

    This component calculates several energy fluxes (e.g., solar radiation,
    long wave radiation, sensible heat, latent heat) to estimate the net total energy
    flux that can be used as the input for the snow energy balance component.

    The net total energy flux (Q_sum) is calculated as:

    Q_sum = Qn_SW + Qn_LW + Qh + Qe + Qa + Qc
    - Qn_SW: net short wave energy flux
    - Qn_LW: net long wave energy flux
    - Qh: sensible heat energy flux
    - Qe: latent heat energy flux
    - Qa: net energy flux advected by moving water
    - Qc: net energy flux via conduction from snow to soil


    Parameters
    ----------
    grid : ModelGrid
        A Landlab model grid object
    start_datetime : string
        start time in format of "yyyy/mm/dd hh:mm:ss"
    rho_H2O : float (default 1000 kg/m3)
        water density
    rho_air : float (default 1.2614 kg/m3)
        air density
    Cp_air  : float (default 1005.7 J kg-1 K-1)
        air heat capacity


    Examples  # TODO add example
    --------
    """

    _name = "Meteorology"

    _unit_agnostic = False

    _info = {
        # TODO: inputs/outputs not sure Ri, De, Dh, Dn, em_air, W_p, e_air,
        #  T_dew, e_sat_air, e_surf, e_sat_surf,
        # input fields
        "atmosphere_bottom_air__temperature": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "deg_C",
            "mapping": "node",
            "doc": "",
        },  # T_air
        "land_surface__temperature": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "deg_C",
            "mapping": "node",
            "doc": "",
        },  # T_surf
        "land_surface__latitude": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "J kg-1 K-1",
            "mapping": "node",
            "doc": "snow heat capacity",
        },  # lat_deg
        "land_surface__longitude": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "J kg-1 K-1",
            "mapping": "node",
            "doc": "snow heat capacity",
        },  # lon_deg
        "land_surface__aspect_angle": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "radians",
            "mapping": "node",
            "doc": "land surface aspect",
        },  # aspect
        "land_surface__slope_angle": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "radians",
            "mapping": "node",
            "doc": "land surface slope",
        },  # slope
        "snowpack__depth": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "snow depth",
        },  # h_snow
        "land_surface__albedo": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "land surface albedo",
        },  # albedo
        "land_surface__emissivity": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "land surface emissivity",
        },  # em_surf
        "atmosphere_aerosol_dust__reduction_of_transmittance": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "aerosol dust reduction of transmittance",
        },  # dust_atten
        "atmosphere_bottom_air__brutsaert_emissivity_canopy_factor": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "brutsaert emissivity canopy factor",
        },  # canopy_factor
        "atmosphere_bottom_air__brutsaert_emissivity_cloud_factor": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "brutsaert emissivity cloud factor",
        },  # cloud_factor
        "atmosphere_bottom_air_water-vapor__relative_saturation": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "air water vapor relative saturation",
        },  # RH
        "atmosphere_bottom_air__pressure": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "mbar",
            "mapping": "node",
            "doc": "air water vapor relative saturation",
        },  # p0
        "atmosphere_bottom_air_flow__log_law_roughness_length": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "meter",
            "mapping": "node",
            "doc": "air flow log law roughness length",
        },  # z0_air
        "atmosphere_bottom_air_flow__speed_reference_height": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "meter",
            "mapping": "node",
            "doc": "air flow speed reference height",
        },  # z
        "atmosphere_bottom_air_flow__reference-height_speed": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "m s-1",
            "mapping": "node",
            "doc": "air flow speed reference height",
        },  # uz
        # output fields (5 var)
        "land_surface_net-total-energy__energy_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W m-2",
            "mapping": "node",
            "doc": "net total energy flux",
        },  # Q_sum
        "land_surface_net-shortwave-radiation__energy_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W m-2",
            "mapping": "node",
            "doc": "net short wave radiation energy flux",
        },  # Qn_SW
        "land_surface_net-longwave-radiation__energy_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W m-2",
            "mapping": "node",
            "doc": "net long wave radiation energy flux",
        },  # Qn_LW
        "atmosphere_bottom_air_land_net-latent-heat__energy_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W m-2",
            "mapping": "node",
            "doc": "net latent heat energy flux",
        },  # Qe
        "atmosphere_bottom_air_land_net-sensible-heat__energy_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W m-2",
            "mapping": "node",
            "doc": "net sensible heat energy flux",
        },  # Qh
    }

    def __init__(
        self,
        grid,
        rho_H2O=1000,
        rho_air=1.2614,
        Cp_air=1005.7,
    ):
        """Initialize Meteorology component"""

        super().__init__(grid)

        # parameters
        self._rho_H2O = rho_H2O  # kg m-3
        self._rho_air = rho_air  # kg m-3
        self._Cp_air = Cp_air  # J kg-1 K-1

        # constants
        self.set_constants()

        # input fields
        self._T_air = grid.at_node["atmosphere_bottom_air__temperature"]
        self._T_surf = grid.at_node["land_surface__temperature"]
        self._lat_deg = grid.at_node["land_surface__latitude"]
        self._lon_deg = grid.at_node["land_surface__longitude"]
        self._aspect = grid.at_node["land_surface__aspect_angle"]
        self._slope = grid.at_node["land_surface__slope_angle"]

        # TODO: check those default values
        if "snowpack__depth" in grid.at_node:
            self._h_swe = grid.at_node["snowpack__depth"]
        else:
            self._h_swe = grid.add_zeros("snowpack__depth", at="node")

        if "land_surface__albedo" in grid.at_node:
            self._albedo = grid.at_node["land_surface__albedo"]
        else:
            self._albedo = grid.add_full("land_surface__albedo", 0.3)

        if "land_surface__emissivity" in grid.at_node:
            self._em_surf = grid.at_node["land_surface__emissivity"]
        else:
            self._em_surf = grid.add_full("land_surface__emissivity", 0.98)

        if "atmosphere_aerosol_dust__reduction_of_transmittance" in grid.at_node:
            self._dust_atten = grid.at_node[
                "atmosphere_aerosol_dust__reduction_of_transmittance"
            ]
        else:
            self._dust_atten = grid.add_full(
                "atmosphere_aerosol_dust__reduction_of_transmittance", 0.08
            )

        if "atmosphere_bottom_air__brutsaert_emissivity_canopy_factor" in grid.at_node:
            self._canopy_factor = grid.at_node[
                "atmosphere_bottom_air__brutsaert_emissivity_canopy_factor"
            ]
        else:
            self._canopy_factor = grid.add_full(
                "atmosphere_bottom_air__brutsaert_emissivity_canopy_factor", 0.0
            )

        if "atmosphere_bottom_air__brutsaert_emissivity_cloud_factor" in grid.at_node:
            self._cloud_factor = grid.at_node[
                "atmosphere_bottom_air__brutsaert_emissivity_cloud_factor"
            ]
        else:
            self._cloud_factor = grid.add_full(
                "atmosphere_bottom_air__brutsaert_emissivity_cloud_factor", 0.0
            )

        if "atmosphere_bottom_air_water-vapor__relative_saturation" in grid.at_node:
            self._RH = grid.at_node[
                "atmosphere_bottom_air_water-vapor__relative_saturation"
            ]
        else:
            self._RH = grid.add_full(
                "atmosphere_bottom_air_water-vapor__relative_saturation", 0.5
            )

        if "atmosphere_bottom_air__pressure" in grid.at_node:
            self._p0 = grid.at_node["atmosphere_bottom_air__pressure"]
        else:
            self._p0 = grid.add_full("atmosphere_bottom_air__pressure", 1000)

        if "atmosphere_bottom_air_flow__log_law_roughness_length" in grid.at_node:
            self._z0_air = grid.at_node[
                "atmosphere_bottom_air_flow__log_law_roughness_length"
            ]
        else:
            self._z0_air = grid.add_full(
                "atmosphere_bottom_air_flow__log_law_roughness_length", 0.02
            )

        if "atmosphere_bottom_air_flow__speed_reference_height" in grid.at_node:
            self._z = grid.at_node["atmosphere_bottom_air_flow__speed_reference_height"]
        else:
            self._z = grid.add_full(
                "atmosphere_bottom_air_flow__speed_reference_height", 2
            )

        if "atmosphere_bottom_air_flow__reference-height_speed" in grid.at_node:
            self._uz = grid.at_node[
                "atmosphere_bottom_air_flow__reference-height_speed"
            ]
        else:
            self._uz = grid.add_full(
                "atmosphere_bottom_air_flow__reference-height_speed", 2
            )

        # output fields
        super().initialize_output_fields()
        self._Q_sum = grid.at_node["land_surface_net-total-energy__energy_flux"]
        self._Qn_SW = grid.at_node["land_surface_net-shortwave-radiation__energy_flux"]
        self._Qn_LW = grid.at_node["land_surface_net-longwave-radiation__energy_flux"]
        self._Qe = grid.at_node[
            "atmosphere_bottom_air_land_net-latent-heat__energy_flux"
        ]
        self._Qh = grid.at_node[
            "atmosphere_bottom_air_land_net-sensible-heat__energy_flux"
        ]

    @property
    def rho_H2O(self):
        return self._rho_H2O

    @rho_H2O.setter
    def rho_H2O(self, rho_H2O):
        assert rho_H2O > 0, "assign rho_H2O with positive value"
        self._rho_H2O = rho_H2O

    @property
    def rho_air(self):
        return self._rho_air

    @rho_air.setter
    def rho_air(self, rho_air):
        assert rho_air > 0, "assign rho_air with positive value"
        self._rho_air = rho_air

    @property
    def Cp_air(self):
        return self._Cp_air

    @Cp_air.setter
    def Cp_air(self, Cp_air):
        assert Cp_air > 0, "assign Cp_air with positive value"
        self._Cp_air = Cp_air

    def set_constants(self):
        # Define physical constants
        self.g = np.float64(9.81)  # [m s-2, gravity]
        self.kappa = np.float64(0.408)  # [1]  (von Karman)
        self.rho_H2O = np.float64(1000)  # [kg m-3]
        self.rho_air = np.float64(1.2614)  # [kg m-3]
        self.Cp_air = np.float64(1005.7)  # [J kg-1 K-1]
        self.Lv = np.float64(2500000)  # [J kg-1] Latent heat of vaporiz.
        self.Lf = np.float64(334000)  # [J kg-1 = W s kg-1], Latent heat of fusion
        self.sigma = np.float64(5.67e-8)  # [W m-2 K-4]  (Stefan-Boltzman constant)
        self.C_to_K = np.float64(273.15)  # (add to convert deg C to K)

        self.twopi = np.float64(2) * np.pi
        self.one_seventh = np.float64(1) / 7
        self.hours_per_day = np.float64(24)
        self.secs_per_day = np.float64(3600) * self.hours_per_day

        self.latent_heat_constant = np.float64(0.662)

    def update_bulk_richardson_number(self):
        pass

    def update_bulk_aero_conductance(self):
        pass

    def update_sensible_heat_flux(self):
        pass

    def update_saturation_vapor_pressure(self):
        pass

    def update_vapor_pressure(self):
        pass

    def update_dew_point(self):
        pass

    def update_precipitable_water_content(self):
        pass

    def update_latent_heat_flux(self):
        pass

    def update_conduction_heat_flux(self):
        pass

    def update_advection_heat_flux(self):
        pass

    def update_julian_day(self):
        pass

    def update_net_shortwave_radiation(self):
        pass

    def update_em_air(self):
        pass

    def update_net_longwave_radiation(self):
        pass

    def update_net_energy_flux(self):
        pass

    def run_one_step(self, dt):
        # update input fields in case there is new input
        self._P = self._grid.at_node["atmosphere_water__precipitation_leq-volume_flux"]
        self._T_air = self._grid.at_node["atmosphere_bottom_air__temperature"]
        self._T_surf = self._grid.at_node["land_surface__temperature"]
        self._Q_sum = self._grid.at_node["land_surface_net-total-energy__energy_flux"]
        self._rho_snow = self._grid.at_node[
            "snowpack__z_mean_of_mass-per-volume_density"
        ]
        self._Cp_snow = self._grid.at_node[
            "snowpack__z_mean_of_mass-specific_isobaric_heat_capacity"
        ]

        # update state variable
        self.update_bulk_richardson_number()
        self.update_bulk_aero_conductance()
        self.update_sensible_heat_flux()
        self.update_saturation_vapor_pressure(MBAR=True)
        self.update_saturation_vapor_pressure(MBAR=True, SURFACE=True)
        self.update_vapor_pressure()
        self.update_dew_point()
        self.update_precipitable_water_content()
        self.update_vapor_pressure(SURFACE=True)  # this function is called twice
        self.update_latent_heat_flux()  # (uses e_air and e_surf)
        self.update_conduction_heat_flux()
        self.update_advection_heat_flux()
        self.update_julian_day()
        self.update_net_shortwave_radiation()
        self.update_em_air()
        self.update_net_longwave_radiation()
        self.update_net_energy_flux()

    def set_aspect_angle(self):
        # ------------------------------------------------------
        # Read aspect grid.  Alpha must be CW from north.
        # NB!  RT aspect grids have NaNs on edges.
        # ---------------------------------------------------------
        # RT files ending in "_mf-angle.rtg" and "fd-aspect.rtg"
        # contain aspect values.  The former are in [0, 2 Pi]
        # while the latter are in [-Pi, Pi] and both measure
        # CCW from due east.
        # ---------------------------------------------------------
        # #aspects = rtg_files.read_grid(self.aspect_grid_file, self.rti, RTG_type="FLOAT")
        # alpha = (np.pi / 2) - aspects
        # alpha = (self.twopi + alpha) % self.twopi
        # # -----------------------------------------------
        # w_nan = np.where(np.logical_not(np.isfinite(alpha)))
        # n_nan = np.size(w_nan[0])
        # if n_nan != 0:
        #     alpha[w_nan] = np.float64(0)
        #
        # self.alpha = alpha
        pass

    def set_slope_angle(self):
        # -------------------------------------------------
        # Read slope grid & convert to slope angle, beta
        # NB!  RT slope grids have NaNs on edges.
        # -------------------------------------------------
        # slopes = rtg_files.read_grid(self.slope_grid_file, self.rti, RTG_type="FLOAT")
        # beta = np.arctan(slopes)
        # beta = (self.twopi + beta) % self.twopi
        # # ---------------------------------------------
        # w_nan = np.where(np.logical_not(np.isfinite(beta)))
        # n_nan = np.size(w_nan[0])
        # if n_nan != 0:
        #     beta[w_nan] = np.float64(0)
        # # ------------------------------------------------------------------
        # w_bad = np.where(np.logical_or((beta < 0), (beta > np.pi / 2)))
        # n_bad = np.size(w_bad[0])
        # if n_bad != 0:
        #     print("ERROR: In met_base.py, some slope angles are out")
        #     print("       of range.  Returning without setting beta.")
        #     print()
        #     return
        #
        # self.beta = beta
        pass
