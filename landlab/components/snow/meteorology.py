"""Landlab component that calculates the energy fluxes

This component calculates several energy fluxes (e.g., solar radiation,
long wave radiation, sensible heat, latent heat) to estimate the net total energy
flux that can be used as the input for the snow energy balance component.
The code is implemented based on the TopoFlow meteorology component
(by Scott D. Packham) with some adjustments.
https://github.com/peckhams/topoflow36/blob/master/topoflow/components/met_base.py

@author: Tian Gan Oct, 2023

"""

import datetime

import numpy as np
from dateutil.relativedelta import relativedelta

from landlab import Component

from . import solar_funcs


class Meteorology(Component):

    """Calculate several energy fluxes and climate variables.

    This component calculates several energy fluxes (solar radiation,
    long wave radiation, sensible heat, latent heat) and related climate variables
    (e.g., air vapor pressure, dew point, and air emissivity).

    The estimated total net energy flux (Q_sum) can be used as the input for
    the SnowEnergyBalance component. Q_sum is calculated as:

    Q_sum = Qn_SW + Qn_LW + Qh + Qe + Qa + Qc
    - Qn_SW: net short wave energy flux
    - Qn_LW: net long wave energy flux
    - Qh: sensible heat energy flux
    - Qe: latent heat energy flux
    - Qa: net energy flux advected by moving water (negligible; Qa=0)
    - Qc: net energy flux via conduction from snow to soil (negligible; Qc=0)


    Parameters
    ----------
    grid : ModelGrid
        A Landlab model grid object
    start_datetime : string
        start time in format of "yyyy-mm-dd hh:mm:ss"
    GMT_offset: int (default 0)
        GMT offset at the locations of interest. It should be entered as an integer
        or an array of integers between 0 and 12 with negative values for
        locations west of the prime meridian.
    rho_H2O : float (default 1000 kg/m3)
        water density
    rho_air : float (default 1.2614 kg/m3)
        air density
    Cp_air  : float (default 1005.7 J kg-1 K-1)
        air heat capacity
    satterlund : bool (default False)
        if true, use Satterlund method for saturation vapor pressure


    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import Meteorology
    >>> grid = RasterModelGrid((2, 2))
    >>> grid.add_full("atmosphere_bottom_air__temperature", 1, at="node")
    array([ 1.,  1.,  1.,  1.])
    >>> grid.add_full("land_surface__temperature", -1, at="node")
    array([-1., -1., -1., -1.])
    >>> grid.add_full("land_surface__latitude", 40, at="node")
    array([ 40.,  40.,  40.,  40.])
    >>> grid.add_full("land_surface__longitude", -105, at="node")
    array([-105., -105., -105., -105.])
    >>> grid.add_full(
    ...     "atmosphere_bottom_air_water-vapor__relative_saturation", 0.4, at="node"
    ... )
    array([ 0.4,  0.4,  0.4,  0.4])
    >>> dt = 60 * 60 * 48
    >>> met = Meteorology(
    ...     grid, start_datetime="2023-01-01 12:00:00", GMT_offset=-7, satterlund=False
    ... )
    >>> met.run_one_step(dt)
    >>> grid.at_node["land_surface_net-total-energy__energy_flux"]
    array([ 498.41796974,  498.41796974,  498.41796974,  498.41796974])
    """

    _name = "Meteorology"

    _unit_agnostic = False

    _info = {
        # input fields (4 req var, 13 opt var)
        "atmosphere_bottom_air__temperature": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "deg_C",
            "mapping": "node",
            "doc": "atmosphere bottom air temperature",
        },  # T_air
        "land_surface__temperature": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "deg_C",
            "mapping": "node",
            "doc": "land surface temperature",
        },  # T_surf (snow pack temp)
        "land_surface__latitude": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "J kg-1 K-1",
            "mapping": "node",
            "doc": "latitude",
        },  # lat_deg
        "land_surface__longitude": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "J kg-1 K-1",
            "mapping": "node",
            "doc": "longitude",
        },  # lon_deg
        "land_surface__aspect_angle": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "radians",
            "mapping": "node",
            "doc": "land surface aspect",
        },  # alpha (aspect) [0, 2pi]
        "land_surface__slope_angle": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "radians",
            "mapping": "node",
            "doc": "land surface slope",
        },  # beta (slope) [0, pi/2]
        # "snowpack__depth": {
        #     "dtype": float,
        #     "intent": "in",
        #     "optional": True,
        #     "units": "m",
        #     "mapping": "node",
        #     "doc": "snow depth",
        # },  # h_snow
        "land_surface__albedo": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "land surface albedo",
        },  # albedo (snow pack)
        "land_surface__emissivity": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "land surface emissivity",
        },  # em_surf (snow pack)
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
            "doc": "bottom air pressure",
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
            "doc": "air flow speed at reference height",
        },  # uz
        # output fields (18 var)
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
            "doc": "net shortwave radiation energy flux",
        },  # Qn_SW
        "land_surface_net-longwave-radiation__energy_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W m-2",
            "mapping": "node",
            "doc": "net longwave radiation energy flux",
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
        "snowpack_land_surface_net-conduction-heat__energy_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W m-2",
            "mapping": "node",
            "doc": "net energy flux via conduction from snow to soil",
        },  # Qc  # name added by me
        "snowpack_land_surface_net-advection-heat__energy_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W m-2",
            "mapping": "node",
            "doc": "net energy flux advected from moving water",
        },  # Qa  # name added by me
        "atmosphere_bottom_air__neutral_bulk_heat_aerodynamic_conductance": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m s-1",
            "mapping": "node",
            "doc": "neutral bulk heat aerodynamic conductance",
        },  # Dn
        "atmosphere_bottom_air__bulk_sensible_heat_aerodynamic_conductance": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m s-1",
            "mapping": "node",
            "doc": "bulk sensible heat aerodynamic conductance",
        },  # Dh
        "atmosphere_bottom_air__bulk_latent_heat_aerodynamic_conductance": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m s-1",
            "mapping": "node",
            "doc": "bulk latent heat aerodynamic conductance",
        },  # De
        "land_surface_air_water-vapor__saturated_partial_pressure": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "mbar",
            "mapping": "node",
            "doc": "surface saturation vapor pressure",
        },  # e_sat_surf
        "land_surface_air_water-vapor__partial_pressure": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "mbar",
            "mapping": "node",
            "doc": "surface vapor pressure",
        },  # e_surf
        "atmosphere_bottom_air_water-vapor__saturated_partial_pressure": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "mbar",
            "mapping": "node",
            "doc": "air saturation vapor pressure",
        },  # e_sat_air
        "atmosphere_bottom_air_water-vapor__partial_pressure": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "mbar",
            "mapping": "node",
            "doc": "air vapor pressure",
        },  # e_air
        "atmosphere_bottom_air_water-vapor__dew_point_temperature": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "deg_C",
            "mapping": "node",
            "doc": "dew point temperature",
        },  # T_dew
        "atmosphere_bottom_air_flow__bulk_richardson_number": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "bulk Richardson number",
        },  # Ri
        "atmosphere_bottom_air__emissivity": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "air emissivity",
        },  # em_air
        "atmosphere_air-column_water-vapor__liquid-equivalent_depth": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "cm",
            "mapping": "node",
            "doc": "precipitable depth",
        },  # W_p
    }

    def __init__(
        self,
        grid,
        start_datetime,
        GMT_offset=0,
        rho_H2O=1000,
        rho_air=1.2614,
        Cp_air=1005.7,
        satterlund=False,
    ):
        """Initialize Meteorology component"""

        super().__init__(grid)

        # physical constants
        self._g = np.float64(9.81)  # [m s-2, gravity]
        self._kappa = np.float64(0.408)  # [1]  (von Karman)
        self._Lv = np.float64(2500000)  # [J kg-1] Latent heat of vaporization.
        self._Lf = np.float64(334000)  # [J kg-1 = W s kg-1], Latent heat of fusion
        self._sigma = np.float64(5.67e-8)  # [W m-2 K-4]  (Stefan-Boltzman constant)
        self._C_to_K = np.float64(273.15)  # (add to convert deg C to K)

        self._one_seventh = np.float64(1) / 7
        self._hours_per_day = np.float64(24)
        self._latent_heat_constant = np.float64(0.622)  # TODO: changed by me

        # parameters
        self._GMT_offset = GMT_offset
        self._rho_H2O = rho_H2O  # kg m-3
        self._rho_air = rho_air  # kg m-3
        self._Cp_air = Cp_air  # J kg-1 K-1
        self._satterlund = satterlund  # bool

        # datetime & Julian day
        try:
            self._datetime_obj = datetime.datetime.strptime(
                start_datetime, "%Y-%m-%d %H:%M:%S"
            )
        except Exception as e:
            raise e

        self._julian_day = solar_funcs.Julian_Day(
            self._datetime_obj.month,
            self._datetime_obj.day,
            self._datetime_obj.hour,
            self._datetime_obj.year,
        )

        # input fields
        self._T_air = grid.at_node["atmosphere_bottom_air__temperature"]
        self._T_surf = grid.at_node["land_surface__temperature"]
        self._lat_deg = grid.at_node["land_surface__latitude"]
        self._lon_deg = grid.at_node["land_surface__longitude"]

        if "land_surface__aspect_angle" in grid.at_node:
            self._alpha = grid.at_node["land_surface__aspect_angle"]
        else:
            self._alpha = grid.add_zeros("land_surface__aspect_angle", at="node")

        if "land_surface__slope_angle" in grid.at_node:
            self._beta = grid.at_node["land_surface__slope_angle"]
        else:
            self._beta = grid.add_zeros("land_surface__slope_angle", at="node")

        # if "snowpack__depth" in grid.at_node:
        #     self._h_snow = grid.at_node["snowpack__depth"]
        # else:
        #     self._h_snow = grid.add_zeros("snowpack__depth", at="node")

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
                "atmosphere_aerosol_dust__reduction_of_transmittance", 0.0
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
                "atmosphere_bottom_air_flow__reference-height_speed", 3
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
        self._Qc = grid.at_node[
            "snowpack_land_surface_net-conduction-heat__energy_flux"
        ]
        self._Qa = grid.at_node["snowpack_land_surface_net-advection-heat__energy_flux"]
        self._Dn = grid.at_node[
            "atmosphere_bottom_air__neutral_bulk_heat_aerodynamic_conductance"
        ]
        self._Dh = grid.at_node[
            "atmosphere_bottom_air__bulk_sensible_heat_aerodynamic_conductance"
        ]
        self._De = grid.at_node[
            "atmosphere_bottom_air__bulk_latent_heat_aerodynamic_conductance"
        ]
        self._e_sat_surf = grid.at_node[
            "land_surface_air_water-vapor__saturated_partial_pressure"
        ]
        self._e_surf = grid.at_node["land_surface_air_water-vapor__partial_pressure"]
        self._e_sat_air = grid.at_node[
            "atmosphere_bottom_air_water-vapor__saturated_partial_pressure"
        ]
        self._e_air = grid.at_node[
            "atmosphere_bottom_air_water-vapor__partial_pressure"
        ]
        self._T_dew = grid.at_node[
            "atmosphere_bottom_air_water-vapor__dew_point_temperature"
        ]
        self._Ri = grid.at_node["atmosphere_bottom_air_flow__bulk_richardson_number"]
        self._em_air = grid.at_node["atmosphere_bottom_air__emissivity"]
        self._W_p = grid.at_node[
            "atmosphere_air-column_water-vapor__liquid-equivalent_depth"
        ]

    @property
    def GMT_offset(self):
        return self._GMT_offset

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

    def update_bulk_richardson_number(self):
        """calculate Ri"""
        # see Price & Dunne 1976
        top = self._g * self._z * (self._T_air - self._T_surf)  # TODO: changed by me
        bot = self._uz**2.0 * (self._T_air + self._C_to_K)
        self._Ri[:] = top / bot

    def update_bulk_aero_conductance(self):
        """calculate Dn, Dh, De"""

        # see Price & Dunne 1976
        # calculate Dn
        arg = self._kappa / np.log(self._z / self._z0_air)  # TODO: changed by me
        Dn = self._uz * arg**2.0

        # check if pixels are neutral
        w = self._T_air != self._T_surf  # (boolean array)
        nw = w.sum()

        # if all pixels are neutral, set Dh = De = Dn
        if nw == 0:
            self._Dn[:] = Dn
            self._Dh[:] = Dn
            self._De[:] = Dn
            return

        # if one or more pixels are not neutral, make correction using Ri
        Dh = Dn.copy()
        ws = self._Ri > 0  # If (Ri > 0) or (T_air > T_surf), then STABLE.
        wu = np.invert(ws)  # where unstable

        Dh[ws] = Dh[ws] / (1.0 + 10 * self._Ri[ws])
        Dh[wu] = Dh[wu] * (1.0 - 10 * self._Ri[wu])

        self._Dn[:] = Dn
        self._Dh[:] = Dh
        self._De[:] = Dh  # assumed equal

    def update_sensible_heat_flux(self):
        """calculate Qh"""

        # see Price & Dunne 1976
        delta_T = self._T_air - self._T_surf
        self._Qh[:] = (self._rho_air * self._Cp_air) * self._Dh * delta_T

    def update_saturation_vapor_pressure(self, MBAR=False, SURFACE=False):
        """calculate e_sat_surface, e_sat_air"""

        if SURFACE:
            T = self._T_surf
        else:
            T = self._T_air

        if not self._satterlund:
            # use Brutsaert method (Dingman 2015 p148)
            term1 = (np.float64(17.3) * T) / (T + np.float64(237.3))
            e_sat = np.float64(0.611) * np.exp(term1)  # [kPa]
        else:
            # use Satterlund method
            term1 = np.float64(2353) / (T + np.float64(273.15))
            e_sat = np.float64(10) ** (np.float64(11.4) - term1)  # [Pa]
            e_sat = e_sat / np.float64(1000)  # [kPa]

        if MBAR:
            e_sat = e_sat * np.float64(10)  # [mbar]

        if SURFACE:
            self._e_sat_surf[:] = e_sat
        else:
            self._e_sat_air[:] = e_sat

    def update_vapor_pressure(self, SURFACE=False):
        """calculate e_surf, e_air"""

        if SURFACE:
            self._e_surf[:] = self._e_sat_surf * self._RH
        else:
            self._e_air[:] = self._e_sat_air * self._RH

    def update_dew_point(self):
        """calculate T_dew"""

        # https: // en.wikipedia.org / wiki / Dew_point
        # e_air in mbar units

        a = 6.1121  # [mbar]
        b = 18.678
        c = 257.14  # [deg C]

        log_term = np.log(self._e_air / a)
        self._T_dew[:] = c * log_term / (b - log_term)  # [deg C]

    def update_precipitable_water_content(self):
        """calculate W_p"""

        arg = np.float64(0.0614 * self._T_dew)
        self._W_p[:] = np.float64(1.12) * np.exp(arg)  # [cm]

    def update_latent_heat_flux(self):
        """calculate Qe"""
        # see Price & Dunne 1976
        const = self._latent_heat_constant
        factor = self._rho_air * self._Lv * self._De
        delta_e = self._e_air - self._e_surf
        self._Qe[:] = factor * delta_e * (const / self._p0)

    def update_conduction_heat_flux(self):
        """calculate Qc"""

        # Qc is set as 0
        pass

    def update_advection_heat_flux(self):
        """calculate Qa"""

        # Qa is set as 0
        pass

    def update_julian_day(self, dt):
        """calculate julian_day, dt in seconds"""

        # Update the datetime_obj
        delta = relativedelta(seconds=dt)
        self._datetime_obj += delta

        # update Julian day
        self._julian_day = solar_funcs.Julian_Day(
            self._datetime_obj.month,
            self._datetime_obj.day,
            self._datetime_obj.hour,
            self._datetime_obj.year,
        )
        # get TSN_offset
        dec_part = self._julian_day - np.int16(self._julian_day)
        clock_hour = dec_part * self._hours_per_day
        solar_noon = solar_funcs.True_Solar_Noon(
            self._julian_day,
            self._lon_deg,
            self._GMT_offset,
            DST_offset=None,
            year=self._datetime_obj.year,
        )
        self._TSN_offset = clock_hour - solar_noon

    def update_net_shortwave_radiation(self):
        """calculate Qn_SW"""

        Qn_SW = solar_funcs.Clear_Sky_Radiation(
            self._lat_deg,
            self._julian_day,
            self._W_p,
            self._TSN_offset,
            self._alpha,
            self._beta,
            self._albedo,
            self._dust_atten,
        )

        # TODO: added by me (Dingman 2015 p227)
        tau_cloud = 0.355 + 0.68 * (1 - self._cloud_factor)
        tau_cloud[tau_cloud > 1] = 1
        tau_canopy = np.exp(-3.91 * self._canopy_factor)
        Qn_SW = tau_cloud * tau_canopy * Qn_SW

        self._Qn_SW[:] = Qn_SW

    def update_em_air(self):
        """calculate em_air"""

        T_air_K = self._T_air + self._C_to_K

        if not self._satterlund:
            # Brutsaert method (Dingman 2002 P196)
            e_air_kPa = self._e_air / np.float64(10)  # [kPa]
            F = self._canopy_factor
            C = self._cloud_factor
            term1 = (1.0 - F) * 1.72 * (e_air_kPa / T_air_K) ** self._one_seventh
            term2 = 1.0 + (0.22 * C**2.0)
            self._em_air[:] = (term1 * term2) + F
        else:
            # satterlund method
            e_air_mbar = self._e_air
            eterm = np.exp(-1 * (e_air_mbar) ** (T_air_K / 2016))
            self._em_air[:] = 1.08 * (1.0 - eterm)

    def update_net_longwave_radiation(self):
        """calculate Qn_LW"""
        # see Dingman 2015 p231
        T_air_K = self._T_air + self._C_to_K
        T_surf_K = self._T_surf + self._C_to_K
        LW_in = self._em_air * self._sigma * (T_air_K) ** 4.0
        LW_out = self._em_surf * self._sigma * (T_surf_K) ** 4.0

        # radiation from the air that is reflected from the surface
        LW_out += (1.0 - self._em_surf) * LW_in

        self._Qn_LW[:] = LW_in - LW_out

    def update_net_energy_flux(self):
        """calculate Q_sum"""

        self._Q_sum[:] = (
            self._Qn_SW + self._Qn_LW + self._Qh + self._Qe + self._Qa + self._Qc
        )

    def run_one_step(self, dt):
        # update input fields in case there is new input
        self._T_air = self._grid.at_node["atmosphere_bottom_air__temperature"]
        self._T_surf = self._grid.at_node["land_surface__temperature"]
        # self._h_snow = self._grid.at_node["snowpack__depth"]
        self._albedo = self._grid.at_node["land_surface__albedo"]
        self._em_surf = self._grid.at_node["land_surface__emissivity"]
        self._dust_atten = self._grid.at_node[
            "atmosphere_aerosol_dust__reduction_of_transmittance"
        ]
        self._canopy_factor = self._grid.at_node[
            "atmosphere_bottom_air__brutsaert_emissivity_canopy_factor"
        ]
        self._cloud_factor = self._grid.at_node[
            "atmosphere_bottom_air__brutsaert_emissivity_cloud_factor"
        ]
        self._RH = self._grid.at_node[
            "atmosphere_bottom_air_water-vapor__relative_saturation"
        ]
        self._p0 = self._grid.at_node["atmosphere_bottom_air__pressure"]
        self._z0_air = self._grid.at_node[
            "atmosphere_bottom_air_flow__log_law_roughness_length"
        ]
        self._z = self._grid.at_node[
            "atmosphere_bottom_air_flow__speed_reference_height"
        ]
        self._uz = self._grid.at_node[
            "atmosphere_bottom_air_flow__reference-height_speed"
        ]

        # update state variable
        self.update_bulk_richardson_number()  # Ri
        self.update_bulk_aero_conductance()  # Dn, Dh, De
        self.update_sensible_heat_flux()  # Qh
        self.update_saturation_vapor_pressure(MBAR=True)  # e_sat_air
        self.update_saturation_vapor_pressure(MBAR=True, SURFACE=True)  # e_sat_surf
        self.update_vapor_pressure()  # e_air
        self.update_dew_point()  # T_dew
        self.update_precipitable_water_content()  # W_p [after update_dew_point()]
        self.update_vapor_pressure(SURFACE=True)  # e_surf
        self.update_latent_heat_flux()  # Qe
        self.update_conduction_heat_flux()  # Qc
        self.update_advection_heat_flux()  # Qa
        self.update_julian_day(dt)  # julian_day
        self.update_net_shortwave_radiation()  # Qn_SW
        self.update_em_air()  # em_air
        self.update_net_longwave_radiation()  # Qn_LW
        self.update_net_energy_flux()  # Q_sum
