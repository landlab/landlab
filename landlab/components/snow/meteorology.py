"""Landlab component that calculates the energy fluxes

This component calculates several energy fluxes (e.g., solar radiation,
long wave radiation, sensible heat, latent heat) to estimate the net total energy
flux that can be used as the input for the SnowEnergyBalance component.
The code is implemented based on the TopoFlow meteorology component
(by Scott D. Packham) with some adjustments.

* https://github.com/peckhams/topoflow36/blob/master/topoflow/components/met_base.py

@author: Tian Gan Sept 2024
"""

import datetime

import numpy as np
from dateutil.relativedelta import relativedelta

from landlab import Component

from . import vendor_solar_funcs

_G = 9.81  # gravitational acceleration [m / s2]
_SIGMA = 5.67e-8  # Stefan-Boltzman constant [W m-2 K-4]
_KAPPA = 0.408  # (von Karman)
_LATENT_HEAT_OF_VAPORIZATION = 2500000  # latent heat of vaporization [J / kg]
_LATENT_HEAT_OF_FUSION = 334000  # latent heat of fusion [J / kg] (to liquid)
_LATENT_HEAT_CONSTANT = 0.622
_C_TO_K = 273.15  # convert deg C to K [deg_K]
_KPA_TO_MBAR = 10.0  # convert kpa to mbar
_ONE_SEVENTH = 1 / 7
_HOURS_PER_DAY = 24


class Meteorology(Component):
    """Calculate several energy fluxes and climate variables.

    This component calculates several energy fluxes (solar radiation,
    long wave radiation, sensible heat, latent heat) and related climate variables
    (e.g., air vapor pressure, dew point, and air emissivity).

    The estimated total net energy flux (q_sum) can be used as the input for
    the SnowEnergyBalance component. q_sum is calculated as:

    q_sum = qn_sw + qn_lw + qh + qe + qa + qc
    - qn_sw: net short wave energy flux
    - qn_lw: net long wave energy flux
    - qh: sensible heat energy flux
    - qe: latent heat energy flux
    - qa: net energy flux advected by moving water (negligible; qa=0)
    - qc: net energy flux via conduction from snow to soil (negligible; qc=0)


    Parameters
    ----------
    grid : ModelGrid
        A Landlab model grid object.
    start_datetime : string
        Start time (local time) in format of "yyyy-mm-dd hh:mm:ss".
    gmt_offset : int or array_like, optional
        GMT offset at the location of interest. It should be entered as an integer or an
        array of integers between 0 and 12 with negative values for locations west of
        the prime meridian.
    rho_air : float, optional
        Air density [kg / m3].
    cp_air  : float, optional
        Air heat capacity [J / kg / K].
    roughness_length : float, optional
        Log law roughness length, which is the height at which the wind speed
        theoretically becomes zero due to surface roughness [m].
    reference_height : float, optional
        A specific height above the ground where wind speed is measured or modeled. [m]
    satterlund : bool, optional
        Indicate the method for calcuating saturation vapor pressure and air emissivity.
        If true, use Satterlund method. If false, use Brutsaert method.
    calc_input: bool, optional
        Indicate whether to calculate land surface temperature as input.
        If false, user provides land surface temperature data as input. If true, use dew
        point temperature to estimate its values.
    clear_sky: bool, optional
        Indicate whether to account for cloud and canopy factor impact.
        If false, the cloud and canopy factors will be used for the calculation of
        the net shortwave energy flux.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.snow import Meteorology
    >>> grid = RasterModelGrid((2, 2))
    >>> grid.add_full("atmosphere_bottom_air__temperature", 8, at="node")
    array([8., 8., 8., 8.])
    >>> grid.add_full("land_surface__temperature", 4, at="node")
    array([4., 4., 4., 4.])
    >>> grid.add_full("land_surface__latitude", 40, at="node")
    array([40., 40., 40., 40.])
    >>> grid.add_full("land_surface__longitude", -75, at="node")
    array([-75., -75., -75., -75.])
    >>> grid.add_full("atmosphere_bottom_air_flow__reference-height_speed", 2.0)
    array([2., 2., 2., 2.])
    >>> grid.add_full("atmosphere_bottom_air__pressure", 1000)
    array([1000., 1000., 1000., 1000.])
    >>> grid.add_full("atmosphere_bottom_air__brutsaert_emissivity_canopy_factor", 0.98)
    array([0.98, 0.98, 0.98, 0.98])
    >>> grid.add_full("atmosphere_bottom_air__brutsaert_emissivity_cloud_factor", 0.62)
    array([0.62, 0.62, 0.62, 0.62])
    >>> grid.add_full("atmosphere_aerosol_dust__reduction_of_transmittance", 0.1)
    array([0.1, 0.1, 0.1, 0.1])
    >>> grid.add_full("land_surface__slope_angle", 0.2)
    array([0.2, 0.2, 0.2, 0.2])
    >>> grid.add_full("land_surface__aspect_angle", 0.5)
    array([0.5, 0.5, 0.5, 0.5])
    >>> grid.add_full("land_surface__albedo", 0.3)
    array([0.3, 0.3, 0.3, 0.3])
    >>> dt = 60 * 60
    >>> met = Meteorology(grid, start_datetime="2023-01-02 12:00:00", gmt_offset=-5)
    >>> met.run_one_step(dt)
    >>> grid.at_node["land_surface_net-total-energy__energy_flux"]
    array([35.68442908, 35.68442908, 35.68442908, 35.68442908])
    """

    _name = "Meteorology"

    _unit_agnostic = False

    _info = {
        # input fields (15 var)
        "atmosphere_bottom_air__temperature": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "deg_C",
            "mapping": "node",
            "doc": "atmosphere bottom air temperature",
        },
        "land_surface__temperature": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "deg_C",
            "mapping": "node",
            "doc": "land surface temperature",
        },
        "land_surface__latitude": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "J kg-1 K-1",
            "mapping": "node",
            "doc": "latitude",
        },
        "land_surface__longitude": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "J kg-1 K-1",
            "mapping": "node",
            "doc": "longitude",
        },
        "land_surface__aspect_angle": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "radians",
            "mapping": "node",
            "doc": "land surface aspect",
        },
        "land_surface__slope_angle": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "radians",
            "mapping": "node",
            "doc": "land surface slope",
        },
        "snowpack__liquid-equivalent_depth": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "snow water equivalent depth",
        },
        "land_surface__albedo": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "land surface albedo",
        },
        "land_surface__emissivity": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "land surface emissivity",
        },
        "atmosphere_aerosol_dust__reduction_of_transmittance": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "aerosol dust reduction of transmittance",
        },
        "atmosphere_bottom_air__brutsaert_emissivity_canopy_factor": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "brutsaert emissivity canopy factor",
        },
        "atmosphere_bottom_air__brutsaert_emissivity_cloud_factor": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "brutsaert emissivity cloud factor",
        },
        "atmosphere_bottom_air_water-vapor__relative_saturation": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "air water vapor relative saturation",
        },
        "atmosphere_bottom_air__pressure": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "mbar",
            "mapping": "node",
            "doc": "bottom air pressure",
        },
        "atmosphere_bottom_air_flow__reference-height_speed": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "m s-1",
            "mapping": "node",
            "doc": "air flow speed at reference height",
        },  # uz
        # output fields (15 var)
        "land_surface_net-total-energy__energy_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W m-2",
            "mapping": "node",
            "doc": "net total energy flux",
        },
        "land_surface_net-shortwave-radiation__energy_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W m-2",
            "mapping": "node",
            "doc": "net shortwave radiation energy flux",
        },
        "land_surface_net-longwave-radiation__energy_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W m-2",
            "mapping": "node",
            "doc": "net longwave radiation energy flux",
        },
        "atmosphere_bottom_air_land_net-latent-heat__energy_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W m-2",
            "mapping": "node",
            "doc": "net latent heat energy flux",
        },
        "atmosphere_bottom_air_land_net-sensible-heat__energy_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "W m-2",
            "mapping": "node",
            "doc": "net sensible heat energy flux",
        },
        "atmosphere_bottom_air__bulk_sensible_heat_aerodynamic_conductance": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m s-1",
            "mapping": "node",
            "doc": "bulk sensible heat aerodynamic conductance",
        },
        "atmosphere_bottom_air__bulk_latent_heat_aerodynamic_conductance": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m s-1",
            "mapping": "node",
            "doc": "bulk latent heat aerodynamic conductance",
        },
        "land_surface_air_water-vapor__saturated_partial_pressure": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "mbar",
            "mapping": "node",
            "doc": "surface saturation vapor pressure",
        },
        "land_surface_air_water-vapor__partial_pressure": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "mbar",
            "mapping": "node",
            "doc": "surface vapor pressure",
        },
        "atmosphere_bottom_air_water-vapor__saturated_partial_pressure": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "mbar",
            "mapping": "node",
            "doc": "air saturation vapor pressure",
        },
        "atmosphere_bottom_air_water-vapor__partial_pressure": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "mbar",
            "mapping": "node",
            "doc": "air vapor pressure",
        },
        "atmosphere_bottom_air_water-vapor__dew_point_temperature": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "deg_C",
            "mapping": "node",
            "doc": "dew point temperature",
        },
        "atmosphere_bottom_air_flow__bulk_richardson_number": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "bulk Richardson number",
        },
        "atmosphere_bottom_air__emissivity": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "air emissivity",
        },
        "atmosphere_air-column_water-vapor__liquid-equivalent_depth": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "cm",
            "mapping": "node",
            "doc": "precipitable depth",
        },
    }

    def __init__(
        self,
        grid,
        start_datetime,
        gmt_offset=0,
        rho_air=1.2614,
        cp_air=1005.7,
        roughness_length=0.02,
        reference_height=10,
        satterlund=False,
        calc_input=False,
        clear_sky=False,
    ):
        """Initialize Meteorology component"""

        super().__init__(grid)
        self._gmt_offset = gmt_offset
        self._rho_air = rho_air  # kg m-3
        self._cp_air = cp_air  # J kg-1 K-1
        self._roughness_length = roughness_length  # m
        self._reference_height = reference_height  # m
        self._satterlund = satterlund  # bool
        self._calc_input = calc_input  # bool
        self._clear_sky = clear_sky

        self._datetime_obj = datetime.datetime.strptime(
            start_datetime, "%Y-%m-%d %H:%M:%S"
        )

        # input fields
        if "land_surface__aspect_angle" not in grid.at_node:
            grid.add_zeros("land_surface__aspect_angle", at="node")

        if "land_surface__slope_angle" not in grid.at_node:
            grid.add_zeros("land_surface__slope_angle", at="node")

        if "snowpack__liquid-equivalent_depth" not in grid.at_node:
            grid.add_zeros("snowpack__liquid-equivalent_depth", at="node")

        if "land_surface__albedo" not in grid.at_node:
            grid.add_full("land_surface__albedo", 0.3, at="node")

        if "land_surface__emissivity" not in grid.at_node:
            grid.add_full("land_surface__emissivity", 0.98, at="node")

        if "atmosphere_aerosol_dust__reduction_of_transmittance" not in grid.at_node:
            grid.add_zeros(
                "atmosphere_aerosol_dust__reduction_of_transmittance", at="node"
            )

        if (
            "atmosphere_bottom_air__brutsaert_emissivity_canopy_factor"
            not in grid.at_node
        ):
            grid.add_zeros(
                "atmosphere_bottom_air__brutsaert_emissivity_canopy_factor", at="node"
            )

        if (
            "atmosphere_bottom_air__brutsaert_emissivity_cloud_factor"
            not in grid.at_node
        ):
            grid.add_zeros(
                "atmosphere_bottom_air__brutsaert_emissivity_cloud_factor", at="node"
            )

        if "atmosphere_bottom_air_water-vapor__relative_saturation" not in grid.at_node:
            grid.add_full(
                "atmosphere_bottom_air_water-vapor__relative_saturation", 0.5, at="node"
            )

        if "atmosphere_bottom_air__pressure" not in grid.at_node:
            grid.add_full("atmosphere_bottom_air__pressure", 1013.25, at="node")

        if "atmosphere_bottom_air_flow__reference-height_speed" not in grid.at_node:
            grid.add_full(
                "atmosphere_bottom_air_flow__reference-height_speed", 3, at="node"
            )

        # output fields
        super().initialize_output_fields()

    @property
    def gmt_offset(self):
        return self._gmt_offset

    @gmt_offset.setter
    def gmt_offset(self, gmt_offset):
        if not isinstance(gmt_offset, int) or gmt_offset < -12 or gmt_offset > 12:
            raise ValueError("GMT offset must be an integer between -12 and 12.")
        self._gmt_offset = gmt_offset

    @property
    def rho_air(self):
        return self._rho_air

    @rho_air.setter
    def rho_air(self, rho_air):
        if rho_air <= 0:
            raise ValueError("air density must be positive")
        self._rho_air = rho_air

    @property
    def cp_air(self):
        return self._cp_air

    @cp_air.setter
    def cp_air(self, cp_air):
        if cp_air <= 0:
            raise ValueError("air heat capacity must be positive")
        self._cp_air = cp_air

    @property
    def roughness_length(self):
        return self._roughness_length

    @roughness_length.setter
    def roughness_length(self, roughness_length):
        if roughness_length <= 0:
            raise ValueError("roughness length must be positive")
        self._roughness_length = roughness_length

    @property
    def reference_height(self):
        return self._reference_height

    @reference_height.setter
    def reference_height(self, reference_height):
        if reference_height <= 0:
            raise ValueError("reference height must be positive")
        self._reference_height = reference_height

    @property
    def satterlund(self):
        return self._satterlund

    @satterlund.setter
    def satterlund(self, satterlund):
        if not isinstance(satterlund, bool):
            raise ValueError("satterlund must be a boolean value")
        self._satterlund = satterlund

    @property
    def calc_input(self):
        return self._calc_input

    @calc_input.setter
    def calc_input(self, calc_input):
        if not isinstance(calc_input, bool):
            raise ValueError("calc_input must be a boolean value")
        self._calc_input = calc_input

    @property
    def clear_sky(self):
        return self._clear_sky

    @clear_sky.setter
    def clear_sky(self, clear_sky):
        if not isinstance(clear_sky, bool):
            raise ValueError("clear_sky must be a boolean value")
        self._clear_sky = clear_sky

    @property
    def julian_day(self):

        return vendor_solar_funcs.Julian_Day(
            self._datetime_obj.month,
            self._datetime_obj.day,
            self._datetime_obj.hour,
            self._datetime_obj.year,
        )

    @staticmethod
    def calc_bulk_richardson_number(
        reference_height,
        air_temp,
        surf_temp,
        wind_speed,
        out=None,
    ):
        """Calculate bulk richardson number (Ri)

        Parameters
        ----------

        reference_height : float
            A specific height above the ground where wind speed
            is measured or modeled [m].
        air_temp : array_like
            Air temperature [deg_C].
        surf_temp : array_like
            Land surface temperature [deg_C].
        wind_speed : array_like
            The wind speed measured or modeled at the reference height [m / s].

        Returns
        -------
        array_like
            Bulk richardson number.
        """

        # see Price & Dunne 1976
        top = _G * reference_height * (np.asarray(air_temp) - np.asarray(surf_temp))
        bot = np.asarray(wind_speed) ** 2 * (np.asarray(air_temp) + _C_TO_K)

        return np.divide(top, bot, out=out)

    @staticmethod
    def calc_bulk_aero_conductance(
        reference_height,
        roughness_length,
        richardson_number,
        wind_speed,
        air_temp,
        surf_temp,
    ):
        """calculate air bulk aerodynamic conductance

        Parameters
        ----------
        reference_height : float
            A specific height above the ground where wind speed
            is measured or modeled [m].
        roughness_length : float
            Log law roughness length, which is the height at which the wind speed
            theoretically becomes zero due to surface roughness [m].
        richardson_number : aray_like
            Bulk richardson number.
        wind_speed : array_like
            The wind speed measured or modeled at the reference height [m / s].
        air_temp : array_like
            Air temperature [deg_C].
        surf_temp : array_like
            Land surface temperature [deg_C].

        Returns
        -------
        array_like
            Bulk sensible or latent heat aerodynamic conductance [m / s].
        """

        # changed for landlab:
        # - changed how arg is calculated, remove h_snow
        # - changed to return only one value to represent Dh and De (Dh=De)

        # see Price & Dunne 1976
        # calculate Dn
        arg = _KAPPA / np.log(reference_height / roughness_length)
        conductance_neutral = np.asarray(wind_speed) * arg**2.0

        # check if pixels are neutral
        w = np.asarray(air_temp) != np.asarray(surf_temp)
        nw = np.asarray(w).sum()

        # if all pixels are neutral, set Dh = De = Dn
        aero_conductance = np.asarray(conductance_neutral)

        # if one or more pixels are not neutral, make correction using Ri, set Dh = De
        if nw != 0:
            ws = (
                np.asarray(richardson_number) > 0
            )  # If (Ri > 0) or (T_air > T_surf), then STABLE.
            wu = np.invert(ws)  # where unstable

            aero_conductance[ws] = aero_conductance[ws] / (
                1.0 + 10 * np.asarray(richardson_number)[ws]
            )
            aero_conductance[wu] = aero_conductance[wu] * (
                1.0 - 10 * np.asarray(richardson_number)[wu]
            )

        return aero_conductance

    @staticmethod
    def calc_sensible_heat_flux(
        rho_air, cp_air, aero_conductance, air_temp, surf_temp, out=None
    ):
        """Calculate net sensible heat energy flux

        Parameters
        ----------
        rho_air : float
            Air density [kg / m3].
        cp_air  : float
            Air heat capacity [J / kg / K].
        aero_conductance:array_like
            Bulk sensible heat aerodynamic conductance [m / s].
        air_temp : array_like
            Air temperature [deg_C].
        surf_temp : array_like
            Land surface temperature [deg_C].

        Returns
        -------
        array_like
            Net sensible heat energy flux [W / m2].
        """

        # see Price & Dunne 1976
        return np.multiply(
            np.asarray(air_temp) - np.asarray(surf_temp),
            rho_air * cp_air * np.asarray(aero_conductance),
            out=out,
        )

    @staticmethod
    def calc_saturation_vapor_pressure(
        temp, satterlund=False, millibar=False, out=None
    ):
        """Calculate saturation vapor pressure

        Parameters
        ----------

        temp : array_like
            Air temperature [deg_C].
        satterlund : bool, optional
            If true, use Satterlund method. If false, use Brutsaert method.
        millibar : bool, optional
            If true, the units for the output is millibar. If false,
            the units is kPa.

        Returns
        -------
        array_like
            Saturation vapor pressure [kPa or millibar].
        """

        if not satterlund:
            # use Brutsaert method (Dingman 2015 p148)
            out = np.multiply(
                0.611,
                np.exp((17.3 * np.asarray(temp)) / (np.asarray(temp) + 237.3)),
                out=out,
            )
        else:
            # use Satterlund method
            out = np.divide(
                10 ** (11.4 - (2353 / (np.asarray(temp) + _C_TO_K))), 1000, out=out
            )

        if millibar:
            out = np.multiply(out, _KPA_TO_MBAR, out=np.asarray(out))

        return out

    @staticmethod
    def calc_vapor_pressure(sat_vapor_pressure, relative_humidity, out=None):
        """Calculate vapor pressure

        Parameters
        ----------

        sat_vapor_pressure : array_like
            Saturation vapor pressure [kPa or millibar].
        relative_humidity : array_like
            Relative humidity.

        Returns
        -------
        array_like
            Vapor pressure [kPa or millibar].
        """

        return np.multiply(sat_vapor_pressure, relative_humidity, out=out)

    @staticmethod
    def calc_dew_point(air_vapor_pressure, out=None):
        """Calculate dew point temperature

        Parameters
        ----------
        air_vapor_pressure : array_like
            Air vapor pressure [millibar].

        Returns
        -------
        array_like
            Dew point temperature [deg_C].
        """

        # https: // en.wikipedia.org / wiki / Dew_point

        a = 6.1121  # [mbar]
        b = 18.678
        c = 257.14  # [deg C]

        log_term = np.log(np.asarray(air_vapor_pressure) / a)

        return np.divide(c * log_term, b - log_term, out=out)

    @staticmethod
    def calc_surf_temp(dew_point, h_swe, out=None):
        """Calculate land surface temperature

        Parameters
        ----------
        dew_point : array_like
            Dew point temperature [deg_C].
        h_swe : array_like
            Snow water equivalent [m].

        Returns
        -------
        array_like
            Land surface temperature [deg_C].
        """

        result = np.where(np.asarray(h_swe) > 0, np.minimum(dew_point, 0), dew_point)

        if out is not None:
            out[:] = result
            return out
        else:
            return result

    @staticmethod
    def calc_precipitable_water_content(dew_point, out=None):
        """Calculate precipitable water depth

        Parameters
        ----------
        dew point : array_like
            Dew point temperature [deg_C].

        Returns
        -------
        array_like
            Precipitable water content [cm].
        """

        return np.multiply(np.exp(0.0614 * np.asarray(dew_point)), 1.12, out=out)

    @staticmethod
    def calc_latent_heat_flux(
        rho_air,
        aero_conductance,
        air_pressure,
        air_vapor_pressure,
        surf_vapor_pressure,
        out=None,
    ):
        """Calculate net latent heat flux

        Parameters
        ----------
        rho_air : float
            Air density [kg / m3].
        aero_conductance :array_like
            Bulk latent heat aerodynamic conductance [m / s].
        air_pressure : array_like
            Bottom air pressure [millibar].
        air_vapor_pressure : array_like
            Air vapor pressure [millibar].
        surf_vapor_pressure : array_like
            Surface vapor pressure [millibar].

        Returns
        -------
        array_like
            Net latent heat energy flux [W / m2].
        """

        # see Price & Dunne 1976
        term1 = rho_air * _LATENT_HEAT_OF_VAPORIZATION * np.asarray(aero_conductance)

        term2 = (np.asarray(air_vapor_pressure) - np.asarray(surf_vapor_pressure)) * (
            _LATENT_HEAT_CONSTANT / np.asarray(air_pressure)
        )

        return np.multiply(term1, term2, out=out)

    @staticmethod
    def calc_tsn_offset(julian_day, year, longitude, gmt_offset=0, dts_offset=None):
        """Calculate true solar noon offset.

        Parameters
        ----------
        julian_day : float
            Julian day for a given datetime.
        year : int
            Year for a given datetime.
        longitude : array_like
            longitude at the location of interest [degree].
        gmt_offset : int or array_like, optional
            GMT offset at the location of interest. It should be entered as an integer
            or an array of integers between 0 and 12 with negative values for locations
            west of the prime meridian.
        dts_offset : int, optional
            Daylight Saving Time offset [hour].

        Returns
        -------
        array_like
            True solar noon offset [hour].
        """
        # Changed for landlab. Removed the code for julian day calculation
        dec_part = julian_day - np.int16(julian_day)
        clock_hour = dec_part * _HOURS_PER_DAY
        solar_noon = vendor_solar_funcs.True_Solar_Noon(
            julian_day,
            np.asarray(longitude),
            np.asarray(gmt_offset),
            DST_offset=dts_offset,
            year=year,
        )
        tsn_offset = clock_hour - solar_noon

        return tsn_offset

    @staticmethod
    def calc_net_shortwave_radiation(
        julian_day,
        tsn_offset,
        latitude,
        prec_water,
        alpha,
        beta,
        albedo,
        dust_atten,
        canopy_factor,
        cloud_factor,
        clear_sky=True,
        out=None,
    ):
        """Calculate net shortwave radiation energy flux

        Parameters
        ----------
        julian_day: float
            Julian day for a given datetime.
        tsn_offset: float
            True solar noon offset. [hour]
        latitude: array_like
            Latitude at the location of interest [degree].
        prec_water: array_like
            Precipitable water content [cm].
        alpha: array_like
            Land surface aspect [radians].
        beta: array_like
            Land surface slope [radians].
        albedo : array_like
            Land surface albedo.
        dust_atten: array_like
            Aerosol dust reduction of transmittance.
        canopy_factor: array_like
            Brutsaert emissivity canopy factor.
        cloud_factor: array_like
            Brutsaert emissivity cloud factor.
        clear_sky: bool, optional
            Indicate whether to account for cloud and canopy factor impact.
            If false, the cloud and canopy factors will be used for calculation.

        Returns
        -------
        array_like
            Net shortwave energy flux [W / m2].
        """

        k_cs = vendor_solar_funcs.Clear_Sky_Radiation(
            np.asarray(latitude),
            julian_day,
            np.asarray(prec_water),
            tsn_offset,
            np.asarray(alpha),
            np.asarray(beta),
            np.asarray(albedo),
            np.asarray(dust_atten),
        )

        qn_sw = k_cs * (1 - np.asarray(albedo))

        # changed for landlab: add cloud and canopy factor impact (Dingman 2015 p227)
        if clear_sky:
            tau_cloud = tau_canopy = 1
        else:
            tau_cloud = np.asarray(0.355 + 0.68 * (1 - np.asarray(cloud_factor)))
            tau_cloud[tau_cloud > 1] = 1
            tau_canopy = np.exp(-3.91 * np.asarray(canopy_factor))

        return np.multiply(tau_cloud * tau_canopy, qn_sw, out=out)

    @staticmethod
    def calc_air_emissivity(
        air_temp,
        air_vapor_pressure,
        canopy_factor,
        cloud_factor,
        satterlund=False,
        out=None,
    ):
        """Calculate air emissivity

        Parameters
        ----------
        air_temp : array_like
            Air temperature [deg_C]
        air_vapor_pressure : array_like
            Air vapor pressure [millibar].
        cloud_factor : array_like
            Brutsaert emissivity cloud factor.
        canopy_factor : array_like
            Brutsaert emissivity canopy factor.
        satterlund : bool, optional
            If true, use Satterlund method. If false, use Brutsaert method.

        Returns
        -------
        array_like
            Air emissivity [-].
        """

        air_temp_k = np.asarray(air_temp) + _C_TO_K

        if not satterlund:
            # Brutsaert method (Dingman 2002 P196)
            air_vapor_pressure_kpa = np.asarray(air_vapor_pressure) / _KPA_TO_MBAR
            term1 = (
                (1.0 - np.asarray(canopy_factor))
                * 1.72
                * (air_vapor_pressure_kpa / air_temp_k) ** _ONE_SEVENTH
            )
            term2 = 1.0 + (0.22 * np.asarray(cloud_factor) ** 2.0)

            out = np.add(term1 * term2, canopy_factor, out=out)

        else:
            # satterlund method
            term = np.exp(-1 * (np.asarray(air_vapor_pressure)) ** (air_temp_k / 2016))
            out = np.multiply(1.08, 1.0 - term, out=out)

        return out

    @staticmethod
    def calc_net_longwave_radiation(
        air_temp, surf_temp, air_emissivity, surf_emissivity, out=None
    ):
        """Calculate net longwave radiation energy flux

        Parameters
        ----------
        air_temp : array_like
            Air temperature [deg_C]
        surf_temp : array_like
            Land surface temperature [deg_C].
        air_emissivity : array_like
            Air emissivity.
        surf_emissivity : array_like
            Land surface emissivity.

        Returns
        -------
        array_like
            Net longwave radiation energy flux [W / m2].
        """
        # see Dingman 2015 p231.
        # This function doesn't account for cloud and canopy factor impact
        air_temp_k = np.asarray(air_temp) + _C_TO_K
        surf_temp_k = np.asarray(surf_temp) + _C_TO_K
        lw_in = np.asarray(air_emissivity) * _SIGMA * air_temp_k**4.0
        lw_out = np.asarray(surf_emissivity) * _SIGMA * surf_temp_k**4.0
        # radiation from the air that is reflected from the surface
        lw_out += (1.0 - np.asarray(surf_emissivity)) * lw_in

        return np.subtract(lw_in, lw_out, out=out)

    @staticmethod
    def calc_net_energy_flux(
        shortwave_energy_flux,
        longwave_energy_flux,
        sensible_heat_flux,
        latent_heat_flux,
        out=None,
    ):
        """Calculate net total energy flux

        Parameters
        ----------
        shortwave_energy_flux: array_like
            Net shortwave radiation energy flux [W / m2].
        longwave_energy_flux: array_like
            Net longwave radiation energy flux [W / m2].
        sensible_heat_flux:array_like
            Net sensible heat energy flux [W / m2].
        latent_heat_flux: array_like
            Net latent heat energy flux [W / m2].

        Returns
        -------
        array_like
            Net total energy flux [W / m2].
        """

        return np.add(
            np.asarray(sensible_heat_flux) + np.asarray(latent_heat_flux),
            np.asarray(shortwave_energy_flux) + np.asarray(longwave_energy_flux),
            out=out,
        )

    def run_one_step(self, dt):
        # update state variables
        Meteorology.calc_saturation_vapor_pressure(
            self.grid.at_node["atmosphere_bottom_air__temperature"],
            satterlund=self._satterlund,
            millibar=True,
            out=self.grid.at_node[
                "atmosphere_bottom_air_water-vapor__saturated_partial_pressure"
            ],
        )

        Meteorology.calc_vapor_pressure(
            self.grid.at_node[
                "atmosphere_bottom_air_water-vapor__saturated_partial_pressure"
            ],
            self.grid.at_node["atmosphere_bottom_air_water-vapor__relative_saturation"],
            out=self.grid.at_node[
                "atmosphere_bottom_air_water-vapor__partial_pressure"
            ],
        )

        Meteorology.calc_dew_point(
            self.grid.at_node["atmosphere_bottom_air_water-vapor__partial_pressure"],
            out=self.grid.at_node[
                "atmosphere_bottom_air_water-vapor__dew_point_temperature"
            ],
        )

        if self._calc_input:  # changed for landlab: add calc_input parameter
            Meteorology.calc_surf_temp(
                self.grid.at_node[
                    "atmosphere_bottom_air_water-vapor__dew_point_temperature"
                ],
                self.grid.at_node["snowpack__liquid-equivalent_depth"],
                out=self.grid.at_node["land_surface__temperature"],
            )

        Meteorology.calc_saturation_vapor_pressure(
            self.grid.at_node["land_surface__temperature"],
            satterlund=self._satterlund,
            millibar=True,
            out=self.grid.at_node[
                "land_surface_air_water-vapor__saturated_partial_pressure"
            ],
        )

        Meteorology.calc_bulk_richardson_number(
            self._reference_height,
            self.grid.at_node["atmosphere_bottom_air__temperature"],
            self.grid.at_node["land_surface__temperature"],
            self.grid.at_node["atmosphere_bottom_air_flow__reference-height_speed"],
            out=self.grid.at_node["atmosphere_bottom_air_flow__bulk_richardson_number"],
        )

        aero_conductance = Meteorology.calc_bulk_aero_conductance(
            self._reference_height,
            self._roughness_length,
            self.grid.at_node["atmosphere_bottom_air_flow__bulk_richardson_number"],
            self.grid.at_node["atmosphere_bottom_air_flow__reference-height_speed"],
            self.grid.at_node["atmosphere_bottom_air__temperature"],
            self.grid.at_node["land_surface__temperature"],
        )
        self.grid.at_node[
            "atmosphere_bottom_air__bulk_sensible_heat_aerodynamic_conductance"
        ][:] = aero_conductance
        self.grid.at_node[
            "atmosphere_bottom_air__bulk_latent_heat_aerodynamic_conductance"
        ][:] = aero_conductance

        Meteorology.calc_sensible_heat_flux(
            self._rho_air,
            self._cp_air,
            self.grid.at_node[
                "atmosphere_bottom_air__bulk_sensible_heat_aerodynamic_conductance"
            ],
            self.grid.at_node["atmosphere_bottom_air__temperature"],
            self.grid.at_node["land_surface__temperature"],
            out=self.grid.at_node[
                "atmosphere_bottom_air_land_net-sensible-heat__energy_flux"
            ],
        )

        Meteorology.calc_precipitable_water_content(
            self.grid.at_node[
                "atmosphere_bottom_air_water-vapor__dew_point_temperature"
            ],
            out=self.grid.at_node[
                "atmosphere_air-column_water-vapor__liquid-equivalent_depth"
            ],
        )

        Meteorology.calc_vapor_pressure(
            self.grid.at_node[
                "land_surface_air_water-vapor__saturated_partial_pressure"
            ],
            self.grid.at_node["atmosphere_bottom_air_water-vapor__relative_saturation"],
            out=self.grid.at_node["land_surface_air_water-vapor__partial_pressure"],
        )

        Meteorology.calc_latent_heat_flux(
            self._rho_air,
            self.grid.at_node[
                "atmosphere_bottom_air__bulk_latent_heat_aerodynamic_conductance"
            ],
            self.grid.at_node["atmosphere_bottom_air__pressure"],
            self.grid.at_node["atmosphere_bottom_air_water-vapor__partial_pressure"],
            self.grid.at_node["land_surface_air_water-vapor__partial_pressure"],
            out=self.grid.at_node[
                "atmosphere_bottom_air_land_net-latent-heat__energy_flux"
            ],
        )

        tsn_offset = Meteorology.calc_tsn_offset(
            self.julian_day,
            self._datetime_obj.year,
            self.grid.at_node["land_surface__longitude"],
            self._gmt_offset,
        )

        Meteorology.calc_net_shortwave_radiation(
            self.julian_day,
            tsn_offset,
            self.grid.at_node["land_surface__latitude"],
            self.grid.at_node[
                "atmosphere_air-column_water-vapor__liquid-equivalent_depth"
            ],
            self.grid.at_node["land_surface__aspect_angle"],
            self.grid.at_node["land_surface__slope_angle"],
            self.grid.at_node["land_surface__albedo"],
            self.grid.at_node["atmosphere_aerosol_dust__reduction_of_transmittance"],
            self.grid.at_node[
                "atmosphere_bottom_air__brutsaert_emissivity_canopy_factor"
            ],
            self.grid.at_node[
                "atmosphere_bottom_air__brutsaert_emissivity_cloud_factor"
            ],
            clear_sky=self._clear_sky,
            out=self.grid.at_node["land_surface_net-shortwave-radiation__energy_flux"],
        )

        Meteorology.calc_air_emissivity(
            self.grid.at_node["atmosphere_bottom_air__temperature"],
            self.grid.at_node["atmosphere_bottom_air_water-vapor__partial_pressure"],
            self.grid.at_node[
                "atmosphere_bottom_air__brutsaert_emissivity_canopy_factor"
            ],
            self.grid.at_node[
                "atmosphere_bottom_air__brutsaert_emissivity_cloud_factor"
            ],
            satterlund=self._satterlund,
            out=self.grid.at_node["atmosphere_bottom_air__emissivity"],
        )

        Meteorology.calc_net_longwave_radiation(
            self.grid.at_node["atmosphere_bottom_air__temperature"],
            self.grid.at_node["land_surface__temperature"],
            self.grid.at_node["atmosphere_bottom_air__emissivity"],
            self.grid.at_node["land_surface__emissivity"],
            out=self.grid.at_node["land_surface_net-longwave-radiation__energy_flux"],
        )

        Meteorology.calc_net_energy_flux(
            self.grid.at_node[
                "atmosphere_bottom_air_land_net-sensible-heat__energy_flux"
            ],
            self.grid.at_node[
                "atmosphere_bottom_air_land_net-latent-heat__energy_flux"
            ],
            self.grid.at_node["land_surface_net-shortwave-radiation__energy_flux"],
            self.grid.at_node["land_surface_net-longwave-radiation__energy_flux"],
            out=self.grid.at_node["land_surface_net-total-energy__energy_flux"],
        )

        time_delta = relativedelta(seconds=dt)
        self._datetime_obj += time_delta
