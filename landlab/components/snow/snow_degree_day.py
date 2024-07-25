"""Landlab component that simulates snowpack dynamics.

This component simulates snowmelt process using the degree-day method.
The code is implemented based on the TopoFlow snow component (by Scott D. Packham)
with some adjustments.

* https://github.com/peckhams/topoflow36/blob/master/topoflow/components/snow_degree_day.py
* https://github.com/peckhams/topoflow36/blob/master/topoflow/components/snow_base.py

@author: Tian Gan June 2024
"""

import numpy as np

from landlab import Component


SECONDS_PER_DAY = 86400.0
MM_PER_M = 1000.0


class SnowDegreeDay(Component):
    """Simulate snowmelt process using degree-day method.

    The degree-day method uses c0 (degree-day coefficient),
    t0 (degree-day threshold temperature) and t_air (air temperature)
    to determine the snow melt rate.
    When t_air is larger than t0, the melt process will start.
    The melt rate is calculated as sm = (t_air - t0) * c0.
    (e.g., units for c0 and sm could be mm day-1 K-1 and mm day-1 respectively.)

    Parameters
    ----------
    grid : ModelGrid
        A Landlab model grid object
    c0: float, optional
        Degree-day coefficient [mm / day / K].
    t0: float, optional
        Degree-day threshold temperature [deg_C].
    rho_water : float, optinal
        Water density [kg / m3].
    t_rain_snow : float, optional
        Temperature threshold for precipitation accurs as rain and snow with
        equal frequency [deg_C].

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.snow import SnowDegreeDay

    >>> grid = RasterModelGrid((2, 2))
    >>> grid.add_full("atmosphere_bottom_air__temperature", 2, at="node")
    array([2., 2., 2., 2.])
    >>> grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0, at="node")
    array([0.,  0.,  0.,  0.])
    >>> grid.add_full("snowpack__liquid-equivalent_depth", 0.005)  # m
    array([0.005,  0.005,  0.005,  0.005])

    >>> sm = SnowDegreeDay(grid, c0=3)
    >>> sm.run_one_step(8.64e4)  # 1 day in sec
    >>> grid.at_node["snowpack__melt_volume_flux"]
    array([5.78703704e-08,   5.78703704e-08,   5.78703704e-08,
         5.78703704e-08])
    >>> grid.at_node["snowpack__liquid-equivalent_depth"]
    array([0.,  0.,  0.,  0.])
    >>> sm.vol_sm
    0.02
    >>> sm.vol_swe
    0.0
    >>> sm.total_p_snow
    array([0.,  0.,  0.,  0.])
    >>> sm.total_sm
    array([0.005,  0.005,  0.005,  0.005])
    """

    _name = "SnowDegreeDay"

    _unit_agnostic = False

    _info = {
        # input fields
        "atmosphere_bottom_air__temperature": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "deg_C",
            "mapping": "node",
            "doc": "atmosphere bottom air temperature",
        },  # t_air
        "atmosphere_water__precipitation_leq-volume_flux": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/s",
            "mapping": "node",
            "doc": "precipitation rate (in water equivalent)",
        },  # p
        "snowpack__liquid-equivalent_depth": {
            "dtype": float,
            "intent": "inout",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "snow water equivalent depth",
        },  # h_swe
        "snowpack__z_mean_of_mass-per-volume_density": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "kg m-3",
            "mapping": "node",
            "doc": "snow density",
        },  # rho_snow
        # output fields
        "snowpack__depth": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "snow depth",
        },  # h_snow
        "snowpack__melt_volume_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/s",
            "mapping": "node",
            "doc": "snow melt volume flux",
        },  # sm
    }

    def __init__(self, grid, c0=2, t0=0, rho_water=1000, t_rain_snow=1):
        """Initialize SnowDegreeDay component."""

        super().__init__(grid)

        self.c0 = c0
        self._t0 = t0
        self._rho_water = rho_water
        self._t_rain_snow = t_rain_snow

        # calculated vars
        self._grid_area = grid.dx * grid.dy
        self._vol_sm = 0  # snowpack__domain_time_integral_of_melt_volume_flux
        self._vol_swe = 0  # snowpack__domain_integral_of_liquid-equivalent_depth
        # change for landlab
        self._total_p_snow = np.zeros(grid.number_of_nodes)  # add new variable
        self._total_sm = np.zeros(grid.number_of_nodes)  # add new variable

        if "snowpack__liquid-equivalent_depth" in grid.at_node:
            self._h_swe = grid.at_node["snowpack__liquid-equivalent_depth"]
            self._vol_swe = self._update_total_snowpack_water_volume(
                self._h_swe, self._grid_area
            )
        else:
            self._h_swe = grid.add_zeros("snowpack__liquid-equivalent_depth", at="node")

        if "snowpack__z_mean_of_mass-per-volume_density" in grid.at_node:
            self._rho_snow = grid.at_node["snowpack__z_mean_of_mass-per-volume_density"]
        else:
            self._rho_snow = grid.add_full(
                "snowpack__z_mean_of_mass-per-volume_density", 300
            )

        # output data fields
        super().initialize_output_fields()
        self._h_snow = grid.at_node["snowpack__depth"]
        self._h_snow[:] = self._h_swe[:] * self.density_ratio
        self._sm = grid.at_node["snowpack__melt_volume_flux"]

    @property
    def c0(self):
        return self._c0

    @c0.setter
    def c0(self, c0):
        if np.any(c0 < 0):
            raise ValueError("degree-day coefficent must be positive")
        self._c0 = c0 / (SECONDS_PER_DAY * MM_PER_M)

    @property
    def t0(self):
        return self._t0

    @t0.setter
    def t0(self, t0):
        self._t0 = t0

    @property
    def rho_water(self):
        return self._rho_water

    @rho_water.setter
    def rho_water(self, rho_water):
        assert rho_water > 0, "assign rho_water with positive value"
        self._rho_water = rho_water

    @property
    def t_rain_snow(self):
        return self._t_rain_snow

    @t_rain_snow.setter
    def t_rain_snow(self, t_rain_snow):
        self._t_rain_snow = t_rain_snow

    @property
    def p_snow(self):
        """precipitation as snow (in water equivalent) at current time step
        (units: m)."""
        return self.p_snow

    @property
    def total_p_snow(self):
        """accumulative precipitation as snow (in water equivalent) over the
        model time (units: m)."""
        return self._total_p_snow

    @property
    def total_sm(self):
        """accumulative snow melt (in water equivalent) over the
        model time (units: m)."""
        return self._total_sm

    @property
    def vol_sm(self):
        """Total volume of snow melt (in water equivalent) for the domain and
        over the model time (units: m^3)."""
        return self._vol_sm

    @property
    def vol_swe(self):
        """Total volume of snow water equivalent for the domain and
        over the model time (units: m^3)."""
        return self._vol_swe

    @property
    def density_ratio(self):
        """Ratio of the densities between water and snow."""
        return self._rho_water / self._rho_snow

    @staticmethod
    def calc_precip_snow(precip_rate, air_temp, t_rain_snow):
        """Calculate snow precipitation.

        Parameters
        ----------
        p : array-like
            Precipitation rate [m/s].
        air_temp : array-like
            Air temperature [degC].
        t_rain_snow : float
            Temperature threshold below which precipitation falls as
            snow [degC].
        """
        return np.asarray(precip_rate) * (np.asarray(air_temp) <= t_rain_snow)

    @staticmethod
    def calc_snow_melt_rate(air_temp, c0, threshold_temp, out=None):
        """Calculate snow melt rate.

        Parameters
        ----------
        air_temp : array-like
            Air temperature [degC].
        c0 : float
            Degree-day coefficient [m / s / degC].
        threshold_temp : float
            Threshold temperatue above which melt occures [degC].

        Returns
        -------
        array-like
            Melt rate [m / s]
        """
        c0 = np.asarray(c0)
        air_temp = np.asarray(air_temp)

        return np.clip(c0 * (air_temp - threshold_temp), a_min=0.0, a_max=None, out=out)

    @staticmethod
    def calc_swe(precip_rate, melt_rate, swe, dt=1.0, out=None):
        """Calculate snow water equivalent.

        Parameters
        ----------
        precip_rate : array-like
            Precipitation rate [m/s].
        melt_rate : array-like
            Melt rate [m/s].
        swe : array-like
            Snow water equivalent [m].
        dt : float, optional
            Time step [s].

        Returns
        -------
        array-like
            New snow water equivalent [m].
        """
        out = np.add(
            swe, (np.asarray(precip_rate) - np.asarray(melt_rate)) * dt, out=out
        )
        return out.clip(min=0.0, max=None, out=out)

    @staticmethod
    def calc_snow_depth(h_swe, density_ratio, out=None):
        """Calculate snow depth from snow water equivalent.

        Parameters
        ----------
        h_swe : array-like
            Snow water equivalent [m].
        density_ratio : array-like
            Density ratio of water to snow [-].

        Returns
        -------
        array-like
            Snow depth [m].
        """
        return np.multiply(h_swe, density_ratio, out=out)

    @staticmethod
    def _update_total_p_snow(dt, p_snow, total_p_snow):  # add for landlab (new method)
        update_total_p_snow = total_p_snow + p_snow * dt

        return update_total_p_snow

    @staticmethod
    def _update_total_snowmelt(dt, sm, total_sm):  # add for landlab (new method)
        update_total_sm = total_sm + sm * dt

        return update_total_sm

    @staticmethod
    def _update_sm_integral(total_sm, grid_area):
        update_vol_sm = np.sum(total_sm * grid_area)

        return update_vol_sm

    @staticmethod
    def _update_total_snowpack_water_volume(h_swe, grid_area):
        volume = h_swe * grid_area
        update_vol_swe = np.sum(volume)

        return update_vol_swe

    def run_one_step(self, dt):
        """Run component for a time step.

        Parameters
        ----------
        dt : float
            Duration to run the component for [s].
        """
        # update input fields in case there is new input
        air_temp = self._grid.at_node["atmosphere_bottom_air__temperature"]

        # update state variables
        precip_snow = SnowDegreeDay.calc_precip_snow(
            self._grid.at_node["atmosphere_water__precipitation_leq-volume_flux"],
            air_temp,
            self._t_rain_snow,
        )

        SnowDegreeDay.calc_snow_melt_rate(
            air_temp,
            self._c0,
            self._t0,
            out=self.grid.at_node["snowpack__melt_volume_flux"],
        )
        np.clip(
            self.grid.at_node["snowpack__melt_volume_flux"],
            a_min=None,
            a_max=self._h_swe / dt,
            out=self.grid.at_node["snowpack__melt_volume_flux"],
        )

        SnowDegreeDay.calc_swe(
            precip_snow,
            self.grid.at_node["snowpack__melt_volume_flux"],
            self._h_swe,
            dt=dt,
            out=self._h_swe,
        )
        SnowDegreeDay.calc_snow_depth(self._h_swe, self.density_ratio, out=self._h_snow)

        # add for landlab (new method)
        self._total_p_snow[:] = self._update_total_p_snow(
            dt, precip_snow, self._total_p_snow
        )
        # add for landlab (new method)
        self._total_sm[:] = self._update_total_snowmelt(dt, self._sm, self._total_sm)
        # after update_total_snowmelt()
        self._vol_sm = self._update_sm_integral(self._total_sm, self._grid_area)
        # after update_swe()
        self._vol_swe = self._update_total_snowpack_water_volume(
            self._h_swe, self._grid_area
        )
