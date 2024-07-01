"""Landlab component that simulates snowpack dynamics.

This component simulates snowmelt process using the degree-day method.
The code is implemented based on the TopoFlow snow component (by Scott D. Packham)
with some adjustments.
https://github.com/peckhams/topoflow36/blob/master/topoflow/components/snow_degree_day.py
https://github.com/peckhams/topoflow36/blob/master/topoflow/components/snow_base.py

@author: Tian Gan June 2024

"""

import numpy as np

from landlab import Component


class SnowDegreeDay(Component):
    """Simulate snowmelt process using degree-day method.

      The degree-day method uses c0 (degree-day coefficient),
      T0 (degree-day threshold temperature) and T_air (air temperature)
      to determine the snow melt rate.
      When T_air is larger than T0, the melt process will start.
      The melt rate is calculated as SM = (T_air -T0) * c0

    Parameters
    ----------
    grid : ModelGrid
        A Landlab model grid object
    rho_H2O : float (default 1000 kg/m3)
        water density
    T_rain_snow : float (default 1 deg_C)
        temperature threshold for precipitation accurs as rain and snow with
        equal frequency

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.snow import SnowDegreeDay
    >>> grid = RasterModelGrid((2, 2), xy_spacing=(1, 1))
    >>> grid.add_full("atmosphere_bottom_air__temperature", 2, at="node")
    array([2., 2., 2., 2.])
    >>> grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0, at="node")
    array([0.,  0.,  0.,  0.])
    >>> grid.add_full("snowpack__liquid-equivalent_depth", 0.005)  # m
    array([0.005,  0.005,  0.005,  0.005])
    >>> grid.add_full("snowpack__degree-day_coefficient", 3)  # mm -day -K
    array([3.,  3.,  3.,  3.])
    >>> sm = SnowDegreeDay(grid)
    >>> sm.run_one_step(np.float64(8.64e4))  # 1 day in sec
    >>> grid.at_node["snowpack__melt_volume_flux"]
    array([5.78703704e-08,   5.78703704e-08,   5.78703704e-08,
         5.78703704e-08])
    >>> grid.at_node["snowpack__time_integral_melt_volume_flux"]
    array([0.005,  0.005,  0.005,  0.005])
    >>> grid.at_node["snowpack__liquid-equivalent_depth"]
    array([0.,  0.,  0.,  0.])
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
        },  # T_air
        "atmosphere_water__precipitation_leq-volume_flux": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/s",
            "mapping": "node",
            "doc": "precipitation rate (in water equivalent)",
        },  # P
        "snowpack__degree-day_coefficient": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "mm day-1 K-1",
            "mapping": "node",
            "doc": "degree-day coefficient",
        },  # c0
        "snowpack__degree-day_threshold_temperature": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "deg_C",
            "mapping": "node",
            "doc": "degree-day threshold temperature to start melt",
        },  # T0
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
        },
        "snowpack__melt_volume_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/s",
            "mapping": "node",
            "doc": "snow melt volume flux",
        },
        "atmosphere_water__time_integral_snowfall_leq-volume_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "total precipitation as snow (in water equivalent) during model time",
        },  # change for landlab (add new variable)
        "snowpack__time_integral_melt_volume_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "total snow melt (in water equivalent) during model time",
        },  # change for landlab (add new variable)
    }

    def __init__(self, grid, rho_H2O=1000, T_rain_snow=1):
        """Initialize SnowDegreeDay component"""

        super().__init__(grid)

        # parameters
        self._rho_H2O = rho_H2O
        self._T_rain_snow = T_rain_snow

        # calculated vars
        self._P_snow = None
        self._vol_SM = 0  # snowpack__domain_time_integral_of_melt_volume_flux
        self._vol_swe = 0  # snowpack__domain_integral_of_liquid-equivalent_depth
        self._grid_area = grid.dx * grid.dy

        # input data fields
        self._T_air = grid.at_node["atmosphere_bottom_air__temperature"]
        self._P = grid.at_node["atmosphere_water__precipitation_leq-volume_flux"]
        if "snowpack__degree-day_coefficient" in grid.at_node:
            self._c0 = grid.at_node["snowpack__degree-day_coefficient"]
        else:
            self._c0 = grid.add_full("snowpack__degree-day_coefficient", 2)

        if "snowpack__degree-day_threshold_temperature" in grid.at_node:
            self._T0 = grid.at_node["snowpack__degree-day_threshold_temperature"]
        else:
            self._T0 = grid.add_full("snowpack__degree-day_threshold_temperature", 0)

        if "snowpack__liquid-equivalent_depth" in grid.at_node:
            self._h_swe = grid.at_node["snowpack__liquid-equivalent_depth"]
            self._update_total_snowpack_water_volume()
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
        self._SM = grid.at_node["snowpack__melt_volume_flux"]
        self._total_P_snow = grid.at_node[
            "atmosphere_water__time_integral_snowfall_leq-volume_flux"
        ]  # add for landlab (new variable)
        self._total_SM = grid.at_node[
            "snowpack__time_integral_melt_volume_flux"
        ]  # add for landlab (new variable)

    @property
    def rho_H2O(self):
        return self._rho_H2O

    @rho_H2O.setter
    def rho_H2O(self, rho_H2O):
        assert rho_H2O > 0, "assign rho_H2O with positive value"
        self._rho_H2O = rho_H2O

    @property
    def T_rain_snow(self):
        return self._T_rain_snow

    @T_rain_snow.setter
    def T_rain_snow(self, T_rain_snow):
        self._T_rain_snow = T_rain_snow

    @property
    def vol_SM(self):
        return self._vol_SM

    @property
    def vol_swe(self):
        return self._vol_swe

    @property
    def density_ratio(self):
        return self._rho_H2O / self._rho_snow

    def _update_snowmelt_rate(self, dt):
        # melt rate based on degree-day coefficient
        sm_c0 = (self._c0 / np.float64(8.64e7)) * (
            self._T_air - self._T0  # convert mm -day -K to  m/s
        )
        sm_c0 = np.maximum(sm_c0, np.float64(0))

        # melt rate based on available swe (enforced max meltrate)
        sm_max = self._h_swe / dt  # original code use self._h_snow (line 526)

        # actual meltrate
        sm_act = np.minimum(sm_c0, sm_max)
        self._SM[:] = sm_act[:]

    def _update_prec_snow(self):
        self._P_snow = self._P * (self._T_air <= self._T_rain_snow)

    def _update_swe(self, dt):
        # increase due to precipitation as snow
        swe_in = self._P_snow * dt
        self._h_swe += swe_in

        # decrease due to melt
        swe_out = self._SM * dt
        self._h_swe -= swe_out
        self._h_swe[self._h_swe < 0] = 0

    def _update_snow_depth(self):
        self._h_snow[:] = self._h_swe[:] * self.density_ratio

    def _update_total_prec_snow(self, dt):  # add for landlab (new method)
        self._total_P_snow += self._P_snow * dt

    def _update_total_snowmelt(self, dt):  # add for landlab (new method)
        self._total_SM += self._SM[:] * dt

    def _update_SM_integral(self):
        self._vol_SM = np.sum(self._total_SM * self._grid_area)

    def _update_total_snowpack_water_volume(self):
        volume = self._h_swe * self._grid_area
        self._vol_swe = np.sum(volume)

    def run_one_step(self, dt):
        # update input fields in case there is new input
        self._T_air = self._grid.at_node["atmosphere_bottom_air__temperature"]
        self._P = self._grid.at_node["atmosphere_water__precipitation_leq-volume_flux"]
        self._c0 = self._grid.at_node["snowpack__degree-day_coefficient"]
        self._T0 = self._grid.at_node["snowpack__degree-day_threshold_temperature"]

        # update state variables
        self._update_prec_snow()
        self._update_snowmelt_rate(dt)
        self._update_swe(dt)
        self._update_snow_depth()

        self._update_total_prec_snow(dt)  # add for landlab (new method)
        self._update_total_snowmelt(dt)  # add for landlab (new method)

        self._update_SM_integral()  # after update_total_snowmelt()
        self._update_total_snowpack_water_volume()  # after update_swe()
