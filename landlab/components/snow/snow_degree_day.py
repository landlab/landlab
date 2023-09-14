"""Landlab component that simulates snowpack dynamics.

This component simulates snowmelt process using the degree-day method.
It is implemented based on the TopoFlow snow component by Scott D. Packham.
https://github.com/peckhams/topoflow36/blob/master/topoflow/components/snow_degree_day.py
https://github.com/peckhams/topoflow36/blob/master/topoflow/components/snow_base.py

@author: Tian Gan Sept, 2023

Check:
- the way to calculate the max meltrate using h_swe value
- the input and output format (as grid or para, two time series var from topoflow are not included here)
- need to store current time ?

"""
import numpy as np
from landlab import Component


class SnowDegreeDay(Component):

    """Simulate snowmelt process using degree-day method.

      The degree-day method uese c0 (degree-day coefficient), T0 (degree-day threshold temperature)
      and T_air (air temperature) to determine the meltrate. When T_air is larger than T0, the melt process will start.
      The melt rate is calculated as SM = (T_air -T0) * c0

    Parameters
    ----------
    grid : ModelGrid
        A Landlab model grid object
    rho_snow : float (default 200 kg/m3)
        snow density
    rho_H2O : float (default 1000 kg/m3)
        water density
    T_rain_snow : float (default 1 deg_C)
        temperature threshold for precipitation accurs as rain and snow with equal frequency

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.snow import SnowDegreeDay

    setup model grid and input
    >>> grid = RasterModelGrid((2, 2))
    >>> grid.add_full('atmosphere_bottom_air__temperature', 2, at='node')
    array([ 2.,  2.,  2.,  2.])
    >>> grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0, at='node')
    array([ 0.,  0.,  0.,  0.])
    >>> grid.add_full('snowpack__liquid-equivalent_depth', 0.005)  # m
    array([ 0.005,  0.005,  0.005,  0.005])
    >>> grid.add_full('snowpack__degree-day_coefficient', 3)  # mm -day -K
    array([ 3.,  3.,  3.,  3.])

    run one step
    >>> sm = SnowDegreeDay(grid)
    >>> sm.run_one_step(np.float64(8.64E4)) # 1 day in sec

    get model output
    >>> grid.at_node['snowpack__melt_volume_flux']
    array([  5.78703704e-08,   5.78703704e-08,   5.78703704e-08,
         5.78703704e-08])
    >>> grid.at_node['snowpack__time_integral_melt_volume_flux']
    array([ 0.005,  0.005,  0.005,  0.005])
    >>> grid.at_node['snowpack__liquid-equivalent_depth']
    array([ 0.,  0.,  0.,  0.])

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
        },
        "atmosphere_water__precipitation_leq-volume_flux": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/s",
            "mapping": "node",
            "doc": "precipitation rate (in water equivalent)",
        },

        "snowpack__degree-day_coefficient": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "mm day-1 K-1",
            "mapping": "node",
            "doc": "degree-day coefficient",
        },

        "snowpack__degree-day_threshold_temperature": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "deg_C",
            "mapping": "node",
            "doc": "degree-day threshold temperature to start melt",
        },

        "snowpack__liquid-equivalent_depth": {
            "dtype": float,
            "intent": "inout",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "snow water equivalent depth",
        },

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
            "doc": "total precipitation as snow (in water equivalent) during model time ",
        },  # variable defined by me

        "snowpack__time_integral_melt_volume_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "total snow melt (in water equivalent) during model time",
        },  # variable defined by me

        # "snowpack__domain_integral_of_liquid-equivalent_depth": {
        #     "dtype": float,
        #     "intent": "out",
        #     "optional": True,
        #     "units": "m3",
        #     "mapping": "node",
        #     "doc": "??",
        # }, # TODO this is time series output, how to represent in landlab component?
    }

    def __init__(
        self,
        grid,
        rho_snow=200,
        rho_H2O=1000,
        T_rain_snow=1
    ):

        """ Initialize SnowDegreeDay component """

        super().__init__(grid)

        # parameters
        self._rho_snow = rho_snow
        self._rho_H2O = rho_H2O
        self._T_rain_snow = T_rain_snow
        self._density_ratio = self._rho_H2O/self._rho_snow  # density_ratio > 1 (h_snow/h_swe)
        self._P_snow = None
        # self._domain_vol_SM = None  # snowpack__domain_time_integral_of_melt_volume_flux
        # self._domain_vol_swe = None  # snowpack__domain_integral_of_liquid-equivalent_depth

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
        else:
            self._h_swe = grid.add_zeros("snowpack__liquid-equivalent_depth", at='node')

        # output data fields
        super().initialize_output_fields()
        self._h_snow = grid.at_node["snowpack__depth"]
        self._h_snow[:] = self._h_swe[:] * self.density_ratio
        self._SM = grid.at_node["snowpack__melt_volume_flux"]
        self._total_P_snow = grid.at_node["atmosphere_water__time_integral_snowfall_leq-volume_flux"]  # defined by me
        self._total_SM = grid.at_node["snowpack__time_integral_melt_volume_flux"]  # defined by me

    @property
    def rho_snow(self):
        return self._rho_snow

    @rho_snow.setter
    def rho_snow(self,rho_snow):
        assert rho_snow > 0, 'assign rho_snow with positive value'
        self._rho_snow = rho_snow

    @property
    def rho_H2O(self):
        return self._rho_H2O

    @rho_H2O.setter
    def rho_H2O(self, rho_H2O):
        assert rho_H2O > 0, 'assign rho_H2O with positive value'
        self._rho_H2O = rho_H2O

    @property
    def T_rain_snow(self):
        return self._T_rain_snow

    @property
    def density_ratio(self):
        return self._rho_H2O/self._rho_snow

    def _update_snowmelt_rate(self, dt):
        # melt rate based on degree-day coefficient
        sm_c0 = (self._c0 / np.float64(8.64E7)) * (self._T_air - self._T0)  # convert mm -day -K to  m/s
        sm_c0 = np.maximum(sm_c0, np.float64(0))

        # melt rate based on available swe (enforced max meltrate)
        sm_max = self._h_swe / dt  # TODO original code use self._h_snow (line 526)

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

    def _update_total_prec_snow(self, dt):
        self._total_P_snow += self._P_snow * dt

    def _update_total_snowmelt(self, dt):
        self._total_SM += self._SM[:] * dt

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

        self._update_total_prec_snow(dt)  # added by me
        self._update_total_snowmelt(dt)  # added by me
