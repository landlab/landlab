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
      t0 (degree-day threshold temperature) and t_air (air temperature)
      to determine the snow melt rate.
      When t_air is larger than t0, the melt process will start.
      The melt rate is calculated as sm = (t_air -t0) * c0.
      (e.g., units for c0 and sm could be mm day-1 K-1 and mm day-1 respectively.)

    Parameters
    ----------
    grid : ModelGrid
        A Landlab model grid object
    c0: float, optional
        degree-day coefficient (units: mm day-1 K-1)
    t0: float, optional
        degree-day threshold temperature (units: deg_C)
    rho_water : float, optinal
        water density (units: kg/m3)
    t_rain_snow : float, optional
        temperature threshold for precipitation accurs as rain and snow with
        equal frequency (units: deg_C)

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

    SECOND_PER_DAY = 86400000.0

    def __init__(self, grid, c0=2, t0=0, rho_water=1000, t_rain_snow=1):
        """Initialize SnowDegreeDay component"""

        super().__init__(grid)

        # parameters
        self.c0 = c0
        self._t0 = t0
        self._rho_water = rho_water
        self._t_rain_snow = t_rain_snow

        # calculated vars
        self._grid_area = grid.dx * grid.dy
        self._vol_sm = 0  # snowpack__domain_time_integral_of_melt_volume_flux
        self._vol_swe = 0  # snowpack__domain_integral_of_liquid-equivalent_depth
        self._p_snow = np.zeros(grid.number_of_nodes)
        # change for landlab
        self._total_p_snow = np.zeros(grid.number_of_nodes)  # add new variable
        self._total_sm = np.zeros(grid.number_of_nodes)  # add new variable

        # input data fields
        self._t_air = grid.at_node["atmosphere_bottom_air__temperature"]
        self._p = grid.at_node["atmosphere_water__precipitation_leq-volume_flux"]

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
        self._sm = grid.at_node["snowpack__melt_volume_flux"]

    @property
    def c0(self):
        return self._c0

    @c0.setter
    def c0(self, c0):
        self._c0 = c0

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
        """ precipitation as snow (in water equivalent) at current time step
        (units: m). """
        return self.p_snow

    @property
    def total_p_snow(self):
        """ accumulative precipitation as snow (in water equivalent) over the
        model time (units: m)."""
        return self._total_p_snow

    @property
    def total_sm(self):
        """ accumulative snow melt (in water equivalent) over the
        model time (units: m)."""
        return self._total_sm

    @property
    def vol_sm(self):
        """Total volume of snow melt (in water equivalent) for the domain and
        over the model time (units: m^3). """
        return self._vol_sm

    @property
    def vol_swe(self):
        """Total volume of snow water equivalent for the domain and
        over the model time (units: m^3)."""
        return self._vol_swe

    @property
    def density_ratio(self):
        """ Ratio of the densities between water and snow. """
        return self._rho_water / self._rho_snow

    def _update_snowmelt_rate(self, dt):
        # melt rate based on degree-day coefficient
        sm_c0 = (self._c0 / self.SECOND_PER_DAY) * (
            self._t_air - self._t0  # convert mm -day -K to  m/s
        )
        sm_c0 = np.maximum(sm_c0, 0)

        # melt rate based on available swe (enforced max meltrate)
        sm_max = self._h_swe / dt  # original code use self._h_snow (line 526)

        # actual meltrate
        sm_act = np.minimum(sm_c0, sm_max)
        self._sm[:] = sm_act[:]

    def _update_prec_snow(self):
        self._p_snow = self._p * (self._t_air <= self._t_rain_snow)

    def _update_swe(self, dt):
        # increase due to precipitation as snow
        swe_in = self._p_snow * dt
        self._h_swe += swe_in

        # decrease due to melt
        swe_out = self._sm * dt
        self._h_swe -= swe_out
        self._h_swe[self._h_swe < 0] = 0

    def _update_snow_depth(self):
        self._h_snow[:] = self._h_swe[:] * self.density_ratio

    def _update_total_prec_snow(self, dt):  # add for landlab (new method)
        self._total_p_snow += self._p_snow * dt

    def _update_total_snowmelt(self, dt):  # add for landlab (new method)
        self._total_sm += self._sm[:] * dt

    def _update_sm_integral(self):
        self._vol_sm = np.sum(self._total_sm * self._grid_area)

    def _update_total_snowpack_water_volume(self):
        volume = self._h_swe * self._grid_area
        self._vol_swe = np.sum(volume)

    def run_one_step(self, dt):
        # update input fields in case there is new input
        self._t_air = self._grid.at_node["atmosphere_bottom_air__temperature"]
        self._p = self._grid.at_node["atmosphere_water__precipitation_leq-volume_flux"]

        # update state variables
        self._update_prec_snow()
        self._update_snowmelt_rate(dt)
        self._update_swe(dt)
        self._update_snow_depth()

        self._update_total_prec_snow(dt)  # add for landlab (new method)
        self._update_total_snowmelt(dt)  # add for landlab (new method)

        self._update_sm_integral()  # after update_total_snowmelt()
        self._update_total_snowpack_water_volume()  # after update_swe()
