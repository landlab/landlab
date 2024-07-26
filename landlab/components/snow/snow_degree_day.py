"""Landlab component that simulates snowpack dynamics.

This component simulates snowmelt process using the degree-day method.
The code is implemented based on the TopoFlow snow component (by Scott
D. Packham) with some adjustments.

* https://github.com/peckhams/topoflow36/blob/master/topoflow/components/snow_degree_day.py
* https://github.com/peckhams/topoflow36/blob/master/topoflow/components/snow_base.py

@author: Tian Gan June 2024
"""

import numpy as np

from landlab import Component

_SECONDS_PER_DAY = 86400.0
_MM_PER_M = 1000.0


class SnowDegreeDay(Component):
    """Simulate snowmelt using the degree-day method.

    The degree-day method uses the degree-day coefficient, `c0`,
    the degree-day threshold temperature, `threshold_temp`, and
    air temperature, to determine the snow melt rate. When the
    air temperature is larger than the threshold temperature,
    the melt process will start. Melt rate calculated as::

        melt_rate = (air_temp - threshold_temp) * c0.

    Units of `c0` are **L / θ / T** (typically, **mm / deg_C / day**),
    so melt rate has units of **L / T**.

    Parameters
    ----------
    grid : ModelGrid
        A Landlab model grid object
    c0: float, optional
        Degree-day coefficient [mm / day / deg_C].
    threshold_temp : float, optional
        Degree-day threshold temperature [deg_C].
    rho_water : float, optinal
        Water density [kg / m3].
    rain_snow_temp : float, optional
        Temperature threshold for precipitation accurs as rain and snow with
        equal frequency [deg_C].

    Notes
    -----

    The *snowpack__melt_volume_flux* field represents the melt rate as
    calculated by the degree-day equation and is not limited by the
    amount of snow that is available to melt. That is, this melt rate
    may be non-zero even in locations that don't have any snow.
    The ``total_snow_melt_at_node`` attribute keeps track of the actual
    amount of snow that was melted at each node.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.snow import SnowDegreeDay

    >>> grid = RasterModelGrid((2, 2))

    >>> grid.add_zeros("atmosphere_water__precipitation_leq-volume_flux", at="node")
    array([0., 0., 0., 0.])
    >>> grid.at_node["atmosphere_bottom_air__temperature"] = [2.0, 2.0, 1.0, -1.0]
    >>> grid.at_node["snowpack__liquid-equivalent_depth"] = [0.005, 0.01, 0.005, 0.005]

    >>> sm = SnowDegreeDay(grid, c0=3.0, threshold_temp=0.0)
    >>> sm.run_one_step(8.64e4)  # 1 day in sec

    >>> grid.at_node["snowpack__melt_volume_flux"] * 86400.0
    array([0.006, 0.006, 0.003, 0.   ])
    >>> grid.at_node["snowpack__liquid-equivalent_depth"]
    array([0.   , 0.004, 0.002, 0.005])

    >>> sm.total_snow_precip_at_node
    array([0., 0., 0., 0.])
    >>> sm.total_snow_melt_at_node
    array([0.005, 0.006, 0.003, 0.   ])
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

    def __init__(
        self, grid, c0=2, threshold_temp=0.0, rho_water=1000.0, rain_snow_temp=1.0
    ):
        """Initialize SnowDegreeDay component."""

        super().__init__(grid)

        self.c0 = c0
        self.t0 = threshold_temp
        self.rho_water = rho_water
        self.t_rain_snow = rain_snow_temp

        self._grid_area = grid.dx * grid.dy

        if "snowpack__z_mean_of_mass-per-volume_density" not in grid.at_node:
            grid.add_full(
                "snowpack__z_mean_of_mass-per-volume_density", 300.0, at="node"
            )

        if "snowpack__liquid-equivalent_depth" not in grid.at_node:
            grid.add_zeros("snowpack__liquid-equivalent_depth", at="node")

        self._total_p_snow = grid.zeros(at="node")
        self._total_sm = grid.zeros(at="node")

        super().initialize_output_fields()

        SnowDegreeDay.calc_snow_depth(
            grid.at_node["snowpack__liquid-equivalent_depth"],
            self.density_ratio,
            out=grid.at_node["snowpack__depth"],
        )

    @property
    def c0(self):
        return self._c0

    @c0.setter
    def c0(self, c0):
        if np.any(c0 < 0):
            raise ValueError("degree-day coefficent must be positive")
        self._c0 = c0 / (_SECONDS_PER_DAY * _MM_PER_M)

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
        if rho_water <= 0.0:
            raise ValueError("water density must be positive")
        self._rho_water = rho_water

    @property
    def t_rain_snow(self):
        return self._t_rain_snow

    @t_rain_snow.setter
    def t_rain_snow(self, t_rain_snow):
        self._t_rain_snow = t_rain_snow

    @property
    def total_snow_precip_at_node(self):
        """Accumulated precipitation as snow over the model time.

        Snow precipitation is given in water equivalent [m].
        """
        return self._total_p_snow

    @property
    def total_snow_melt_at_node(self):
        """Accumulated snow melt over the model time.

        Snow melt is given in water equivalent [m].
        """
        return self._total_sm

    @property
    def density_ratio(self):
        """Ratio of the densities between water and snow."""
        return (
            self._rho_water
            / self.grid.at_node["snowpack__z_mean_of_mass-per-volume_density"]
        )

    @staticmethod
    def calc_precip_snow(precip_rate, air_temp, rain_snow_temp):
        """Calculate snow precipitation.

        Parameters
        ----------
        p : array_like
            Precipitation rate [m/s].
        air_temp : array_like
            Air temperature [deg_C].
        rain_snow_temp : float
            Temperature threshold below which precipitation falls as
            snow [deg_C].
        """
        return np.asarray(precip_rate) * (np.asarray(air_temp) <= rain_snow_temp)

    @staticmethod
    def calc_snow_melt_rate(air_temp, c0, threshold_temp, out=None):
        """Calculate snow melt rate.

        Parameters
        ----------
        air_temp : array_like
            Air temperature [deg_C].
        c0 : float
            Degree-day coefficient [m / s / deg_C].
        threshold_temp : float
            Threshold temperatue above which melt occures [deg_C].

        Returns
        -------
        array_like
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
        precip_rate : array_like
            Precipitation rate [m/s].
        melt_rate : array_like
            Melt rate [m/s].
        swe : array_like
            Snow water equivalent [m].
        dt : float, optional
            Time step [s].

        Returns
        -------
        array_like
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
        h_swe : array_like
            Snow water equivalent [m].
        density_ratio : array_like
            Density ratio of water to snow [-].

        Returns
        -------
        array_like
            Snow depth [m].
        """
        return np.multiply(h_swe, density_ratio, out=out)

    def run_one_step(self, dt):
        """Run component for a time step.

        Parameters
        ----------
        dt : float
            Duration to run the component for [s].
        """
        initial_swe = self.grid.at_node["snowpack__liquid-equivalent_depth"].copy()

        # update state variables
        precip_rate = SnowDegreeDay.calc_precip_snow(
            self._grid.at_node["atmosphere_water__precipitation_leq-volume_flux"],
            self._grid.at_node["atmosphere_bottom_air__temperature"],
            self._t_rain_snow,
        )

        SnowDegreeDay.calc_snow_melt_rate(
            self._grid.at_node["atmosphere_bottom_air__temperature"],
            self._c0,
            self._t0,
            out=self.grid.at_node["snowpack__melt_volume_flux"],
        )

        SnowDegreeDay.calc_swe(
            precip_rate,
            self.grid.at_node["snowpack__melt_volume_flux"],
            self.grid.at_node["snowpack__liquid-equivalent_depth"],
            dt=dt,
            out=self.grid.at_node["snowpack__liquid-equivalent_depth"],
        )
        SnowDegreeDay.calc_snow_depth(
            self.grid.at_node["snowpack__liquid-equivalent_depth"],
            self.density_ratio,
            out=self.grid.at_node["snowpack__depth"],
        )

        self._total_p_snow += precip_rate * dt
        self._total_sm += np.minimum(
            self.grid.at_node["snowpack__melt_volume_flux"] * dt,
            initial_swe + precip_rate * dt,
        )
