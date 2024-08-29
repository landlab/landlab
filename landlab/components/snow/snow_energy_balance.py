"""Landlab component that simulates snowpack dynamics.

This component simulates snowmelt process using the snow energy balance method.
The code is implemented based on the TopoFlow snow component (by Scott D. Packham)
with some adjustments.

* https://github.com/peckhams/topoflow36/blob/master/topoflow/components/snow_energy_balance.py
* https://github.com/peckhams/topoflow36/blob/master/topoflow/components/snow_base.py

@author: Tian Gan July 2024
"""

import numpy as np

from landlab import Component

_LATENT_HEAT_OF_FUSION = 334000  # latent heat of fusion [J/kg] (to liquid)


class SnowEnergyBalance(Component):
    """Simulate snowmelt process using snow energy balance method.

    Snow energy balance method accounts for various energy fluxes that interact with
    the snow layer to simulate snowpack dynamics.
    This component uses the net total energy flux (q_sum) and the snowpack cold content
    (cold_content) to determine the snow melt rate:
    - If the net total energy (q_sum * time) is larger than cold content, the melt process
    starts using the remaining energy (q_sum * time - cold_content).
    - If the net total energy is positive but less than cold content, the snowpack is
    warming and the cold content decreases.
    - If the total energy is negative, the snow is cooling and the cold content
    increases.

    q_sum is the sum of the various net energy fluxes from net short
    wave radiation, long wave radition, sensible heat, and latent heat, etc.
    q_sum can be estimated using the meterology comopnent in landlab or
    provided as known input data.

    cold_content can be estimated as:
    cold_content = rho_snow * cp_snow * (melting_point - snow_temp) * h_snow
    - rho_snow: snow density
    - cp_snow: snow heat capacity
    - melting_point: snow melting point temperature
    - snow_temp: snowpack average temperature
    - h_snow: snow depth

    * This component uses land surface temperature (surf_temp) as an approximation to
      represent snowpack average temperature (snow_temp).

    Parameters
    ----------
    grid : ModelGrid
        A Landlab model grid object
    rho_water : float, optional
        water density [kg / m3].
    rho_air : float, optional
        air density [kg / m3].
    cp_air  : float, optional
        air heat capacity [J / kg / K]
    cp_snow: float, optional
        snow heat capacity [J / kg / K]
    rain_snow_temp : float, optional
        Temperature threshold for precipitation accurs as rain and snow with
        equal frequency [deg_C].
    melting_point : float, optional
        snow melting-point temperature [deg_C].

    Notes
    -----

    The *snowpack__melt_volume_flux* field represents the melt rate as
    calculated by the snow energy balance equation and is not limited by the
    amount of snow that is available to melt. That is, this melt rate
    may be non-zero even in locations that don't have any snow.
    The ``total_snow_melt_at_node`` attribute keeps track of the actual
    amount of snow that was melted at each node.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.snow import SnowEnergyBalance

    >>> grid = RasterModelGrid((2, 2))

    >>> grid.add_zeros("atmosphere_water__precipitation_leq-volume_flux", at="node")
    array([0., 0., 0., 0.])
    >>> grid.at_node["atmosphere_bottom_air__temperature"] = [2.0, 2.0, 1.0, -1.0]
    >>> grid.at_node["snowpack__liquid-equivalent_depth"] = [0.005, 0.01, 0.005, 0.005]
    >>> grid.at_node["land_surface__temperature"] = [0.0, -1, 1, -1]
    >>> grid.add_full("land_surface_net-total-energy__energy_flux", 334, at="node")
    array([334., 334., 334., 334.])
    >>> grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 200, at="node")
    array([200., 200., 200., 200.])

    >>> sm = SnowEnergyBalance(grid)
    >>> sm.run_one_step(1000)

    >>> grid.at_node["snowpack__melt_volume_flux"]
    array([1.00000000e-06, 9.37425150e-07, 1.00000000e-06, 9.68712575e-07])
    >>> grid.at_node["snowpack__liquid-equivalent_depth"]
    array([0.004     , 0.00906257, 0.004     , 0.00403129])
    >>> grid.at_node["snowpack__energy-per-area_cold_content"]
    array([0., 0., 0., 0.])

    >>> sm.total_snow_precip_at_node
    array([0., 0., 0., 0.])
    >>> sm.total_snow_melt_at_node
    array([0.001     , 0.00093743, 0.001     , 0.00096871])
    """

    _name = "SnowEnergyBalance"

    _unit_agnostic = False

    _info = {
        # input fields
        "atmosphere_water__precipitation_leq-volume_flux": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/s",
            "mapping": "node",
            "doc": "precipitation rate (in water equivalent)",
        },
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
            "intent": "in",
            "optional": False,
            "units": "deg_C",
            "mapping": "node",
            "doc": "land surface temperature",
        },
        "land_surface_net-total-energy__energy_flux": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "W m-2",
            "mapping": "node",
            "doc": "net total energy flux",
        },
        "snowpack__liquid-equivalent_depth": {
            "dtype": float,
            "intent": "inout",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "snow water equivalent depth",
        },
        "snowpack__z_mean_of_mass-per-volume_density": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "kg m-3",
            "mapping": "node",
            "doc": "snow density",
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
        "snowpack__energy-per-area_cold_content": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "J m-2",
            "mapping": "node",
            "doc": "snowpack cold content",
        },
    }

    def __init__(
        self,
        grid,
        rho_water=1000,
        rho_air=1.2614,
        cp_air=1005.7,
        cp_snow=2090.0,
        rain_snow_temp=1,
        melting_point=0,
    ):
        """Initialize SnowEnergyBalance component"""

        super().__init__(grid)

        self.rho_water = rho_water
        self.rho_air = rho_air
        self.cp_air = cp_air
        self.cp_snow = cp_snow
        self.rain_snow_temp = rain_snow_temp
        self.melting_point = melting_point
        self._total_p_snow = grid.zeros(at="node")
        self._total_sm = grid.zeros(at="node")

        if "snowpack__z_mean_of_mass-per-volume_density" not in grid.at_node:
            grid.add_full(
                "snowpack__z_mean_of_mass-per-volume_density", 300.0, at="node"
            )

        if "snowpack__liquid-equivalent_depth" not in grid.at_node:
            grid.add_zeros("snowpack__liquid-equivalent_depth", at="node")

        if "snowpack__melt_volume_flux" not in grid.at_node:
            grid.add_empty("snowpack__melt_volume_flux", at="node")
        grid.at_node["snowpack__melt_volume_flux"].fill(0.0)

        if "snowpack__depth" not in grid.at_node:
            grid.add_empty("snowpack__depth", at="node")

        SnowEnergyBalance.calc_snow_depth(
            grid.at_node["snowpack__liquid-equivalent_depth"],
            self.density_ratio,
            out=grid.at_node["snowpack__depth"],
        )

        if "snowpack__energy-per-area_cold_content" not in grid.at_node:
            grid.add_empty("snowpack__energy-per-area_cold_content", at="node")

        SnowEnergyBalance.initialize_cold_content(
            grid.at_node["snowpack__z_mean_of_mass-per-volume_density"],
            grid.at_node["snowpack__depth"],
            grid.at_node["land_surface__temperature"],
            2090,  # self.cp_snow,
            0,  # self.melting_point,
            out=grid.at_node["snowpack__energy-per-area_cold_content"],
        )

    @property
    def rho_water(self):
        return self._rho_water

    @rho_water.setter
    def rho_water(self, rho_water):
        if rho_water <= 0.0:
            raise ValueError("water density must be positive")
        self._rho_water = rho_water

    @property
    def rho_air(self):
        return self._rho_air

    @rho_air.setter
    def rho_air(self, rho_air):
        if rho_air <= 0:
            raise ValueError("air temperature density must be positive")
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
    def cp_snow(self):
        return self._cp_snow

    @cp_snow.setter
    def cp_snow(self, cp_snow):
        if cp_snow <= 0:
            raise ValueError("snow heat capacity must be positive")
        self._cp_snow = cp_snow

    @property
    def rain_snow_temp(self):
        return self._rain_snow_temp

    @rain_snow_temp.setter
    def rain_snow_temp(self, rain_snow_temp):
        self._rain_snow_temp = rain_snow_temp

    @property
    def melting_point(self):
        return self._melting_point

    @melting_point.setter
    def melting_point(self, melting_point):
        self._melting_point = melting_point

    @property
    def density_ratio(self):
        """Ratio of the densities between water and snow."""
        return (
            self._rho_water
            / self.grid.at_node["snowpack__z_mean_of_mass-per-volume_density"]
        )

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

    @staticmethod
    def initialize_cold_content(
        rho_snow, h_snow, surf_temp, cp_snow, melting_point, out=None
    ):
        """Initialize cold content

        Parameters
        ----------

        rho_snow: array_like
            Snow density [kg / m3].
        h_snow : array_like
            Snow depth [m].
        surf_temp : array_like
            Land surface temperature [deg_C].
        cp_snow: float
            Snow heat capacity [J / kg / K].
        melting_point: float
            Snow melting-point temperature [deg_C].

        Returns
        -------
        array_like
            Snowpack cold content [J / m^2].
        """
        cold_content = (
            np.asarray(rho_snow)
            * cp_snow
            * (melting_point - np.asarray(surf_temp))
            * np.asarray(h_snow)
        )

        return np.clip(cold_content, a_min=0, a_max=None, out=out)

    @staticmethod
    def calc_precip_snow(precip_rate, air_temp, rain_snow_temp):
        """Calculate snow precipitation.

        Parameters
        ----------
        precip_rate : array_like
            Precipitation rate [m/s].
        air_temp : array_like
            Air temperature [deg_C].
        rain_snow_temp : float
            Temperature threshold below which precipitation falls as
            snow [deg_C].
        """
        return np.asarray(precip_rate) * (np.asarray(air_temp) <= rain_snow_temp)

    @staticmethod
    def calc_snow_melt_rate(q_sum, cold_content, rho_water, dt, out=None):
        """Calculate snow melt rate.

        Parameters
        ----------

        q_sum: array_like
            net total energy flux [W / m2].
        cold_content: array_like
            snowpack cold content [J / m2].
        rho_water: float
            water density [kg / m3].
        dt : float
            Duration to run the component for [s].

        Returns
        -------
        array_like
            Melt rate [m / s]
        """
        q_rem = (np.asarray(q_sum) * dt - np.asarray(cold_content)) / dt
        melt_rate = q_rem / (rho_water * _LATENT_HEAT_OF_FUSION)

        return np.clip(melt_rate, a_min=0.0, a_max=None, out=out)

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

    @staticmethod
    def calc_cold_content(q_sum, cold_content, dt, out=None):
        """Calculate snowpack cold content

        Parameters
        ----------
        q_sum: array_like
            net total energy flux [W / m2].
        cold_content: array_like
            snowpack cold content [J / m2].
        dt : float, optional
            Time step [s].

        Returns
        -------
        array_like
            New snowpack cold content [J / m2].
        """
        out = np.subtract(cold_content, q_sum * dt, out)

        return out.clip(min=0, max=None, out=out)

    def run_one_step(self, dt):
        """Run component for a time step.

        Parameters
        ----------
        dt : float
            Duration to run the component for [s].
        """
        initial_swe = self.grid.at_node["snowpack__liquid-equivalent_depth"].copy()

        # update state variables
        precip_rate = SnowEnergyBalance.calc_precip_snow(
            self._grid.at_node["atmosphere_water__precipitation_leq-volume_flux"],
            self._grid.at_node["atmosphere_bottom_air__temperature"],
            self.rain_snow_temp,
        )

        SnowEnergyBalance.calc_snow_melt_rate(
            self._grid.at_node["land_surface_net-total-energy__energy_flux"],
            self._grid.at_node["snowpack__energy-per-area_cold_content"],
            self.rho_water,
            dt=dt,
            out=self.grid.at_node["snowpack__melt_volume_flux"],
        )

        SnowEnergyBalance.calc_swe(
            precip_rate,
            self.grid.at_node["snowpack__melt_volume_flux"],
            self.grid.at_node["snowpack__liquid-equivalent_depth"],
            dt=dt,
            out=self.grid.at_node["snowpack__liquid-equivalent_depth"],
        )

        SnowEnergyBalance.calc_snow_depth(
            self.grid.at_node["snowpack__liquid-equivalent_depth"],
            self.density_ratio,
            out=self.grid.at_node["snowpack__depth"],
        )

        SnowEnergyBalance.calc_cold_content(
            self._grid.at_node["land_surface_net-total-energy__energy_flux"],
            self._grid.at_node["snowpack__energy-per-area_cold_content"],
            dt=dt,
            out=self._grid.at_node["snowpack__energy-per-area_cold_content"],
        )

        self._total_p_snow += precip_rate * dt
        self._total_sm += np.minimum(
            self.grid.at_node["snowpack__melt_volume_flux"] * dt,
            initial_swe + precip_rate * dt,
        )
