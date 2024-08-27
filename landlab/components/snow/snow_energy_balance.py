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


class SnowEnergyBalance(Component):

    """Simulate snowmelt process using snow energy balance method.

    This component accounts for energy fluxes to simulate snowpack dynamics
    (e.g., solar radiation, long wave radiation, sensible heat, latent heat).
    It uses the net total energy flux (Q_sum) and the snowpack cold content (Ecc)
    to determine the snow melt rate (SM):
    - If the net total energy (Q_sum * dt) is larger than Ecc, the melt process starts
    using the remaining energy (Q_sum*dt - Ecc).
    - If the net total energy is positive but less than Ecc, the snowpack is warming
    and the Ecc decrease.
    - If the total energy is negative, the snow is cooling and the Ecc increase.

    Q_sum = Qn_SW + Qn_LW + Qh + Qe + Qa + Qc
    - Qn_SW: net short wave energy flux
    - Qn_LW: net long wave energy flux
    - Qh: sensible heat energy flux
    - Qe: latent heat energy flux
    - Qa: net energy flux advected by moving water
    - Qc: net energy flux via conduction from snow to soil

    Ecc = rho_snow * Cp_snow * (T_air - T_surf) * h_snow
    - rho_snow: snow density
    - Cp_snow: snow heat capacity
    - T_air: air temperature
    - T_surf: land surface temperature (~ snow surface temperature)
    - h_snow: snow depth

    Parameters
    ----------
    grid : ModelGrid
        A Landlab model grid object
    rho_H2O : float (default 1000 kg/m3)
        water density
    rho_air : float (default 1.2614 kg/m3)
        air density
    Cp_air  : float (default 1005.7 J kg-1 K-1)
        air heat capacity
    T_rain_snow : float (default 1 deg_C)
        temperature threshold for precipitation accurs as rain and snow with
        equal frequency
    T0_cc : float (default 0 deg_C)
        melting-point temperature
    grid_area: float (default 0 m2)
        space area represented by a single grid cell.
        if value = 0, grid_area=grid.dy * grid.dx

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.snow import SnowEnergyBalance
    >>> grid = RasterModelGrid((2, 2))
    >>> grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0, at="node")
    array([ 0.,  0.,  0.,  0.])
    >>> grid.add_full("atmosphere_bottom_air__temperature", 1, at="node")
    array([ 1.,  1.,  1.,  1.])
    >>> grid.add_full("land_surface__temperature", -1, at="node")
    array([-1., -1., -1., -1.])
    >>> grid.add_full(
    ...     "land_surface_net-total-energy__energy_flux", 2e3 + 334, at="node"
    ... )
    array([ 2334.,  2334.,  2334.,  2334.])
    >>> grid.add_full("snowpack__liquid-equivalent_depth", 1, at="node")
    array([ 1.,  1.,  1.,  1.])
    >>> grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 200, at="node")
    array([ 200.,  200.,  200.,  200.])
    >>> grid.add_full(
    ...     "snowpack__z_mean_of_mass-specific_isobaric_heat_capacity", 2000, at="node"
    ... )
    array([ 2000.,  2000.,  2000.,  2000.])
    >>> sm = SnowEnergyBalance(grid, grid_area=100)
    >>> grid.at_node["snowpack__depth"]
    array([ 5.,  5.,  5.,  5.])
    >>> grid.at_node["snowpack__energy-per-area_cold_content"]
    array([ 2000000.,  2000000.,  2000000.,  2000000.])
    >>> sm.run_one_step(1000)
    >>> grid.at_node["snowpack__melt_volume_flux"]
    array([  1.00000000e-06,   1.00000000e-06,   1.00000000e-06,   1.00000000e-06])
    >>> grid.at_node["snowpack__liquid-equivalent_depth"]
    array([ 0.999,  0.999,  0.999,  0.999])
    >>> grid.at_node["snowpack__energy-per-area_cold_content"]
    array([ 0.,  0.,  0.,  0.])
    """

    _name = "SnowEnergyBalance"

    _unit_agnostic = False

    _info = {
        # input fields (4 req + 3 opt)
        "atmosphere_water__precipitation_leq-volume_flux": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/s",
            "mapping": "node",
            "doc": "precipitation rate (in water equivalent)",
        },  # P
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
        },  # T_surf (this model uses snow surface temp/land surface temp
            # to represent snowpack average temp)
        "land_surface_net-total-energy__energy_flux": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "W m-2",
            "mapping": "node",
            "doc": "net total energy flux",
        },  # Q_sum
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
        "snowpack__z_mean_of_mass-specific_isobaric_heat_capacity": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "J kg-1 K-1",
            "mapping": "node",
            "doc": "snow heat capacity",
        },  # Cp_snow
        # output fields (5 var)
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
        },  # SM
        "snowpack__energy-per-area_cold_content": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "J m-2",
            "mapping": "node",
            "doc": "snowpack cold content",
        },  # Ecc
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

    def __init__(
        self,
        grid,
        rho_H2O=1000,
        rho_air=1.2614,
        Cp_air=1005.7,
        T_rain_snow=1,
        T0_cc=0,
        grid_area=0,
    ):
        """Initialize SnowEnergyBalance component"""

        super().__init__(grid)

        # parameters
        self._rho_H2O = rho_H2O  # kg m-3
        self._rho_air = rho_air  # kg m-3
        self._Cp_air = Cp_air  # J kg-1 K-1
        self._T_rain_snow = T_rain_snow  # deg_C
        self._T0_cc = T0_cc  # deg_C
        self._grid_area = grid_area if grid_area > 0 else self._grid.dy * self._grid.dx

        # constant
        self.Lv = np.float64(
            2500000
        )  # latent heat of vaporization [J/kg] (to gas) TODO may not be used
        self.Lf = np.float64(334000)  # latent heat of fusion [J/kg] (to liquid)

        # calculated vars
        self._P_snow = None  # precipitation as snow (m/s)
        self._vol_SM = 0  # snowpack__domain_time_integral_of_melt_volume_flux (m3)
        self._vol_swe = 0  # snowpack__domain_integral_of_liquid-equivalent_depth (m3)

        # input fields
        self._P = grid.at_node["atmosphere_water__precipitation_leq-volume_flux"]
        self._T_air = grid.at_node["atmosphere_bottom_air__temperature"]
        self._T_surf = grid.at_node["land_surface__temperature"]
        self._Q_sum = grid.at_node["land_surface_net-total-energy__energy_flux"]

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

        if "snowpack__z_mean_of_mass-specific_isobaric_heat_capacity" in grid.at_node:
            self._Cp_snow = grid.at_node[
                "snowpack__z_mean_of_mass-specific_isobaric_heat_capacity"
            ]
        else:
            self._Cp_snow = grid.add_full(
                "snowpack__z_mean_of_mass-specific_isobaric_heat_capacity", 2090.0
            )

        # output fields
        super().initialize_output_fields()
        self._h_snow = grid.at_node["snowpack__depth"]
        self._h_snow[:] = self._h_swe[:] * self.density_ratio
        self._SM = grid.at_node["snowpack__melt_volume_flux"]
        self._Ecc = grid.at_node["snowpack__energy-per-area_cold_content"]
        self._initialize_cold_content()
        self._total_P_snow = grid.at_node[
            "atmosphere_water__time_integral_snowfall_leq-volume_flux"
        ]  # defined by me
        self._total_SM = grid.at_node[
            "snowpack__time_integral_melt_volume_flux"
        ]  # defined by me

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

    @property
    def T_rain_snow(self):
        return self._T_rain_snow

    @T_rain_snow.setter
    def T_rain_snow(self, T_rain_snow):
        self._T_rain_snow = T_rain_snow

    @property
    def T0_cc(self):
        return self._T0_cc

    @T0_cc.setter
    def T0_cc(self, T0_cc):
        self._T0_cc = T0_cc

    @property
    def grid_area(self):
        return self._grid_area

    @property
    def vol_SM(self):
        return self._vol_SM

    @property
    def vol_swe(self):
        return self._vol_swe

    @property
    def density_ratio(self):
        return self._rho_H2O / self._rho_snow

    def _initialize_cold_content(self):
        del_T = self._T0_cc - self._T_surf
        self._Ecc[:] = self._rho_snow * self._Cp_snow * del_T * self._h_snow
        np.maximum(self._Ecc, np.float64(0), out=self._Ecc)

    def _update_cold_content(self, dt):
        # This model estimates Ecc mainly based on the net energy input
        # if Ecc is positive, there is no melt
        E_in = self._Q_sum * dt
        Ecc = np.maximum(self._Ecc - E_in, np.float64(0))
        self._Ecc[:] = Ecc[:]

    def _update_snowmelt_rate(self, dt):
        # melt rate based on available energy
        E_in = self._Q_sum * dt
        E_rem = np.maximum(E_in - self._Ecc, np.float64(0))
        Qm = E_rem / dt  # energy flux for melting W m-2
        sm_en = Qm / (self._rho_H2O * self.Lf)

        # melt rate based on available swe (enforced max meltrate)
        sm_max = self._h_swe / dt

        # actual melt rate
        sm_act = np.minimum(sm_en, sm_max)
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
        self._h_snow[:] = self._h_swe * self.density_ratio

    def _update_total_prec_snow(self, dt):  # defined by me
        self._total_P_snow += self._P_snow * dt

    def _update_total_snowmelt(self, dt):  # defined by me
        self._total_SM += self._SM[:] * dt

    def _update_SM_integral(self):
        self._vol_SM = np.sum(self._total_SM * self._grid_area)

    def _update_total_snowpack_water_volume(self):
        volume = self._h_swe * self._grid_area
        self._vol_swe = np.sum(volume)

    def run_one_step(self, dt):
        # update input fields in case there is new input
        self._P = self._grid.at_node["atmosphere_water__precipitation_leq-volume_flux"]
        self._T_air = self._grid.at_node["atmosphere_bottom_air__temperature"]
        self._Q_sum = self._grid.at_node["land_surface_net-total-energy__energy_flux"]
        self._rho_snow = self._grid.at_node[
            "snowpack__z_mean_of_mass-per-volume_density"
        ]
        self._Cp_snow = self._grid.at_node[
            "snowpack__z_mean_of_mass-specific_isobaric_heat_capacity"
        ]

        # update state variables
        self._update_prec_snow()
        self._update_snowmelt_rate(dt)
        self._update_swe(dt)
        self._update_snow_depth()  # after update_swe()
        self._update_cold_content(dt)  # after update_snow_depth()

        self._update_total_prec_snow(dt)  # change for landlab (new method)
        self._update_total_snowmelt(dt)  # change for landlab (new method)

        self._update_SM_integral()  # after update_total_snowmelt()
        self._update_total_snowpack_water_volume()  # after update_swe()
