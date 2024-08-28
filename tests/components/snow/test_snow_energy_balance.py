"""
Unit tests for landlab.components.snow.snow_energy_balance

@author Tian Gan  Sept 2023
"""

import pytest
from numpy.testing import assert_almost_equal

from landlab import RasterModelGrid
from landlab.components import SnowEnergyBalance


def test_create_fields():
    """check to make sure the right fields are created"""

    grid = RasterModelGrid((2, 2))
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 1, at="node")
    grid.add_full("atmosphere_bottom_air__temperature", 1, at="node")
    grid.add_full("land_surface__temperature", -1, at="node")
    grid.add_full("land_surface_net-total-energy__energy_flux", 1, at="node")
    grid.add_full("snowpack__liquid-equivalent_depth", 1, at="node")

    SnowEnergyBalance(grid)
    assert_almost_equal(
        grid.at_node["snowpack__z_mean_of_mass-per-volume_density"], 300
    )
    assert_almost_equal(
        grid.at_node["snowpack__z_mean_of_mass-specific_isobaric_heat_capacity"], 2090.0
    )
    assert_almost_equal(grid.at_node["snowpack__depth"], 10 / 3)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 0)
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 2090000)
    assert_almost_equal(
        grid.at_node["atmosphere_water__time_integral_snowfall_leq-volume_flux"], 0
    )
    assert_almost_equal(grid.at_node["snowpack__time_integral_melt_volume_flux"], 0)


def test_assign_parameters():
    """Test when parameters are updated"""

    grid = RasterModelGrid((2, 2))
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 1, at="node")
    grid.add_full("atmosphere_bottom_air__temperature", 1, at="node")
    grid.add_full("land_surface__temperature", -1, at="node")
    grid.add_full("land_surface_net-total-energy__energy_flux", 1, at="node")

    sm = SnowEnergyBalance(
        grid,
        rho_H2O=1003,
        rho_air=1.2,
        Cp_air=1004,
        T_rain_snow=1.2,
        T0_cc=0.1,
        grid_area=1000,
    )

    # constants
    assert sm.Lv == 2500000, "wrong Lv value"
    assert sm.Lf == 334000, "wrong Lf value"

    # parameters
    assert sm.rho_H2O == 1003, "wrong rho_snow value"
    assert sm.rho_air == 1.2, "wrong rho_H2O value"
    assert sm.Cp_air == 1004, "wrong Cp_air value"
    assert sm.T_rain_snow == 1.2, "wrong T_rain_snow value"
    assert sm.T0_cc == 0.1, "wrong T0_cc value"
    assert sm.grid_area == 1000, "wrong grid_area value"

    with pytest.raises(AssertionError):
        sm.rho_H2O = -1

    with pytest.raises(AssertionError):
        sm.rho_air = -1

    with pytest.raises(AssertionError):
        sm.Cp_air = -1


def test_snow_accumulation():
    """Test when there is only snow accumulation"""

    grid = RasterModelGrid((2, 2))
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0.001, at="node")
    grid.add_full("atmosphere_bottom_air__temperature", -1, at="node")
    grid.add_full("land_surface__temperature", -1, at="node")
    grid.add_full("land_surface_net-total-energy__energy_flux", 0, at="node")
    grid.add_full("snowpack__liquid-equivalent_depth", 1, at="node")
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 200, at="node")
    grid.add_full(
        "snowpack__z_mean_of_mass-specific_isobaric_heat_capacity", 2000, at="node"
    )

    sm = SnowEnergyBalance(grid, grid_area=100)
    assert sm.vol_SM == 0, f"wrong vol_SM value is {sm.vol_SM}"
    assert sm.vol_swe == 400, f"wrong vol_swe value is {sm.vol_swe}"
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 2e6)

    sm.run_one_step(1000)  # dt = 1000 sec
    assert_almost_equal(grid.at_node["snowpack__depth"], 10)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 0)
    # TODO: may need to change Ecc using new code for update_code_content
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 2e6)
    assert_almost_equal(
        grid.at_node["atmosphere_water__time_integral_snowfall_leq-volume_flux"], 1
    )
    assert_almost_equal(grid.at_node["snowpack__time_integral_melt_volume_flux"], 0)

    assert sm.vol_SM == 0, f"wrong vol_SM value is {sm.vol_SM}"
    assert sm.vol_swe == 800, f"wrong vol_swe value is {sm.vol_swe}"


def test_snow_melt():
    """Test when there is only snow melt"""

    grid = RasterModelGrid((2, 2))
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0, at="node")
    grid.add_full("atmosphere_bottom_air__temperature", 14, at="node")
    grid.add_full("land_surface__temperature", 0.009, at="node")
    grid.add_full(
        "land_surface_net-total-energy__energy_flux", 156, at="node"
    )  # 334000 energy to melt 0.001 m/s
    grid.add_full("snowpack__liquid-equivalent_depth", 0.20, at="node")

    sm = SnowEnergyBalance(grid, grid_area=100)
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 0)

    sm.run_one_step(1000)  # dt = 1000 sec
    assert_almost_equal(grid.at_node["snowpack__depth"], 0.66510978043912183)
    assert_almost_equal(
        grid.at_node["snowpack__melt_volume_flux"], 4.6706586826347305e-07
    )
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 0)
    assert_almost_equal(
        grid.at_node["atmosphere_water__time_integral_snowfall_leq-volume_flux"], 0
    )
    assert_almost_equal(
        grid.at_node["snowpack__time_integral_melt_volume_flux"], 4.6706586826347305e-04
    )

    assert sm.vol_SM == 0.18682634730538922, f"wrong vol_SM value is {sm.vol_SM}"
    assert sm.vol_swe == 79.813173652694616, f"wrong vol_swe value is {sm.vol_swe}"


def test_snow_melt_accumulation():
    """Test when there is both snow melt and accumulation"""

    grid = RasterModelGrid((2, 2))
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0.002, at="node")
    grid.add_full("atmosphere_bottom_air__temperature", 0.5, at="node")
    grid.add_full("land_surface__temperature", -1, at="node")
    grid.add_full(
        "land_surface_net-total-energy__energy_flux", 2e3 + 334, at="node"
    )  # 334000 energy to melt 0.001 m/s
    grid.add_full("snowpack__liquid-equivalent_depth", 1, at="node")
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 200, at="node")
    grid.add_full(
        "snowpack__z_mean_of_mass-specific_isobaric_heat_capacity", 2000, at="node"
    )

    sm = SnowEnergyBalance(grid, grid_area=100)
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 2e6)

    sm.run_one_step(1000)  # dt = 1000 sec
    assert_almost_equal(grid.at_node["snowpack__liquid-equivalent_depth"], 2.999)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 1e-6)
    # TODO may need to update Ecc > 0 swe !=0 using new code for update_code_content
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 0)
    assert_almost_equal(
        grid.at_node["atmosphere_water__time_integral_snowfall_leq-volume_flux"], 2
    )
    assert_almost_equal(grid.at_node["snowpack__time_integral_melt_volume_flux"], 1e-3)

    assert sm.vol_SM == 0.4, f"wrong vol_SM value is {sm.vol_SM}"
    assert sm.vol_swe == 1199.6000000000001, f"wrong vol_swe value is {sm.vol_swe}"


def test_multiple_steps():
    """Test multiple runs with changing inputs at each time step"""

    # step1: snow melt
    grid = RasterModelGrid((2, 2))
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0, at="node")
    grid.add_full("atmosphere_bottom_air__temperature", 1, at="node")
    grid.add_full("land_surface__temperature", -1, at="node")
    grid.add_full(
        "land_surface_net-total-energy__energy_flux", 2e3 + 334, at="node"
    )  # 334000*1000 energy to melt 1 m/s
    grid.add_full("snowpack__liquid-equivalent_depth", 1, at="node")
    init_swe = grid.at_node["snowpack__liquid-equivalent_depth"].copy()
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 200, at="node")
    grid.add_full(
        "snowpack__z_mean_of_mass-specific_isobaric_heat_capacity", 2000, at="node"
    )

    sm = SnowEnergyBalance(grid, grid_area=100)
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 2e6)

    sm.run_one_step(1000)  # dt = 1000 sec
    assert_almost_equal(grid.at_node["snowpack__liquid-equivalent_depth"], 0.999)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 1e-6)
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 0)
    assert_almost_equal(
        grid.at_node["atmosphere_water__time_integral_snowfall_leq-volume_flux"], 0
    )
    assert_almost_equal(grid.at_node["snowpack__time_integral_melt_volume_flux"], 1e-3)

    assert sm.vol_SM == 0.4, f"wrong vol_SM value is {sm.vol_SM}"
    assert sm.vol_swe == 399.6, f"wrong vol_swe value is {sm.vol_swe}"

    # step2: change P, T_air, Q_sum for snow accumulation
    grid.at_node["atmosphere_bottom_air__temperature"].fill(-1)
    grid.at_node["atmosphere_water__precipitation_leq-volume_flux"].fill(0.003)
    grid.at_node["land_surface_net-total-energy__energy_flux"].fill(0)

    sm.run_one_step(1000)
    assert_almost_equal(grid.at_node["snowpack__depth"], 19.995)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 0)
    # TODO: may need to change Ecc because swe != 0 using new code
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 0)
    assert_almost_equal(
        grid.at_node["atmosphere_water__time_integral_snowfall_leq-volume_flux"], 3
    )
    assert_almost_equal(grid.at_node["snowpack__time_integral_melt_volume_flux"], 1e-3)

    assert sm.vol_SM == 0.4, f"wrong vol_SM value is {sm.vol_SM}"
    assert sm.vol_swe == 1599.6000000000001, f"wrong vol_swe value is {sm.vol_swe}"
    #
    # step3: change T_air and Q_sum for snow melt & accumulation
    grid.at_node["atmosphere_bottom_air__temperature"].fill(1)
    grid.at_node["atmosphere_water__precipitation_leq-volume_flux"].fill(0.003)
    # TODO: may need to change Q_sum to account for Ecc>0 for step2 using new code
    grid.at_node["land_surface_net-total-energy__energy_flux"].fill(334)

    sm.run_one_step(1000)
    assert_almost_equal(grid.at_node["snowpack__depth"], 34.99)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 1e-6)
    # TODO: may need to change Ecc>0 because swe != 0 using new code
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 0)
    assert_almost_equal(
        grid.at_node["atmosphere_water__time_integral_snowfall_leq-volume_flux"], 6
    )
    assert_almost_equal(grid.at_node["snowpack__time_integral_melt_volume_flux"], 2e-3)

    assert sm.vol_SM == 0.8, f"wrong vol_SM value is {sm.vol_SM}"
    assert sm.vol_swe == 2799.2000000000003, f"wrong vol_swe value is {sm.vol_swe}"

    # mass balance check
    in_out = (
        grid.at_node["atmosphere_water__time_integral_snowfall_leq-volume_flux"]
        - grid.at_node["snowpack__time_integral_melt_volume_flux"]
    )
    store = grid.at_node["snowpack__liquid-equivalent_depth"] - init_swe
    assert_almost_equal(in_out, store, err_msg="Error for mass balance", decimal=11)
