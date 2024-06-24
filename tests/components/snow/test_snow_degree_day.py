"""
Unit tests for landlab.components.snow.snow_degree_day

@author Tian Gan  Sept 2023
"""

import numpy as np
import pytest
from numpy.testing import assert_almost_equal

from landlab import RasterModelGrid
from landlab.components.snow import SnowDegreeDay


def test_create_fields():
    """check to make sure the right fields are created"""

    grid = RasterModelGrid((2, 2))
    grid.add_full("atmosphere_bottom_air__temperature", -2, at="node")
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 1, at="node")
    grid.add_full("snowpack__liquid-equivalent_depth", 1, at="node")
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 500, at="node")

    SnowDegreeDay(grid, rho_H2O=1000)
    assert_almost_equal(grid.at_node["snowpack__degree-day_coefficient"], 2)
    assert_almost_equal(grid.at_node["snowpack__degree-day_threshold_temperature"], 0)
    assert_almost_equal(grid.at_node["snowpack__depth"], 2)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 0)
    assert_almost_equal(
        grid.at_node["atmosphere_water__time_integral_snowfall_leq-volume_flux"], 0
    )
    assert_almost_equal(grid.at_node["snowpack__time_integral_melt_volume_flux"], 0)


def test_assign_parameters():
    """Test when parameters are updated"""
    grid = RasterModelGrid((2, 2))
    grid.add_full("atmosphere_bottom_air__temperature", 1, at="node")
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0, at="node")
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 300, at="node")
    sm = SnowDegreeDay(grid, rho_H2O=1001, T_rain_snow=1.1, grid_area=100)

    assert sm.rho_H2O == 1001, "wrong rho_H2O value"
    assert sm.T_rain_snow == 1.1, "wrong T_rain_snow value"
    assert sm.grid_area == 100, "wrong grid_area value"

    with pytest.raises(AssertionError):
        sm.rho_H2O = -1


def test_snow_accumulation():
    """Test when there is only snow accumulation"""

    grid = RasterModelGrid((2, 2), xy_spacing=(2, 2))
    grid.add_full("atmosphere_bottom_air__temperature", -2, at="node")
    grid.add_full(
        "atmosphere_water__precipitation_leq-volume_flux",
        0.001 / np.float64(8.64e4),
        at="node",
    )
    grid.add_full("snowpack__liquid-equivalent_depth", 0.001)
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 200, at="node")

    sm = SnowDegreeDay(grid, grid_area=100)
    assert sm.vol_SM == 0, f"wrong vol_SM value is {sm.vol_SM}"
    assert sm.vol_swe == 0.4, f"wrong vol_swe value is {sm.vol_swe}"

    sm.run_one_step(np.float64(8.64e4))  # 1 day in sec
    assert_almost_equal(grid.at_node["snowpack__depth"], 0.01)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 0)
    assert_almost_equal(
        grid.at_node["atmosphere_water__time_integral_snowfall_leq-volume_flux"], 0.001
    )
    assert_almost_equal(grid.at_node["snowpack__time_integral_melt_volume_flux"], 0)
    assert sm.vol_swe == 0.8, f"wrong vol_swe value is {sm.vol_swe}"
    assert sm.vol_SM == 0, f"wrong vol_SM value is  {sm.vol_SM}"

    sm = SnowDegreeDay(grid)
    assert sm.grid_area == 4, f"wrong grid_area value is {sm.grid_area}"


def test_snow_melt():
    """Test when there is only snow melt"""

    grid = RasterModelGrid((2, 2))
    grid.add_full("atmosphere_bottom_air__temperature", 2, at="node")
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0, at="node")
    grid.add_full("snowpack__liquid-equivalent_depth", 0.005)  # m
    grid.add_full("snowpack__degree-day_coefficient", 3)  # mm -day -K
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 200, at="node")

    sm = SnowDegreeDay(grid, grid_area=100)
    sm.run_one_step(np.float64(8.64e4))  # 1 day in sec
    assert_almost_equal(grid.at_node["snowpack__degree-day_threshold_temperature"], 0)
    assert_almost_equal(grid.at_node["snowpack__liquid-equivalent_depth"], 0)
    assert_almost_equal(grid.at_node["snowpack__depth"], 0)
    assert_almost_equal(
        grid.at_node["snowpack__melt_volume_flux"],
        0.005 / np.float64(8.64e4),
        decimal=11,
    )
    assert_almost_equal(
        grid.at_node["atmosphere_water__time_integral_snowfall_leq-volume_flux"], 0
    )
    assert_almost_equal(grid.at_node["snowpack__time_integral_melt_volume_flux"], 0.005)
    assert sm.vol_swe == 0, f" wrong vol_swe value is {sm.vol_swe}"
    assert sm.vol_SM == 2, f"wrong vol_SM value is  {sm.vol_SM}"


def test_snow_melt_accumulation():
    """Test when there is both snow melt and accumulation"""

    grid = RasterModelGrid((2, 2))
    grid.add_full("atmosphere_bottom_air__temperature", 0.5, at="node")
    grid.add_full(
        "atmosphere_water__precipitation_leq-volume_flux",
        0.002 / np.float64(8.64e4),
        at="node",
    )  # m
    grid.add_full("snowpack__liquid-equivalent_depth", 0.005)  # m
    grid.add_full("snowpack__degree-day_coefficient", 3)  # mm -K -Day
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 200, at="node")

    sm = SnowDegreeDay(grid, grid_area=100)
    sm.run_one_step(np.float64(8.64e4))  # 1 day in sec
    assert_almost_equal(grid.at_node["snowpack__degree-day_threshold_temperature"], 0)
    assert_almost_equal(grid.at_node["snowpack__liquid-equivalent_depth"], 5.5 / 1000)
    assert_almost_equal(grid.at_node["snowpack__depth"], 5.5 * 5 / 1000)
    assert_almost_equal(
        grid.at_node["snowpack__melt_volume_flux"],
        (1.5 / 1000) / np.float64(8.64e4),
        decimal=11,
    )
    assert_almost_equal(
        grid.at_node["atmosphere_water__time_integral_snowfall_leq-volume_flux"],
        2 / 1000,
    )
    assert_almost_equal(
        grid.at_node["snowpack__time_integral_melt_volume_flux"], 1.5 / 1000
    )

    assert sm.vol_swe == 2.2, f"wrong vol_swe value is {sm.vol_swe}"
    assert sm.vol_SM == 0.6, f"wrong vol_SM value is {sm.vol_SM}"


def test_multiple_runs():
    """Test multiple runs with changing inputs at each time step"""

    # step1: snow melt
    grid = RasterModelGrid((2, 2))
    grid.add_full("atmosphere_bottom_air__temperature", 1, at="node")
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0, at="node")
    grid.add_full("snowpack__liquid-equivalent_depth", 0.005)  # m
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 200, at="node")
    init_swe = grid.at_node["snowpack__liquid-equivalent_depth"].copy()
    sm = SnowDegreeDay(grid, grid_area=100)

    sm.run_one_step(np.float64(8.64e4))  # 1 day in sec
    assert_almost_equal(grid.at_node["snowpack__degree-day_threshold_temperature"], 0)
    assert_almost_equal(grid.at_node["snowpack__liquid-equivalent_depth"], 0.003)
    assert_almost_equal(grid.at_node["snowpack__depth"], 0.015)
    assert_almost_equal(
        grid.at_node["snowpack__melt_volume_flux"],
        0.002 / np.float64(8.64e4),
        decimal=11,
    )
    assert_almost_equal(
        grid.at_node["atmosphere_water__time_integral_snowfall_leq-volume_flux"], 0
    )
    assert_almost_equal(grid.at_node["snowpack__time_integral_melt_volume_flux"], 0.002)

    # step2: change T_air and P for snow accumulation
    grid.at_node["atmosphere_bottom_air__temperature"].fill(-3)
    grid.at_node["atmosphere_water__precipitation_leq-volume_flux"].fill(
        0.002 / np.float64(8.64e4)
    )
    grid.at_node["snowpack__z_mean_of_mass-per-volume_density"].fill(500)

    sm.run_one_step(np.float64(8.64e4))  # 1 day in sec
    assert_almost_equal(grid.at_node["snowpack__degree-day_threshold_temperature"], 0)
    assert_almost_equal(grid.at_node["snowpack__liquid-equivalent_depth"], 0.005)
    assert_almost_equal(grid.at_node["snowpack__depth"], 0.01)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 0)
    assert_almost_equal(
        grid.at_node["atmosphere_water__time_integral_snowfall_leq-volume_flux"], 0.002
    )
    assert_almost_equal(grid.at_node["snowpack__time_integral_melt_volume_flux"], 0.002)

    # step3: change c0 and T0 for snow melt
    grid.at_node["atmosphere_bottom_air__temperature"].fill(1)
    grid.at_node["snowpack__degree-day_coefficient"].fill(1)
    grid.at_node["snowpack__degree-day_threshold_temperature"].fill(0.5)

    sm.run_one_step(np.float64(8.64e4))  # 1 day in sec
    assert_almost_equal(grid.at_node["snowpack__liquid-equivalent_depth"], 0.0065)
    assert_almost_equal(grid.at_node["snowpack__depth"], 0.013)
    assert_almost_equal(
        grid.at_node["snowpack__melt_volume_flux"],
        (0.5 / 1000) / np.float64(8.64e4),
        decimal=11,
    )
    assert_almost_equal(
        grid.at_node["atmosphere_water__time_integral_snowfall_leq-volume_flux"], 0.004
    )
    assert_almost_equal(
        grid.at_node["snowpack__time_integral_melt_volume_flux"], 0.0025
    )
    assert sm.vol_swe == 2.6, f"wrong vol_swe value is {sm.vol_swe}"
    assert sm.vol_SM == 1, f"wrong vol_SM value is {sm.vol_SM}"

    # check mass balance
    in_out = (
        grid.at_node["atmosphere_water__time_integral_snowfall_leq-volume_flux"]
        - grid.at_node["snowpack__time_integral_melt_volume_flux"]
    )

    store = grid.at_node["snowpack__liquid-equivalent_depth"] - init_swe

    assert_almost_equal(in_out, store, err_msg="Error for mass balance", decimal=11)
