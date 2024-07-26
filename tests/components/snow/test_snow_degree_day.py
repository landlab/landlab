"""
Unit tests for landlab.components.snow.snow_degree_day

@author Tian Gan  Sept 2023
"""

import numpy as np
import pytest
from numpy.testing import assert_almost_equal

from landlab import RasterModelGrid
from landlab.components.snow import SnowDegreeDay

SECONDS_PER_DAY = 86400.0


def test_create_fields():
    """check to make sure the right fields are created"""

    grid = RasterModelGrid((2, 2))
    grid.add_full("atmosphere_bottom_air__temperature", -2, at="node")
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 1, at="node")
    grid.add_full("snowpack__liquid-equivalent_depth", 1, at="node")
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 500, at="node")

    SnowDegreeDay(grid, rho_water=1000)
    assert_almost_equal(grid.at_node["snowpack__depth"], 2)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 0)


def test_assign_parameters():
    """Test when parameters are not default value"""
    grid = RasterModelGrid((2, 2))
    grid.add_full("atmosphere_bottom_air__temperature", 1, at="node")
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0, at="node")
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 300, at="node")
    sm = SnowDegreeDay(
        grid, c0=1.4, threshold_temp=2, rho_water=1001, rain_snow_temp=1.1
    )

    assert sm.rho_water == 1001
    assert sm.rain_snow_temp == 1.1

    with pytest.raises(ValueError):
        sm.c0 = -1
    with pytest.raises(ValueError):
        sm.rho_water = -1


def test_snow_accumulation():
    """Test when there is only snow accumulation"""

    grid = RasterModelGrid((2, 2), xy_spacing=(10, 10))
    grid.add_full("atmosphere_bottom_air__temperature", -2, at="node")
    grid.add_full(
        "atmosphere_water__precipitation_leq-volume_flux",
        0.001 / SECONDS_PER_DAY,
        at="node",
    )
    grid.add_full("snowpack__liquid-equivalent_depth", 0.001)
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 200, at="node")

    sm = SnowDegreeDay(grid)
    assert np.sum(sm.total_snow_melt_at_node * grid.dx * grid.dy) == pytest.approx(0.0)
    assert np.sum(
        grid.at_node["snowpack__liquid-equivalent_depth"] * grid.dx * grid.dy
    ) == pytest.approx(0.4)

    sm.run_one_step(SECONDS_PER_DAY)  # 1 day in sec
    assert_almost_equal(grid.at_node["snowpack__depth"], 0.01)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 0)
    assert_almost_equal(sm.total_snow_precip_at_node, 0.001)
    assert_almost_equal(sm.total_snow_melt_at_node, 0)

    assert np.sum(
        grid.at_node["snowpack__liquid-equivalent_depth"] * grid.dx * grid.dy
    ) == pytest.approx(0.8)
    assert np.sum(sm.total_snow_melt_at_node * grid.dx * grid.dy) == pytest.approx(0.0)


def test_snow_melt():
    """Test when there is only snow melt"""

    grid = RasterModelGrid((2, 2), xy_spacing=(10, 10))
    grid.add_full("atmosphere_bottom_air__temperature", 2, at="node")
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0, at="node")
    grid.add_full("snowpack__liquid-equivalent_depth", 0.005)  # m
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 200, at="node")

    sm = SnowDegreeDay(grid, c0=3)
    sm.run_one_step(dt=SECONDS_PER_DAY)  # 1 day in sec
    assert_almost_equal(grid.at_node["snowpack__liquid-equivalent_depth"], 0)
    assert_almost_equal(grid.at_node["snowpack__depth"], 0)
    assert_almost_equal(
        grid.at_node["snowpack__melt_volume_flux"] * SECONDS_PER_DAY, 0.006
    )
    assert_almost_equal(sm.total_snow_precip_at_node, 0)
    assert_almost_equal(sm.total_snow_melt_at_node, 0.005)

    assert np.sum(
        grid.at_node["snowpack__liquid-equivalent_depth"] * grid.dx * grid.dy
    ) == pytest.approx(0.0)
    assert np.sum(sm.total_snow_melt_at_node * grid.dx * grid.dy) == pytest.approx(2.0)


def test_snow_melt_accumulation():
    """Test when there is snow melt and accumulation at the same time"""

    grid = RasterModelGrid((2, 2), xy_spacing=(10, 10))
    grid.add_full("atmosphere_bottom_air__temperature", 0.5, at="node")
    grid.add_full(
        "atmosphere_water__precipitation_leq-volume_flux",
        0.002 / SECONDS_PER_DAY,
        at="node",
    )
    grid.add_full("snowpack__liquid-equivalent_depth", 0.005)
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 200, at="node")

    sm = SnowDegreeDay(grid, c0=3)
    sm.run_one_step(SECONDS_PER_DAY)  # 1 day in sec
    assert_almost_equal(grid.at_node["snowpack__liquid-equivalent_depth"], 0.0055)
    assert_almost_equal(grid.at_node["snowpack__depth"], 0.0055 * 5)
    assert_almost_equal(
        grid.at_node["snowpack__melt_volume_flux"] * SECONDS_PER_DAY, 0.0015
    )
    assert_almost_equal(sm.total_snow_precip_at_node, 0.002)
    assert_almost_equal(sm.total_snow_melt_at_node, 0.0015)

    assert np.sum(
        grid.at_node["snowpack__liquid-equivalent_depth"] * grid.dx * grid.dy
    ) == pytest.approx(2.2)
    assert np.sum(sm.total_snow_melt_at_node * grid.dx * grid.dy) == pytest.approx(0.6)


def test_multiple_runs():
    """Test multiple runs with changing inputs at each time step"""

    # step1: snow melt
    grid = RasterModelGrid((2, 2), xy_spacing=(10, 10))
    grid.add_full("atmosphere_bottom_air__temperature", 1, at="node")
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0, at="node")
    grid.add_full("snowpack__liquid-equivalent_depth", 0.005)  # m
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 200, at="node")
    init_swe = grid.at_node["snowpack__liquid-equivalent_depth"].copy()
    sm = SnowDegreeDay(grid)

    sm.run_one_step(SECONDS_PER_DAY)  # 1 day in sec
    assert_almost_equal(grid.at_node["snowpack__liquid-equivalent_depth"], 0.003)
    assert_almost_equal(grid.at_node["snowpack__depth"], 0.015)
    assert_almost_equal(
        grid.at_node["snowpack__melt_volume_flux"] * SECONDS_PER_DAY, 0.002
    )
    assert_almost_equal(sm.total_snow_precip_at_node, 0)
    assert_almost_equal(sm.total_snow_melt_at_node, 0.002)

    # step2: change t_air and p for snow accumulation
    grid.at_node["atmosphere_bottom_air__temperature"].fill(-3)
    grid.at_node["atmosphere_water__precipitation_leq-volume_flux"].fill(
        0.002 / SECONDS_PER_DAY
    )
    grid.at_node["snowpack__z_mean_of_mass-per-volume_density"].fill(500)

    sm.run_one_step(SECONDS_PER_DAY)  # 1 day in sec
    assert_almost_equal(grid.at_node["snowpack__liquid-equivalent_depth"], 0.005)
    assert_almost_equal(grid.at_node["snowpack__depth"], 0.01)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 0)
    assert_almost_equal(sm.total_snow_precip_at_node, 0.002)
    assert_almost_equal(sm.total_snow_melt_at_node, 0.002)

    # step3: change c0 and threshold_temp for snow melt
    grid.at_node["atmosphere_bottom_air__temperature"].fill(1)
    sm.c0 = 1
    sm.threshold_temp = 0.5

    sm.run_one_step(SECONDS_PER_DAY)  # 1 day in sec
    assert_almost_equal(grid.at_node["snowpack__liquid-equivalent_depth"], 0.0065)
    assert_almost_equal(grid.at_node["snowpack__depth"], 0.013)
    assert_almost_equal(
        grid.at_node["snowpack__melt_volume_flux"] * SECONDS_PER_DAY, 0.0005
    )
    assert_almost_equal(sm.total_snow_precip_at_node, 0.004)
    assert_almost_equal(sm.total_snow_melt_at_node, 0.0025)

    assert np.sum(
        grid.at_node["snowpack__liquid-equivalent_depth"] * grid.dx * grid.dy
    ) == pytest.approx(2.6)
    assert np.sum(sm.total_snow_melt_at_node * grid.dx * grid.dy) == pytest.approx(1.0)

    # check mass balance
    in_out = sm.total_snow_precip_at_node - sm.total_snow_melt_at_node

    store = grid.at_node["snowpack__liquid-equivalent_depth"] - init_swe

    assert_almost_equal(in_out, store, err_msg="Error for mass balance", decimal=11)


@pytest.mark.parametrize(
    "precip", (1.0, [2.0], [1.0, 2.0, 3.0], [[1.0, 2.0], [4.0, 5.0]])
)
def test_calc_p_snow(precip):
    assert_almost_equal(SnowDegreeDay.calc_precip_snow(precip, -3.0, 1.0), precip)
    assert_almost_equal(SnowDegreeDay.calc_precip_snow(precip, 3.0, 1.0), 0.0)


@pytest.mark.parametrize("c0", (3.0, [3.0], [3.0, 3.0, 3.0], [[3.0, 3.0], [3.0, 3.0]]))
def test_calc_snow_melt_rate(c0):
    snow_melt_rate = SnowDegreeDay.calc_snow_melt_rate(2.0, c0, 0.0)
    assert np.allclose(snow_melt_rate, 6.0)

    snow_melt_rate = SnowDegreeDay.calc_snow_melt_rate(-2.0, c0, 0.0)
    assert np.allclose(snow_melt_rate, 0.0)


def test_update_swe():
    h_swe = np.full(4, 0.01)
    sm_rate = np.full(4, 0.001)
    p_snow = np.full(4, 0.003)

    SnowDegreeDay.calc_swe(p_snow, sm_rate, h_swe, dt=5.0, out=h_swe)

    assert np.all(h_swe == 0.02)


def test_update_snow_depth():
    h_snow = SnowDegreeDay.calc_snow_depth([0.01, 0.01, 0.01], 2.0)
    assert np.all(h_snow == 0.02)
