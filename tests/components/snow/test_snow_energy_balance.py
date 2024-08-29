"""
Unit tests for landlab.components.snow.snow_energy_balance

@author Tian Gan  Sept 2023
"""

import numpy as np
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
    assert_almost_equal(grid.at_node["snowpack__depth"], 10 / 3)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 0)
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 2090000)


def test_assign_parameters():
    """Test when parameters are not default value"""

    grid = RasterModelGrid((2, 2))
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 1, at="node")
    grid.add_full("atmosphere_bottom_air__temperature", 1, at="node")
    grid.add_full("land_surface__temperature", -1, at="node")
    grid.add_full("land_surface_net-total-energy__energy_flux", 1, at="node")

    sm = SnowEnergyBalance(
        grid,
        rho_water=1001,
        rho_air=1.2,
        cp_air=1005,
        cp_snow=2080.0,
        rain_snow_temp=1.5,
        melting_point=0.1,
    )

    # parameters
    assert sm.rho_water == 1001
    assert sm.rho_air == 1.2
    assert sm.cp_air == 1005
    assert sm.cp_snow == 2080.0
    assert sm.rain_snow_temp == 1.5
    assert sm.melting_point == 0.1

    with pytest.raises(ValueError):
        sm.rho_water = -1

    with pytest.raises(ValueError):
        sm.rho_air = -1

    with pytest.raises(ValueError):
        sm.cp_air = -1

    with pytest.raises(ValueError):
        sm.cp_snow = -1


def test_snow_accumulation():
    """Test when there is only snow accumulation"""

    grid = RasterModelGrid((2, 2), xy_spacing=(10, 10))
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0.001, at="node")
    grid.add_full("atmosphere_bottom_air__temperature", -1, at="node")
    grid.add_full("land_surface__temperature", -1, at="node")
    grid.add_full("land_surface_net-total-energy__energy_flux", 0, at="node")
    grid.add_full("snowpack__liquid-equivalent_depth", 1, at="node")
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 200, at="node")

    sm = SnowEnergyBalance(grid)
    assert np.sum(sm.total_snow_melt_at_node * grid.dx * grid.dy) == pytest.approx(0.0)
    assert np.sum(
        grid.at_node["snowpack__liquid-equivalent_depth"] * grid.dx * grid.dy
    ) == pytest.approx(400)

    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 2090000)
    assert_almost_equal(grid.at_node["snowpack__depth"], 5)

    sm.run_one_step(1000)
    assert_almost_equal(grid.at_node["snowpack__depth"], 10)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 0)
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 2090000)
    assert np.sum(sm.total_snow_melt_at_node * grid.dx * grid.dy) == pytest.approx(0.0)
    assert np.sum(
        grid.at_node["snowpack__liquid-equivalent_depth"] * grid.dx * grid.dy
    ) == pytest.approx(800)


def test_snow_melt():
    """Test when there is only snow melt"""

    grid = RasterModelGrid((2, 2), xy_spacing=(10, 10))
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0, at="node")
    grid.add_full("atmosphere_bottom_air__temperature", 14, at="node")
    grid.add_full("land_surface__temperature", 1, at="node")
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 200, at="node")
    grid.add_full(
        "land_surface_net-total-energy__energy_flux", 334000, at="node"
    )  # energy used to melt 1m with 1000s, cold content=0
    grid.add_full("snowpack__liquid-equivalent_depth", 2, at="node")

    sm = SnowEnergyBalance(grid)
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 0)
    assert_almost_equal(grid.at_node["snowpack__depth"], 10)

    sm.run_one_step(1000)
    assert_almost_equal(grid.at_node["snowpack__depth"], 5)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 0.001)
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 0)

    assert np.sum(sm.total_snow_melt_at_node * grid.dx * grid.dy) == pytest.approx(400)
    assert np.sum(
        grid.at_node["snowpack__liquid-equivalent_depth"] * grid.dx * grid.dy
    ) == pytest.approx(400)


def test_snow_melt_accumulation():
    """Test when there is both snow melt and accumulation"""

    grid = RasterModelGrid((2, 2), xy_spacing=(10, 10))
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0.002, at="node")
    grid.add_full("atmosphere_bottom_air__temperature", 0.5, at="node")
    grid.add_full("land_surface__temperature", -1, at="node")
    grid.add_full(
        "land_surface_net-total-energy__energy_flux", 334000 + 2090, at="node"
    )  # energy to melt 1m in 1000 sec (consider existing cold content)
    grid.add_full("snowpack__liquid-equivalent_depth", 1, at="node")
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 200, at="node")

    sm = SnowEnergyBalance(grid)
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 2090000)

    sm.run_one_step(1000)
    assert_almost_equal(grid.at_node["snowpack__liquid-equivalent_depth"], 2)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 0.001)
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 0)

    assert np.sum(sm.total_snow_melt_at_node * grid.dx * grid.dy) == pytest.approx(400)
    assert np.sum(
        grid.at_node["snowpack__liquid-equivalent_depth"] * grid.dx * grid.dy
    ) == pytest.approx(800)


def test_multiple_steps():
    """Test multiple runs with changing inputs at each time step"""

    # step1: snow melt
    grid = RasterModelGrid((2, 2), xy_spacing=(10, 10))
    grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0, at="node")
    grid.add_full("atmosphere_bottom_air__temperature", 1, at="node")
    grid.add_full("land_surface__temperature", -1, at="node")
    grid.add_full(
        "land_surface_net-total-energy__energy_flux", 2090 + 33400, at="node"
    )  # energy to melt 0.1m in 1000sec
    grid.add_full("snowpack__liquid-equivalent_depth", 1, at="node")
    grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 200, at="node")
    init_swe = grid.at_node["snowpack__liquid-equivalent_depth"].copy()

    sm = SnowEnergyBalance(grid)
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 2090000)

    sm.run_one_step(1000)
    assert_almost_equal(grid.at_node["snowpack__liquid-equivalent_depth"], 0.9)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 0.0001)
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 0)
    assert_almost_equal(sm.total_snow_precip_at_node, 0)
    assert_almost_equal(sm.total_snow_melt_at_node, 0.1)

    # step2: change prec_rate, air_temp, q_sum for snow accumulation
    grid.at_node["atmosphere_bottom_air__temperature"].fill(-1)
    grid.at_node["atmosphere_water__precipitation_leq-volume_flux"].fill(0.0001)
    grid.at_node["land_surface_net-total-energy__energy_flux"].fill(0)

    sm.run_one_step(1000)
    assert_almost_equal(grid.at_node["snowpack__depth"], 5)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 0)
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 0)
    assert_almost_equal(sm.total_snow_precip_at_node, 0.1)
    assert_almost_equal(sm.total_snow_melt_at_node, 0.1)

    # step3: change air_temp and q_sum for snow melt & accumulation
    grid.at_node["atmosphere_bottom_air__temperature"].fill(0.5)
    grid.at_node["atmosphere_water__precipitation_leq-volume_flux"].fill(0.002)
    grid.at_node["land_surface_net-total-energy__energy_flux"].fill(334000)

    sm.run_one_step(1000)
    assert_almost_equal(grid.at_node["snowpack__depth"], 10)
    assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 0.001)
    assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 0)
    assert_almost_equal(sm.total_snow_precip_at_node, 2.1)
    assert_almost_equal(sm.total_snow_melt_at_node, 1.1)

    # mass balance check
    in_out = sm.total_snow_precip_at_node - sm.total_snow_melt_at_node

    store = grid.at_node["snowpack__liquid-equivalent_depth"] - init_swe

    assert_almost_equal(in_out, store, err_msg="Error for mass balance", decimal=11)


@pytest.mark.parametrize(
    "precip", (1.0, [2.0], [1.0, 2.0, 3.0], [[1.0, 2.0], [4.0, 5.0]])
)
def test_calc_p_snow(precip):
    assert_almost_equal(SnowEnergyBalance.calc_precip_snow(precip, -3.0, 1.0), precip)
    assert_almost_equal(SnowEnergyBalance.calc_precip_snow(precip, 3.0, 1.0), 0.0)


def test_calc_snow_melt_rate():
    q_sum = np.full(4, 334000)
    cold_content = np.full(4, 0)
    rho_water = 1000
    dt = 1000
    snow_melt_rate = SnowEnergyBalance.calc_snow_melt_rate(
        q_sum, cold_content, rho_water, dt
    )
    assert np.allclose(snow_melt_rate, 0.001)


def test_update_swe():
    h_swe = np.full(4, 0.01)
    sm_rate = np.full(4, 0.001)
    p_snow = np.full(4, 0.003)

    SnowEnergyBalance.calc_swe(p_snow, sm_rate, h_swe, dt=5.0, out=h_swe)

    assert np.all(h_swe == 0.02)


def test_update_snow_depth():
    h_snow = SnowEnergyBalance.calc_snow_depth([0.01, 0.01, 0.01], 2.0)
    assert np.all(h_snow == 0.02)


def test_initialize_cold_content():
    rho_snow = np.full(4, 300)
    h_snow = np.full(4, 1)
    surf_temp = np.full(4, -1)
    cold_content = SnowEnergyBalance.initialize_cold_content(
        rho_snow, h_snow, surf_temp, 2090, 0
    )
    assert np.all(cold_content == 627000.0)


def test_update_cold_content():
    q_sum = np.full(4, 1000)
    cold_content = np.array([1000, 2000, 3000, 4000])

    SnowEnergyBalance.calc_cold_content(q_sum, cold_content, 2, out=cold_content)
    assert np.allclose(cold_content, np.array([0, 0, 1000, 2000]))
