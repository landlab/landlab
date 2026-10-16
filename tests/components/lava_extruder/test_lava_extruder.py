#! /usr/bin/env python
"""
Unit tests for landlab.components.lava_extruder.lava_extruder
"""

import numpy as np
import pytest
import requireit
from numpy.testing import assert_array_equal

import landlab.components.lava_extruder.lava_helper_functions as lhf
from landlab import RasterModelGrid
from landlab.components import LavaExtruder
from landlab.core.component_utils import iter_adaptive_time_steps


def test_calc_lava_speed():
    velocity = lhf.calc_lava_speed(
        np.array([1.0, 2.0, 3.0, 4.0]),
        np.array([-0.1, 0.5, 0.0, -1.0]),
        lava_density=3.0,
        gravity=2.0,
        viscosity=3.0,
        yield_stress=1.0,
    )

    assert_array_equal(velocity, np.array([0.0, 1.0, 0.0, 10.0]))


def test_lava_flow():
    max_dt = 30.0
    time_to_run = 600.0

    mg = RasterModelGrid((3, 10), xy_spacing=(100.0, 100.0))

    # create the fields in the grid
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("lava__thickness", at="node")
    mg.at_node["topographic__elevation"][:] += mg.x_of_node

    mg.set_fixed_value_boundaries_at_grid_edges(True, True, True, True)

    # instantiate:
    lf = LavaExtruder(
        mg,
        vent_coordinates=[[800, 100]],
        eruptive_volume_flux=10.0,
        viscosity=100.0,
        specific_heat_capacity=300,
    )

    # perform the loop:
    for dt in iter_adaptive_time_steps(time_to_run, calc_dt=lambda: max_dt):
        lf.run_one_step(dt)

    # lava_thickness_target = np.array(
    #     [
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.00010652427957657601,
    #         8.58421868745749e-02,
    #         1.26872326e-01,
    #         1.26951209e-01,
    #         1.26951660e-01,
    #         1.26951662e-01,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #     ]
    # )

    # lid_thickness_target = np.array(
    #     [
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.015916797630756093,
    #         0.01040743512055451,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #         0.01,
    #     ]
    # )

    # lava_energy_target = np.array(
    #     [
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         18197649456.94168,
    #         147377488026.09647,
    #         361968200695.9663,
    #         692051597250.8744,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #     ]
    # )

    # lava_temperature_target = np.array(
    #     [
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         800.0,
    #         817.7077224998175,
    #         943.3208178680401,
    #         1152.003526587297,
    #         1473.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #         300.0,
    #     ]
    # )

    # freezing_rate_target = np.array(
    #     [
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.00010652427957657601,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #         0.0,
    #     ]
    # )

    # assert_array_almost_equal(mg.at_node["lid__thickness"], lid_thickness_target)
    # assert_array_almost_equal(mg.at_node["lava__thickness"], lava_thickness_target)
    # assert_array_almost_equal(mg.at_node["lava__thermal_energy"], lava_energy_target)
    # assert_array_almost_equal(mg.at_node["lava__temperature"], lava_temperature_target)
    # assert_array_almost_equal(mg.at_node["freezing__rate"], freezing_rate_target)


@pytest.mark.parametrize(
    "name, value",
    (
        ("lava_density", -1),
        ("emissivity", 1.5),
        ("specific_heat_capacity", -1),
        ("thermal_conductivity", -1),
        ("viscosity", -1),
        ("latent_heat_of_freezing", -1),
        ("eruptive_volume_flux", -1),
        ("gravity", -1),
        ("yield_stress", -1),
        ("eruption_temperature", -10),
        ("freezing_temperature", -10),
        ("ambient_temperature", -1),
        ("max_lava_thickness", -1),
        ("basal_thickness", -1),
        ("max_lid_thickness", -1),
        ("min_lid_thickness", -1),
        ("courant", -1),
        ("min_dt", -1),
    ),
)
def test_bad_param_raises(name, value):
    grid = RasterModelGrid((3, 3))
    grid.add_zeros("topographic__elevation", at="node")
    grid.add_zeros("lava__thickness", at="node")

    with pytest.raises(requireit.ValidationError, match=f"^{name} must be"):
        LavaExtruder(grid, **{name: value})


def test_vent_outside_domain():
    grid = RasterModelGrid((3, 3))
    grid.add_zeros("topographic__elevation", at="node")
    grid.add_zeros("lava__thickness", at="node")

    with pytest.raises(ValueError, match="Vent coordinates must be within the domain"):
        LavaExtruder(grid, vent_coordinates=[(10, 10)])
