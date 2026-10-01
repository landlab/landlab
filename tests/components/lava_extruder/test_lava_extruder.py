#! /usr/bin/env python
"""
Unit tests for landlab.components.lava_extruder.lava_extruder
"""

import numpy as np
import pytest
import requireit
from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_equal

import landlab.components.lava_extruder.lava_helper_functions as lhf
from landlab import RasterModelGrid
from landlab.components import LavaExtruder
from landlab.core.component_utils import iter_adaptive_time_steps


def test_thermal_energy_to_temperature():
    temperature = lhf.thermal_energy_to_temperature(
        energy=1e6,
        lava_thickness=1.0,
        cell_area=1.0,
        lava_density=1000,
        specific_heat_capacity=1.0,
    )

    assert_almost_equal(temperature, 1000.0)


def test_vent_flux_calculation():
    vent_flux = lhf.vent_flux(
        eruptive_volume_flux=10.0,
        cell_area=2.0,
    )
    assert_almost_equal(vent_flux, 5.0)


def test_vent_thermal_energy():
    energy = lhf.vent_thermal_energy(
        temperature=1200.0,
        lava_thickness=1.0,
        cell_area=1.0,
        lava_density=1000,
        specific_heat_capacity=1.0,
    )
    assert_almost_equal(energy, 1.2e6)


def test_calc_conductive_heat_flux():
    conductive_heat_flux = lhf.conductive_heat_flux(
        lava_temperature=1500,
        area=1.0,
        conductive_thickness=1.0,
        thermal_conductivity=1.0,
        ambient_temperature=500.0,
    )
    assert_almost_equal(conductive_heat_flux, 1000.0)


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


def test_volume_conservation():
    max_dt = 30.0
    time_to_run = 300.0

    mg = RasterModelGrid((3, 10), xy_spacing=(100.0, 100.0))

    # create the fields in the grid
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("lava__thickness", at="node")
    mg.at_node["topographic__elevation"][:] += mg.x_of_node

    initial_volume = np.sum(mg.at_node["topographic__elevation"] * mg.cell_area_at_node)
    mg.set_fixed_value_boundaries_at_grid_edges(True, True, True, True)

    eruptive_volume_flux = 10.0

    le = LavaExtruder(
        mg,
        vent_coordinates=[[500, 100]],
        eruptive_volume_flux=eruptive_volume_flux,
        viscosity=1.0e10,
        specific_heat_capacity=300,
        latent_heat_of_freezing=1,
    )

    # perform the loop:
    for dt in iter_adaptive_time_steps(time_to_run, calc_dt=lambda: max_dt):
        le.run_one_step(dt)

    actual_added_volume = eruptive_volume_flux * time_to_run
    model_added_volume = np.sum(
        (mg.at_node["topographic__elevation"] + mg.at_node["lava__thickness"])
        * mg.cell_area_at_node
    )

    assert_almost_equal(actual_added_volume, model_added_volume - initial_volume)


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


def test_advection_mass_and_energy_conservation():
    max_dt = 100.0
    time_to_run = 600.0

    lava_density = 2700.0
    specific_heat_capacity = 1000.0
    freezing_temperature = 400.0
    eruption_temperature = 600.0
    ambient_temperature = 150.0

    mg = RasterModelGrid((3, 20), xy_spacing=(100.0, 100.0))

    # create the fields in the grid
    mg.add_zeros("topographic__elevation", at="node")
    mg.at_node["topographic__elevation"][:] += mg.x_of_node

    mg.add_zeros("lava__thickness", at="node")
    mg.at_node["lava__thickness"][mg.find_nearest_node([1500.0, 100.0])] += 10.0

    mg.add_zeros("lava__temperature", at="node")
    mg.at_node["lava__temperature"][
        mg.at_node["lava__thickness"] > 0
    ] = eruption_temperature

    initial_lava_volume = np.sum(mg.at_node["lava__thickness"] * mg.cell_area_at_node)
    initial_lava_energy = np.sum(
        mg.at_node["lava__thickness"]
        * (mg.at_node["lava__temperature"])
        * mg.cell_area_at_node
        * specific_heat_capacity
        * lava_density
    )

    mg.set_closed_boundaries_at_grid_edges(True, True, True, True)

    le = LavaExtruder(
        mg,
        vent_coordinates=None,
        eruptive_volume_flux=0.0,
        viscosity=8.0e3,
        specific_heat_capacity=specific_heat_capacity,
        lava_density=lava_density,
        thermal_conductivity=1e-50,
        freezing_temperature=freezing_temperature,
        eruption_temperature=eruption_temperature,
        emissivity=0,
        ambient_temperature=ambient_temperature,
    )

    # perform the loop:
    for dt in iter_adaptive_time_steps(time_to_run, calc_dt=lambda: max_dt):
        le.run_one_step(dt)

    final_lava_volume = np.sum(mg.at_node["lava__thickness"] * mg.cell_area_at_node)
    final_lava_energy = np.sum(
        mg.at_node["lava__thickness"]
        * (mg.at_node["lava__temperature"])
        * mg.cell_area_at_node
        * specific_heat_capacity
        * lava_density
    )

    assert_almost_equal(final_lava_energy / 1e10, initial_lava_energy / 1e10)
    assert_almost_equal(final_lava_volume, initial_lava_volume)


def test_freezing_mass_conservation():
    max_dt = 100.0
    time_to_run = 600.0

    T_freeze = 800

    mg = RasterModelGrid((3, 20), xy_spacing=(100.0, 100.0))

    # create the fields in the grid
    mg.add_zeros("topographic__elevation", at="node")

    mg.add_zeros("lava__thickness", at="node")
    mg.at_node["lava__thickness"] += 10.0

    mg.add_zeros("lava__temperature", at="node")
    mg.at_node["lava__temperature"][mg.at_node["lava__thickness"] > 0] = T_freeze

    initial_lava_volume = np.sum(mg.at_node["lava__thickness"] * mg.cell_area_at_node)

    mg.set_closed_boundaries_at_grid_edges(True, True, True, True)

    le = LavaExtruder(
        mg,
        vent_coordinates=None,
        eruptive_volume_flux=0.0,
        viscosity=1.0e50,
        freezing_temperature=T_freeze,
        specific_heat_capacity=1.0,
        latent_heat_of_freezing=1.0,
    )

    # perform the loop:
    for dt in iter_adaptive_time_steps(time_to_run, calc_dt=lambda: max_dt):
        le.run_one_step(dt)

    final_topo_volume = np.sum(
        mg.at_node["topographic__elevation"] * mg.cell_area_at_node
    )
    assert_almost_equal(final_topo_volume, initial_lava_volume)


def test_vent_outside_domain():
    grid = RasterModelGrid((3, 3))
    grid.add_zeros("topographic__elevation", at="node")
    grid.add_zeros("lava__thickness", at="node")

    with pytest.raises(ValueError, match="Vent coordinates must be within the domain"):
        LavaExtruder(grid, vent_coordinates=[(10, 10)])


def test_bad_initial_lava_with_temperature_below_freezing():
    T_freeze = 800

    mg = RasterModelGrid((3, 3), xy_spacing=(100.0, 100.0))
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("lava__thickness", at="node")
    mg.at_node["lava__thickness"] += 10.0

    mg.add_zeros("lava__temperature", at="node")
    mg.at_node["lava__temperature"][mg.at_node["lava__thickness"] > 0] = T_freeze - 100

    with pytest.raises(
        ValueError, match="initial temperature that is at or above the lava freezing"
    ):
        LavaExtruder(
            mg,
            vent_coordinates=None,
            eruptive_volume_flux=0.0,
            viscosity=1.0e50,
            freezing_temperature=T_freeze,
            specific_heat_capacity=1.0,
            latent_heat_of_freezing=1.0,
        )
