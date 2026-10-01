#! /usr/bin/env python
"""
Unit tests for landlab.components.lava_flow.lava_flow
"""

import numpy as np
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_raises

from landlab import RasterModelGrid
from landlab.components import LavaFlow


def test_lava_flow():
    dt = 30.0
    time_to_run = 600.0

    mg = RasterModelGrid((3, 10), xy_spacing=(100.0, 100.0))

    # create the fields in the grid
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("lava__thickness", at="node")
    mg.at_node["topographic__elevation"][:] += mg.x_of_node

    mg.set_fixed_value_boundaries_at_grid_edges(True, True, True, True)

    # instantiate:
    lf = LavaFlow(
        mg, vent_coordinates=[800, 100], volume_flux=10.0, viscosity=100.0, Cp=300
    )

    # perform the loop:
    elapsed_time = 0.0  # total time in simulation
    while elapsed_time < time_to_run:
        if elapsed_time + dt > time_to_run:
            dt = time_to_run - elapsed_time
        lf.run_one_step(dt)
        elapsed_time += dt

    lava_thickness_target = np.array(
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            5.61827142e-05,
            8.57112135e-02,
            1.26872326e-01,
            1.26951209e-01,
            1.26951660e-01,
            1.26951662e-01,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    lid_thickness_target = np.array(
        [
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.0151681,
            0.01098781,
            0.01034984,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
            0.01,
        ]
    )

    lava_energy_target = np.array(
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            4.194765050074071e10,
            1.8516629389592276e11,
            3.614779190884509e11,
            6.920515972460206e11,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    print(mg.at_node["lava__thermal_energy"])

    lava_temperature_target = np.array(
        [
            300.0,
            300.0,
            300.0,
            300.0,
            300.0,
            300.0,
            300.0,
            300.0,
            300.0,
            300.0,
            300.0,
            300.0,
            300.0,
            300.0,
            800.0,
            840.91751298,
            980.07219682,
            1151.52675705,
            1473.0,
            300.0,
            300.0,
            300.0,
            300.0,
            300.0,
            300.0,
            300.0,
            300.0,
            300.0,
            300.0,
            300.0,
        ]
    )

    freezing_rate_target = np.array(
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            3.96145227e-05,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    assert_array_almost_equal(mg.at_node["lid__thickness"], lid_thickness_target)
    assert_array_almost_equal(mg.at_node["lava__thickness"], lava_thickness_target)
    assert_array_almost_equal(mg.at_node["lava__thermal_energy"], lava_energy_target)
    assert_array_almost_equal(mg.at_node["lava__temperature"], lava_temperature_target)
    assert_array_almost_equal(mg.at_node["freezing__rate"], freezing_rate_target)


def test_exception_handling():
    grid = RasterModelGrid((3, 3))
    grid.add_zeros("topographic__elevation", at="node")
    grid.add_zeros("lava__thickness", at="node")
    assert_raises(ValueError, LavaFlow, grid, rho_lava=-1)
    assert_raises(ValueError, LavaFlow, grid, Cp=-1)
    assert_raises(ValueError, LavaFlow, grid, k_lava=-1)
    assert_raises(ValueError, LavaFlow, grid, viscosity=-1)
    assert_raises(ValueError, LavaFlow, grid, L_freeze=-1)
    assert_raises(ValueError, LavaFlow, grid, volume_flux=-1)
    assert_raises(ValueError, LavaFlow, grid, gravity=-1)
    assert_raises(ValueError, LavaFlow, grid, yield_stress=-1)
    assert_raises(ValueError, LavaFlow, grid, T_erupt=-10)
    assert_raises(ValueError, LavaFlow, grid, T_freeze=-10)
    assert_raises(ValueError, LavaFlow, grid, T_surf=-1)
    assert_raises(ValueError, LavaFlow, grid, max_lava_thickness=-1)
    assert_raises(ValueError, LavaFlow, grid, basal_thickness=-1)
    assert_raises(ValueError, LavaFlow, grid, max_lid_thickness=-1)
    assert_raises(ValueError, LavaFlow, grid, min_lid_thickness=-1)
    assert_raises(ValueError, LavaFlow, grid, courant=-1)
    assert_raises(ValueError, LavaFlow, grid, min_dt=-1)
