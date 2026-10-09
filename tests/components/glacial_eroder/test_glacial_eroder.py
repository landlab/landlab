#! /usr/bin/env python
"""
Unit tests for landlab.components.glacial_eroder.glacial_eroder
"""

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_raises

from landlab import RasterModelGrid
from landlab.components import GlacialEroder


@pytest.mark.parametrize(
    "fixed,dynamic,error",
    [
        (False, False, "Must have either fixed or dynamic ice thickness"),
        (True, False, None),
        (False, True, None),
        (True, True, "Cannot have both fixed and dynamic ice thickness"),
    ],
)
def test_requires_exactly_one_ice_thickness_field(fixed, dynamic, error):
    grid = RasterModelGrid((3, 3))
    grid.add_zeros("topographic__elevation", at="node")
    if fixed:
        grid.add_zeros("fixed_ice__thickness", at="node")
    if dynamic:
        grid.add_zeros("dynamic_ice__thickness", at="node")

    if error:
        with pytest.raises(ValueError, match=error):
            GlacialEroder(grid)
    else:
        component = GlacialEroder(grid)
        assert component._evolve_ice_thickness is dynamic
        assert component._ice_string == (
            "dynamic_ice__thickness" if dynamic else "fixed_ice__thickness"
        )


@pytest.mark.parametrize("at", ("node",))
def test_ela_height_as_string(at):
    grid = RasterModelGrid((5, 5))
    grid.add_ones("topographic__elevation", at="node")
    grid.add_zeros("fixed_ice__thickness", at="node")
    ela_height = grid.add_ones("ela_height", at=at)
    ela_height[0] = 999

    glacial_comp = GlacialEroder(grid, ela_height="ela_height")

    assert_array_almost_equal(glacial_comp._ela_height, ela_height)


@pytest.mark.parametrize("at", ("node",))
def test_ela_height_as_number(at):
    grid = RasterModelGrid((5, 5))
    grid.add_ones("topographic__elevation", at="node")
    grid.add_zeros("fixed_ice__thickness", at="node")
    ela_height = grid.ones(at=at) * 0.01

    glacial_comp = GlacialEroder(grid, ela_height=ela_height)

    assert_array_almost_equal(glacial_comp._ela_height, ela_height)


@pytest.mark.parametrize("at", ("patch", "corner", "face", "cell", "link"))
def test_diffusion_as_bad_string(at):
    grid = RasterModelGrid((5, 5))
    grid.add_ones("topographic__elevation", at="node")
    grid.add_zeros("fixed_ice__thickness", at="node")
    grid.add_ones("ela_height", at=at)

    with pytest.raises(ValueError):
        GlacialEroder(grid, ela_height="ela_height")


@pytest.mark.parametrize("at", ("patch", "corner", "face", "cell", "link"))
def test_diffusion_as_bad_number(at):
    grid = RasterModelGrid((5, 5))
    grid.add_ones("topographic__elevation", at="node")
    ela_height = grid.ones(at=at) * 0.01

    with pytest.raises(ValueError):
        GlacialEroder(grid, ela_height=ela_height)


def test_glacier_dynamic_advection():
    dt = 10 * 60 * 60 * 24 * 365.25  # one year in seconds
    time_to_run = 10000 * 60 * 60 * 24 * 365.25  # one thousand year in seconds

    mg = RasterModelGrid(
        (4, 11), xy_spacing=(1000.0, 1000.0), xy_of_lower_left=(-5000.0, 0.0)
    )

    # create the fields in the grid
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("dynamic_ice__thickness", at="node")

    mg.set_closed_boundaries_at_grid_edges(True, True, True, True)

    slope = 10
    mg.at_node["topographic__elevation"] += mg.x_of_node * np.sin(np.deg2rad(slope))

    positive_topo = mg.at_node["topographic__elevation"] > 0.0
    mg.at_node["dynamic_ice__thickness"][positive_topo] = 2000.0
    mg.at_node["dynamic_ice__thickness"][mg.boundary_nodes] = 0.0

    # instantiate:
    ge = GlacialEroder(
        mg,
        ela_height=0.0,
        ela_thickness_band=1.0,
        characteristic_u=20.0 / 60 / 60 / 24 / 365.25,
        glacial_erosion_coefficient=1e-20,
        Q_geo=0,
        Q_accum=0,
        Q_ablat=0,
        latent_heat=1e10,
    )

    # perform the loop:
    elapsed_time = 0.0
    while elapsed_time < time_to_run:
        if elapsed_time + dt > time_to_run:
            dt = time_to_run - elapsed_time
        ge.run_one_step(dt)
        elapsed_time += dt

    ice_thickness_target = np.array(
        [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            101.46802122,
            84.25360602,
            128.74497872,
            271.88224825,
            960.52386498,
            1407.17248361,
            1700.70971928,
            1865.42226128,
            1479.12305695,
            0,
            0,
            101.46802122,
            84.25360602,
            128.74497872,
            271.88224825,
            960.52386498,
            1407.17248361,
            1700.70971928,
            1865.42226128,
            1479.12305695,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
    )

    assert_array_almost_equal(
        mg.at_node["dynamic_ice__thickness"], ice_thickness_target
    )


def test_fixed_glacier():
    dt = 10 * 60 * 60 * 24 * 365.25  # one year in seconds
    time_to_run = 1000 * 60 * 60 * 24 * 365.25  # one thousand year in seconds

    mg = RasterModelGrid(
        (4, 11), xy_spacing=(1000.0, 1000.0), xy_of_lower_left=(-5000.0, 0.0)
    )

    # create the fields in the grid
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("fixed_ice__thickness", at="node")

    mg.set_closed_boundaries_at_grid_edges(True, True, True, True)

    slope = 10
    mg.at_node["topographic__elevation"] += mg.x_of_node * np.sin(np.deg2rad(slope))

    positive_topo = mg.at_node["topographic__elevation"] > 0.0
    mg.at_node["fixed_ice__thickness"][positive_topo] = 2000.0
    mg.at_node["fixed_ice__thickness"][mg.boundary_nodes] = 0.0

    # instantiate:
    ge = GlacialEroder(
        mg,
        ela_height=0.0,
        ela_thickness_band=1.0,
        characteristic_u=20.0 / 60 / 60 / 24 / 365.25,
        glacial_erosion_coefficient=1e-20,
        Q_geo=0,
        Q_accum=0,
        Q_ablat=0,
        latent_heat=1e10,
    )

    # perform the loop:
    elapsed_time = 0.0
    while elapsed_time < time_to_run:
        if elapsed_time + dt > time_to_run:
            dt = time_to_run - elapsed_time
        ge.run_one_step(dt)
        elapsed_time += dt

    ice_thickness_target = np.array(
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
            0.0,
            0.0,
            2000.0,
            2000.0,
            2000.0,
            2000.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            2000.0,
            2000.0,
            2000.0,
            2000.0,
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

    assert_array_almost_equal(mg.at_node["fixed_ice__thickness"], ice_thickness_target)


def test_exception_handling():
    grid = RasterModelGrid((3, 3))
    grid.add_zeros("topographic__elevation", at="node")
    grid.add_zeros("fixed_ice__thickness", at="node")
    assert_raises(ValueError, GlacialEroder, grid, rho_ice=-1)
    assert_raises(ValueError, GlacialEroder, grid, latent_heat=-1)
    assert_raises(ValueError, GlacialEroder, grid, ice_flow_coefficient=-1)
    assert_raises(ValueError, GlacialEroder, grid, gravity=-1)
    assert_raises(ValueError, GlacialEroder, grid, glacial_erosion_coefficient=-1)
    assert_raises(ValueError, GlacialEroder, grid, courant=-1)
