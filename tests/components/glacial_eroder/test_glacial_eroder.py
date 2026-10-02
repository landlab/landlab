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
    # In this test, the ice thickness should evolve by flowing down the imposed
    # slope.
    assert_array_almost_equal(
        mg.at_node["dynamic_ice__thickness"], ice_thickness_target
    )


def test_fixed_glacier_advection():
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
    # Check that the ice thickness after advection matches the target. In this test,
    # the ice thickness should not actually advect.
    assert_array_almost_equal(mg.at_node["fixed_ice__thickness"], ice_thickness_target)


def test_dynamic_glacier_erosion():
    dt = 50.0 * 60 * 60 * 24 * 365.25  # one year in seconds
    time_to_run = 500.0 * 60 * 60 * 24 * 365.25  # five hundred years in seconds

    mg = RasterModelGrid((11, 11), xy_spacing=(50.0, 50.0))

    # create the fields in the grid
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("dynamic_ice__thickness", at="node")

    uplift_rate = (
        10.0
        * np.maximum(
            1.0 - np.hypot(mg.x_of_node - 250.0, mg.y_of_node - 250.0) / 250.0,
            0.0,
        )
        / 60
        / 60
        / 24
        / 365.25
    )  # convert to meters per year

    mg.set_fixed_value_boundaries_at_grid_edges(True, True, True, True)

    # instantiate:
    ge = GlacialEroder(
        mg,
        ela_height=500.0,
        ela_thickness_band=0.0,
        characteristic_u=20.0 / 60 / 60 / 24 / 365.25,
        glacial_erosion_coefficient=1,
        max_Q_accum=1 / 60 / 60 / 24 / 365.25,
    )

    # perform the loop:
    elapsed_time = 0.0  # total time in simulation
    while elapsed_time < time_to_run:
        if elapsed_time + dt > time_to_run:
            dt = time_to_run - elapsed_time
        mg.at_node["topographic__elevation"] += uplift_rate * dt
        ge.run_one_step(dt)
        elapsed_time += dt
    # Check that the ice thickness after advection matches the target. In this test,
    # the ice thickness should not actually advect.

    elevation_target = np.array(
        [
            527.85997328,
            876.87300532,
            999.97992185,
            876.87433079,
            527.84313301,
            757.34169801,
            1394.41787353,
            1837.68791559,
            1999.96720345,
            1837.69032937,
            1394.41545868,
            757.32915759,
            527.85997328,
            1394.41787353,
            2171.53763553,
            2763.89352016,
            2999.96265339,
            2763.89614707,
            2171.53574177,
            1394.41362014,
            527.859824,
            876.87300532,
            1837.68791559,
            2763.89352016,
            3585.74756562,
            3999.96130881,
            3585.75082454,
            2763.89516471,
            1837.68244052,
            876.87260319,
            999.97992185,
            1999.96720345,
            2999.96265339,
            3999.96130881,
            4999.9594712,
            3999.9594712,
            2999.96226441,
            1999.95941779,
            999.97982921,
            876.87433079,
            1837.69032937,
            2763.89614707,
            3585.75082454,
            3999.9594712,
            3585.78643763,
            2763.89288727,
            1837.6825396,
            876.87420361,
            527.84313301,
            1394.41545868,
            2171.53574177,
            2763.89516471,
            2999.96226441,
            2763.89288727,
            2171.53401391,
            1394.44872454,
            527.83275852,
            757.32915759,
            1394.41362014,
            1837.68244052,
            1999.95941779,
            1837.6825396,
            1394.44872454,
            757.32765189,
            527.859824,
            876.87260319,
            999.97982921,
            876.87420361,
            527.83275852,
        ]
    )

    topo_mask = mg.at_node["topographic__elevation"] > 0.0
    assert_array_almost_equal(
        mg.at_node["topographic__elevation"][topo_mask], elevation_target
    )


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
