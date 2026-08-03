"""Tests for the GeoEnthalpyDelta component.

These tests exercise the fixed-grid, enthalpy-style finite-volume update
described in Lorenzo-Trueba et al., *GeoEnthalpy-Delta v1.0*: sediment
mass conservation, the positivity-preserving thickness update, the
feeder-mask construction, and the CFL-limited stable time step.
"""

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import HexModelGrid
from landlab import RasterModelGrid
from landlab.components import GeoEnthalpyDelta


def make_grid(nrows=20, ncols=30, dx=0.5, dy=0.5):
    grid = RasterModelGrid((nrows, ncols), xy_spacing=(dx, dy))
    x = grid.x_of_node.reshape(grid.shape)
    basement = grid.add_zeros("basement__elevation", at="node")
    basement[:] = (-x).reshape(-1)
    return grid


def test_requires_raster_model_grid():
    grid = HexModelGrid((5, 5))
    grid.add_zeros("basement__elevation", at="node")
    with pytest.raises(TypeError):
        GeoEnthalpyDelta(grid)


def test_feeder_width_smaller_than_spacing_raises():
    grid = make_grid()
    with pytest.raises(ValueError):
        GeoEnthalpyDelta(grid, feeder_width=1.0e-6)


def test_fields_created_and_initial_topography():
    grid = make_grid()
    esd = GeoEnthalpyDelta(grid)

    assert "sediment__thickness" in grid.at_node
    assert "topographic__elevation" in grid.at_node
    assert_array_equal(grid.at_node["sediment__thickness"], 0.0)
    assert_array_equal(
        grid.at_node["topographic__elevation"], grid.at_node["basement__elevation"]
    )
    assert esd.time_elapsed == 0.0
    assert esd.sea_level == esd.sea_level_start


def test_feeder_mask_matches_requested_width():
    dy = 0.5
    grid = make_grid(dy=dy)
    feeder_width = 2.0
    esd = GeoEnthalpyDelta(grid, feeder_width=feeder_width)

    n_feeder = np.count_nonzero(esd._feeder_mask)
    # The feeder should span close to feeder_width, not several times wider.
    assert n_feeder * dy == pytest.approx(feeder_width + dy, abs=dy)


def test_calc_stable_time_step_matches_cfl_formula():
    grid = make_grid(dx=0.2, dy=0.2)
    esd = GeoEnthalpyDelta(
        grid,
        topset_diffusivity_x=1.0,
        topset_diffusivity_y=1.0,
        foreset_diffusivity_x=3.0,
        foreset_diffusivity_y=2.0,
        cfl=0.4,
    )
    d_max = 3.0
    expected = 0.4 / (2.0 * d_max * (1.0 / 0.2**2 + 1.0 / 0.2**2))
    assert esd.calc_stable_time_step() == pytest.approx(expected)


def test_mass_is_conserved():
    grid = make_grid(nrows=20, ncols=40, dx=0.2, dy=0.2)
    sediment_flux = 0.5
    esd = GeoEnthalpyDelta(
        grid,
        sediment_flux=sediment_flux,
        feeder_width=1.0,
        sea_level_start=-5.0,
        topset_threshold_x=0.25,
        topset_threshold_y=0.1,
        foreset_threshold_x=2.0,
        foreset_threshold_y=2.0,
    )
    dt = esd.calc_stable_time_step()
    nsteps = 200
    for _ in range(nsteps):
        esd.run_one_step(dt)

    modeled_volume = np.sum(grid.at_node["sediment__thickness"]) * grid.dx * grid.dy
    expected_volume = sediment_flux * nsteps * dt
    assert modeled_volume == pytest.approx(expected_volume, rel=1e-6)


def test_thickness_stays_nonnegative():
    grid = make_grid(nrows=20, ncols=40, dx=0.2, dy=0.2)
    esd = GeoEnthalpyDelta(
        grid,
        sediment_flux=0.5,
        feeder_width=1.0,
        sea_level_start=-5.0,
        topset_threshold_x=0.25,
        topset_threshold_y=0.1,
        foreset_threshold_x=2.0,
        foreset_threshold_y=2.0,
    )
    dt = esd.calc_stable_time_step()
    for _ in range(200):
        esd.run_one_step(dt)
        assert np.all(grid.at_node["sediment__thickness"] >= 0.0)


def test_sea_level_rises_with_time():
    grid = make_grid()
    esd = GeoEnthalpyDelta(grid, sea_level_start=-5.0, sea_level_rise_rate=0.1)
    dt = esd.calc_stable_time_step()
    for _ in range(10):
        esd.run_one_step(dt)
    assert esd.sea_level == pytest.approx(-5.0 + 0.1 * esd.time_elapsed)


def test_run_one_step_feeder_mask_override_persists():
    grid = make_grid(nrows=20, ncols=30)
    esd = GeoEnthalpyDelta(grid, sediment_flux=0.5, feeder_width=1.0)
    dt = esd.calc_stable_time_step()

    default_mask = esd._feeder_mask.copy()
    custom_mask = np.zeros_like(default_mask)
    custom_mask[5:8] = True

    esd.run_one_step(dt, feeder_mask=custom_mask)
    assert_array_equal(esd._feeder_mask, custom_mask)

    # The override should stick for later calls that don't pass a mask.
    esd.run_one_step(dt)
    assert_array_equal(esd._feeder_mask, custom_mask)
    assert not np.array_equal(esd._feeder_mask, default_mask)
