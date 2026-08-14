"""
Unit tests for the GeoEnthalpyDelta component.

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

    assert "sediment__thickness" not in grid.at_node
    assert "topographic__elevation" in grid.at_node
    assert_array_equal(esd.sediment_thickness, 0.0)
    assert_array_equal(
        grid.at_node["topographic__elevation"], grid.at_node["basement__elevation"]
    )
    assert esd.time_elapsed == 0.0
    assert esd.sea_level == esd._sea_level_start


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
        topset_diffusivity=(1.0, 1.0),
        foreset_diffusivity=(3.0, 2.0),
        cfl=0.4,
    )
    d_max = 3.0
    expected = 0.4 / (2.0 * d_max * (1.0 / 0.2**2 + 1.0 / 0.2**2))
    assert esd._calc_stable_time_step() == pytest.approx(expected)


def test_mass_is_conserved():
    grid = make_grid(nrows=20, ncols=40, dx=0.2, dy=0.2)
    sediment_flux = 0.5
    esd = GeoEnthalpyDelta(
        grid,
        sediment_flux=sediment_flux,
        feeder_width=1.0,
        sea_level_start=-5.0,
        topset_threshold=(0.25, 0.1),
        foreset_threshold=(2.0, 2.0),
    )
    nsteps = 200
    for _ in range(nsteps):
        esd.run_one_step()

    modeled_volume = np.sum(esd.sediment_thickness) * grid.dx * grid.dy
    expected_volume = sediment_flux * esd.time_elapsed
    assert modeled_volume == pytest.approx(expected_volume, rel=1e-6)


def test_thickness_stays_nonnegative():
    grid = make_grid(nrows=20, ncols=40, dx=0.2, dy=0.2)
    esd = GeoEnthalpyDelta(
        grid,
        sediment_flux=0.5,
        feeder_width=1.0,
        sea_level_start=-5.0,
        topset_threshold=(0.25, 0.1),
        foreset_threshold=(2.0, 2.0),
    )
    for _ in range(200):
        esd.run_one_step()
        assert np.all(esd.sediment_thickness >= 0.0)


def test_sea_level_rises_with_time():
    grid = make_grid()
    esd = GeoEnthalpyDelta(grid, sea_level_start=-5.0, sea_level_rise_rate=0.1)
    for _ in range(10):
        esd.run_one_step()
    assert esd.sea_level == pytest.approx(-5.0 + 0.1 * esd.time_elapsed)


def test_run_one_step_dt_defaults_to_stable_time_step():
    grid = make_grid()
    esd = GeoEnthalpyDelta(grid)
    esd.run_one_step()
    assert esd.time_elapsed == pytest.approx(esd._calc_stable_time_step())


def test_run_one_step_subcycles_dt_larger_than_stable_limit():
    n_substeps = 5

    grid_direct = make_grid()
    esd_direct = GeoEnthalpyDelta(grid_direct, sediment_flux=0.5, feeder_width=1.0)
    stable_dt = esd_direct._calc_stable_time_step()
    esd_direct.run_one_step(dt=n_substeps * stable_dt)

    grid_manual = make_grid()
    esd_manual = GeoEnthalpyDelta(grid_manual, sediment_flux=0.5, feeder_width=1.0)
    for _ in range(n_substeps):
        esd_manual.run_one_step()

    assert esd_direct.time_elapsed == pytest.approx(esd_manual.time_elapsed)
    np.testing.assert_allclose(
        esd_direct.sediment_thickness, esd_manual.sediment_thickness, atol=1e-12
    )
