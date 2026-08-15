"""
Unit tests for the GeoEnthalpyDelta component.

These tests exercise the fixed-grid, enthalpy-style finite-volume update
described in Lorenzo-Trueba et al., *GeoEnthalpy-Delta v1.0*: sediment
mass conservation, the positivity-preserving thickness update, the
sediment__influx feeder pattern, the fixed internal basement, and the
CFL-limited stable time step.
"""

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from requireit import ValidationError

from landlab import HexModelGrid
from landlab import RasterModelGrid
from landlab.components import GeoEnthalpyDelta
from landlab.field.errors import FieldError


def make_grid(nrows=20, ncols=30, dx=0.5, dy=0.5, sea_level=0.0):
    grid = RasterModelGrid((nrows, ncols), xy_spacing=(dx))
    x = grid.x_of_node.reshape(grid.shape)
    topo = grid.add_zeros("topographic__elevation", at="node")
    topo[:] = (-x).reshape(-1)
    grid.add_zeros("sediment__influx", at="node")
    grid.add_field("sea_level__elevation", sea_level, at="grid")
    return grid


def set_west_boundary_influx(grid, total_flux, width, y_center=None):
    """Populate sediment__influx at west-boundary nodes.

    Mimics the component's old feeder_width/feeder_y_center semantics:
    total_flux is spread evenly over the west-boundary nodes within
    width of y_center (grid's y-midpoint by default).
    """
    y_of_row = grid.y_of_node.reshape(grid.shape)[:, 0]
    if y_center is None:
        y_center = 0.5 * (y_of_row.min() + y_of_row.max())
    mask = np.abs(y_of_row - y_center) <= 0.5 * width + 1.0e-12
    influx_xy = grid.at_node["sediment__influx"].reshape(grid.shape)
    influx_xy[mask, 0] = total_flux / np.count_nonzero(mask)


def test_requires_raster_model_grid():
    grid = HexModelGrid((5, 5))
    grid.add_zeros("topographic__elevation", at="node")
    with pytest.raises(TypeError):
        GeoEnthalpyDelta(grid)


def test_missing_topographic_elevation_field_raises():
    grid = RasterModelGrid((20, 30), xy_spacing=0.5)
    with pytest.raises(FieldError):
        GeoEnthalpyDelta(grid)


def test_thresholds_must_be_nonnegative():
    grid = make_grid()
    with pytest.raises(ValidationError):
        GeoEnthalpyDelta(grid, topset_threshold=(-0.1, 0.0))
    with pytest.raises(ValidationError):
        GeoEnthalpyDelta(grid, foreset_threshold=(0.0, -2.0))


def test_diffusivities_must_be_positive():
    grid = make_grid()
    with pytest.raises(ValidationError):
        GeoEnthalpyDelta(grid, topset_diffusivity=(0.0, 1.0))
    with pytest.raises(ValidationError):
        GeoEnthalpyDelta(grid, foreset_diffusivity=(1.0, 0.0))


def test_cfl_must_be_in_zero_one_interval():
    grid = make_grid()
    with pytest.raises(ValidationError):
        GeoEnthalpyDelta(grid, cfl=0.0)
    with pytest.raises(ValidationError):
        GeoEnthalpyDelta(grid, cfl=1.1)


def test_fields_created_and_initial_topography():
    grid = make_grid()
    esd = GeoEnthalpyDelta(grid)

    assert "sediment__thickness" in grid.at_node
    assert_array_equal(grid.at_node["sediment__thickness"], 0.0)
    assert esd.time_elapsed == 0.0
    assert esd.sea_level == grid.at_grid["sea_level__elevation"]


def test_basement_is_computed_from_initial_topo_and_thickness():
    grid = make_grid()
    grid.at_node["topographic__elevation"][:] += 0.3
    grid.add_field("sediment__thickness", np.full(grid.number_of_nodes, 0.3), at="node")

    esd = GeoEnthalpyDelta(grid)

    expected_basement = (
        grid.at_node["topographic__elevation"] - grid.at_node["sediment__thickness"]
    )
    assert_array_equal(esd._basement, expected_basement)


def test_basement_stays_fixed_across_steps():
    grid = make_grid(nrows=20, ncols=40, dx=0.2, dy=0.2, sea_level=-5.0)
    set_west_boundary_influx(grid, total_flux=0.5, width=1.0)
    esd = GeoEnthalpyDelta(
        grid, topset_threshold=(0.25, 0.1), foreset_threshold=(2.0, 2.0)
    )
    basement_before = esd._basement.copy()

    for _ in range(20):
        esd.run_one_step()

    assert_array_equal(esd._basement, basement_before)
    np.testing.assert_allclose(
        grid.at_node["topographic__elevation"] - grid.at_node["sediment__thickness"],
        basement_before,
        atol=1e-12,
    )


def test_influx_only_enters_at_specified_west_boundary_nodes():
    grid = make_grid(nrows=10, ncols=10, dx=0.2, dy=0.2, sea_level=-5.0)
    influx_xy = grid.at_node["sediment__influx"].reshape(grid.shape)
    influx_xy[3, 0] = 0.5  # only row 3 gets sediment supply

    esd = GeoEnthalpyDelta(
        grid, topset_threshold=(0.25, 0.1), foreset_threshold=(2.0, 2.0)
    )
    esd.run_one_step()

    thickness_xy = grid.at_node["sediment__thickness"].reshape(grid.shape)
    assert thickness_xy[3, 0] > 0.0
    other_rows = np.arange(grid.shape[0]) != 3
    assert np.all(thickness_xy[other_rows, 0] == 0.0)


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
    grid = make_grid(nrows=20, ncols=40, dx=0.2, dy=0.2, sea_level=-5.0)
    set_west_boundary_influx(grid, total_flux=0.5, width=1.0)
    esd = GeoEnthalpyDelta(
        grid,
        topset_threshold=(0.25, 0.1),
        foreset_threshold=(2.0, 2.0),
    )
    nsteps = 200
    for _ in range(nsteps):
        esd.run_one_step()

    modeled_volume = np.sum(grid.at_node["sediment__thickness"]) * grid.dx * grid.dy
    expected_volume = np.sum(grid.at_node["sediment__influx"]) * esd.time_elapsed
    assert modeled_volume == pytest.approx(expected_volume, rel=1e-6)


def test_thickness_stays_nonnegative():
    grid = make_grid(nrows=20, ncols=40, dx=0.2, dy=0.2, sea_level=-5.0)
    set_west_boundary_influx(grid, total_flux=0.5, width=1.0)
    esd = GeoEnthalpyDelta(
        grid,
        topset_threshold=(0.25, 0.1),
        foreset_threshold=(2.0, 2.0),
    )
    for _ in range(200):
        esd.run_one_step()
        assert np.all(grid.at_node["sediment__thickness"] >= 0.0)


def test_sea_level_is_read_from_grid_field():
    grid = make_grid(sea_level=-5.0)
    esd = GeoEnthalpyDelta(grid)
    assert esd.sea_level == -5.0

    grid.at_grid["sea_level__elevation"] = -3.0
    assert esd.sea_level == -3.0


def test_sea_level_setter_updates_grid_field():
    grid = make_grid(sea_level=-5.0)
    esd = GeoEnthalpyDelta(grid)

    esd.sea_level = -4.0
    assert esd.sea_level == -4.0
    assert grid.at_grid["sea_level__elevation"] == -4.0


def test_sea_level_can_be_varied_externally_between_steps():
    grid = make_grid(sea_level=-5.0)
    esd = GeoEnthalpyDelta(grid)
    for _ in range(10):
        esd.run_one_step()
        esd.sea_level += 0.1
    assert esd.sea_level == pytest.approx(-5.0 + 0.1 * 10)


def test_run_one_step_dt_defaults_to_stable_time_step():
    grid = make_grid()
    esd = GeoEnthalpyDelta(grid)
    esd.run_one_step()
    assert esd.time_elapsed == pytest.approx(esd._calc_stable_time_step())


def test_run_one_step_subcycles_dt_larger_than_stable_limit():
    n_substeps = 5

    grid_direct = make_grid()
    set_west_boundary_influx(grid_direct, total_flux=0.5, width=1.0)
    esd_direct = GeoEnthalpyDelta(grid_direct)
    stable_dt = esd_direct._calc_stable_time_step()
    esd_direct.run_one_step(dt=n_substeps * stable_dt)

    grid_manual = make_grid()
    set_west_boundary_influx(grid_manual, total_flux=0.5, width=1.0)
    esd_manual = GeoEnthalpyDelta(grid_manual)
    for _ in range(n_substeps):
        esd_manual.run_one_step()

    assert esd_direct.time_elapsed == pytest.approx(esd_manual.time_elapsed)
    np.testing.assert_allclose(
        grid_direct.at_node["sediment__thickness"],
        grid_manual.at_node["sediment__thickness"],
        atol=1e-12,
    )
