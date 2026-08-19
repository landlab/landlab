"""
Unit tests for the GeoEnthalpyDelta component.

These tests exercise the fixed-grid, enthalpy-style finite-volume update
described in Lorenzo-Trueba et al., *GeoEnthalpy-Delta v1.0*: sediment
mass conservation, the positivity-preserving thickness update, the
sediment__influx feeder pattern, the internally re-derived basement, and
the CFL-limited stable time step.
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
    grid = RasterModelGrid((nrows, ncols), xy_spacing=(dx, dy))
    grid.status_at_node[:] = grid.BC_NODE_IS_CORE
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


def test_closed_boundary_nodes_raise():
    grid = make_grid()
    grid.status_at_node[0] = grid.BC_NODE_IS_CLOSED
    esd = GeoEnthalpyDelta(grid)
    with pytest.raises(ValueError):
        esd.run_one_step()


def test_default_fixed_value_boundary_raises():
    # RasterModelGrid's default perimeter status is FIXED_VALUE, not
    # CORE, and every node must be core for run_one_step to proceed.
    grid = RasterModelGrid((20, 30), xy_spacing=0.5)
    topo = grid.add_zeros("topographic__elevation", at="node")
    topo[:] = -grid.x_of_node
    grid.add_zeros("sediment__influx", at="node")
    grid.add_field("sea_level__elevation", 0.0, at="grid")
    assert np.any(grid.status_at_node == grid.BC_NODE_IS_FIXED_VALUE)

    esd = GeoEnthalpyDelta(grid)
    with pytest.raises(ValueError):
        esd.run_one_step()


def test_all_core_boundary_is_accepted():
    grid = make_grid()
    assert np.all(grid.status_at_node == grid.BC_NODE_IS_CORE)
    esd = GeoEnthalpyDelta(grid)
    esd.run_one_step()  # should not raise


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

    assert "sediment__thickness" not in grid.at_node
    assert_array_equal(esd.sediment_thickness, 0.0)
    assert esd.time_elapsed == 0.0
    assert esd.sea_level == grid.at_grid["sea_level__elevation"]


def test_basement_is_computed_from_initial_topo_and_thickness():
    grid = make_grid()
    grid.at_node["topographic__elevation"][:] += 0.3

    esd = GeoEnthalpyDelta(grid, sediment_thickness=0.3)

    expected_basement = grid.at_node["topographic__elevation"] - 0.3
    esd.run_one_step()
    basement_after = grid.at_node["topographic__elevation"] - esd.sediment_thickness
    np.testing.assert_allclose(basement_after, expected_basement, atol=1e-12)


def test_basement_stays_invariant_across_steps():
    grid = make_grid(nrows=20, ncols=40, dx=0.2, dy=0.2, sea_level=-5.0)
    set_west_boundary_influx(grid, total_flux=0.5, width=1.0)
    esd = GeoEnthalpyDelta(
        grid, topset_threshold=(0.25, 0.1), foreset_threshold=(2.0, 2.0)
    )
    basement_before = (
        grid.at_node["topographic__elevation"] - esd.sediment_thickness
    ).copy()

    for _ in range(20):
        esd.run_one_step()

    np.testing.assert_allclose(
        grid.at_node["topographic__elevation"] - esd.sediment_thickness,
        basement_before,
        atol=1e-12,
    )


def test_external_elevation_change_is_absorbed_into_basement():
    grid = make_grid(nrows=20, ncols=40, dx=0.2, dy=0.2, sea_level=-5.0)
    set_west_boundary_influx(grid, total_flux=0.5, width=1.0)
    esd = GeoEnthalpyDelta(
        grid, topset_threshold=(0.25, 0.1), foreset_threshold=(2.0, 2.0)
    )
    esd.run_one_step()
    thickness_before = esd.sediment_thickness.copy()

    # Simulate an external process (e.g. GIA) shifting the elevation.
    grid.at_node["topographic__elevation"] += 1.5

    esd.run_one_step()

    # Our tracked thickness keeps evolving from where it was; the shift
    # was absorbed into basement, not thickness.
    assert not np.allclose(esd.sediment_thickness, thickness_before)
    assert np.all(esd.sediment_thickness >= 0.0)


def test_influx_only_enters_at_specified_west_boundary_nodes():
    grid = make_grid(nrows=10, ncols=10, dx=0.2, dy=0.2, sea_level=-5.0)
    influx_xy = grid.at_node["sediment__influx"].reshape(grid.shape)
    influx_xy[3, 0] = 0.5  # only row 3 gets sediment supply

    esd = GeoEnthalpyDelta(
        grid, topset_threshold=(0.25, 0.1), foreset_threshold=(2.0, 2.0)
    )
    esd.run_one_step()

    thickness_xy = esd.sediment_thickness.reshape(grid.shape)
    assert thickness_xy[3, 0] > 0.0
    other_rows = np.arange(grid.shape[0]) != 3
    assert np.all(thickness_xy[other_rows, 0] == 0.0)


def test_calc_lateral_flux_below_threshold_is_zero():
    grid = make_grid(dx=1.0, dy=1.0, sea_level=100.0)
    esd = GeoEnthalpyDelta(grid, foreset_threshold=0.5, foreset_diffusivity=3.0)
    eta = np.array([[0.2], [0.0]])  # |slope| = 0.2 < threshold 0.5
    qy = esd._calc_lateral_flux(eta)
    assert_array_equal(qy, 0.0)


def test_calc_lateral_flux_above_threshold_matches_formula():
    grid = make_grid(dx=1.0, dy=1.0, sea_level=100.0)
    esd = GeoEnthalpyDelta(grid, foreset_threshold=0.5, foreset_diffusivity=3.0)
    eta = np.array([[2.0], [0.0]])  # slope = 2.0, above threshold
    qy = esd._calc_lateral_flux(eta)
    assert qy[0, 0] == 0.0  # outer y-faces carry no flux
    assert qy[2, 0] == 0.0
    assert qy[1, 0] == pytest.approx(3.0 * (2.0 - 0.5))


def test_limit_lateral_flux_caps_outflow_to_available_thickness():
    grid = make_grid(dx=1.0, dy=1.0)
    esd = GeoEnthalpyDelta(grid)
    thickness = np.array([[0.2], [5.0]])  # node 0 can only supply 0.2 in dt=1.0
    qy = np.array([[0.0], [5.0], [0.0]])  # candidate flux leaving node 0 is 5.0
    limited = esd._limit_lateral_flux(qy, thickness, dt=1.0)
    assert limited[1, 0] == pytest.approx(0.2)


def test_calc_downstream_flux_limited_by_diffusive_capacity():
    grid = make_grid(dx=1.0, dy=1.0, sea_level=100.0)
    esd = GeoEnthalpyDelta(grid, foreset_threshold=0.5, foreset_diffusivity=2.0)
    eta = np.array([[3.0, 0.0]])
    thickness = np.array([[10.0, 10.0]])  # plenty of sediment, not the constraint
    qy = np.zeros((2, 2))
    qx = esd._calc_downstream_flux(eta, thickness, qy, dt=1.0)
    assert qx[0, 0] == 0.0  # no flux enters across the west domain edge
    assert qx[0, 1] == pytest.approx(2.0 * (3.0 - 0.5))  # diffusive capacity
    assert qx[0, 2] == 0.0  # east domain boundary


def test_advance_substep_mixes_influx_into_thickness_at_any_node():
    grid = make_grid(dx=1.0, dy=1.0, sea_level=100.0)
    esd = GeoEnthalpyDelta(grid, foreset_threshold=100.0, foreset_diffusivity=1.0)
    basement = np.zeros((3, 3))
    thickness = np.zeros((3, 3))
    influx = np.zeros((3, 3))
    influx[1, 1] = 2.0  # interior node, not on the west boundary

    new_thickness = esd._advance_substep(basement, thickness, influx, dt=1.0)

    # Threshold is high enough that nothing transports away this step, so
    # the node should hold exactly what was injected.
    expected = influx[1, 1] * 1.0 / (grid.dx * grid.dy)
    assert new_thickness[1, 1] == pytest.approx(expected)
    assert new_thickness[0, 0] == 0.0


def test_apply_conservative_update_matches_flux_divergence():
    grid = make_grid(dx=1.0, dy=1.0)
    esd = GeoEnthalpyDelta(grid)
    thickness = np.array([[1.0]])
    qx = np.array([[2.0, 0.5]])
    qy = np.array([[0.3], [0.1]])
    new_thickness = esd._apply_conservative_update(thickness, qx, qy, dt=1.0)
    assert new_thickness[0, 0] == pytest.approx(1.0 + (2.0 - 0.5) + (0.3 - 0.1))


def test_apply_conservative_update_clips_at_zero():
    grid = make_grid(dx=1.0, dy=1.0)
    esd = GeoEnthalpyDelta(grid)
    thickness = np.array([[0.1]])
    qx = np.array([[0.0, 5.0]])  # large outflow, more than available
    qy = np.array([[0.0], [0.0]])
    new_thickness = esd._apply_conservative_update(thickness, qx, qy, dt=1.0)
    assert new_thickness[0, 0] == 0.0


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


def test_calc_stable_time_step_matches_cfl_formula_with_dx_ne_dy():
    grid = make_grid(dx=0.2, dy=0.4)  # non-even spacing
    assert grid.dx != grid.dy
    esd = GeoEnthalpyDelta(
        grid,
        topset_diffusivity=(1.0, 1.0),
        foreset_diffusivity=(3.0, 2.0),
        cfl=0.4,
    )
    d_max = 3.0
    expected = 0.4 / (2.0 * d_max * (1.0 / 0.2**2 + 1.0 / 0.4**2))
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

    modeled_volume = np.sum(esd.sediment_thickness) * grid.dx * grid.dy
    expected_volume = np.sum(grid.at_node["sediment__influx"]) * esd.time_elapsed
    assert modeled_volume == pytest.approx(expected_volume, rel=1e-6)


def test_mass_is_conserved_with_dx_ne_dy():
    grid = make_grid(nrows=20, ncols=40, dx=0.2, dy=0.4, sea_level=-5.0)
    assert grid.dx != grid.dy
    set_west_boundary_influx(grid, total_flux=0.5, width=1.0)
    esd = GeoEnthalpyDelta(
        grid,
        topset_threshold=(0.25, 0.1),
        foreset_threshold=(2.0, 2.0),
    )
    nsteps = 200
    for _ in range(nsteps):
        esd.run_one_step()

    modeled_volume = np.sum(esd.sediment_thickness) * grid.dx * grid.dy
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
        assert np.all(esd.sediment_thickness >= 0.0)


def test_sea_level_is_read_from_grid_field():
    grid = make_grid(sea_level=-5.0)
    esd = GeoEnthalpyDelta(grid)
    assert esd.sea_level == -5.0

    grid.at_grid["sea_level__elevation"] = -3.0
    assert esd.sea_level == -3.0


def test_sea_level_has_no_setter():
    grid = make_grid(sea_level=-5.0)
    esd = GeoEnthalpyDelta(grid)
    with pytest.raises(AttributeError):
        esd.sea_level = -4.0


def test_sea_level_can_be_varied_externally_between_steps():
    grid = make_grid(sea_level=-5.0)
    esd = GeoEnthalpyDelta(grid)
    for _ in range(10):
        esd.run_one_step()
        grid.at_grid["sea_level__elevation"] += 0.1
    assert esd.sea_level == pytest.approx(-5.0 + 0.1 * 10)


def test_run_one_step_dt_defaults_to_stable_time_step():
    grid = make_grid()
    esd = GeoEnthalpyDelta(grid)
    esd.run_one_step()
    assert esd.time_elapsed == pytest.approx(esd._calc_stable_time_step())


def test_run_one_step_dt_must_be_positive_and_finite():
    grid = make_grid()
    esd = GeoEnthalpyDelta(grid)
    with pytest.raises(ValidationError):
        esd.run_one_step(dt=0.0)
    with pytest.raises(ValidationError):
        esd.run_one_step(dt=-1.0)
    with pytest.raises(ValidationError):
        esd.run_one_step(dt=np.inf)


def test_run_one_step_influx_must_be_nonnegative():
    grid = make_grid()
    esd = GeoEnthalpyDelta(grid)
    grid.at_node["sediment__influx"][0] = -1.0
    with pytest.raises(ValidationError):
        esd.run_one_step()


def test_run_one_step_influx_must_be_finite():
    grid = make_grid()
    esd = GeoEnthalpyDelta(grid)
    grid.at_node["sediment__influx"][0] = np.inf
    with pytest.raises(ValidationError):
        esd.run_one_step()
    grid.at_node["sediment__influx"][0] = np.nan
    with pytest.raises(ValidationError):
        esd.run_one_step()


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
        esd_direct.sediment_thickness,
        esd_manual.sediment_thickness,
        atol=1e-12,
    )
