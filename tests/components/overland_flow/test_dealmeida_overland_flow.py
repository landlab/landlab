"""
Unit tests for landlab.components.overland_flow.OverlandFlow

last updated: 3/14/16
"""

import numpy as np
import pytest
from numpy.testing import assert_almost_equal
from numpy.testing import assert_equal

from landlab import RasterModelGrid
from landlab.components.overland_flow import OverlandFlow
from landlab.components.overland_flow._links import horizontal_link_ids
from landlab.components.overland_flow.generate_overland_flow_deAlmeida import (
    NoWaterError,
)
from landlab.core.errors import ValidationError
from landlab.grid.nodestatus import NodeStatus

(_SHAPE, _SPACING, _ORIGIN) = ((32, 240), (25, 25), (0.0, 0.0))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


@pytest.fixture
def deAlm():
    grid = RasterModelGrid((32, 240), xy_spacing=25)
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("topographic__elevation", at="node")
    grid.add_zeros("surface_water__discharge", at="link")
    return OverlandFlow(grid, mannings_n=0.01, h_init=0.001)


def _left_edge_horizontal_ids(shape):
    return horizontal_link_ids(shape)[:, 0]


@pytest.mark.parametrize(
    "kwds",
    (
        {"h_init": -1.0},
        {"alpha": -1.0},
        {"alpha": 0.0},
        {"g": 0.0},
        {"g": -1.0},
        {"theta": -0.1},
        {"theta": 1.1},
        {"rainfall_intensity": -1.0},
    ),
)
def test_dealmeida_invalid_parameter(kwds):
    grid = RasterModelGrid((4, 5))
    grid.add_empty("surface_water__depth", at="node")
    grid.add_empty("topographic__elevation", at="node")
    grid.add_empty("surface_water__discharge", at="link")
    with pytest.raises(ValidationError):
        OverlandFlow(grid, **kwds)


def test_update_active_nodes_cached():
    grid = RasterModelGrid((4, 5))
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("topographic__elevation", at="node")

    overland_flow = OverlandFlow(grid, h_init=0.0)

    active_nodes = overland_flow._update_active_nodes()
    assert_equal(
        active_nodes,
        [
            *[0, 1, 1, 1, 0],
            *[1, 1, 1, 1, 1],
            *[1, 1, 1, 1, 1],
            *[0, 1, 1, 1, 0],
        ],
    )
    assert overland_flow._update_active_nodes() is active_nodes


def test_update_active_nodes_clear_cache():
    grid = RasterModelGrid((4, 5))
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("topographic__elevation", at="node")
    grid.status_at_node[:] = NodeStatus.CORE

    overland_flow = OverlandFlow(grid, h_init=0.0)

    old_values = overland_flow._update_active_nodes()
    new_values = overland_flow._update_active_nodes(clear_cache=True)
    assert new_values is not old_values
    assert np.all(new_values == old_values)


def test_calc_time_step_no_water():
    grid = RasterModelGrid((4, 5))
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("topographic__elevation", at="node")

    overland_flow = OverlandFlow(grid, h_init=0.0)
    with pytest.raises(NoWaterError, match="no water on landscape"):
        overland_flow.calc_time_step()

    with pytest.raises(ValueError, match="no water on landscape and dt not provided"):
        overland_flow.run_one_step(dt=None)

    assert overland_flow.run_one_step(1.0) == 1.0

    grid.status_at_node[:] = NodeStatus.CLOSED
    grid.at_node["surface_water__depth"][:] = 1.0

    overland_flow._update_active_nodes(clear_cache=True)
    with pytest.raises(NoWaterError, match="no active links"):
        overland_flow.calc_time_step()


def test_deAlm_name(deAlm):
    assert deAlm.name == "OverlandFlow"


def test_deAlm_input_var_names(deAlm):
    assert deAlm.input_var_names == ("surface_water__depth", "topographic__elevation")


def test_deAlm_output_var_names(deAlm):
    assert deAlm.output_var_names == ("surface_water__depth",)


def test_deAlm_var_units(deAlm):
    assert set(deAlm.input_var_names) | set(deAlm.output_var_names) | set(
        deAlm.optional_var_names
    ) == set(dict(deAlm.units).keys())

    assert deAlm.var_units("surface_water__depth") == "m"
    assert deAlm.var_units("surface_water__discharge") == "m3/s"
    assert deAlm.var_units("water_surface__gradient") == "-"
    assert deAlm.var_units("topographic__elevation") == "m"


def test_grid_shape(deAlm):
    assert deAlm.grid.number_of_node_rows == _SHAPE[0]
    assert deAlm.grid.number_of_node_columns == _SHAPE[1]


@pytest.mark.parametrize("name", ("surface_water__depth", "water_surface__gradient"))
def test_with_optional_at_link_field(name):
    grid = RasterModelGrid((8, 32), xy_spacing=25.0)

    z_at_node = grid.zeros(at="node")
    h_at_node = -0.3 * grid.x_of_node

    grid.at_node["surface_water__depth"] = z_at_node.copy()
    grid.at_node["topographic__elevation"] = h_at_node.copy()
    grid.set_closed_boundaries_at_grid_edges(True, True, True, True)

    OverlandFlow(
        grid, mannings_n=0.01, h_init=1e-9, rainfall_intensity=0.001, steep_slopes=True
    ).run_one_step(900.0)
    expected = grid.at_node["surface_water__depth"]

    grid = RasterModelGrid((8, 32), xy_spacing=25.0)

    grid.at_node["surface_water__depth"] = z_at_node.copy()
    grid.at_node["topographic__elevation"] = h_at_node.copy()
    grid.set_closed_boundaries_at_grid_edges(True, True, True, True)

    grid.add_empty(name, at="link")
    OverlandFlow(
        grid, mannings_n=0.01, h_init=1e-9, rainfall_intensity=0.001, steep_slopes=True
    ).run_one_step(900.0)
    actual = grid.at_node["surface_water__depth"]

    assert_almost_equal(actual, expected)


def test_deAlm_analytical():
    grid = RasterModelGrid((32, 240), xy_spacing=25)
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("topographic__elevation", at="node")
    grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
    left_inactive_ids = _left_edge_horizontal_ids(grid.shape)
    deAlm = OverlandFlow(grid, mannings_n=0.01, h_init=0.001)
    time = 0.0

    while time < 500.0:
        grid.at_link["surface_water__discharge"][left_inactive_ids] = grid.at_link[
            "surface_water__discharge"
        ][left_inactive_ids + 1]
        dt = deAlm.calc_time_step()
        assert deAlm.overland_flow(dt) == dt
        h_boundary = ((7.0 / 3.0) * (0.01**2) * (0.4**3) * time) ** (3.0 / 7.0)
        grid.at_node["surface_water__depth"][grid.nodes[1:-1, 1]] = h_boundary
        time += dt

    x = np.arange(0, ((grid.shape[1]) * grid.dx), grid.dx)
    h_analytical = -(7.0 / 3.0) * (0.01**2) * (0.4**2) * (x - (0.4 * 500))

    h_analytical[np.where(h_analytical > 0)] = h_analytical[
        np.where(h_analytical > 0)
    ] ** (3.0 / 7.0)
    h_analytical[np.where(h_analytical < 0)] = 0.0

    hdeAlm = deAlm.h.reshape(grid.shape)
    hdeAlm = hdeAlm[1][1:]
    hdeAlm = np.append(hdeAlm, [0])
    np.testing.assert_almost_equal(h_analytical, hdeAlm, decimal=1)


def test_deAlm_analytical_imposed_dt_short():
    grid = RasterModelGrid((32, 240), xy_spacing=25)
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("topographic__elevation", at="node")
    grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
    left_inactive_ids = _left_edge_horizontal_ids(grid.shape)
    deAlm = OverlandFlow(grid, mannings_n=0.01, h_init=0.001)
    time = 0.0

    while time < 500.0:
        grid.at_link["surface_water__discharge"][left_inactive_ids] = grid.at_link[
            "surface_water__discharge"
        ][left_inactive_ids + 1]
        dt = 10.0
        assert deAlm.overland_flow(dt) == dt
        h_boundary = ((7.0 / 3.0) * (0.01**2) * (0.4**3) * time) ** (3.0 / 7.0)
        grid.at_node["surface_water__depth"][grid.nodes[1:-1, 1]] = h_boundary
        time += dt

    x = np.arange(0, ((grid.shape[1]) * grid.dx), grid.dx)
    h_analytical = -(7.0 / 3.0) * (0.01**2) * (0.4**2) * (x - (0.4 * 500))

    h_analytical[np.where(h_analytical > 0)] = h_analytical[
        np.where(h_analytical > 0)
    ] ** (3.0 / 7.0)
    h_analytical[np.where(h_analytical < 0)] = 0.0

    hdeAlm = deAlm.h.reshape(grid.shape)
    hdeAlm = hdeAlm[1][1:]
    hdeAlm = np.append(hdeAlm, [0])
    np.testing.assert_almost_equal(h_analytical, hdeAlm, decimal=1)


def test_deAlm_analytical_imposed_dt_long():
    grid = RasterModelGrid((32, 240), xy_spacing=25)
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("topographic__elevation", at="node")
    grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
    left_inactive_ids = _left_edge_horizontal_ids(grid.shape)
    deAlm = OverlandFlow(grid, mannings_n=0.01, h_init=0.001)
    time = 0.0

    while time < 500.0:
        grid.at_link["surface_water__discharge"][left_inactive_ids] = grid.at_link[
            "surface_water__discharge"
        ][left_inactive_ids + 1]
        dt = 100.0
        elapsed = deAlm.run_one_step(dt)

        assert elapsed == dt
        assert not np.any(np.isnan(grid.at_node["surface_water__depth"]))
        assert np.all(grid.at_node["surface_water__depth"] >= 0.0)

        h_boundary = ((7.0 / 3.0) * (0.01**2) * (0.4**3) * time) ** (3.0 / 7.0)
        grid.at_node["surface_water__depth"][grid.nodes[1:-1, 1]] = h_boundary
        time += dt

    x = np.arange(0, ((grid.shape[1]) * grid.dx), grid.dx)
    h_analytical = -(7.0 / 3.0) * (0.01**2) * (0.4**2) * (x - (0.4 * 500))

    h_analytical[np.where(h_analytical > 0)] = h_analytical[
        np.where(h_analytical > 0)
    ] ** (3.0 / 7.0)
    h_analytical[np.where(h_analytical < 0)] = 0.0

    hdeAlm = deAlm.h.reshape(grid.shape)
    hdeAlm = hdeAlm[1][1:]
    hdeAlm = np.append(hdeAlm, [0])
    np.testing.assert_almost_equal(h_analytical, hdeAlm, decimal=1)


@pytest.mark.parametrize("slope", (0.0, 0.1, 0.2, 0.3, 0.4, 0.5))
def test_mass_balance(slope):
    grid = RasterModelGrid((8, 32), xy_spacing=25.0)

    grid.add_zeros("surface_water__depth", at="node")
    grid.add_field("topographic__elevation", -slope * grid.x_of_node, at="node")
    grid.set_closed_boundaries_at_grid_edges(True, True, True, True)

    overland_flow = OverlandFlow(
        grid, mannings_n=0.01, h_init=1e-9, rainfall_intensity=0.001, steep_slopes=True
    )
    assert overland_flow.run_one_step(900.0) == 900.0

    expected = (900.0 * 0.001 + 1e-9) * len(grid.core_nodes)
    actual = grid.at_node["surface_water__depth"][grid.core_nodes].sum()

    assert actual == pytest.approx(expected)


def test_deAlm_rainfall_array():
    """
    Make sure that rainfall_intensity can be set with an array, and confirm
    that this returns the same result as setting with a float of the same magnitude.
    """

    mg1 = RasterModelGrid((10, 10), xy_spacing=25)
    mg1.add_zeros("surface_water__depth", at="node")
    mg1.add_zeros("topographic__elevation", at="node")
    mg1.set_closed_boundaries_at_grid_edges(True, True, True, True)
    r1 = 1e-6
    deAlm1 = OverlandFlow(mg1, mannings_n=0.01, h_init=0.001, rainfall_intensity=r1)

    mg2 = RasterModelGrid((10, 10), xy_spacing=25)
    mg2.add_zeros("surface_water__depth", at="node")
    mg2.add_zeros("topographic__elevation", at="node")
    mg2.set_closed_boundaries_at_grid_edges(True, True, True, True)
    r2 = 1e-6 * np.ones(100)
    deAlm2 = OverlandFlow(mg2, mannings_n=0.01, h_init=0.001, rainfall_intensity=r2)

    assert deAlm1.run_one_step(100) == 100
    assert deAlm2.run_one_step(100) == 100
    np.testing.assert_equal(deAlm1.h, deAlm2.h)


@pytest.mark.parametrize("mannings_n", ("foo", 0.03))
def test_mannings_n_as_a_field(mannings_n):
    initial_elevation = [
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0, 1.0, 1.0],
        [2.0, 2.0, 2.0, 2.0, 2.0],
        [3.0, 3.0, 3.0, 3.0, 3.0],
    ]
    initial_water_depth = [
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [0.1, 0.1, 0.1, 0.1, 0.1],
    ]

    grid = RasterModelGrid((4, 5))
    grid.at_node["topographic__elevation"] = initial_elevation
    grid.at_node["surface_water__depth"] = initial_water_depth

    expected = OverlandFlow(grid, mannings_n=0.03, steep_slopes=True)
    expected.run_one_step()

    grid = RasterModelGrid((4, 5))
    grid.at_node["topographic__elevation"] = initial_elevation
    grid.at_node["surface_water__depth"] = initial_water_depth

    if isinstance(mannings_n, str):
        grid.add_empty(mannings_n, at="link").fill(0.03)
    else:
        mannings_n = grid.empty(at="link")
        mannings_n.fill(0.03)

    actual = OverlandFlow(grid, mannings_n=mannings_n, steep_slopes=True)
    actual.run_one_step()

    assert_almost_equal(actual.h, expected.h)
