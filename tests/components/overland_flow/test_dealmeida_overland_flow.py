"""
Unit tests for landlab.components.overland_flow.OverlandFlow

last updated: 3/14/16
"""

import numpy as np
import pytest

from landlab import HexModelGrid
from landlab import RasterModelGrid
from landlab.components.overland_flow import OverlandFlow
from landlab.graph.structured_quad.structured_quad import StructuredQuadGraphTopology
from landlab.grid.nodestatus import NodeStatus

(_SHAPE, _SPACING, _ORIGIN) = ((32, 240), (25, 25), (0.0, 0.0))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def _left_edge_horizontal_ids(shape):
    layout = StructuredQuadGraphTopology(shape)
    return layout.horizontal_links.reshape((shape[0], shape[1] - 1))[:, 0]


@pytest.fixture
def grid():
    grid = RasterModelGrid((3, 4), xy_spacing=1.0)
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("topographic__elevation", at="node")

    return grid


def test_not_a_raster_grid():
    grid = HexModelGrid((4, 5))
    with pytest.raises(NotImplementedError):
        OverlandFlow(grid)


def test_unequal_spacing():
    grid = RasterModelGrid((4, 5), xy_spacing=(3.0, 4.0))
    with pytest.raises(ValueError):
        OverlandFlow(grid)


@pytest.mark.parametrize(
    "key,value",
    (
        ("g", 0.0),
        ("g", -1.0),
        ("theta", -1e-6),
        ("theta", 1.00001),
    ),
)
def test_invalid_keyword_values(grid, key, value):
    with pytest.raises(ValueError):
        OverlandFlow(grid, **{key: value})


def test_mannings_n_as_a_field():
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
    grid.add_empty("foo", at="link")
    grid.at_link["foo"].fill(0.03)

    actual = OverlandFlow(grid, mannings_n="foo", steep_slopes=True)
    actual.run_one_step()

    np.testing.assert_almost_equal(actual.h, expected.h)


def test_deAlm_name(deAlm):
    assert deAlm.name == "OverlandFlow"


def test_deAlm_input_var_names(deAlm):
    assert deAlm.input_var_names == ("surface_water__depth", "topographic__elevation")


def test_deAlm_output_var_names(deAlm):
    assert deAlm.output_var_names == (
        "surface_water__depth",
        "surface_water__discharge",
        "water_surface__gradient",
    )


def test_deAlm_var_units(deAlm):
    assert set(deAlm.input_var_names) | set(deAlm.output_var_names) == set(
        dict(deAlm.units).keys()
    )

    assert deAlm.var_units("surface_water__depth") == "m"
    assert deAlm.var_units("surface_water__discharge") == "m3/s"
    assert deAlm.var_units("water_surface__gradient") == "-"
    assert deAlm.var_units("topographic__elevation") == "m"


def test_grid_shape(deAlm):
    assert deAlm.grid.number_of_node_rows == _SHAPE[0]
    assert deAlm.grid.number_of_node_columns == _SHAPE[1]


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
        deAlm.overland_flow(dt)
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
        deAlm.overland_flow(dt)
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
        deAlm.run_one_step(dt)
        h_boundary = ((7.0 / 3.0) * (0.01**2) * (0.4**3) * time) ** (3.0 / 7.0)
        grid.at_node["surface_water__depth"][grid.nodes[1:-1, 1]] = h_boundary
        time += dt

        assert np.all(grid.at_node["surface_water__depth"] >= 0.0)
        assert not np.any(np.isnan(grid.at_node["surface_water__depth"]))

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
    overland_flow.run_one_step(900.0)

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

    deAlm1.run_one_step(100)
    deAlm2.run_one_step(100)
    np.testing.assert_equal(deAlm1.h, deAlm2.h)


def test_rainfall_intensity_setter(grid):
    overland_flow = OverlandFlow(grid)

    assert overland_flow.rainfall_intensity == pytest.approx(0.0)

    overland_flow.rainfall_intensity = 12.0
    np.testing.assert_almost_equal(
        overland_flow.rainfall_intensity, grid.ones(at="node") * 12.0
    )

    overland_flow.rainfall_intensity = np.arange(grid.number_of_nodes)
    np.testing.assert_almost_equal(
        overland_flow.rainfall_intensity, np.arange(grid.number_of_nodes)
    )

    overland_flow.rainfall_intensity = 0.0
    np.testing.assert_almost_equal(
        overland_flow.rainfall_intensity, grid.zeros(at="node")
    )


def test_bad_rainfall_intensity(grid):
    overland_flow = OverlandFlow(grid)

    with pytest.raises(ValueError):
        overland_flow.rainfall_intensity = -1.0

    rainfall_intensity = grid.ones(at="node")
    rainfall_intensity[-1] = -1e-6

    with pytest.raises(ValueError):
        overland_flow.rainfall_intensity = rainfall_intensity


def test_discharge_mapper():
    grid = RasterModelGrid((3, 4), xy_spacing=1.0)
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("topographic__elevation", at="node")

    overland_flow = OverlandFlow(grid)

    discharge_at_link = grid.ones(at="link", dtype=float)
    discharge_at_node = overland_flow.discharge_mapper(discharge_at_link)

    np.testing.assert_almost_equal(
        discharge_at_node,
        [
            *[0.0, 1.0, 1.0, 1.0],
            *[1.0, 2.0, 2.0, 2.0],
            *[1.0, 2.0, 2.0, 2.0],
        ],
    )


def test_discharge_mapper_out_keyword(grid):
    discharge_at_link = grid.ones(at="link", dtype=float)

    overland_flow = OverlandFlow(grid)
    discharge_at_node = grid.empty(at="node")
    out = overland_flow.discharge_mapper(discharge_at_link, out=discharge_at_node)

    assert out is discharge_at_node
    np.testing.assert_almost_equal(
        discharge_at_node, overland_flow.discharge_mapper(discharge_at_link)
    )


@pytest.mark.parametrize("spacing", (1.0, 2.0))
def test_discharge_mapper_convert_to_volume(spacing):
    grid = RasterModelGrid((3, 4), xy_spacing=spacing)
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("topographic__elevation", at="node")
    discharge_at_link = grid.ones(at="link", dtype=float)

    overland_flow = OverlandFlow(grid)
    discharge_at_node = overland_flow.discharge_mapper(
        discharge_at_link, convert_to_volume=True
    )

    np.testing.assert_almost_equal(
        discharge_at_node,
        overland_flow.discharge_mapper(discharge_at_link, convert_to_volume=False)
        * spacing,
    )


def test_no_water(grid):
    overland_flow = OverlandFlow(grid)
    overland_flow.run_one_step(1.0)

    assert np.allclose(overland_flow.h, 0.0)


def test_no_water_no_time_step(grid):
    overland_flow = OverlandFlow(grid)
    with pytest.raises(ValueError):
        overland_flow.run_one_step()


def test_update_active_node(grid):
    overland_flow = OverlandFlow(grid)

    expected = grid.ones(at="node", dtype=bool)
    expected[grid.nodes_at_corners_of_grid,] = False

    active_nodes = overland_flow.update_active_nodes().copy()
    assert np.all(active_nodes == expected)

    grid.status_at_node[grid.shape[1] + 1] = NodeStatus.CLOSED

    assert np.all(overland_flow.update_active_nodes(clear_cache=False) == active_nodes)
    assert np.any(overland_flow.update_active_nodes(clear_cache=True) != active_nodes)
