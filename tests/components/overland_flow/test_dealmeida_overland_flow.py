# -*- coding: utf-8 -*-
"""
Unit tests for landlab.components.overland_flow.OverlandFlow

last updated: 3/14/16
"""
import numpy as np

from landlab import RasterModelGrid
from landlab.components.overland_flow import OverlandFlow
from landlab.components.overland_flow._links import horizontal_link_ids

(_SHAPE, _SPACING, _ORIGIN) = ((32, 240), (25, 25), (0.0, 0.0))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def _left_edge_horizontal_ids(shape):
    return horizontal_link_ids(shape)[:, 0]


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
        h_boundary = ((7.0 / 3.0) * (0.01 ** 2) * (0.4 ** 3) * time) ** (3.0 / 7.0)
        grid.at_node["surface_water__depth"][grid.nodes[1:-1, 1]] = h_boundary
        time += dt

    x = np.arange(0, ((grid.shape[1]) * grid.dx), grid.dx)
    h_analytical = -(7.0 / 3.0) * (0.01 ** 2) * (0.4 ** 2) * (x - (0.4 * 500))

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
        h_boundary = ((7.0 / 3.0) * (0.01 ** 2) * (0.4 ** 3) * time) ** (3.0 / 7.0)
        grid.at_node["surface_water__depth"][grid.nodes[1:-1, 1]] = h_boundary
        time += dt

    x = np.arange(0, ((grid.shape[1]) * grid.dx), grid.dx)
    h_analytical = -(7.0 / 3.0) * (0.01 ** 2) * (0.4 ** 2) * (x - (0.4 * 500))

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
        h_boundary = ((7.0 / 3.0) * (0.01 ** 2) * (0.4 ** 3) * time) ** (3.0 / 7.0)
        grid.at_node["surface_water__depth"][grid.nodes[1:-1, 1]] = h_boundary
        time += dt

    x = np.arange(0, ((grid.shape[1]) * grid.dx), grid.dx)
    h_analytical = -(7.0 / 3.0) * (0.01 ** 2) * (0.4 ** 2) * (x - (0.4 * 500))

    h_analytical[np.where(h_analytical > 0)] = h_analytical[
        np.where(h_analytical > 0)
    ] ** (3.0 / 7.0)
    h_analytical[np.where(h_analytical < 0)] = 0.0

    hdeAlm = deAlm.h.reshape(grid.shape)
    hdeAlm = hdeAlm[1][1:]
    hdeAlm = np.append(hdeAlm, [0])
    np.testing.assert_almost_equal(h_analytical, hdeAlm, decimal=1)
