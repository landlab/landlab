#! /usr/bin/env python
"""
Unit tests for landlab.components.flexure.flexure
"""
import numpy as np
import pytest

from landlab import RasterModelGrid
from landlab.components import Flexure

(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e3, 10e3), (0.0, 0.0))


def test_method_names():
    grid = RasterModelGrid((20, 20), xy_spacing=10e3)
    assert Flexure(grid, method="airy").method == "airy"
    assert Flexure(grid, method="flexure").method == "flexure"
    with pytest.raises(ValueError):
        Flexure(grid, method="bad-name")


def test_eet_attribute():
    grid = RasterModelGrid((20, 20), xy_spacing=10e3)
    for val in (10e3, 1e3):
        assert Flexure(grid, eet=val).eet == pytest.approx(val)
    with pytest.raises(ValueError):
        assert Flexure(grid, eet=-10e3)


def test_youngs_attribute():
    grid = RasterModelGrid((20, 20), xy_spacing=10e3)
    for val in (10e3, 1e3):
        assert Flexure(grid, youngs=val).youngs == pytest.approx(val)


def test_gravity_attribute():
    grid = RasterModelGrid((20, 20), xy_spacing=10e3)
    for val in (10e3, 1e3):
        assert Flexure(grid, gravity=val).gravity == pytest.approx(val)


def test_rho_mantle_attribute():
    grid = RasterModelGrid((20, 20), xy_spacing=10e3)
    for val in (10e3, 1e3):
        assert Flexure(grid, rho_mantle=val).rho_mantle == pytest.approx(val)


def test_name(flex):
    assert flex.name == "Flexure"


def test_input_var_names(flex):
    assert flex.input_var_names == ("lithosphere__overlying_pressure_increment",)


def test_output_var_names(flex):
    assert flex.output_var_names == ("lithosphere_surface__elevation_increment",)


def test_var_units(flex):
    assert set(flex.input_var_names) | set(flex.output_var_names) == set(
        dict(flex.units).keys()
    )

    assert flex.var_units("lithosphere_surface__elevation_increment") == "m"
    assert flex.var_units("lithosphere__overlying_pressure_increment") == "Pa"


def test_grid_shape(flex):
    assert flex.grid.number_of_node_rows == _SHAPE[0]
    assert flex.grid.number_of_node_columns == _SHAPE[1]


def test_grid_x_extent(flex):
    assert flex.grid.extent[1] == (_SHAPE[1] - 1) * _SPACING[1]


def test_grid_y_extent(flex):
    assert flex.grid.extent[0] == (_SHAPE[0] - 1) * _SPACING[0]


def test_field_getters(flex):
    for name in flex.grid["node"]:
        field = flex.grid["node"][name]
        assert isinstance(field, np.ndarray)
        assert field.shape == (
            flex.grid.number_of_node_rows * flex.grid.number_of_node_columns,
        )

    with pytest.raises(KeyError):
        flex.grid["not_a_var_name"]


def test_field_initialized_to_zero(flex):
    for name in flex.grid["node"]:
        field = flex.grid["node"][name]
        assert np.all(field == 0.0)


def test_update():
    n = 11
    n_mid = (n - 1) // 2
    i_mid = np.ravel_multi_index((n_mid, n_mid), (n, n))
    load_0 = 1e9

    grid = RasterModelGrid((n, n), xy_spacing=1e3)
    flex = Flexure(grid, method="flexure")

    load = grid.at_node["lithosphere__overlying_pressure_increment"]
    load[i_mid] = load_0

    flex.update()
    dz = flex.grid.at_node["lithosphere_surface__elevation_increment"].reshape((n, n))

    assert np.argmax(dz) == i_mid
    assert dz[n_mid, n_mid] > 0.0
    assert np.all(dz[:, n_mid::-1] == dz[:, n_mid:])
    assert np.all(dz[n_mid::-1, :] == dz[n_mid:, :])


def test_subside_loads():
    n, load_0 = 11, 1e9

    grid = RasterModelGrid((n, n), xy_spacing=1e3)
    flex = Flexure(grid, method="flexure")

    grid.at_node["lithosphere__overlying_pressure_increment"][0] = load_0
    flex.update()
    dz_expected = flex.grid.at_node["lithosphere_surface__elevation_increment"]

    load = np.zeros((n, n))
    load[0, 0] = load_0

    dz = flex.subside_loads(load)
    assert np.all(dz.flatten() == pytest.approx(dz_expected))

    out = np.zeros((n, n))
    dz = flex.subside_loads(load, out=out)
    assert dz is out
