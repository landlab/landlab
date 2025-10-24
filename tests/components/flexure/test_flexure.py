#! /usr/bin/env python
"""
Unit tests for landlab.components.flexure.flexure
"""
import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.components import Flexure
from landlab.components.flexure._ext.flexure2d import subside_loads
from landlab.components.flexure.flexure import _find_nonzero_loads
from landlab.components.flexure.flexure import _validate_range


@pytest.mark.parametrize(
    "loads,expected",
    (
        ([[0, 0], [1.5, 0]], ([1.5], ([1], [0]))),
        ([[0, 0], [-1.5, 0]], ([-1.5], ([1], [0]))),
        ([[0, 0], [0, 0]], ([], ([], []))),
        ([[]], ([], ([], []))),
        ([[1, 2], [-1, -2]], ([1, 2, -1, -2], ([0, 0, 1, 1], [0, 1, 0, 1]))),
    ),
)
def test_find_nonzero_loads(loads, expected):
    loads, (rows, cols) = _find_nonzero_loads(loads)
    assert_array_equal(loads, expected[0])
    assert_array_equal(rows, expected[1][0])
    assert_array_equal(cols, expected[1][1])


@pytest.mark.parametrize("loads", ([], [0, 1, 0], [[[0, 1, 0]]]))
def test_find_nonzero_loads_wrong_ndim(loads):
    with pytest.raises(ValueError, match="'loads' must be 2D"):
        _find_nonzero_loads(loads)


@pytest.mark.parametrize("array", ([1, 2, 3, 4], (1, 2, 3, 4), np.array([1, 2, 3, 4])))
def test_validate_range_with_sequences(array):
    actual = _validate_range(array)
    assert np.allclose(actual, array)


def test_validate_range_out_is_in():
    array = np.arange(12)
    actual = _validate_range(array)
    assert np.shares_memory(actual, array)
    assert np.allclose(actual, array)


@pytest.mark.parametrize("low,high", ((2, None), (None, 4), (-1, 2), (3, 9)))
def test_validate_range_out_of_range(low, high):
    with pytest.raises(ValueError, match="array has values"):
        _validate_range([1, 2, 3, 4], low=low, high=high)


@pytest.mark.parametrize("low,high", ((None, None), (0, 9)))
@pytest.mark.parametrize("array", ([1, 2, 3, 4], [], [1], [1, 1, 1, 1], [0, 8.999999]))
def test_validate_range_in_range(array, low, high):
    actual = _validate_range(array, low=low, high=high)
    assert np.allclose(actual, array)


@pytest.mark.parametrize("n", [4, 5, 6, 7, 8, 9, 10])
@pytest.mark.parametrize("method", ["old", "new", "cython"])
def test_one_load_bench(n, method):
    load_0 = 1e9

    size = 2**n + 1

    grid = RasterModelGrid((size, size), xy_spacing=10.0)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    flex = Flexure(grid, method="flexure")

    loads = np.array([load_0])
    row_col_of_load = np.array([2 ** (n - 1)]), np.array([2 ** (n - 1)])

    w = grid.zeros(at="node").reshape(grid.shape)

    if method == "old":
        load_grid = grid.zeros(at="node").reshape(grid.shape)
        load_grid[row_col_of_load] = load_0

        func = flex.subside_loads_slow
        args = (load_grid,)
        kwds = {"out": w}
    elif method == "new":
        func = flex.subside_loads
        args = (loads, row_col_of_load)
        kwds = {"out": w}
    else:
        func = subside_loads
        args = (
            w,
            flex._r,
            loads,
            row_col_of_load[0],
            row_col_of_load[1],
            flex.alpha,
            flex.gamma_mantle,
        )
        kwds = {}

    func(*args, **kwds)

    assert_array_almost_equal(w, w[:, ::-1])
    assert_array_almost_equal(w, w[::-1, :])

    assert_array_equal(
        np.unravel_index(np.argmax(w, axis=None), w.shape), (2 ** (n - 1), 2 ** (n - 1))
    )


@pytest.mark.parametrize("method", ["subside_loads", "subside_loads_slow"])
def test_flexure_deflection_is_proportional_to_load(method):
    n, load_0 = 8, 1e9

    size = 2**n + 1

    grid = RasterModelGrid((size, size), xy_spacing=10.0)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    flex = Flexure(grid, method="flexure")

    loads = grid.zeros(at="node").reshape(grid.shape)
    loads[: size // 2, : size // 2] = load_0

    expected = np.zeros((size, size))
    getattr(flex, method)(loads, out=expected)

    actual = np.zeros((size, size))
    getattr(flex, method)(loads * 10.0, out=actual)

    assert_array_almost_equal(actual, expected * 10.0)


def test_update_bench():
    n, load_0 = 10, 1e9

    n_rows = 2**n + 1
    n_cols = 2**n + 1

    grid = RasterModelGrid((n_rows, n_cols), xy_spacing=10.0)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    flex = Flexure(grid, method="flexure")

    grid.at_node["lithosphere__overlying_pressure_increment"].reshape(grid.shape)[
        2 ** (n - 1), 2 ** (n - 1)
    ] = load_0
    flex.update()
    dz = flex.grid.at_node["lithosphere_surface__elevation_increment"].reshape(
        grid.shape
    )

    assert_array_almost_equal(dz, dz[:, ::-1])
    assert_array_almost_equal(dz, dz[::-1, :])


def test_flexure_cmp():
    n, load_0 = 10, 1e9

    n_rows = 2**n + 1
    n_cols = 2**n + 1

    grid = RasterModelGrid((n_rows, n_cols), xy_spacing=10.0)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    flex = Flexure(grid, method="flexure")

    grid.at_node["lithosphere__overlying_pressure_increment"].reshape(grid.shape)[
        2 ** (n - 1), 2 ** (n - 1)
    ] = load_0
    flex.update()
    expected = flex.grid.at_node["lithosphere_surface__elevation_increment"].reshape(
        grid.shape
    )

    grid = RasterModelGrid((n_rows, n_cols), xy_spacing=10.0)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    flex = Flexure(grid, method="flexure")

    loads = grid.zeros(at="node").reshape(grid.shape)
    loads[2 ** (n - 1), 2 ** (n - 1)] = load_0

    actual = np.zeros((n_rows, n_cols))
    flex.subside_loads_slow(loads, out=actual)

    assert_array_almost_equal(actual, expected)


@pytest.mark.parametrize("n", (10,))
def test_subside_zero_load(n):
    n_rows = n_cols = 2**n + 1

    grid = RasterModelGrid((n_rows, n_cols), xy_spacing=10.0)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    flex = Flexure(grid, method="flexure")

    loads = grid.zeros(at="node").reshape(grid.shape)

    actual = flex.subside_loads(loads)
    assert_array_almost_equal(actual, 0.0)

    actual = grid.empty(at="node")
    flex.subside_loads(loads, out=actual)
    assert_array_almost_equal(actual, 0.0)

    actual = flex.subside_loads(
        (0.0,), row_col_of_load=((0, n_rows - 1), (n_cols - 1, 0))
    )
    assert_array_almost_equal(actual, 0.0)


def test_subside_is_symmetric():
    n, load_0 = 10, 1e9

    n_rows = 2**n + 1
    n_cols = 2**n + 1

    grid = RasterModelGrid((n_rows, n_cols), xy_spacing=10.0)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    flex = Flexure(grid, method="flexure")

    loads = [load_0] * 5
    row_col_of_load = [
        [2 ** (n - 1), 0, 0, n_rows - 1, n_rows - 1],
        [2 ** (n - 1), 0, n_cols - 1, 0, n_cols - 1],
    ]

    dz = np.empty((n_rows, n_cols))
    flex.subside_loads(loads, row_col_of_load=row_col_of_load, out=dz)

    assert_array_almost_equal(dz, dz[:, ::-1])
    assert_array_almost_equal(dz, dz[::-1, :])

    assert not np.allclose(dz, 0.0)


def test_subside_negative_load():
    n, load_0 = 10, 1e9

    n_rows = 2**n + 1
    n_cols = 2**n + 1

    grid = RasterModelGrid((n_rows, n_cols), xy_spacing=10.0)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    flex = Flexure(grid, method="flexure")

    row_col_of_load = [
        [2 ** (n - 1), 0, 0, n_rows - 1, n_rows - 1],
        [2 ** (n - 1), 0, n_cols - 1, 0, n_cols - 1],
    ]

    actual = np.empty((n_rows, n_cols))
    flex.subside_loads([load_0] * 5, row_col_of_load=row_col_of_load, out=actual)

    expected = np.empty((n_rows, n_cols))
    flex.subside_loads([-load_0] * 5, row_col_of_load=row_col_of_load, out=expected)

    assert_array_almost_equal(actual, -expected)


def test_max_deflection_under_load():
    n, load_0 = 10, 1e9

    n_rows = 2**n + 1
    n_cols = 2**n + 1

    grid = RasterModelGrid((n_rows, n_cols), xy_spacing=10.0)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    flex = Flexure(grid, method="flexure")

    loads = [load_0]
    row_col_of_load = [[2 ** (n - 1)], [2 ** (n - 1)]]

    dz = np.zeros((n_rows, n_cols))
    flex.subside_loads(loads, row_col_of_load=row_col_of_load, out=dz)

    assert_array_equal(
        np.unravel_index(np.argmax(dz, axis=None), (n_rows, n_cols)),
        (row_col_of_load[0][0], row_col_of_load[1][0]),
    )


@pytest.mark.parametrize("method", ["airy", "flexure"])
def test_method_names(method):
    grid = RasterModelGrid((20, 20), xy_spacing=1e3)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    assert Flexure(grid, method=method).method == method


@pytest.mark.parametrize("method", ["Airy", "FlExUrE", "", " airy", "airy ", "foo"])
def test_method_bad_name(method):
    grid = RasterModelGrid((20, 20), xy_spacing=1e3)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    with pytest.raises(ValueError, match=f"{method}: method not understood"):
        Flexure(grid, method=method)


@pytest.mark.parametrize("size", (399, 401))
def test_out_is_wrong_size(size):
    grid = RasterModelGrid((20, 20), xy_spacing=1e3)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    flexure = Flexure(grid)
    with pytest.raises(ValueError, match="'out' array has incorrect size"):
        flexure.subside_loads(grid.zeros(at="node"), out=np.empty(size))


@pytest.mark.parametrize("size", (399, 401))
def test_loads_is_wrong_size(size):
    grid = RasterModelGrid((20, 20), xy_spacing=1e3)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    flexure = Flexure(grid)
    with pytest.raises(ValueError, match="'loads' array has incorrect size"):
        flexure.subside_loads(np.zeros(size), row_col_of_load=None)


def test_row_col_of_load_is_wrong_size():
    grid = RasterModelGrid((20, 20), xy_spacing=1e3)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    flexure = Flexure(grid)
    with pytest.raises(ValueError, match="'rows' and 'cols' must have the same length"):
        flexure.subside_loads([0, 0, 0], row_col_of_load=([0, 1, 2], [2, 3]))


@pytest.mark.parametrize("eet", [1e4, 1e3])
def test_eet_attribute(eet):
    grid = RasterModelGrid((20, 20), xy_spacing=1e3)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    assert Flexure(grid, eet=eet).eet == pytest.approx(eet)


@pytest.mark.parametrize("eet", [0.0, -1e3])
def test_eet_bad_value(eet):
    grid = RasterModelGrid((20, 20), xy_spacing=1e3)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    with pytest.raises(ValueError):
        assert Flexure(grid, eet=eet)


def test_youngs_attribute():
    grid = RasterModelGrid((20, 20), xy_spacing=10e3)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    for val in (10e3, 1e3):
        assert Flexure(grid, youngs=val).youngs == pytest.approx(val)


def test_gravity_attribute():
    grid = RasterModelGrid((20, 20), xy_spacing=10e3)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    for val in (10e3, 1e3):
        assert Flexure(grid, gravity=val).gravity == pytest.approx(val)


def test_rho_mantle_attribute():
    grid = RasterModelGrid((20, 20), xy_spacing=10e3)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
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
    n = 5
    n_rows = 2**n + 1
    n_cols = 2**n + 1

    load_0 = 1e9

    grid = RasterModelGrid((n_rows, n_cols), xy_spacing=1e3)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    flex = Flexure(grid, method="flexure")

    load = grid.at_node["lithosphere__overlying_pressure_increment"].reshape(grid.shape)
    load[2 ** (n - 1), 2 ** (n - 1)] = load_0

    flex.update()
    dz = flex.grid.at_node["lithosphere_surface__elevation_increment"].reshape(
        (n_rows, n_cols)
    )

    assert_array_almost_equal(dz, dz[:, ::-1])
    assert_array_almost_equal(dz, dz[::-1, :])

    assert_array_equal(np.argmax(dz), np.argmax(load))
    assert dz[2 ** (n - 1), 2 ** (n - 1)] > 0.0


def test_subside_loads():
    n, load_0 = 11, 1e9

    grid = RasterModelGrid((n, n), xy_spacing=1e3)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
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
