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


@pytest.mark.benchmark(group="grid-size")
@pytest.mark.parametrize("n", [4, 5, 6, 7, 8, 9, 10])
@pytest.mark.parametrize("method", ["old", "new", "cython"])
def test_one_load_bench(benchmark, n, method):
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

    benchmark(func, *args, **kwds)

    assert_array_almost_equal(w, w[:, ::-1])
    assert_array_almost_equal(w, w[::-1, :])

    assert_array_equal(
        np.unravel_index(np.argmax(w, axis=None), w.shape), (2 ** (n - 1), 2 ** (n - 1))
    )


@pytest.mark.benchmark(group="number-of-loads")
@pytest.mark.parametrize("n_loads", [0, 1, 2, 3, 4, 5, 6, 7, 8])
@pytest.mark.parametrize("method", ["old", "new", "cython"])
def test_number_of_loads_bench(benchmark, n_loads, method):
    load_0 = 1e9

    n = 8
    size = 2**n + 1

    grid = RasterModelGrid((size, size), xy_spacing=10.0)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    flex = Flexure(grid, method="flexure")

    w = np.zeros((size, size))

    nodes = np.arange(0, size * size, 2 ** (n - n_loads))
    n_loads = len(nodes)

    loads = np.full(n_loads, load_0)
    row_col_of_load = np.unravel_index(nodes, w.shape)

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

    benchmark(func, *args, **kwds)


@pytest.mark.benchmark(group="row-col-of-grid")
@pytest.mark.parametrize("n_loads", [0, 1, 2, 3, 4, 5, 6, 7, 8])
def test_subside_loads_with_row_col_bench(benchmark, n_loads):
    n, load_0 = 8, 1e9

    size = 2**n + 1

    grid = RasterModelGrid((size, size), xy_spacing=10.0)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    flex = Flexure(grid, method="flexure")

    nodes = np.arange(0, size * size, 2 ** (n - n_loads))
    n_loads = len(nodes)

    loads = np.full(n_loads, load_0)
    row_col_of_load = np.unravel_index(nodes, grid.shape)

    dz = grid.zeros(at="node")
    benchmark(flex.subside_loads, loads, row_col_of_load=row_col_of_load, out=dz)


@pytest.mark.benchmark(group="row-col-of-grid")
@pytest.mark.parametrize("n_loads", [0, 1, 2, 3, 4, 5, 6, 7, 8])
def test_subside_loads_without_row_col_bench(benchmark, n_loads):
    n, load_0 = 8, 1e9

    size = 2**n + 1

    grid = RasterModelGrid((size, size), xy_spacing=10.0)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    flex = Flexure(grid, method="flexure")

    loads = grid.zeros(at="node")

    nodes = np.arange(0, size * size, 2 ** (n - n_loads))
    n_loads = len(nodes)
    loads[nodes] = load_0

    actual = grid.zeros(at="node")
    benchmark(
        flex.subside_loads, loads.reshape(grid.shape), row_col_of_load=None, out=actual
    )


@pytest.mark.benchmark(group="speedup")
@pytest.mark.parametrize("method", ["subside_loads", "subside_loads_slow"])
@pytest.mark.parametrize("n_loads", [0, 1, 2, 3, 4, 5, 6, 7, 8])
def test_flexure_loads_everywhere(benchmark, method, n_loads):
    n, load_0 = 8, 1e9

    size = 2**n + 1

    grid = RasterModelGrid((size, size), xy_spacing=10.0)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    flex = Flexure(grid, method="flexure")

    loads = grid.zeros(at="node")
    nodes = np.arange(0, size * size, 2 ** (n - n_loads))
    n_loads = len(nodes)
    loads[nodes] = load_0

    dz = np.zeros((size, size))
    benchmark(getattr(flex, method), loads, out=dz)

    # assert_array_almost_equal(dz, dz[:, ::-1])
    # assert_array_almost_equal(dz, dz[::-1, :])


# @pytest.mark.benchmark(group="speedup")
# @pytest.mark.parametrize("method", ["subside_loads", "subside_loads_slow"])
# def test_flexure_loads_somewhere(benchmark, method):
#     n, load_0 = 8, 1e9

#     size = 2**n + 1

#     grid = RasterModelGrid((size, size), xy_spacing=10.0)
#     grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
#     flex = Flexure(grid, method="flexure")

#     loads = grid.zeros(at="node").reshape(grid.shape)
#     loads[:size//2, :size//2] = load_0
#     # loads.fill(load_0)

#     dz = np.zeros((size, size))
#     benchmark(getattr(flex, method), loads, out=dz)


@pytest.mark.benchmark(group="speedup")
@pytest.mark.parametrize("method", ["subside_loads", "subside_loads_slow"])
def test_flexure_deflection_is_proportional_to_load(benchmark, method):
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


def test_update_bench(benchmark):
    n, load_0 = 10, 1e9

    n_rows = 2**n + 1
    n_cols = 2**n + 1

    grid = RasterModelGrid((n_rows, n_cols), xy_spacing=10.0)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    flex = Flexure(grid, method="flexure")

    grid.at_node["lithosphere__overlying_pressure_increment"].reshape(grid.shape)[
        2 ** (n - 1), 2 ** (n - 1)
    ] = load_0
    benchmark(flex.update)
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


def test_subside_is_symetric():
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

    dz = np.zeros((n_rows, n_cols))
    flex.subside_loads(loads, row_col_of_load=row_col_of_load, out=dz)

    assert_array_almost_equal(dz, dz[:, ::-1])
    assert_array_almost_equal(dz, dz[::-1, :])


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

    actual = np.zeros((n_rows, n_cols))
    flex.subside_loads([load_0] * 5, row_col_of_load=row_col_of_load, out=actual)

    expected = np.zeros((n_rows, n_cols))
    flex.subside_loads([-load_0] * 5, row_col_of_load=row_col_of_load, out=actual)

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


def test_method_bad_name():
    grid = RasterModelGrid((20, 20), xy_spacing=1e3)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    with pytest.raises(ValueError):
        Flexure(grid, method="bad-name")


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
