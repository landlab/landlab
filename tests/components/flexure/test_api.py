import numpy as np
from pytest import approx
from pytest import raises

from landlab.components.flexure import get_flexure_parameter
from landlab.components.flexure import subside_point_load


def test_calc_flexure_parameter():
    eet, youngs = 65000.0, 7e10
    alpha_1d = get_flexure_parameter(eet, youngs, 1)
    alpha_2d = get_flexure_parameter(eet, youngs, 2)
    assert alpha_2d == approx(84828.72)
    assert alpha_1d == approx(np.power(4.0, 0.25) * alpha_2d)


def test_calc_flexure_parameter_bad_ndim():
    with raises(ValueError):
        get_flexure_parameter(65000.0, 7e10, 0)
    with raises(ValueError):
        get_flexure_parameter(65000.0, 7e10, 3)


def test_subside_point_load():
    params = {"eet": 65000.0, "youngs": 7e10}
    load = 1e9
    loc = (5000.0, 2500.0)

    x = np.arange(0, 10000, 1000.0)
    y = np.arange(0, 5000, 1000.0)
    (x, y) = np.meshgrid(x, y)
    x.shape = (x.size,)
    y.shape = (y.size,)

    dz_one_load = subside_point_load(load, loc, (x, y), params=params)

    n_loads = 16

    dz = subside_point_load(
        np.full(n_loads, load / n_loads),
        np.full((n_loads, 2), loc).T,
        (x, y),
        params=params,
    )

    assert dz.mean() == approx(dz_one_load.mean())
    assert dz.min() == approx(dz_one_load.min())
    assert dz.max() == approx(dz_one_load.max())


def test_point_load_1d_with_scalar_args():
    params = {"eet": 65000.0, "youngs": 7e10}
    load = 1e9

    x = np.arange(0, 10000, 1000.0)

    dz = subside_point_load(load, 5000.0, x, params=params)

    assert dz.shape == x.shape


def test_point_load_1d_is_symetric():
    params = {"eet": 65000.0, "youngs": 7e10}
    load = 1e9

    x = np.arange(0, 11000, 1000.0)
    i_mid = (len(x) - 1) // 2

    dz = subside_point_load(load, x[i_mid], x, params=params)

    assert np.argmax(dz) == i_mid
    assert dz[i_mid] > 0.0
    assert np.all(dz[i_mid::-1] == dz[i_mid:])


def test_point_load_2d_is_symetric():
    params = {"eet": 65000.0, "youngs": 7e10}
    load = 1e9

    n = 11
    n_mid = (n - 1) // 2

    (x, y) = np.meshgrid(np.arange(0.0, n), np.arange(0.0, n))

    dz = subside_point_load(
        load,
        (x[n_mid, n_mid], y[n_mid, n_mid]),
        (x.flatten(), y.flatten()),
        params=params,
    )
    dz.shape = (n, n)

    assert np.argmax(dz) == np.ravel_multi_index((n_mid, n_mid), (n, n))
    assert dz[n_mid, n_mid] > 0.0
    assert np.all(dz[:, n_mid::-1] == dz[:, n_mid:])
    assert np.all(dz[n_mid::-1, :] == dz[n_mid:, :])


def test_subside_point_load_1d():
    params = {"eet": 65000.0, "youngs": 7e10}
    load = 1e9
    loc = (5000.0,)

    x = np.arange(0, 10000, 1000.0)

    dz_one_load = subside_point_load(load, loc, (x,), params=params)

    n_loads = 32
    dz = subside_point_load(
        np.full(n_loads, load / n_loads), (np.full(n_loads, loc),), (x,), params=params
    )

    assert np.all(dz == approx(dz_one_load))


def test_out_keyword():
    params = {"eet": 65000.0, "youngs": 7e10}
    load = 1e9
    loc = (5000.0,)

    x = np.arange(0, 10000, 1000.0)
    out = np.empty_like(x)

    dz = subside_point_load(load, loc, (x,), params=params, out=out)

    assert dz is out


def test_dimension_mismatch():
    params = {"eet": 65000.0, "youngs": 7e10}
    load = 1e9
    loc = (5000.0,)

    x = np.arange(0, 10000, 1000.0)
    y = np.zeros_like(x)
    with raises(ValueError):
        subside_point_load(load, loc, (x, y), params=params)


def test_ndim_too_big():
    params = {"eet": 65000.0, "youngs": 7e10}

    load = 1e9
    loc = ((5000.0, 1000.0, 500.0),)

    x = np.arange(0, 10000, 1000.0)
    y = np.zeros_like(x)
    z = np.zeros_like(x)
    with raises(ValueError):
        subside_point_load(load, loc, (x, y, z), params=params)
