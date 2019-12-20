import numpy as np
from pytest import approx

from landlab.utils import add_halo


def test_add_halo_with_defaults():
    values = np.arange(6, dtype=int).reshape((2, 3))
    with_halo = add_halo(values)
    assert np.all(with_halo[1:-1, 1:-1] == values)
    assert values.dtype == with_halo.dtype


def test_add_halo_to_float_array():
    values = np.arange(6, dtype=float).reshape((2, 3))
    with_halo = add_halo(values)
    assert np.all(with_halo[1:-1, 1:-1] == approx(values))
    assert values.dtype == with_halo.dtype


def test_add_halo_to_bool_array():
    values = np.full((2, 3), True, dtype=bool)
    with_halo = add_halo(values)
    assert np.all(with_halo[1:-1, 1:-1] == values)
    assert values.dtype == with_halo.dtype


def test_add_halo_with_zero_halo():
    values = np.ones((2, 3), dtype=int)
    with_halo = add_halo(values, halo=0)
    assert np.all(with_halo == values)


def test_add_halo_with_bigger_halo():
    values = np.ones((2, 4), dtype=int)
    with_halo = add_halo(values, halo=3, halo_value=0)
    assert np.all(
        with_halo
        == [
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 1, 1, 1, 0, 0, 0],
            [0, 0, 0, 1, 1, 1, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ]
    )


def test_add_halo_where_halo_value_is_zero():
    values = np.ones((3, 4), dtype=int)
    with_halo = add_halo(values, halo_value=0)
    assert np.all(
        with_halo
        == [
            [0, 0, 0, 0, 0, 0],
            [0, 1, 1, 1, 1, 0],
            [0, 1, 1, 1, 1, 0],
            [0, 1, 1, 1, 1, 0],
            [0, 0, 0, 0, 0, 0],
        ]
    )


def test_add_halo_where_halo_value_is_nan():
    values = np.ones((2, 3), dtype=float)
    with_halo = add_halo(values, halo=1, halo_value=np.nan)
    assert np.all(with_halo[1:-1, 1:-1] == approx(values))
    assert np.all(
        np.isnan(with_halo)
        == [
            [True, True, True, True, True],
            [True, False, False, False, True],
            [True, False, False, False, True],
            [True, True, True, True, True],
        ]
    )
