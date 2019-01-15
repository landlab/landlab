import os

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import create_grid

_TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def test_no_grid_value():
    dict_like = {"foo": "bar"}
    with pytest.raises(ValueError):
        create_grid(dict_like)


def test_bad_grid_name():
    dict_like = {"grid": "MagicModelGrid"}
    with pytest.raises(ValueError):
        create_grid(dict_like)


def test_two_grid_types():
    dict_like = {"grid": {"RasterModelGrid": [(4, 5)], "HexModelGrid": [6, 7]}}
    with pytest.raises(ValueError):
        create_grid(dict_like)


def test_bad_field_location():
    dict_like = {
        "grid": {"RasterModelGrid": [(4, 5)]},
        "fields": {"at_bad_location": {"foo": "bar"}},
    }
    with pytest.raises(ValueError):
        create_grid(dict_like)


def test_bad_field_function():
    dict_like = {
        "grid": {"RasterModelGrid": [(4, 5)]},
        "fields": {"at_node": {"new_field_name": {"not_a_function": ["bar", "spam"]}}},
    }
    with pytest.raises(ValueError):
        create_grid(dict_like)


def test_simple_create():
    filename = os.path.join(_TEST_DATA_DIR, "simple_create.yaml")
    mg = create_grid(filename)

    assert mg.number_of_nodes == 20
    assert "topographic__elevation" in mg.at_node

    x_of_node = np.array(
        [
            0.,
            3.,
            6.,
            9.,
            12.,
            0.,
            3.,
            6.,
            9.,
            12.,
            0.,
            3.,
            6.,
            9.,
            12.,
            0.,
            3.,
            6.,
            9.,
            12.,
        ]
    )

    status_at_node = np.array(
        [4, 4, 4, 4, 4, 4, 0, 0, 0, 4, 4, 0, 0, 0, 4, 4, 4, 4, 4, 4], dtype=np.uint8
    )
    topographic__elevation = np.array(
        [
            [-2., 4., 10., 16., 22.],
            [2., 8., 14., 20., 26.],
            [6., 12., 18., 24., 30.],
            [10., 16., 22., 28., 34.],
        ]
    )

    assert_array_equal(mg.x_of_node, x_of_node)
    assert_array_equal(status_at_node, mg.status_at_node)
    assert_array_equal(
        topographic__elevation,
        np.round(mg.at_node["topographic__elevation"].reshape(mg.shape), decimals=2),
    )


def test_esri_ascii_create():
    filename = os.path.join(_TEST_DATA_DIR, "4_x_3_no_nodata_value.asc")
    dict_like = {
        "grid": {"RasterModelGrid": [(4, 3)]},
        "fields": {
            "at_node": {"topographic__elevation": {"read_esri_ascii": [filename]}}
        },
    }
    mg = create_grid(dict_like)

    x_of_node = np.array([2])
    y_of_node = np.array([2])

    topographic__elevation = np.array(
        [9., 10., 11., 6., 7., 8., 3., 4., 5., 0., 1., 2.]
    )

    assert_array_equal(topographic__elevation, mg.at_node["topographic__elevation"])


def test_read_netcdf_create():
    pass
