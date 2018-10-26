#! /usr/bin/env python
"""
Unit tests for landlab.components.flexure.flexure
"""
import pytest

import numpy as np


(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e3, 10e3), (0., 0.))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


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
        assert np.all(field == 0.)
