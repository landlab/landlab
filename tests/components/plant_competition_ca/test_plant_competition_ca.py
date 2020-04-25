"""
Unit tests for landlab.components.plant_competition_ca.plant_competition_ca
"""
import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e0, 10e0), (0.0, 0.0))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def test_name(ca_veg):
    assert ca_veg.name == "Cellular Automata Plant Competition"


def test_input_var_names(ca_veg):
    assert sorted(ca_veg.input_var_names) == [
        "vegetation__cumulative_water_stress",
        "vegetation__plant_functional_type",
    ]


def test_output_var_names(ca_veg):
    assert sorted(ca_veg.output_var_names) == ["plant__age", "plant__live_index"]


def test_var_units(ca_veg):
    assert set(ca_veg.input_var_names) | set(ca_veg.output_var_names) == set(
        dict(ca_veg.units).keys()
    )

    assert ca_veg.var_units("vegetation__cumulative_water_stress") == "None"
    assert ca_veg.var_units("vegetation__plant_functional_type") == "None"
    assert ca_veg.var_units("plant__live_index") == "None"
    assert ca_veg.var_units("plant__age") == "Years"


def test_grid_shape(ca_veg):
    assert ca_veg.grid.number_of_node_rows == _SHAPE[0]
    assert ca_veg.grid.number_of_node_columns == _SHAPE[1]


def test_grid_x_extent(ca_veg):
    assert ca_veg.grid.extent[1] == (_SHAPE[1] - 1) * _SPACING[1]


def test_grid_y_extent(ca_veg):
    assert ca_veg.grid.extent[0] == (_SHAPE[0] - 1) * _SPACING[0]


def test_field_getters(ca_veg):
    for name in ca_veg.grid["node"]:
        field = ca_veg.grid["node"][name]
        assert isinstance(field, np.ndarray)
        assert field.shape == (
            ca_veg.grid.number_of_node_rows * ca_veg.grid.number_of_node_columns,
        )

    for name in ca_veg.grid["cell"]:
        field = ca_veg.grid["cell"][name]
        assert isinstance(field, np.ndarray)
        assert field.shape, (
            ca_veg.grid.number_of_cell_rows * ca_veg.grid.number_of_cell_columns,
        )

    with pytest.raises(KeyError):
        ca_veg.grid["not_a_var_name"]


def test_field_initialized_to_zero(ca_veg):
    for name in ca_veg.grid["node"]:
        field = ca_veg.grid["node"][name]
        assert_array_almost_equal(field, np.zeros(ca_veg.grid.number_of_nodes))
