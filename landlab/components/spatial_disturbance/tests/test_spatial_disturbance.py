"""
Unit tests for \
    landlab.components.resource_redistribution.resource_redistribution
"""
import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e0, 10e0), (0.0, 0.0))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def test_name(rr):
    assert rr.name == 'Resource Redistribution'


def test_input_var_names(rr):
    assert sorted(rr.input_var_names) == [
        "soil__resources",
        "vegetation__plant_functional_type",
    ]


def test_output_var_names(rr):
    assert sorted(rr.output_var_names) == [
        "soil__resources",
        "vegetation__plant_functional_type",
    ]


def test_var_units(rr):
    assert set(rr.input_var_names) | set(rr.output_var_names), set(
        dict(rr.units).keys()
    )

    assert rr.var_units("soil__resources") == "None"
    assert rr.var_units("vegetation__plant_functional_type") == "None"


def test_grid_shape(rr):
    assert rr.grid.number_of_node_rows == _SHAPE[0]
    assert rr.grid.number_of_node_columns == _SHAPE[1]


def test_grid_x_extent(rr):
    assert rr.grid.extent[1] == (_SHAPE[1] - 1) * _SPACING[1]


def test_grid_y_extent(rr):
    assert rr.grid.extent[0] == (_SHAPE[0] - 1) * _SPACING[0]


def test_field_getters(rr):
    for name in rr.grid["node"]:
        field = rr.grid["node"][name]
        assert isinstance(field, np.ndarray)
        assert field.shape == (
            rr.grid.number_of_node_rows * rr.grid.number_of_node_columns,
        )

    for name in rr.grid["cell"]:
        field = rr.grid["cell"][name]
        assert isinstance(field, np.ndarray)
        assert field.shape == (
            rr.grid.number_of_cell_rows * rr.grid.number_of_cell_columns,
        )

    with pytest.raises(KeyError):
        rr.grid["not_a_var_name"]


def test_field_initialized_to_zero(rr):
    for name in rr.grid["node"]:
        field = rr.grid["node"][name]
        assert_array_almost_equal(field, np.zeros(rr.grid.number_of_nodes))
    for name in rr.grid["cell"]:
        field = rr.grid["cell"][name]
        assert_array_almost_equal(field, np.zeros(rr.grid.number_of_cells))
