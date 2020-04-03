"""
Unit tests for \
    landlab.components.soil_moisture.soil_moisture_dynamics
"""
import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e0, 10e0), (0.0, 0.0))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def test_name(sm):
    assert sm.name == "Soil Moisture"


def test_input_var_names(sm):
    assert sorted(sm.input_var_names) == [
        "rainfall__daily_depth",
        "soil_moisture__initial_saturation_fraction",
        "surface__potential_evapotranspiration_rate",
        "vegetation__cover_fraction",
        "vegetation__live_leaf_area_index",
        "vegetation__plant_functional_type",
    ]


def test_output_var_names(sm):
    assert sorted(sm.output_var_names) == [
        "soil_moisture__root_zone_leakage",
        "soil_moisture__saturation_fraction",
        "surface__evapotranspiration",
        "surface__runoff",
        "vegetation__water_stress",
    ]


def test_var_units(sm):
    assert set(sm.input_var_names) | set(sm.output_var_names), set(
        dict(sm.units).keys()
    )

    assert sm.var_units("vegetation__cover_fraction") == "None"
    assert sm.var_units("vegetation__live_leaf_area_index") == "None"
    assert sm.var_units("surface__potential_evapotranspiration_rate") == "mm"
    assert sm.var_units("vegetation__plant_functional_type") == "None"
    assert sm.var_units("vegetation__water_stress") == "None"
    assert sm.var_units("soil_moisture__saturation_fraction") == "None"
    assert sm.var_units("soil_moisture__initial_saturation_fraction") == "None"
    assert sm.var_units("soil_moisture__root_zone_leakage") == "mm"
    assert sm.var_units("surface__runoff") == "mm"
    assert sm.var_units("surface__evapotranspiration") == "mm"
    assert sm.var_units("rainfall__daily_depth") == "mm"


def test_grid_shape(sm):
    assert sm.grid.number_of_node_rows == _SHAPE[0]
    assert sm.grid.number_of_node_columns == _SHAPE[1]


def test_grid_x_extent(sm):
    assert sm.grid.extent[1] == (_SHAPE[1] - 1) * _SPACING[1]


def test_grid_y_extent(sm):
    assert sm.grid.extent[0] == (_SHAPE[0] - 1) * _SPACING[0]


def test_field_getters(sm):
    for name in sm.grid["node"]:
        field = sm.grid["node"][name]
        assert isinstance(field, np.ndarray)
        assert field.shape == (
            sm.grid.number_of_node_rows * sm.grid.number_of_node_columns,
        )

    for name in sm.grid["cell"]:
        field = sm.grid["cell"][name]
        assert isinstance(field, np.ndarray)
        assert field.shape == (
            sm.grid.number_of_cell_rows * sm.grid.number_of_cell_columns,
        )

    with pytest.raises(KeyError):
        sm.grid["not_a_var_name"]


def test_field_initialized_to_zero(sm):
    for name in sm.grid["node"]:
        field = sm.grid["node"][name]
        assert_array_almost_equal(field, np.zeros(sm.grid.number_of_nodes))
    for name in sm.grid["cell"]:
        field = sm.grid["cell"][name]
        assert_array_almost_equal(field, np.zeros(sm.grid.number_of_cells))
