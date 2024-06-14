"""
Unit tests for \
    landlab.components.vegetation_dynamics.vegetation_dynamics
"""

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e0, 10e0), (0.0, 0.0))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def test_name(veg):
    assert veg.name == "Vegetation"


def test_input_var_names(veg):
    assert sorted(veg.input_var_names) == [
        "surface__evapotranspiration",
        "surface__potential_evapotranspiration_30day_mean",
        "surface__potential_evapotranspiration_rate",
        "vegetation__plant_functional_type",
        "vegetation__water_stress",
    ]


def test_output_var_names(veg):
    assert sorted(veg.output_var_names) == [
        "vegetation__cover_fraction",
        "vegetation__dead_biomass",
        "vegetation__dead_leaf_area_index",
        "vegetation__live_biomass",
        "vegetation__live_leaf_area_index",
    ]


def test_var_units(veg):
    assert set(veg.input_var_names) | set(veg.output_var_names) == set(
        dict(veg.units).keys()
    )

    assert veg.var_units("vegetation__live_leaf_area_index") == "None"
    assert veg.var_units("vegetation__dead_leaf_area_index") == "None"
    assert veg.var_units("vegetation__cover_fraction") == "None"
    assert veg.var_units("surface__evapotranspiration") == "mm"
    assert veg.var_units("surface__potential_evapotranspiration_rate") == "mm"
    assert veg.var_units("surface__potential_evapotranspiration_30day_mean") == "mm"
    assert veg.var_units("vegetation__water_stress") == "None"
    assert veg.var_units("vegetation__live_biomass") == "g m^-2 d^-1"
    assert veg.var_units("vegetation__dead_biomass") == "g m^-2 d^-1"
    assert veg.var_units("vegetation__plant_functional_type") == "None"


def test_grid_shape(veg):
    assert veg.grid.number_of_node_rows == _SHAPE[0]
    assert veg.grid.number_of_node_columns == _SHAPE[1]


def test_grid_x_extent(veg):
    assert veg.grid.extent[1] == (_SHAPE[1] - 1) * _SPACING[1]


def test_grid_y_extent(veg):
    assert veg.grid.extent[0] == (_SHAPE[0] - 1) * _SPACING[0]


def test_field_getters(veg):
    for name in veg.grid["node"]:
        field = veg.grid["node"][name]
        assert isinstance(field, np.ndarray)
        assert field.shape == (
            veg.grid.number_of_node_rows * veg.grid.number_of_node_columns,
        )

    for name in veg.grid["cell"]:
        field = veg.grid["cell"][name]
        assert isinstance(field, np.ndarray)
        assert field.shape == (
            veg.grid.number_of_cell_rows * veg.grid.number_of_cell_columns,
        )

    with pytest.raises(KeyError):
        veg.grid["not_a_var_name"]


def test_field_initialized_to_zero(veg):
    for name in veg.grid["node"]:
        field = veg.grid["node"][name]
        assert_array_almost_equal(field, np.zeros(veg.grid.number_of_nodes))
    for name in veg.grid["cell"]:
        field = veg.grid["cell"][name]
        assert_array_almost_equal(field, np.zeros(veg.grid.number_of_cells))
