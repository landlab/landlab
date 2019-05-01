"""
Unit tests for landlab.components.pet.potential_evapotranspiration_field
"""
import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e0, 10e0), (0.0, 0.0))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def test_name(pet):
    assert pet.name == "Potential Evapotranspiration"


def test_input_var_names(pet):
    assert pet.input_var_names == ("radiation__ratio_to_flat_surface",)


def test_output_var_names(pet):
    assert sorted(pet.output_var_names) == [
        "radiation__incoming_shortwave_flux",
        "radiation__net_flux",
        "radiation__net_longwave_flux",
        "radiation__net_shortwave_flux",
        "surface__potential_evapotranspiration_rate",
    ]


def test_var_units(pet):
    assert set(pet.input_var_names) | set(pet.output_var_names) == set(
        dict(pet.units).keys()
    )

    assert pet.var_units("radiation__incoming_shortwave_flux") == "W/m^2"
    assert pet.var_units("radiation__net_flux") == "W/m^2"
    assert pet.var_units("radiation__net_longwave_flux") == "W/m^2"
    assert pet.var_units("radiation__net_shortwave_flux") == "W/m^2"
    assert pet.var_units("radiation__ratio_to_flat_surface") == "None"
    assert pet.var_units("surface__potential_evapotranspiration_rate") == "mm"


def test_grid_shape(pet):
    assert pet.grid.number_of_node_rows == _SHAPE[0]
    assert pet.grid.number_of_node_columns == _SHAPE[1]


def test_grid_x_extent(pet):
    assert pet.grid.extent[1] == (_SHAPE[1] - 1) * _SPACING[1]


def test_grid_y_extent(pet):
    assert pet.grid.extent[0] == (_SHAPE[0] - 1) * _SPACING[0]


def test_field_getters(pet):
    for name in pet.grid["node"]:
        field = pet.grid["node"][name]
        assert isinstance(field, np.ndarray)
        assert field.shape == (
            pet.grid.number_of_node_rows * pet.grid.number_of_node_columns,
        )

    for name in pet.grid["cell"]:
        field = pet.grid["cell"][name]
        assert isinstance(field, np.ndarray)
        assert field.shape == (
            pet.grid.number_of_cell_rows * pet.grid.number_of_cell_columns,
        )

    with pytest.raises(KeyError):
        pet.grid["not_a_var_name"]


def test_field_initialized_to_zero(pet):
    for name in pet.grid["node"]:
        field = pet.grid["node"][name]
        assert_array_almost_equal(field, np.zeros(pet.grid.number_of_nodes))
    for name in pet.grid["cell"]:
        field = pet.grid["cell"][name]
        assert_array_almost_equal(field, np.zeros(pet.grid.number_of_cells))
