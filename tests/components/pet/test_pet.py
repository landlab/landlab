"""
Unit tests for landlab.components.pet.potential_evapotranspiration_field
"""

import numpy as np
import pytest
from landlab import RasterModelGrid
from numpy.testing import assert_array_almost_equal

(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e0, 10e0), (0.0, 0.0))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def test_name(pet):
    assert pet.name == "PotentialEvapotranspiration"


def test_output_var_names(pet):
    assert sorted(pet.output_var_names) == [
        "surface__potential_evapotranspiration_rate",
    ]


def test_var_units(pet):
    assert set(pet.input_var_names) | set(pet.output_var_names) == set(
        dict(pet.units).keys()
    )
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


def test_temperature_validation(pet):
    with pytest.raises(ValueError):
        pet._validate_temperature_range(None, None)

    with pytest.raises(ValueError):
        pet._validate_temperature_range(10.0, 0.0)


def test_value_fixation(pet):
    error_value = 0.0
    fixed_value = 1.0

    for name in pet.grid["cell"]:
        field = pet.grid["cell"][name]
        pet._fix_values(field, error_value, fixed_value)
        assert not np.any(field == error_value)
        assert_array_almost_equal(field, np.ones(pet.grid.number_of_cells))


def test_priestley_taylor_method(pet):
    sample_grid = RasterModelGrid((5, 4), xy_spacing=(0.2, 0.2))
    sample_grid.at_node["topographic__elevation"] = [
        [0.0, 0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0, 1.0],
        [2.0, 2.0, 2.0, 2.0],
        [3.0, 4.0, 4.0, 3.0],
        [4.0, 4.0, 4.0, 4.0],
    ]

    pet._method = "PriestleyTaylor"
    pet._grid = sample_grid
    pet._latitude = 40.0

    pet.update()
    assert not np.allclose(pet._PET_value, 0.0)

def test_penman_method(pet):
    sample_grid = RasterModelGrid((5, 4), xy_spacing=(0.2, 0.2))
    sample_grid.at_node["topographic__elevation"] = [
        [0.0, 0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0, 1.0],
        [2.0, 2.0, 2.0, 2.0],
        [3.0, 4.0, 4.0, 3.0],
        [4.0, 4.0, 4.0, 4.0],
    ]

    pet._method = "PenmanMonteith"
    pet._grid = sample_grid
    pet._latitude = 40.0

    pet.update()
    assert not np.allclose(pet._PET_value, 0.0)
