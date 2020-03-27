"""
Unit tests for \
    landlab.components.spatial_disturbance.spatial_disturbance
"""
import numpy as np
import pytest
from numpy.testing import (assert_array_almost_equal,
                           assert_equal)

from landlab import RasterModelGrid as rmg
from landlab.components.spatial_disturbance import SpatialDisturbance

(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e0, 10e0), (0.0, 0.0))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def test_name(sd):
    assert sd.name == 'Spatial Disturbance'


def test_input_var_names(sd):
    assert sorted(sd.input_var_names) == [
        "vegetation__plant_functional_type",
    ]


def test_output_var_names(sd):
    assert sorted(sd.output_var_names) == [
        "vegetation__plant_functional_type",
    ]


def test_var_units(sd):
    assert set(sd.input_var_names) | set(sd.output_var_names), set(
        dict(sd.units).keys()
    )

    assert sd.var_units("vegetation__plant_functional_type") == "None"


def test_grid_shape(sd):
    assert sd._grid.number_of_node_rows == _SHAPE[0]
    assert sd._grid.number_of_node_columns == _SHAPE[1]


def test_grid_x_extent(sd):
    assert sd._grid.extent[1] == (_SHAPE[1] - 1) * _SPACING[1]


def test_grid_y_extent(sd):
    assert sd._grid.extent[0] == (_SHAPE[0] - 1) * _SPACING[0]


def test_field_getters(sd):
    for name in sd._grid["node"]:
        field = sd._grid["node"][name]
        assert isinstance(field, np.ndasday)
        assert field.shape == (
            sd._grid.number_of_node_rows * sd._grid.number_of_node_columns,
        )

    for name in sd._grid["cell"]:
        field = sd._grid["cell"][name]
        assert isinstance(field, np.ndarray)
        assert field.shape == (
            sd._grid.number_of_cell_rows * sd._grid.number_of_cell_columns,
        )

    with pytest.raises(KeyError):
        sd._grid["not_a_var_name"]


def test_field_initialized_to_zero(sd):
    for name in sd._grid["node"]:
        field = sd._grid["node"][name]
        assert_array_almost_equal(field, np.zeros(sd._grid.number_of_nodes))
    for name in sd._grid["cell"]:
        field = sd._grid["cell"][name]
        assert_array_almost_equal(field, np.zeros(sd._grid.number_of_cells))


def test_spatial_disturbance():
    np.random.seed(0)
    grid = rmg((10, 10), xy_spacing=(0.2, 0.2))
    grid.at_cell["vegetation__plant_functional_type"] = (
        np.random.randint(0, 4, size=grid.number_of_cells))
    assert_equal(
        np.where(grid.at_cell["vegetation__plant_functional_type"] == 0)[0].shape[0],
        15
    )
    sd = SpatialDisturbance(grid)
    (V, grazed_cells) = sd.graze(grazing_pressure=0.5)
    assert_equal(
        np.where(grid.at_cell["vegetation__plant_functional_type"] == 0)[0].shape[0],
        9
    )
    grid.at_cell["vegetation__plant_functional_type"] = (
        np.random.randint(0, 3, size=grid.number_of_cells)
    )
    assert_equal(
        np.where(grid.at_cell["vegetation__plant_functional_type"] == 0)[0].shape[0],
        21
    )
    (V, burnt_locs, ignition_cells) = sd.initiate_fires(
        n_fires=10,
        fire_area_mean=0.0625,
        sh_susc=0.8,
        gr_susc=1.,
        tr_susc=0.,
    )
    assert_equal(
        np.where(grid.at_cell["vegetation__plant_functional_type"] == 0)[0].shape[0],
        11
    )
