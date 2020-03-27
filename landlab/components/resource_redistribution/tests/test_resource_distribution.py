"""
Unit tests for \
    landlab.components.resource_redistribution.resource_redistribution
"""
import numpy as np
import pytest
from numpy.testing import (assert_array_almost_equal,
                           assert_equal,
                           assert_almost_equal,
                           assert_array_equal)

from landlab import RasterModelGrid as rmg
from landlab.components.resource_redistribution import ResourceRedistribution

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


def test_resource_redistribution():
    np.random.seed(0)
    grid = rmg((5, 4), xy_spacing=(0.2, 0.2))
    grid.at_cell["vegetation__plant_functional_type"] = (
        np.random.randint(0, 5, size=grid.number_of_cells))
    assert_equal(
        grid.at_cell["vegetation__plant_functional_type"],
        np.array([4, 0, 3, 3, 3, 1])
    )
    grid.at_cell["soil__resources"] = (
        np.ones(grid.number_of_cells, dtype=float))
    rr = ResourceRedistribution(grid)
    (eroded_soil,
     eroded_soil_shrub,
     burnt_shrub,
     burnt_grass,
     bare_cells) = rr.erode()
    assert_almost_equal(eroded_soil, 0.16)
    (burnt_shrubs_neigh,
     exclusive,
     shrub_exclusive,
     grass_exclusive,
     bare_exclusive,
     eroded_soil_part) = rr.deposit(eroded_soil, eroded_soil_shrub)
    assert_array_almost_equal(burnt_shrubs_neigh, np.array([1, 2, 3, 4, 5]))
    (resource_adjusted,
     eligible_locs_to_adj_neigh,
     elig_locs,
     sed_to_borrow) = rr.re_adjust_resource()
    assert_almost_equal(resource_adjusted, 0.)
    veg_age = np.zeros(rr.grid.number_of_cells, dtype=int)
    veg_age = rr.initialize_Veg_age(V_age=veg_age)
    assert_array_equal(
        veg_age,
        np.zeros(rr.grid.number_of_cells, dtype=int))
    (veg_age, est_1, est_2, est_3, est_4, est_5) = rr.establish(veg_age)
    assert_array_equal(
        grid.at_cell["vegetation__plant_functional_type"],
        np.array([4, 1, 3, 3, 3, 1])
    )
    (veg_age, Pmor_age, Pmor_age_ws) = rr.mortality(veg_age)
    assert_array_equal(
        grid.at_cell["vegetation__plant_functional_type"],
        np.array([4, 1, 0, 0, 0, 1])
    )
    veg_age = rr.update_Veg_age(veg_age)
    assert_array_equal(
        veg_age, np.array([1, 1, 0, 0, 0, 1]))
