#!/usr/bin/env python3
"""
Created on Mon Apr 30 09:17:36 2018

@author: barnhark
"""
import numpy as np
import pytest
from numpy.testing import assert_array_equal  # , assert_array_almost_equal

from landlab import RasterModelGrid
from landlab.bmi import wrap_as_bmi
from landlab.components import LithoLayers
from landlab.components import Lithology


def test_lithology_as_bmi():
    """Test Lithology can be wrapped with a BMI."""
    wrap_as_bmi(Lithology)


def test_bad_layer_method():
    """Test passing a bad name for the layer method."""
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1, 2, 4, 1]
    ids = [1, 2, 1, 2]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    with pytest.raises(ValueError):
        Lithology(mg, thicknesses, ids, attrs, layer_type="spam")


def test_no_topographic__elevation():
    """Test init with no topo__elevation."""
    mg = RasterModelGrid((3, 3))
    thicknesses = [1, 2, 4, 1]
    ids = [1, 2, 1, 2]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    with pytest.raises(ValueError):
        Lithology(mg, thicknesses, ids, attrs)


def test_thickness_ids_wrong_shape():
    """Test wrong size thickness and id shapes."""
    # first with thicknesses and IDs both as ndim = 1 arrays
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1, 2, 4, 1, 5]
    ids = [1, 2, 1, 2]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    with pytest.raises(ValueError):
        Lithology(mg, thicknesses, ids, attrs)

    # next as both as ndim = 2 arrays
    ones = np.ones(mg.number_of_nodes)
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1 * ones, 2 * ones, 4 * ones, 1 * ones, 5 * ones]
    ids = [1 * ones, 2 * ones, 1 * ones, 2 * ones]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    with pytest.raises(ValueError):
        Lithology(mg, thicknesses, ids, attrs)

    # now with thickness as ndim 2 and id as ndim 1
    ones = np.ones(mg.number_of_nodes)
    mg = RasterModelGrid((3, 3))
    thicknesses = [1 * ones, 2 * ones, 4 * ones, 1 * ones, 5 * ones]
    ids = [1, 2, 1, 2]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    with pytest.raises(ValueError):
        Lithology(mg, thicknesses, ids, attrs)


def test_thickness_ndim3():
    """Test too many ndim for thickness."""
    # next as both as ndim = 3 arrays
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    mg = RasterModelGrid((3, 3))
    ones = np.ones((mg.number_of_nodes, 2))
    mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1 * ones, 2 * ones, 4 * ones, 1 * ones, 5 * ones]
    ids = [1, 2, 1, 2]
    with pytest.raises(ValueError):
        Lithology(mg, thicknesses, ids, attrs)


def test_id_ndim3():
    """Test too many ndim for ids."""
    # next as both as ndim = 3 arrays
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    mg = RasterModelGrid((3, 3))
    ones = np.ones(mg.number_of_nodes)

    extra_ones = np.ones((mg.number_of_nodes, 2))
    mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1 * ones, 2 * ones, 4 * ones, 1 * ones, 5 * ones]
    ids = [1 * extra_ones, 2 * extra_ones, 1 * extra_ones, 2 * extra_ones]
    with pytest.raises(ValueError):
        Lithology(mg, thicknesses, ids, attrs)


def test_thickness_nodes_wrong_shape():
    """Test wrong size thickness and id shapes."""
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("topographic__elevation", at="node")
    ones = np.ones(mg.number_of_nodes + 1)
    thicknesses = [1 * ones, 2 * ones, 4 * ones, 1 * ones, 5 * ones]
    ids = [1 * ones, 2 * ones, 1 * ones, 2 * ones, 1 * ones]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    with pytest.raises(ValueError):
        Lithology(mg, thicknesses, ids, attrs)


def test_atts_lack_ids():
    """Test Lithology missing ID."""
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1, 2, 4, 1, 5]
    ids = [1, 2, 1, 2, 1]
    attrs = {"K_sp": {2: 0.0001}, "age": {1: 100, 2: 300}}
    with pytest.raises(ValueError):
        Lithology(mg, thicknesses, ids, attrs)


def test_erode_to_zero_thickness():
    """Test that eroding Lithology to zero thickness raises an error."""
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1, 2, 4, 1, 5]
    ids = [1, 2, 1, 2, 1]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    lith = Lithology(mg, thicknesses, ids, attrs)
    with pytest.raises(ValueError):
        lith.add_layer(-100)


def test_deposit_with_no_rock_id():
    """Test that adding a deposit to Lithology with no id raises an error."""
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1, 2, 4, 1, 5]
    ids = [1, 2, 1, 2, 1]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    lith = Lithology(mg, thicknesses, ids, attrs)
    with pytest.raises(ValueError):
        lith.add_layer(100)


def test_deposit_with_bad_rock_id():
    """Test that adding a deposit to Lithology with no id raises an error."""
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1, 2, 4, 1, 5]
    ids = [1, 2, 1, 2, 1]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    lith = Lithology(mg, thicknesses, ids, attrs)
    with pytest.raises(ValueError):
        lith.add_layer(100, rock_id=3)

    ones = np.ones(mg.number_of_nodes)
    new_ids = [0, 1, 3, 4, 0, 1, 0, 1, 5]
    with pytest.raises(ValueError):
        lith.add_layer(ones, rock_id=new_ids)


def test_adding_existing_attribute():
    """Test adding an existing attribute."""
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1, 2, 4, 1, 5]
    ids = [1, 2, 1, 2, 1]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    lith = Lithology(mg, thicknesses, ids, attrs)

    new_attr = {"K_sp": {1: 0.001, 2: 0.0001}}

    with pytest.raises(ValueError):
        lith.add_property(new_attr)


def test_adding_new_attribute_missing_rock_id():
    """Test adding an new attribute missing an existing rock id."""
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1, 2, 4, 1, 5]
    ids = [1, 2, 1, 2, 1]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    lith = Lithology(mg, thicknesses, ids, attrs)

    new_attr = {"D": {2: 0.0001}}

    with pytest.raises(ValueError):
        lith.add_property(new_attr)


def test_adding_new_attribute_extra_rock_id():
    """Test adding an new attribute with an extra rock id."""
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1, 2, 4, 1, 5]
    ids = [1, 2, 1, 2, 1]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    lith = Lithology(mg, thicknesses, ids, attrs)

    new_attr = {"D": {1: 0.001, 2: 0.0001, 3: 5.3}}

    with pytest.raises(ValueError):
        lith.add_property(new_attr)


def test_adding_new_id_existing_rock_type():
    """Test adding an rock type that already exists."""
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1, 2, 4, 1, 5]
    ids = [1, 2, 1, 2, 1]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    lith = Lithology(mg, thicknesses, ids, attrs)

    new_attr = {"K_sp": {1: 0.001, 5: 0.0001}}

    with pytest.raises(ValueError):
        lith.add_rock_type(new_attr)


def test_adding_new_id_extra_attribute():
    """Test adding an new rock type with an extra attribute."""
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1, 2, 4, 1, 5]
    ids = [1, 2, 1, 2, 1]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    lith = Lithology(mg, thicknesses, ids, attrs)

    new_attr = {"K_sp": {4: 0.001, 5: 0.0001}, "D": {4: 0.001, 5: 0.0001}}

    with pytest.raises(ValueError):
        lith.add_rock_type(new_attr)


def test_adding_new_id_missing_attribute():
    """Test adding an new rock type with an extra attribute."""
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1, 2, 4, 1, 5]
    ids = [1, 2, 1, 2, 1]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    lith = Lithology(mg, thicknesses, ids, attrs)
    new_attr = {"D": {4: 0.001, 5: 0.0001}}
    with pytest.raises(ValueError):
        lith.add_rock_type(new_attr)


def test_updating_attribute_that_doesnt_exist():
    """Test updating an attribute that doesn't exist."""
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1, 2, 4, 1, 5]
    ids = [1, 2, 1, 2, 1]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    lith = Lithology(mg, thicknesses, ids, attrs)
    with pytest.raises(ValueError):
        lith.update_rock_properties("spam", 1, 4)


def test_updating_rock_type_that_doesnt_exist():
    """Test adding an new rock type with an extra attribute."""
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1, 2, 4, 1, 5]
    ids = [1, 2, 1, 2, 1]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    lith = Lithology(mg, thicknesses, ids, attrs)
    with pytest.raises(ValueError):
        lith.update_rock_properties("K_sp", 3, 4)


def test_run_one_step_deposit_no_id_raises_error():
    """Test that giving the run one step method a deposit with no id raises an error."""
    mg = RasterModelGrid((3, 3))
    z = mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1, 2, 4, 1, 5]
    ids = [1, 2, 1, 2, 1]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    lith = Lithology(mg, thicknesses, ids, attrs)
    z += 1
    with pytest.raises(ValueError):
        lith.run_one_step()


def test_run_one_step_erodes_all_raises_error():
    """Test that eroding all material with the run one step method raises an error."""
    mg = RasterModelGrid((3, 3))
    z = mg.add_zeros("topographic__elevation", at="node")
    thicknesses = [1, 2, 4, 1, 5]
    ids = [1, 2, 1, 2, 1]
    attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
    lith = Lithology(mg, thicknesses, ids, attrs)
    z -= 30
    with pytest.raises(ValueError):
        lith.run_one_step()


def test_rock_block_xarray():
    """Test that the xarray method works as expected."""
    sample_depths = np.arange(0, 10, 1)

    mg = RasterModelGrid((3, 3))
    mg.add_zeros("topographic__elevation", at="node")
    layer_ids = np.tile([0, 1, 2, 3], 5)
    layer_elevations = 3.0 * np.arange(-10, 10)
    layer_elevations[-1] = layer_elevations[-2] + 100
    attrs = {"K_sp": {0: 0.0003, 1: 0.0001, 2: 0.0002, 3: 0.0004}}

    lith = LithoLayers(
        mg, layer_elevations, layer_ids, function=lambda x, y: x + y, attrs=attrs
    )
    ds = lith.rock_cube_to_xarray(sample_depths)
    expected_array = np.array(
        [
            [[3.0, 2.0, 2.0], [2.0, 2.0, 2.0], [2.0, 2.0, 1.0]],
            [[3.0, 3.0, 2.0], [3.0, 2.0, 2.0], [2.0, 2.0, 2.0]],
            [[3.0, 3.0, 3.0], [3.0, 3.0, 2.0], [3.0, 2.0, 2.0]],
            [[0.0, 3.0, 3.0], [3.0, 3.0, 3.0], [3.0, 3.0, 2.0]],
            [[0.0, 0.0, 3.0], [0.0, 3.0, 3.0], [3.0, 3.0, 3.0]],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 3.0], [0.0, 3.0, 3.0]],
            [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 3.0]],
            [[1.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            [[1.0, 1.0, 1.0], [1.0, 1.0, 0.0], [1.0, 0.0, 0.0]],
            [[2.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 0.0]],
        ]
    )

    assert_array_equal(ds.rock_type__id.values, expected_array)
