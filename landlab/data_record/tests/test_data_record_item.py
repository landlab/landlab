# -*- coding: utf-8 -*-
"""
Unit tests for landlab.data_record.data_record.DataRecord
Dimension = item_id

Last updated 8/24/2018

"""

import pytest
import numpy as np
from landlab import RasterModelGrid

grid = RasterModelGrid((3,3))
shape = (3,3)
my_items2 = {'grid_element': np.array(('node', 'link'), dtype=str),
             'element_id': np.array([1, 3])}

def test_dr_time_name(dr_item):
    assert dr_item._name == 'DataRecord'

def test_grid_shape(dr_item):
    assert dr_item._grid.number_of_node_rows == shape[0]
    assert dr_item._grid.number_of_node_columns == shape[1]

def test_permitted_locations(dr_item):
    assert dr_item.permitted_locations == grid.groups

def test_coordinates(dr_item):
    assert len(dr_item.dims) == 1
    assert list(dr_item.item_id.values) == [0, 1]
    assert list(dr_item.item_coordinates) == [0, 1]
    assert dr_item.number_of_items == len(my_items2['element_id'])
    with pytest.raises(AttributeError):
        dr_item.time
        dr_item.time_coordinates
        dr_item.number_of_timesteps
        dr_item.earliest_time
        dr_item.latest_time
        dr_item.prior_time

def test_variable_names(dr_item):
    assert dr_item.variable_names == ['grid_element', 'element_id']




#        assert np.isnan()