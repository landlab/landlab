# -*- coding: utf-8 -*-
"""
Unit tests for landlab.data_record.data_record.DataRecord
Dimension = time and item_id

Last updated 8/24/2018

"""

import pytest
import numpy as np
from landlab import RasterModelGrid

grid = RasterModelGrid((3,3))
shape = (3,3)
time=[0.]
my_items3 = {'grid_element':np.array([['node'], ['link']]),
             'element_id': np.array([[1],[3]])}
my_data_vars = {'mean_elevation' : (['time'], [110.]),
                'item_size' : (['item_id', 'time'],
                               np.array([[0.3], [0.4]]))}
#    return DataRecord(grid, time=time, items=my_items3, data_vars=my_data_vars )


def test_dr_time_name(dr_2dim):
    assert dr_2dim._name == 'DataRecord'

def test_grid_shape(dr_2dim):
    assert dr_2dim._grid.number_of_node_rows == shape[0]
    assert dr_2dim._grid.number_of_node_columns == shape[1]

def test_permitted_locations(dr_2dim):
    assert dr_2dim.permitted_locations == grid.groups

def test_coordinates(dr_2dim):
    assert len(dr_2dim.dims) == 2
    assert list(dr_2dim.time.values) == time
    assert list(dr_2dim.time_coordinates) == time
    assert list(dr_2dim.item_id.values) == [0,1]
    assert list(dr_2dim.item_coordinates) == [0,1]
    # properties:
    assert dr_2dim.number_of_timesteps == 1
    assert dr_2dim.number_of_items == 2
    assert dr_2dim.earliest_time == 0.
    assert dr_2dim.latest_time == 0.
    assert np.isnan(dr_2dim.prior_time)

def test_variable_names(dr_2dim):
    assert dr_2dim.variable_names == ['grid_element', 'element_id',
                                      'mean_elevation', 'item_size']