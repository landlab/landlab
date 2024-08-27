"""
Unit tests for landlab.data_record.data_record.DataRecord
Dimension = time

Last updated 8/24/2018

"""

import numpy as np
import pytest

from landlab import RasterModelGrid

grid = RasterModelGrid((3, 3))
shape = (3, 3)
time = [0.0]
data_vars = {"mean_elevation": (["time"], np.array([100]))}
attrs = {"time_units": "y"}


def test_dr_time_name(dr_time):
    assert dr_time._name == "DataRecord"


def test_grid_shape(dr_time):
    assert dr_time._grid.number_of_node_rows == shape[0]
    assert dr_time._grid.number_of_node_columns == shape[1]


def test_permitted_locations(dr_time):
    assert dr_time._permitted_locations == grid.groups


def test_coordinates(dr_time):
    assert len(dr_time.dataset.sizes) == 1
    assert list(dr_time.dataset.time.values) == list(np.array(time))
    assert list(dr_time.time_coordinates) == list(np.array(time))
    # properties:
    assert dr_time.number_of_timesteps == 1
    assert dr_time.earliest_time == 0.0
    assert dr_time.latest_time == 0.0
    assert np.isnan(dr_time.prior_time)
    # no item_id coord:
    with pytest.raises(AttributeError):
        dr_time.item_id
    with pytest.raises(AttributeError):
        dr_time.item_coordinates
    with pytest.raises(AttributeError):
        dr_time.number_of_items


def test_variable_names(dr_time):
    assert dr_time.variable_names == ["mean_elevation"]


def test_add_record(dr_time):
    dr_time.add_record(
        time=[50.0], new_record={"mean_elevation": (["time"], np.array([120]))}
    )
    dr_time.add_record(
        time=[100.0], new_record={"new_variable": (["time"], ["new_data"])}
    )
    assert np.isnan(dr_time.dataset["mean_elevation"].values[2])


def test_get_data(dr_time):
    assert dr_time.get_data(time=[0.0], data_variable="mean_elevation") == 100.0
    assert dr_time.get_data(data_variable="mean_elevation") == [100]


def test_set_data(dr_time):
    dr_time.set_data(time=[0.0], data_variable="mean_elevation", new_value=105.0)
    assert dr_time.dataset["mean_elevation"].values[0] == 105.0
