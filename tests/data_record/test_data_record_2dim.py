"""
Unit tests for landlab.data_record.data_record.DataRecord
Dimension = time and item_id

Last updated 8/24/2018

"""

import numpy as np

from landlab import RasterModelGrid

grid = RasterModelGrid((3, 3))
shape = (3, 3)
time = [0.0]
my_items3 = {
    "grid_element": np.array([["node"], ["link"]]),
    "element_id": np.array([[1], [3]]),
}
my_data_vars = {
    "mean_elevation": (["time"], [110.0]),
    "item_size": (["item_id", "time"], np.array([[0.3], [0.4]])),
}
#    return DataRecord(grid, time=time, items=my_items3, data_vars=my_data_vars )


def test_dr_name(dr_2dim):
    assert dr_2dim._name == "DataRecord"


def test_grid_shape(dr_2dim):
    assert dr_2dim._grid.number_of_node_rows == shape[0]
    assert dr_2dim._grid.number_of_node_columns == shape[1]


def test_permitted_locations(dr_2dim):
    assert dr_2dim._permitted_locations == grid.groups


def test_coordinates(dr_2dim):
    assert len(dr_2dim.dataset.sizes) == 2
    assert list(dr_2dim.dataset.time.values) == time
    assert list(dr_2dim.time_coordinates) == time
    assert list(dr_2dim.dataset.item_id.values) == [0, 1]
    assert list(dr_2dim.item_coordinates) == [0, 1]
    # properties:
    assert dr_2dim.number_of_timesteps == 1
    assert dr_2dim.number_of_items == 2
    assert dr_2dim.earliest_time == 0.0
    assert dr_2dim.latest_time == 0.0
    assert np.isnan(dr_2dim.prior_time)
    dr_2dim.add_record(time=[10.0, 20.0], new_record={"new_data": ((), (120.0, 130.0))})
    assert dr_2dim.prior_time == 10.0


def test_variable_names(dr_2dim):
    assert sorted(dr_2dim.variable_names) == sorted(
        ["element_id", "grid_element", "item_size", "mean_elevation"]
    )


def test_add_record(dr_2dim):
    dr_2dim.add_record(
        time=[10.0],
        item_id=[1],
        new_item_loc={
            "grid_element": np.array([["cell"]]),
            "element_id": np.array([[0]]),
        },
    )
    dr_2dim.add_record(
        time=[20.0], new_record={"mean_elevation": (["time"], np.array([130.0]))}
    )
    assert (
        dr_2dim.dataset["grid_element"].values[1, 1],
        dr_2dim.dataset["mean_elevation"].values[2],
    ) == ("cell", 130.0)
    assert np.isnan(dr_2dim.dataset["element_id"].values[1, 2])
    dr_2dim.add_record(
        time=[10.0],
        item_id=[1],
        new_record={"size": (["item_id", "time"], np.array([[0.3]]))},
    )
    assert np.isnan(dr_2dim.dataset["size"].values[0, 0])


def test_add_item(dr_2dim):
    dr_2dim.add_item(
        time=[10.0],
        new_item={
            "grid_element": np.array([["node"], ["cell"]]),
            "element_id": np.array([[2], [0]]),
        },
        new_item_spec={"size": (["item_id"], [10, 5])},
    )
    assert (
        dr_2dim.dataset["grid_element"].values[3, 1],
        dr_2dim.dataset["element_id"].values[2, 1],
        dr_2dim.dataset["size"].values[3],
    ) == ("cell", 2.0, 5.0)


def test_get_data(dr_2dim):
    assert (
        dr_2dim.get_data(time=[0.0], item_id=[1], data_variable="grid_element")
        == "link"
    )
    assert dr_2dim.get_data(data_variable="mean_elevation") == [110.0]


def test_set_data(dr_2dim):
    dr_2dim.set_data(
        time=[0.0], item_id=[1], data_variable="grid_element", new_value="node"
    )
    dr_2dim.set_data(time=[0.0], data_variable="mean_elevation", new_value=150.0)
    assert all(dr_2dim.dataset["grid_element"].values == "node")
    assert dr_2dim.dataset["mean_elevation"].values[0] == 150.0


def test_ffill_grid_element_and_id(dr_2dim):
    dr_2dim.add_record(
        time=[20.0], new_record={"mean_elevation": (["time"], np.array([130.0]))}
    )
    dr_2dim.ffill_grid_element_and_id()
    assert dr_2dim.dataset["grid_element"].values[0, 0] == (
        dr_2dim.dataset["grid_element"].values[0, 1]
    )
    assert dr_2dim.dataset["element_id"].values[0, 0] == (
        dr_2dim.dataset["element_id"].values[0, 1]
    )
