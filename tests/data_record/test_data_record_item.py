"""
Unit tests for landlab.data_record.data_record.DataRecord
Dimension = item_id

Last updated 8/24/2018

"""

import numpy as np
import pytest

from landlab import RasterModelGrid

grid = RasterModelGrid((3, 3))
shape = (3, 3)
my_items2 = {
    "grid_element": np.array(("node", "link"), dtype=str),
    "element_id": np.array([1, 3]),
}


def test_dr_item_name(dr_item):
    assert dr_item._name == "DataRecord"


def test_grid_shape(dr_item):
    assert dr_item._grid.number_of_node_rows == shape[0]
    assert dr_item._grid.number_of_node_columns == shape[1]


def test_permitted_locations(dr_item):
    assert dr_item._permitted_locations == grid.groups


def test_coordinates(dr_item):
    assert len(dr_item.dataset.sizes) == 1
    assert list(dr_item.dataset.item_id.values) == [0, 1]
    assert list(dr_item.item_coordinates) == [0, 1]
    assert dr_item.number_of_items == len(my_items2["element_id"])
    with pytest.raises(AttributeError):
        dr_item.dataset.time
    with pytest.raises(AttributeError):
        dr_item.time_coordinates
    with pytest.raises(AttributeError):
        dr_item.number_of_timesteps
    with pytest.raises(AttributeError):
        dr_item.earliest_time
    with pytest.raises(AttributeError):
        dr_item.latest_time
    with pytest.raises(AttributeError):
        dr_item.prior_time


def test_variable_names(dr_item):
    assert sorted(dr_item.variable_names) == sorted(["grid_element", "element_id"])


def test_add_item(dr_item):
    dr_item.add_item(
        new_item={
            "grid_element": np.array(["node", "node"]),
            "element_id": np.array([4, 4]),
        },
        new_item_spec={"size": (["item_id"], [10, 5])},
    )
    assert (
        dr_item.dataset["grid_element"].values[3],
        dr_item.dataset["element_id"].values[3],
        dr_item.dataset["size"].values[3],
    ) == ("node", 4.0, 5.0)


def test_get_data(dr_item):
    assert dr_item.get_data(item_id=[1], data_variable="grid_element") == "link"
    assert dr_item.get_data(data_variable="element_id")[1] == 3


def test_set_data(dr_item):
    dr_item.set_data(item_id=[1], data_variable="element_id", new_value=2)
    assert dr_item.dataset["element_id"].values[1] == 2
