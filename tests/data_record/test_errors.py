"""

Unit tests for landlab.data_record.data_record.DataRecord
Test errors/exceptions

Last updated 8/27/2018


"""

import numpy as np
import pytest

from landlab import HexModelGrid
from landlab import RadialModelGrid
from landlab import RasterModelGrid
from landlab import VoronoiDelaunayGrid
from landlab.data_record import DataRecord


def test_misc():
    grid = RasterModelGrid((3, 3))

    # test bad dimension:
    with pytest.raises(ValueError):
        DataRecord(
            grid,
            time=[0.0],
            data_vars={"mean_elev": (["time"], [100.0]), "test": (["bad_dim"], [12])},
        )
    # should return ValueError('Data variable dimensions must be time and/or'
    #                              'item_id')
    # test bad time format:
    with pytest.raises(TypeError):
        DataRecord(grid=grid, time="bad_time")
    # should return TypeError: Time must be a list or an array of length 1

    # test bad datavars format:
    with pytest.raises(TypeError):
        DataRecord(grid, time=[0.0], data_vars=["not a dict"])
    # should return TypeError(('Data variables (data_vars) passed to'
    #                                 ' DataRecord must be a dictionary (see '
    #                                 'documentation for valid structure)'))

    # test bad items format:
    with pytest.raises(TypeError):
        DataRecord(
            grid,
            time=[0.0],
            items=["not a dict"],
            data_vars={"mean_elevation": (["time"], np.array([100]))},
            attrs={"time_units": "y"},
        )
    # Should return TypeError(('You must provide an ''items'' dictionary '
    #                                 '(see documentation for required format)'))
    #
    with pytest.raises(TypeError):
        # Test bad items keys
        DataRecord(
            grid,
            time=[0.0],
            items={
                "grid_element": np.array([["node"], ["link"]]),
                "bad_key": np.array([[1], [3]]),
            },
            data_vars={"mean_elevation": (["time"], np.array([100]))},
            attrs={"time_units": "y"},
        )
    # Should return TypeError(('You must provide an ''items'' dictionary '
    #                                  '(see documentation for required format)'))
    with pytest.raises(TypeError):
        # Test bad attrs
        DataRecord(
            grid,
            time=[0.0],
            items={
                "grid_element": np.array([["node"], ["link"]]),
                "element_id": np.array([[1], [3]]),
            },
            data_vars={"mean_elevation": (["time"], np.array([100]))},
            attrs=["not a dict"],
        )
    # Should return except AttributeError:
    #                 raise TypeError(('Attributes (attrs) passed to DataRecord'
    #                                 'must be a dictionary'))
    #
    # test bad loc and id:
    rmg = RasterModelGrid((3, 3))
    hmg = HexModelGrid((3, 2), spacing=1.0)
    radmg = RadialModelGrid(n_rings=1, nodes_in_first_ring=5, xy_of_center=(0.0, 0.0))
    vdmg = VoronoiDelaunayGrid(np.random.rand(25), np.random.rand(25))
    my_items_bad_loc = {
        "grid_element": np.array(["node", "bad_loc"]),
        "element_id": np.array([1, 3]),
    }
    my_items_bad_loc2 = {"grid_element": "bad_loc", "element_id": np.array([1, 3])}
    my_items_bad_loc3 = {
        "grid_element": np.array(["node", "node", "node"]),
        "element_id": np.array([1, 3]),
    }
    my_items_bad_id = {
        "grid_element": np.array(["node", "link"]),
        "element_id": np.array([1, 300]),
    }
    my_items_bad_id2 = {
        "grid_element": np.array(["node", "link"]),
        "element_id": np.array([1, -300]),
    }
    my_items_bad_id3 = {
        "grid_element": np.array(["node", "link"]),
        "element_id": np.array([1, 2.0]),
    }
    my_items_bad_id4 = {
        "grid_element": np.array(["node", "link"]),
        "element_id": np.array([1, 2, 3]),
    }
    with pytest.raises(ValueError):
        DataRecord(rmg, items=my_items_bad_loc)
    with pytest.raises(ValueError):
        DataRecord(hmg, items=my_items_bad_loc)
    with pytest.raises(ValueError):
        DataRecord(radmg, items=my_items_bad_loc)
    with pytest.raises(ValueError):
        DataRecord(vdmg, items=my_items_bad_loc)
    # should return ValueError(('One or more of the grid elements provided is/are'
    #                             ' not permitted location for this grid type.'))
    with pytest.raises(ValueError):
        DataRecord(rmg, items=my_items_bad_loc2)
    # should return ValueError: Location provided: bad_loc is not a permitted
    #  location for this grid type.
    with pytest.raises(ValueError):
        DataRecord(rmg, items=my_items_bad_loc3)
    # should return ValueError(('grid_element passed to DataRecord must be '
    #  ' the same length as the number of items or 1.'))
    with pytest.raises(ValueError):
        DataRecord(rmg, items=my_items_bad_id)
    with pytest.raises(ValueError):
        DataRecord(hmg, items=my_items_bad_id)
    with pytest.raises(ValueError):
        DataRecord(radmg, items=my_items_bad_id)
    with pytest.raises(ValueError):
        DataRecord(vdmg, items=my_items_bad_id)
    # should return ValueError(('An item residing at ' + at + ' has an '
    #  'element_id larger than the number of'+ at + 'on the grid.'))
    with pytest.raises(ValueError):
        DataRecord(rmg, items=my_items_bad_id2)
    #  should return ValueError(('An item residing at ' + at + ' has '
    # 'an element id below zero. This is not permitted.'))
    with pytest.raises(ValueError):
        DataRecord(rmg, items=my_items_bad_id3)
    # should return ValueError(('You have passed a non-integer element_id to '
    #                              'DataRecord, this is not permitted.'))
    with pytest.raises(ValueError):
        DataRecord(rmg, items=my_items_bad_id4)
    # should return ValueError(('The number of grid_element passed '
    #  ' to Datarecord must be 1 or equal  to the number of element_id '))


# TIME ONLY
def test_dr_time_bad_add_record(dr_time):
    with pytest.raises(TypeError):
        dr_time.add_record(time="bad_time", new_record=None)
    #  should return TypeError: You have passed a time that is not permitted,
    #  must be list or array
    with pytest.raises(KeyError):
        dr_time.add_record(time=[300.0], item_id=[0], new_record=None)


#  should return KeyError: 'This Datarecord does not hold items'
def test_dr_time_bad_get_data_(dr_time):
    with pytest.raises(KeyError):
        dr_time.get_data(time=[0.0], data_variable="bad_variable")
    # should return KeyError: "the variable 'bad_variable' is not in the Datarecord"
    with pytest.raises(KeyError):
        dr_time.get_data(time=[0.0], item_id=[0], data_variable="mean_elevation")
    # should return KeyError: 'This Datarecord does not hold items'
    with pytest.raises(TypeError):
        dr_time.get_data(time=0.0, data_variable="mean_elevation")


# TypeError('time must be a list or a 1-D array')
def test_dr_time_bad_set_data_(dr_time):
    with pytest.raises(KeyError):
        dr_time.set_data(time=[0.0], data_variable="bad_variable", new_value=105.0)
    # should return KeyError: "the variable 'bad_variable' is not in the Datarecord"
    with pytest.raises(IndexError):
        dr_time.set_data(time=[10.0], data_variable="mean_elevation", new_value=105.0)
    # shoulde return IndexError: 'The time you passed is not currently in the
    # Datarecord, you must the value you pass or first create the new time
    # coordinate using the add_record method.'
    with pytest.raises(KeyError):
        dr_time.set_data(
            time=[0.0], item_id=[0.0], data_variable="mean_elevation", new_value=105.0
        )
    # should return KeyError: 'This datarecord does not hold items'
    with pytest.raises(TypeError):
        dr_time.set_data(time=0.0, data_variable="mean_elevation", new_value=105.0)


# TypeError('time must be a list or a 1-d array')
def test_dr_time_bad_properties(dr_time):
    with pytest.raises(AttributeError):
        dr_time.number_of_items
    with pytest.raises(AttributeError):
        dr_time.item_coordinates


# ITEM ONLY
def test_dr_item_bad_add_record(dr_item):
    with pytest.raises(KeyError):
        dr_item.add_record(
            time=[0.0], item_id=[0], new_record={"new_var": ["new_data"]}
        )
    # should return KeyError: 'This Datarecord does not record time'
    with pytest.raises(ValueError):
        dr_item.add_record(item_id=[10], new_record={"new_var": ["new_data"]})
    # should return ValueError: There is no item with item_id 14, modify
    # the value(s) you pass as item_id or create a new item using the method
    # add_item.
    with pytest.raises(ValueError):
        dr_item.add_record(
            item_id=[0],
            new_item_loc={
                "grid_element": np.array(["cell"]),
                "element_id": np.array([0]),
            },
        )


# should raise ValueError('Use set_data to change the location of'
#                                     ' an item in this DataRecord')
def test_dr_item_bad_add_item(dr_item):
    with pytest.raises(KeyError):
        dr_item.add_item(
            time=[0.0],
            new_item={
                "grid_element": np.array(["node", "link"]),
                "element_id": np.array([1, 5]),
            },
            new_item_spec={"size": (["item_id"], np.array([0.2, 0.3]))},
        )


# should return KeyError: This Datarecord does not record time
def test_dr_item_bad_get_data(dr_item):
    with pytest.raises(KeyError):
        dr_item.get_data(time=[0.0], item_id=[0], data_variable="element_id")
    # should return KeyError: 'This Datarecord does not record time.'
    with pytest.raises(TypeError):
        dr_item.get_data(item_id=0, data_variable="element_id")
    # TypeError('item_id must be a list or a 1-D array')
    with pytest.raises(KeyError):
        dr_item.get_data(item_id=[0], data_variable="bad_variable")


# should return KeyError: "the variable 'bad_variable' is not in the Datarecord"
def test_dr_item_bad_set_data(dr_item):
    with pytest.raises(KeyError):
        dr_item.set_data(
            time=[0.0], item_id=[1], data_variable="element_id", new_value=2
        )
    # should return KeyError: 'This Datarecord does not record time.'
    with pytest.raises(KeyError):
        dr_item.get_data(time=[0.0], item_id=[0], data_variable="bad_variable")
    # should return KeyError: "the variable 'bad_variable' is not in the Datarecord"
    with pytest.raises(IndexError):
        dr_item.set_data(item_id=[3], data_variable="element_id", new_value=2)
    # should return IndexError: The item_id you passed does not exist in this
    # Datarecord
    with pytest.raises(TypeError):
        dr_item.set_data(item_id=3.0, data_variable="element_id", new_value=2)
    # should return TypeError: item_id must be a list or a 1-D array
    with pytest.raises(ValueError):
        dr_item.set_data(item_id=[1], data_variable="element_id", new_value=2.0)
    # ValueError('You have passed a non-integer element_id to DataRecord,
    # this is not permitted')
    with pytest.raises(ValueError):
        dr_item.set_data(item_id=[1], data_variable="element_id", new_value=-2)


# ValueError('You have passed an element id below zero. This is not permitted')
def test_dr_item_bad_properties(dr_item):
    with pytest.raises(AttributeError):
        dr_item.number_of_timesteps
    with pytest.raises(AttributeError):
        dr_item.time_coordinates
    with pytest.raises(AttributeError):
        dr_item.earliest_time
    with pytest.raises(AttributeError):
        dr_item.latest_time
    with pytest.raises(AttributeError):
        dr_item.prior_time


# ITEM AND TIME
def test_dr_2dim_bad_add_record(dr_2dim):
    with pytest.raises(TypeError):
        dr_2dim.add_record(
            time=10.0,
            item_id=[1],
            new_item_loc={
                "grid_element": np.array([["cell"]]),
                "element_id": np.array([[0]]),
            },
        )
    # should return TypeError: You have passed a time that is not permitted,
    # must be list or array.
    with pytest.raises(ValueError):
        dr_2dim.add_record(
            time=[10.0],
            item_id=[10],
            new_item_loc={
                "grid_element": np.array([["cell"]]),
                "element_id": np.array([[0]]),
            },
        )
    # should return ValueError: There is no item with item_id 10, modify the
    # value(s) you pass as item_id or create a new item using the method add_item.
    with pytest.raises(KeyError):
        dr_2dim.add_record(
            time=[10.0],
            item_id=[0],
            new_item_loc={"grid_element": np.array([["cell"]])},
        )
    # should return KeyError: 'You must provide a new_item_loc dictionnary with
    # both grid_element and element_id'
    with pytest.raises(TypeError):
        dr_2dim.add_record(
            time=[10.0],
            item_id=1,
            new_item_loc={
                "grid_element": np.array([["cell"]]),
                "element_id": np.array([[0]]),
            },
        )


# should return except TypeError:
# raise TypeError('item_id must be a list or a 1D array')
def test_dr_2dim_bad_add_item(dr_2dim):
    with pytest.raises(ValueError):
        dr_2dim.add_item(
            new_item={
                "grid_element": np.array([["node"], ["cell"]]),
                "element_id": np.array([[2], [0]]),
            },
            new_item_spec={"size": (["item_id"], [10, 5])},
        )
    # should return ValueError: The items previously defined in this Datarecord have
    # dimensions "time" and "item_id", please provide a "time" for the
    # new item(s)
    with pytest.raises(TypeError):
        dr_2dim.add_item(
            time=[10.0],
            new_item=("not a dict"),
            new_item_spec={"size": (["item_id"], [10, 5])},
        )
    # should return AttributeError raise TypeError: You must provide an new_item
    # dictionary (see documentation for required format)
    with pytest.raises(KeyError):
        dr_2dim.add_item(
            time=[10.0],
            new_item={"grid_element": np.array([["node"], ["cell"]])},
            new_item_spec={"size": (["item_id"], [10, 5])},
        )
    # should return AttributeError: You must provide an new_item dictionnary (see
    # documentation for required format)
    with pytest.raises(TypeError):
        dr_2dim.add_item(
            time=10.0,
            new_item={
                "grid_element": np.array([["node"], ["cell"]]),
                "element_id": np.array([[2], [0]]),
            },
            new_item_spec={"size": (["item_id"], [10, 5])},
        )


# should return TypeError: You have passed a time that is not permitted,
# must be list or array.
def test_dr_2dim_bad_get_data(dr_2dim):
    with pytest.raises(KeyError):
        dr_2dim.get_data(time=[0.0], item_id=1, data_variable="bad_variable")
    with pytest.raises(KeyError):
        dr_2dim.get_data(data_variable="bad_variable")
    # should return KeyError: "the variable 'bad_variable' is not in the Datarecord"
    with pytest.raises(IndexError):
        dr_2dim.get_data(time=[0.0], item_id=[10], data_variable="grid_element")
    # should return IndexError: The item_id you passed does not exist in this
    # Datarecord
    with pytest.raises(TypeError):
        dr_2dim.get_data(time=[0.0], item_id=0, data_variable="element_id")


# TypeError('item_id must be a list or a 1-D array')


def test_dr_2dim_bad_set_data(dr_2dim):
    with pytest.raises(KeyError):
        dr_2dim.set_data(
            time=[0.0], item_id=[1], data_variable="bad_variable", new_value="node"
        )
    with pytest.raises(ValueError):
        dr_2dim.set_data(
            time=[0.0], item_id=[1], data_variable="grid_element", new_value="cell"
        )
    with pytest.raises(IndexError):
        dr_2dim.set_data(
            time=[0.0], item_id=[1.0], data_variable="grid_element", new_value="node"
        )
    with pytest.raises(IndexError):
        dr_2dim.set_data(
            time=[130.0], item_id=[1], data_variable="grid_element", new_value="node"
        )
    with pytest.raises(TypeError):
        dr_2dim.set_data(
            time=[0.0], item_id=1, data_variable="mean_elevation", new_value=105.0
        )


# TypeError('item_id must be a list or a 1-d array')


# NO DIM
def test_dr_nodim_bad_get_data(dr_nodim):
    dr_nodim.add_record(new_record={"mean_elev": (100.0)})
    with pytest.raises(KeyError):
        dr_nodim.get_data(item_id=[0], data_variable="mean_elev")


# should return KeyError: 'This Datarecord does not hold items'
