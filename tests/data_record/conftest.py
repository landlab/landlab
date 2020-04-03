import numpy as np
import pytest

from landlab import RasterModelGrid
from landlab.data_record import DataRecord

grid = RasterModelGrid((3, 3))


@pytest.fixture
def dr_time():
    time = [0.0]
    data_vars = {"mean_elevation": (["time"], np.array([100]))}
    attrs = {"time_units": "y"}
    return DataRecord(grid=grid, time=time, data_vars=data_vars, attrs=attrs)


@pytest.fixture
def dr_item():
    my_items2 = {
        "grid_element": np.array(("node", "link"), dtype=str),
        "element_id": np.array([1, 3]),
    }
    return DataRecord(grid=grid, items=my_items2)


@pytest.fixture
def dr_2dim():
    time = [0.0]
    my_items3 = {
        "grid_element": np.array([["node"], ["link"]]),
        "element_id": np.array([[1], [3]]),
    }
    my_data_vars = {
        "mean_elevation": (["time"], [110.0]),
        "item_size": (["item_id", "time"], np.array([[0.3], [0.4]])),
    }
    return DataRecord(grid=grid, time=time, items=my_items3, data_vars=my_data_vars)


@pytest.fixture
def dr_nodim():
    return DataRecord(grid=grid)
