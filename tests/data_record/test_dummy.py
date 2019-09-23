import numpy as np
import pytest

from landlab import RasterModelGrid
from landlab.data_record import DataRecord


def test_ok_dummy():
    grid = RasterModelGrid((3, 3))
    element_id = [0, 0, 0, -9999, 1, 2, 3, 4, 5]
    volumes = [4, 5, 1, 2, 3, 4, 5, 6, 7]
    dr = DataRecord(
        grid,
        dummy_elements={"link": [9999, 1234, -9999], "node": [9999, 1234, -9999]},
        items={"grid_element": "node", "element_id": np.array(element_id)},
        data_vars={"volumes": (["item_id"], np.array(volumes))},
    )

    dr.add_item(
        new_item={"grid_element": np.array(["node"]), "element_id": np.array([9999])},
        new_item_spec={"volumes": (["item_id"], [5])},
    )


@pytest.mark.parametrize("dmmy", [0, 8])
def test_bad_dummy_init(dmmy):
    grid = RasterModelGrid((3, 3))
    element_id = [0, 0, 0, -9999, 1, 2, 3, 4, 5]
    volumes = [4, 5, 1, 2, 3, 4, 5, 6, 7]
    with pytest.raises(ValueError):
        DataRecord(
            grid,
            dummy_elements={"node": [dmmy]},
            items={"grid_element": "node", "element_id": np.array(element_id)},
            data_vars={"volumes": (["item_id"], np.array(volumes))},
        )


@pytest.mark.parametrize("dmmy", [-100, 9, 999])
def test_add_bad_dummy(dmmy):
    grid = RasterModelGrid((3, 3))
    element_id = [0, 0, 0, -9999, 1, 2, 3, 4, 5]
    volumes = [4, 5, 1, 2, 3, 4, 5, 6, 7]
    dr = DataRecord(
        grid,
        dummy_elements={"link": [9999, 1234, -9999], "node": [9999, 1234, -9999]},
        items={"grid_element": "node", "element_id": np.array(element_id)},
        data_vars={"volumes": (["item_id"], np.array(volumes))},
    )
    with pytest.raises(ValueError):
        dr.add_item(
            new_item={
                "grid_element": np.array(["node"]),
                "element_id": np.array([dmmy]),
            },
            new_item_spec={"volumes": (["item_id"], [5])},
        )
