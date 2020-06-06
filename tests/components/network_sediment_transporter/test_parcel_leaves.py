import numpy as np
import pytest

from landlab.components import NetworkSedimentTransporter
from landlab.data_record import DataRecord


def test_parcel_leaves(example_nmg, example_flow_director):

    example_nmg.at_link["reach_length"] = ([10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0],)

    time = [0.0]

    items = {"grid_element": "link", "element_id": np.array([[6]])}

    initial_volume = np.array([[1]])
    abrasion_rate = np.array([0.0001])

    variables = {
        "starting_link": (["item_id"], np.array([6])),
        "abrasion_rate": (["item_id"], abrasion_rate),
        "density": (["item_id"], np.array([2650])),
        "time_arrival_in_link": (["item_id", "time"], np.array([[0.71518937]])),
        "active_layer": (["item_id", "time"], np.array([[1]])),
        "location_in_link": (["item_id", "time"], np.array([[0]])),
        "D": (["item_id", "time"], np.array([[0.05]])),
        "volume": (["item_id", "time"], initial_volume),
    }

    one_parcel = DataRecord(
        example_nmg,
        items=items,
        time=time,
        data_vars=variables,
        dummy_elements={"link": [NetworkSedimentTransporter.OUT_OF_NETWORK]},
    )

    timesteps = 5

    example_nmg.at_link["flow_depth"] = example_nmg.at_link["flow_depth"] * 20

    nst = NetworkSedimentTransporter(
        example_nmg,
        one_parcel,
        example_flow_director,
        bed_porosity=0.03,
        g=9.81,
        fluid_density=1000,
        transport_method="WilcockCrowe",
    )

    dt = 60 * 60 * 24  # (seconds) daily timestep

    # Parcel should exit, then runtime error is raised
    with pytest.raises(RuntimeError):
        for t in range(timesteps):
            nst.run_one_step(dt)
