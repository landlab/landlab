import numpy as np
from numpy.testing import assert_array_almost_equal

from landlab.components import NetworkSedimentTransporter
from landlab.data_record import DataRecord


def test_abrasion(example_nmg, example_parcels, example_flow_director):
    time = [0.0]  # probably not the sensible way to do this...

    items = {"grid_element": "link", "element_id": np.array([[6], [6]])}

    initial_volume = np.array([[1], [1]])
    abrasion_rate = np.array([0.0001, 0])

    variables = {
        "starting_link": (["item_id"], np.array([6, 6])),
        "abrasion_rate": (["item_id"], abrasion_rate),
        "density": (["item_id"], np.array([2650, 90000])),
        "time_arrival_in_link": (["item_id", "time"], np.array([[0.5], [0]])),
        "active_layer": (["item_id", "time"], np.array([[1], [1]])),
        "location_in_link": (["item_id", "time"], np.array([[0], [0]])),
        "D": (["item_id", "time"], np.array([[0.05], [0.05]])),
        "volume": (["item_id", "time"], initial_volume),
    }

    two_parcels = DataRecord(
        example_nmg,
        items=items,
        time=time,
        data_vars=variables,
        dummy_elements={"link": [NetworkSedimentTransporter.OUT_OF_NETWORK]},
    )

    timesteps = 8

    example_nmg.at_link["flow_depth"] = (
        example_nmg.at_link["flow_depth"] * 5
    )  # high transport rate

    nst = NetworkSedimentTransporter(
        example_nmg,
        two_parcels,
        example_flow_director,
        bed_porosity=0.03,
        g=9.81,
        fluid_density=1000,
        transport_method="WilcockCrowe",
    )

    dt = 60 * 60 * 24  # (seconds) daily timestep

    for t in range(0, (timesteps * dt), dt):
        nst.run_one_step(dt)

    # Parcel volume should decrease according to abrasion rate
    volume_after_transport = np.squeeze(np.transpose(initial_volume)) * np.exp(
        nst._distance_traveled_cumulative * -abrasion_rate
    )

    assert_array_almost_equal(
        volume_after_transport, two_parcels.dataset.volume[0:2, -1]
    )
