import numpy as np
from numpy.testing import assert_array_equal

from landlab.components import FlowDirectorSteepest
from landlab.components import NetworkSedimentTransporter
from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid


def test_first_in_last_out():
    y_of_node = (0, 0, 0, 0, 0, 0)
    x_of_node = (0, 100, 200, 300, 400, 500)
    nodes_at_link = ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5))

    nmg_constant_slope = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)

    # add variables to nmg
    nmg_constant_slope.at_node["topographic__elevation"] = [
        5.0,
        4.0,
        3.0,
        2.0,
        1.0,
        0.0,
    ]
    nmg_constant_slope.at_node["bedrock__elevation"] = [5.0, 4.0, 3.0, 2.0, 1.0, 0.0]
    nmg_constant_slope.at_link["channel_slope"] = [0.001, 0.001, 0.001, 0.001, 0.001]
    nmg_constant_slope.at_link["flow_depth"] = [1.75, 1.75, 1.75, 1.75, 1.75]
    nmg_constant_slope.at_link["reach_length"] = [
        100.0,
        100.0,
        100.0,
        100.0,
        100.0,
    ]  # m
    nmg_constant_slope.at_link["channel_width"] = 15 * np.ones(
        nmg_constant_slope.size("link")
    )

    flow_director = FlowDirectorSteepest(nmg_constant_slope)
    flow_director.run_one_step()

    timesteps = 20

    time = [0.0]

    element_id = np.zeros(100, dtype=int)
    element_id = np.expand_dims(element_id, axis=1)

    items = {"grid_element": "link", "element_id": element_id}

    abrasion_rate = np.zeros(np.size(element_id))
    initial_volume = np.ones(np.shape(element_id))
    time_arrival_in_link = np.arange(0, 0.1, 0.001)
    time_arrival_in_link = np.expand_dims(time_arrival_in_link, axis=1)

    variables = {
        "starting_link": (["item_id"], np.squeeze(element_id)),
        "abrasion_rate": (["item_id"], abrasion_rate),
        "density": (["item_id"], 2650 * np.ones(np.size(element_id))),
        "time_arrival_in_link": (["item_id", "time"], time_arrival_in_link),
        "active_layer": (["item_id", "time"], initial_volume),
        "location_in_link": (["item_id", "time"], element_id),
        "D": (["item_id", "time"], initial_volume * 0.05),
        "volume": (["item_id", "time"], initial_volume),
    }

    hundred_boring_parcels = DataRecord(
        nmg_constant_slope,
        items=items,
        time=time,
        data_vars=variables,
        dummy_elements={"link": [NetworkSedimentTransporter.OUT_OF_NETWORK]},
    )

    nst = NetworkSedimentTransporter(
        nmg_constant_slope,
        hundred_boring_parcels,
        flow_director,
        bed_porosity=0.03,
        g=9.81,
        fluid_density=1000,
        transport_method="WilcockCrowe",
    )

    dt = 60 * 60  # (seconds) 1 hour timestep

    for _ in range(0, (timesteps * dt), dt):
        nst.run_one_step(dt)

    first_in_parcel = np.argmin(time_arrival_in_link)
    last_in_parcel = np.argmax(time_arrival_in_link)
    link_of_first_in_parcel = hundred_boring_parcels.dataset.element_id.values[
        first_in_parcel, :
    ]
    link_of_last_in_parcel = hundred_boring_parcels.dataset.element_id.values[
        last_in_parcel, :
    ]

    First_in_lags_behind = np.greater_equal(
        link_of_last_in_parcel, link_of_first_in_parcel
    )
    SO_TRUE = np.ones(np.shape(First_in_lags_behind), dtype=bool)

    assert_array_equal(SO_TRUE, First_in_lags_behind)
    # Asserts that the last-in parcel is consistently in either the same
    # link, or a farther downstream link than the first in parcel
