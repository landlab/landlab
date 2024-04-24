import numpy as np
from numpy.testing import assert_array_almost_equal

from landlab.components import FlowDirectorSteepest
from landlab.components import NetworkSedimentTransporter
from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid


def test_recycling():
    y_of_node = (0, 0, 0, 0)
    x_of_node = (0, 100, 200, 300)
    nodes_at_link = ((0, 1), (1, 2), (2, 3))

    nmg_constant_slope = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)

    # add variables to nmg
    nmg_constant_slope.at_node["topographic__elevation"] = [3.0, 2.0, 1.0, 0.0]
    nmg_constant_slope.at_node["bedrock__elevation"] = [3.0, 2.0, 1.0, 0.0]
    nmg_constant_slope.at_link["channel_slope"] = [0.001, 0.001, 0.001]
    nmg_constant_slope.at_link["reach_length"] = [100.0, 100.0, 100.0]  # m
    nmg_constant_slope.at_link["channel_width"] = 15 * np.ones(
        nmg_constant_slope.size("link")
    )
    nmg_constant_slope.at_link["flow_depth"] = 3 * np.ones(
        nmg_constant_slope.size("link")
    )

    flow_director = FlowDirectorSteepest(nmg_constant_slope)
    flow_director.run_one_step()

    time = [0.0]

    items = {"grid_element": "link", "element_id": np.array([[0]])}

    initial_volume = np.array([[1]])
    abrasion_rate = np.array([0])

    variables = {
        "starting_link": (["item_id"], np.array([0])),
        "abrasion_rate": (["item_id"], abrasion_rate),
        "density": (["item_id"], np.array([2650])),
        "time_arrival_in_link": (["item_id", "time"], np.array([[0.71518937]])),
        "active_layer": (["item_id", "time"], np.array([[1]])),
        "location_in_link": (["item_id", "time"], np.array([[0]])),
        "D": (["item_id", "time"], np.array([[0.05]])),
        "volume": (["item_id", "time"], initial_volume),
    }

    parcels = DataRecord(
        nmg_constant_slope,
        items=items,
        time=time,
        data_vars=variables,
        dummy_elements={"link": [NetworkSedimentTransporter.OUT_OF_NETWORK]},
    )

    nst = NetworkSedimentTransporter(
        nmg_constant_slope,
        parcels,
        flow_director,
        bed_porosity=0.03,
        g=9.81,
        fluid_density=1000,
        transport_method="WilcockCrowe",
    )

    dt = 60  # (seconds) 1 min timestep

    timesteps = 8

    for _ in range(0, (timesteps * dt), dt):
        # RECYCLE sediment: what left the network gets added back in at top.
        parcels.dataset.location_in_link.values[
            parcels.dataset.element_id.values
            == NetworkSedimentTransporter.OUT_OF_NETWORK
        ] = 0
        parcels.dataset.element_id.values[
            parcels.dataset.element_id.values
            == NetworkSedimentTransporter.OUT_OF_NETWORK
        ] = 0
        nst.run_one_step(dt)

        print("done with another timestep")
        print(parcels.dataset.element_id.values)

    Parcel_element_id = parcels.dataset.element_id.values

    Parcel_element_id_Should_Be = np.array(
        [[0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 0.0, 0.0]]
    )

    assert_array_almost_equal(
        Parcel_element_id_Should_Be, Parcel_element_id, decimal=-1
    )
