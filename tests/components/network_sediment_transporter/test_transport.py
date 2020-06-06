import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal

from landlab.components import FlowDirectorSteepest, NetworkSedimentTransporter
from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid


def test_no_flow_no_transport(example_nmg, example_parcels, example_flow_director):

    timesteps = 3
    example_nmg.at_link["flow_depth"] = 0 * np.ones(example_nmg.size("link"))

    nst = NetworkSedimentTransporter(
        example_nmg,
        example_parcels,
        example_flow_director,
        bed_porosity=0.03,
        g=9.81,
        fluid_density=1000,
    )

    dt = 60 * 60 * 24  # (seconds) daily timestep

    for t in range(0, (timesteps * dt), dt):
        nst.run_one_step(dt)

        # Need to define original_node_elev after a timestep has passed.
        if t / (60 * 60 * 24) == 1:
            original_node_elev = example_nmg.at_node["topographic__elevation"]

    # No flow? Parcels should stay in same locations in links
    assert_array_almost_equal(
        example_parcels.dataset.location_in_link[:, 0],
        example_parcels.dataset.location_in_link[:, 1],
    )

    # No flow? Total parcel volume should stay the same
    assert_array_equal(
        np.sum(example_parcels.dataset["volume"].values, axis=0)[0],
        np.sum(example_parcels.dataset["volume"].values, axis=0)[1],
    )

    # No flow? Elevations should stay the same
    assert_array_equal(
        original_node_elev, example_nmg.at_node["topographic__elevation"]
    )


def test_defined_parcel_transport():

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
    nmg_constant_slope.at_link["flow_depth"] = 2 * np.ones(
        nmg_constant_slope.size("link")
    )

    flow_director = FlowDirectorSteepest(nmg_constant_slope)
    flow_director.run_one_step()

    timesteps = 11

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

    one_parcel = DataRecord(
        nmg_constant_slope,
        items=items,
        time=time,
        data_vars=variables,
        dummy_elements={"link": [NetworkSedimentTransporter.OUT_OF_NETWORK]},
    )

    nst = NetworkSedimentTransporter(
        nmg_constant_slope,
        one_parcel,
        flow_director,
        bed_porosity=0.03,
        g=9.81,
        fluid_density=1000,
        transport_method="WilcockCrowe",
    )

    dt = 60  # (seconds) 1 min timestep

    distance_traveled = np.arange(0.0, timesteps)
    # distance_traveled = np.arange(0.,timesteps)

    for t in range(0, (timesteps * dt), dt):
        nst.run_one_step(dt)
        distance_traveled[np.int(t / dt)] = nst._distance_traveled_cumulative
    # NEED TO CALCULATE THINGS HERE.
    # Transport distance should match?
    Distance_Traveled_Should_Be = [
        21.61998527,
        43.23997054,
        64.85995581,
        86.47994109,
        108.09992636,
        129.71991163,
        151.3398969,
        172.95988217,
        194.57986744,
        216.19985272,
        237.81983799,
    ]

    assert_array_almost_equal(
        Distance_Traveled_Should_Be, distance_traveled, decimal=-1
    )
