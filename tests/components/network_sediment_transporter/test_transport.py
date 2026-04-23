import numpy as np
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

from landlab.components import FlowDirectorSteepest
from landlab.components import NetworkSedimentTransporter
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

    reach_length = 100.0

    # add variables to nmg
    nmg_constant_slope.at_node["topographic__elevation"] = [3.0, 2.0, 1.0, 0.0]
    nmg_constant_slope.at_node["bedrock__elevation"] = [3.0, 2.0, 1.0, 0.0]
    nmg_constant_slope.at_link["channel_slope"] = [0.01, 0.01, 0.01]
    nmg_constant_slope.at_link["reach_length"] = reach_length * np.ones(
        nmg_constant_slope.size("link")
    )  # m
    nmg_constant_slope.at_link["channel_width"] = 15 * np.ones(
        nmg_constant_slope.size("link")
    )
    nmg_constant_slope.at_link["flow_depth"] = 0.5 * np.ones(
        nmg_constant_slope.size("link")
    )

    flow_director = FlowDirectorSteepest(nmg_constant_slope)
    flow_director.run_one_step()

    timesteps = 3

    time = [0.0]

    initial_volume = np.array([[1]])
    abrasion_rate = np.array([0])
    D = 0.05
    rhos = 2650
    parcel_grid_element = 1  # place parcel on link 1

    items = {"grid_element": "link", "element_id": np.array([[parcel_grid_element]])}

    variables = {
        "starting_link": (["item_id"], np.array([0])),
        "abrasion_rate": (["item_id"], abrasion_rate),
        "density": (["item_id"], np.array([rhos])),
        "time_arrival_in_link": (["item_id", "time"], np.array([[0.71518937]])),
        "active_layer": (["item_id", "time"], np.array([[1]])),
        "location_in_link": (["item_id", "time"], np.array([[0]])),
        "D": (["item_id", "time"], np.array([[D]])),
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
        assert len(nst._distance_traveled_cumulative) == 1
        distance_traveled[int(t / dt)] = nst._distance_traveled_cumulative[0]

    # TEST A: "by hand" calculation of W&C transport for link 1 (where parcel resides)
    S = np.abs(
        (
            nmg_constant_slope.at_node["bedrock__elevation"][2]
            - nmg_constant_slope.at_node["bedrock__elevation"][1]
        )
        / (x_of_node[2] - x_of_node[1])
    )
    h = nmg_constant_slope.at_link["flow_depth"][parcel_grid_element]
    w = nmg_constant_slope.at_link["channel_width"][parcel_grid_element]
    g = 9.81
    rho = 1000

    tau = rho * g * h * S
    # tau_star_rm = rho/

    tau_star_rm = 0.036  # no sand reference Shields

    tau_rm = tau_star_rm * (rhos - rho) * g * D  # reference *shear*

    tau_ri = tau_rm  # uniform D hiding function goes to 1

    phi = tau / tau_ri  #

    # phi > 1.35
    Wstar = 14 * (1 - 0.894 / np.sqrt(phi)) ** 4.5

    ustar = np.sqrt(tau / rho)

    # qb calc per unit width, should match nst._qb
    qb_calc = (Wstar * (ustar**3)) / (
        (rhos / rho - 1) * g
    )  # single parcel, so qb = qbi

    # TEST B: transport via sum of parcel motion
    # nst._qb should match sum(parcel velocities * parcel volumes) / link length

    Qb_parcelmotion = (
        np.sum(nst._pvelocity * initial_volume)
        / nmg_constant_slope.at_link["reach_length"][parcel_grid_element]
    )

    qb_parcelmotion = Qb_parcelmotion / w

    # ASSERT these three are equal

    assert_array_almost_equal(nst._qb[parcel_grid_element], qb_calc, decimal=7)

    assert_array_almost_equal(nst._qb[parcel_grid_element], qb_parcelmotion, decimal=7)


def test_defined_parcel_transport_WCd50():
    y_of_node = (0, 0, 0, 0)
    x_of_node = (0, 100, 200, 300)
    nodes_at_link = ((0, 1), (1, 2), (2, 3))

    nmg_constant_slope = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)

    reach_length = 100.0

    # add variables to nmg
    nmg_constant_slope.at_node["topographic__elevation"] = [3.0, 2.0, 1.0, 0.0]
    nmg_constant_slope.at_node["bedrock__elevation"] = [3.0, 2.0, 1.0, 0.0]
    nmg_constant_slope.at_link["channel_slope"] = [0.01, 0.01, 0.01]
    nmg_constant_slope.at_link["reach_length"] = reach_length * np.ones(
        nmg_constant_slope.size("link")
    )  # m
    nmg_constant_slope.at_link["channel_width"] = 15 * np.ones(
        nmg_constant_slope.size("link")
    )
    nmg_constant_slope.at_link["flow_depth"] = 0.5 * np.ones(
        nmg_constant_slope.size("link")
    )

    flow_director = FlowDirectorSteepest(nmg_constant_slope)
    flow_director.run_one_step()

    timesteps = 3

    time = [0.0]

    initial_volume = np.array([[1]])
    abrasion_rate = np.array([0])
    D = 0.05
    rhos = 2650
    parcel_grid_element = 1  # place parcel on link 1

    items = {"grid_element": "link", "element_id": np.array([[parcel_grid_element]])}

    variables = {
        "starting_link": (["item_id"], np.array([0])),
        "abrasion_rate": (["item_id"], abrasion_rate),
        "density": (["item_id"], np.array([rhos])),
        "time_arrival_in_link": (["item_id", "time"], np.array([[0.71518937]])),
        "active_layer": (["item_id", "time"], np.array([[1]])),
        "location_in_link": (["item_id", "time"], np.array([[0]])),
        "D": (["item_id", "time"], np.array([[D]])),
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
        transport_method="WilcockCroweD50",
    )

    dt = 60  # (seconds) 1 min timestep

    distance_traveled = np.arange(0.0, timesteps)
    # distance_traveled = np.arange(0.,timesteps)

    for t in range(0, (timesteps * dt), dt):
        nst.run_one_step(dt)
        assert len(nst._distance_traveled_cumulative) == 1
        distance_traveled[int(t / dt)] = nst._distance_traveled_cumulative[0]

    # TEST A: "by hand" calculation of W&C transport for link 1 (where parcel resides)
    S = np.abs(
        (
            nmg_constant_slope.at_node["bedrock__elevation"][2]
            - nmg_constant_slope.at_node["bedrock__elevation"][1]
        )
        / (x_of_node[2] - x_of_node[1])
    )
    h = nmg_constant_slope.at_link["flow_depth"][parcel_grid_element]
    w = nmg_constant_slope.at_link["channel_width"][parcel_grid_element]
    g = 9.81
    rho = 1000

    tau = rho * g * h * S
    # tau_star_rm = rho/

    tau_star_rm = 0.036  # no sand reference Shields

    tau_rm = tau_star_rm * (rhos - rho) * g * D  # reference *shear*

    tau_ri = tau_rm  # uniform D hiding function goes to 1

    phi = tau / tau_ri  #

    # phi > 1.35
    Wstar = 14 * (1 - 0.894 / np.sqrt(phi)) ** 4.5

    ustar = np.sqrt(tau / rho)

    # qb calc per unit width, should match nst._qb
    qb_calc = (Wstar * (ustar**3)) / (
        (rhos / rho - 1) * g
    )  # single parcel, so qb = qbi

    # TEST B: transport via sum of parcel motion
    # nst._qb should match sum(parcel velocities * parcel volumes) / link length

    Qb_parcelmotion = (
        np.sum(nst._pvelocity * initial_volume)
        / nmg_constant_slope.at_link["reach_length"][parcel_grid_element]
    )

    qb_parcelmotion = Qb_parcelmotion / w

    # ASSERT these three are equal

    assert_array_almost_equal(nst._qb[parcel_grid_element], qb_calc, decimal=7)

    assert_array_almost_equal(nst._qb[parcel_grid_element], qb_parcelmotion, decimal=7)
