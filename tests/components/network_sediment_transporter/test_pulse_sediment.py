# -*- coding: utf-8 -*-
import numpy as np
from numpy.testing import assert_array_almost_equal

from landlab.components import FlowDirectorSteepest, NetworkSedimentTransporter
from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid


def test_add_pulse():

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

    nst.run_one_step(dt)

    # ONE TIMESTEP BEFORE PULSE
    # TIMESTEP 1 should have NANS.

    num_pulse_parcels = 2

    newpar_element_id = np.zeros(num_pulse_parcels, dtype=int)
    newpar_element_id = np.expand_dims(newpar_element_id, axis=1)

    new_starting_link = np.squeeze(newpar_element_id)

    np.random.seed(0)

    new_time_arrival_in_link = nst._time * np.ones(np.shape(newpar_element_id))

    new_volume = 0.5 * np.ones(
        np.shape(newpar_element_id)
    )  # (m3) the volume of each parcel

    new_lithology = ["pulse_material"] * np.size(
        newpar_element_id
    )  # a lithology descriptor for each parcel

    new_active_layer = np.ones(
        np.shape(newpar_element_id)
    )  # 1 = active/surface layer; 0 = subsurface layer

    new_density = 2650 * np.ones(np.size(newpar_element_id))  # (kg/m3)

    new_location_in_link = np.random.rand(np.size(newpar_element_id), 1)

    new_abrasion_rate = 0 * np.ones(np.size(newpar_element_id))

    new_D = 0.03 * np.ones(np.shape(newpar_element_id))

    newpar_grid_elements = np.array(
        np.empty((np.shape(newpar_element_id)), dtype=object)
    )  # BUG: should be able to pass ["link"], but datarecord fills it into an incorrect array shape-- the length of parcels (NOT new parcels)
    newpar_grid_elements.fill("link")

    new_parcels = {
        "grid_element": newpar_grid_elements,
        "element_id": newpar_element_id,
    }

    new_variables = {
        "starting_link": (["item_id"], new_starting_link),
        "abrasion_rate": (["item_id"], new_abrasion_rate),
        "density": (["item_id"], new_density),
        "lithology": (["item_id"], new_lithology),
        "time_arrival_in_link": (["item_id", "time"], new_time_arrival_in_link),
        "active_layer": (["item_id", "time"], new_active_layer),
        "location_in_link": (["item_id", "time"], new_location_in_link),
        "D": (["item_id", "time"], new_D),
        "volume": (["item_id", "time"], new_volume),
    }

    parcels.add_item(
        time=[nst._time], new_item=new_parcels, new_item_spec=new_variables
    )

    nst.run_one_step(dt)

    print(parcels.dataset.element_id.values)
    Parcel_element_id = parcels.dataset.element_id.values

    Parcel_element_id_Should_Be = np.array(
        [[0.0, 0.0, 0.0], [np.nan, 0.0, 0.0], [np.nan, 0.0, 0.0]]
    )

    assert_array_almost_equal(
        Parcel_element_id_Should_Be, Parcel_element_id, decimal=-1
    )
