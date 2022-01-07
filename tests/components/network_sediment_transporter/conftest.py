import numpy as np
import pytest

from landlab.components import FlowDirectorSteepest, NetworkSedimentTransporter
from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid


@pytest.fixture()
def example_nmg():
    y_of_node = (0, 100, 200, 200, 300, 400, 400, 125)
    x_of_node = (0, 0, 100, -50, -100, 50, -150, -100)
    nodes_at_link = ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6))

    grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)

    # add variables to grid
    grid.add_field(
        "topographic__elevation", [0.0, 0.1, 0.3, 0.2, 0.35, 0.45, 0.5, 0.6], at="node"
    )
    grid.add_field(
        "bedrock__elevation", [0.0, 0.1, 0.3, 0.2, 0.35, 0.45, 0.5, 0.6], at="node"
    )
    grid.add_field(
        "reach_length",
        [10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0],
        at="link",
    )  # m

    grid.add_field("channel_width", 15 * np.ones(grid.size("link")), at="link")

    grid.add_field(
        "flow_depth", 0.01121871 * np.ones(grid.size("link")), at="link"
    )  # Why such an odd, low flow depth? Well, it doesn't matter... and oops.-AP

    return grid


@pytest.fixture()
def example_parcels(example_nmg):
    element_id = np.array(
        [0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6], dtype=int,
    )  # current link for each parcel

    element_id = np.expand_dims(element_id, axis=1)
    starting_link = np.squeeze(element_id)  # starting link for each parcel

    np.random.seed(0)

    time_arrival_in_link = np.random.rand(
        np.size(element_id), 1
    )  # time of arrival in each link -- larger numbers are younger
    volume = np.ones(np.shape(element_id))  # (m3) the volume of each parcel
    D = 0.05 * np.ones(
        np.shape(element_id)
    )  # (m) the diameter of grains in each parcel
    abrasion_rate = 0.0001 * np.ones(
        np.size(element_id)
    )  # 0 = no abrasion; abrasion rates are positive mass loss coefficients (mass loss / METER)
    active_layer = np.ones(
        np.shape(element_id)
    )  # 1 = active/surface layer; 0 = subsurface layer

    density = 2650 * np.ones(np.size(element_id))  # (kg/m3)

    location_in_link = np.zeros(
        np.shape(element_id)
    )  # [0 1], 0 is upstream end of link, 1 is downstream end

    D[0] = 0.075
    D[5] = 0.0001  # make one of them sand

    volume[2] = 0.3

    time = [0.0]  # probably not the sensible way to do this...

    items = {"grid_element": "link", "element_id": element_id}

    variables = {
        "starting_link": (["item_id"], starting_link),
        "abrasion_rate": (["item_id"], abrasion_rate),
        "density": (["item_id"], density),
        "time_arrival_in_link": (["item_id", "time"], time_arrival_in_link),
        "active_layer": (["item_id", "time"], active_layer),
        "location_in_link": (["item_id", "time"], location_in_link),
        "D": (["item_id", "time"], D),
        "volume": (["item_id", "time"], volume),
    }

    parcels = DataRecord(
        example_nmg,
        items=items,
        time=time,
        data_vars=variables,
        dummy_elements={"link": [NetworkSedimentTransporter.OUT_OF_NETWORK]},
    )
    return parcels


@pytest.fixture()
def example_flow_director(example_nmg):

    fd = FlowDirectorSteepest(example_nmg)
    fd.run_one_step()
    return fd
