import numpy as np
import pytest

from landlab import ExampleData
from landlab.components import FlowDirectorSteepest, NetworkSedimentTransporter
from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid
from landlab.io import read_shapefile


@pytest.fixture()
def synthetic():
    y_of_node = (0, 100, 200, 200, 300, 400, 400, 125)
    x_of_node = (0, 0, 100, -50, -100, 50, -150, -100)

    nodes_at_link = ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6))

    grid1 = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
    grid1.at_node["topographic__elevation"] = [
        0.0,
        0.08,
        0.25,
        0.15,
        0.25,
        0.4,
        0.8,
        0.8,
    ]
    grid1.at_node["bedrock__elevation"] = [0.0, 0.08, 0.25, 0.15, 0.25, 0.4, 0.8, 0.8]
    grid1.at_link["flow_depth"] = 2.5 * np.ones(grid1.number_of_links)
    grid1.at_link["channel_slope"] = np.ones(
        grid1.number_of_links
    )  # AP ?? Does it get overwritten?? check
    grid1.at_link["reach_length"] = 200 * np.ones(grid1.number_of_links)  # m
    grid1.at_link["channel_width"] = 1 * np.ones(grid1.number_of_links)

    # element_id is the link on which the parcel begins.
    element_id = np.repeat(np.arange(grid1.number_of_links), 30)
    element_id = np.expand_dims(element_id, axis=1)

    volume = 0.05 * np.ones(np.shape(element_id))  # (m3)
    active_layer = np.ones(np.shape(element_id))  # 1= active, 0 = inactive
    density = 2650 * np.ones(np.size(element_id))  # (kg/m3)
    abrasion_rate = 0 * np.ones(np.size(element_id))  # (mass loss /m)

    # Lognormal GSD
    medianD = 0.085  # m
    mu = np.log(medianD)
    sigma = np.log(2)  # assume that D84 = sigma*D50
    np.random.seed(0)
    D = np.random.lognormal(
        mu, sigma, np.shape(element_id)
    )  # (m) the diameter of grains in each parcel

    time_arrival_in_link = np.random.rand(np.size(element_id), 1)
    location_in_link = np.random.rand(np.size(element_id), 1)

    lithology = ["quartzite"] * np.size(element_id)

    variables = {
        "abrasion_rate": (["item_id"], abrasion_rate),
        "density": (["item_id"], density),
        "lithology": (["item_id"], lithology),
        "time_arrival_in_link": (["item_id", "time"], time_arrival_in_link),
        "active_layer": (["item_id", "time"], active_layer),
        "location_in_link": (["item_id", "time"], location_in_link),
        "D": (["item_id", "time"], D),
        "volume": (["item_id", "time"], volume),
    }

    items = {"grid_element": "link", "element_id": element_id}

    parcels1 = DataRecord(
        grid1,
        items=items,
        time=[0.0],
        data_vars=variables,
        dummy_elements={"link": [NetworkSedimentTransporter.OUT_OF_NETWORK]},
    )

    fd1 = FlowDirectorSteepest(grid1, "topographic__elevation")
    fd1.run_one_step()

    nst1 = NetworkSedimentTransporter(
        grid1,
        parcels1,
        fd1,
        bed_porosity=0.3,
        g=9.81,
        fluid_density=1000,
        transport_method="WilcockCrowe",
    )
    timesteps = 2  # total number of timesteps
    dt = 60 * 60 * 24 * 1  # length of timestep (seconds)
    for t in range(0, (timesteps * dt), dt):
        nst1.run_one_step(dt)

    return nst1


@pytest.fixture()
def methow():
    example_data_dir = ExampleData("io/shapefile", case="methow").base
    shp_file = example_data_dir / "MethowSubBasin.shp"
    points_shapefile = example_data_dir / "MethowSubBasin_Nodes_4.shp"

    grid2 = read_shapefile(
        shp_file,
        points_shapefile=points_shapefile,
        node_fields=["usarea_km2", "Elev_m"],
        link_fields=["usarea_km2", "Length_m"],
        link_field_conversion={
            "usarea_km2": "drainage_area",
            "Slope": "channel_slope",
            "Length_m": "reach_length",
        },
        node_field_conversion={
            "usarea_km2": "drainage_area",
            "Elev_m": "topographic__elevation",
        },
        threshold=0.01,
    )
    grid2.at_node["bedrock__elevation"] = grid2.at_node["topographic__elevation"].copy()
    grid2.at_link["channel_width"] = 1 * np.ones(grid2.number_of_links)
    grid2.at_link["flow_depth"] = 0.5 * np.ones(grid2.number_of_links)

    # element_id is the link on which the parcel begins.
    element_id = np.repeat(np.arange(grid2.number_of_links), 50)
    element_id = np.expand_dims(element_id, axis=1)

    volume = 1 * np.ones(np.shape(element_id))  # (m3)
    active_layer = np.ones(np.shape(element_id))  # 1= active, 0 = inactive
    density = 2650 * np.ones(np.size(element_id))  # (kg/m3)
    abrasion_rate = 0 * np.ones(np.size(element_id))  # (mass loss /m)

    # Lognormal GSD
    medianD = 0.15  # m
    mu = np.log(medianD)
    sigma = np.log(2)  # assume that D84 = sigma*D50
    np.random.seed(0)
    D = np.random.lognormal(
        mu, sigma, np.shape(element_id)
    )  # (m) the diameter of grains in each parcel

    time_arrival_in_link = np.random.rand(np.size(element_id), 1)
    location_in_link = np.random.rand(np.size(element_id), 1)

    lithology = ["quartzite"] * np.size(element_id)

    variables = {
        "abrasion_rate": (["item_id"], abrasion_rate),
        "density": (["item_id"], density),
        "lithology": (["item_id"], lithology),
        "time_arrival_in_link": (["item_id", "time"], time_arrival_in_link),
        "active_layer": (["item_id", "time"], active_layer),
        "location_in_link": (["item_id", "time"], location_in_link),
        "D": (["item_id", "time"], D),
        "volume": (["item_id", "time"], volume),
    }

    items = {"grid_element": "link", "element_id": element_id}

    parcels2 = DataRecord(
        grid2,
        items=items,
        time=[0.0],
        data_vars=variables,
        dummy_elements={"link": [NetworkSedimentTransporter.OUT_OF_NETWORK]},
    )

    fd2 = FlowDirectorSteepest(grid2, "topographic__elevation")
    fd2.run_one_step()

    nst2 = NetworkSedimentTransporter(
        grid2,
        parcels2,
        fd2,
        bed_porosity=0.3,
        g=9.81,
        fluid_density=1000,
        transport_method="WilcockCrowe",
    )
    timesteps = 2  # total number of timesteps
    dt = 60 * 60 * 24 * 1  # length of timestep (seconds)
    for t in range(0, (timesteps * dt), dt):
        nst2.run_one_step(dt)

    return nst2
