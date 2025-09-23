import numpy as np
import pytest

from landlab import RasterModelGrid
from landlab.components import DepressionFinderAndRouter
from landlab.components import FlowAccumulator
from landlab.grid.create_network import AtMostNodes
from landlab.grid.create_network import network_grid_from_raster


@pytest.fixture
def grid():
    """small watershed used for tests"""
    mg = RasterModelGrid((8, 7), 10)
    mg.at_node["topographic__elevation"] = np.array(
        [
            5,
            3,
            2,
            1,
            2,
            3,
            5,
            7,
            4,
            1.5,
            3,
            3,
            4,
            7,
            7.5,
            5,
            4,
            2,
            4,
            5,
            7.5,
            7.5,
            5,
            3,
            5,
            6,
            6.5,
            9,
            9,
            5.5,
            4,
            4.1,
            5.5,
            8,
            10,
            10.5,
            6,
            8,
            7,
            7,
            7,
            11,
            12,
            7,
            9,
            8,
            9,
            10,
            12,
            14,
            12,
            11.5,
            11,
            11.5,
            12,
            14,
        ]
    )
    mg.at_node["id"] = mg.nodes
    # add flow direction and drainage area and correct for any pits
    fr = FlowAccumulator(mg, "topographic__elevation", flow_director="D8")  # add
    fr.run_one_step()
    df_4 = DepressionFinderAndRouter(mg)
    df_4.map_depressions()
    return mg


@pytest.fixture
def nmgrid_f(grid):
    """fine-scale network model grid representation of watershed"""
    nmg_f = network_grid_from_raster(
        grid,
        reducer=AtMostNodes(count=4),
        minimum_channel_threshold=200,  # upstream drainage area to truncate network, in m^2
        include=["drainage_area", "topographic__elevation"],
    )
    nmg_f.at_link["drainage_area"] = nmg_f.map_mean_of_link_nodes_to_link(
        "drainage_area"
    )
    return nmg_f


@pytest.fixture
def nmgrid_c(grid):
    """coarse-scale network model grid representation of watershed"""
    nmg_c = network_grid_from_raster(
        grid,
        reducer=AtMostNodes(count=2),
        minimum_channel_threshold=300,
        include=["drainage_area", "topographic__elevation"],
    )
    nmg_c.at_link["drainage_area"] = nmg_c.map_mean_of_link_nodes_to_link(
        "drainage_area"
    )
    return nmg_c
