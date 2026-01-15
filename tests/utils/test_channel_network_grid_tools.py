import numpy as np
import pandas as pd
import pytest

import landlab.utils.channel_network_grid_tools as gt
from landlab import RasterModelGrid
from landlab.components import DepressionFinderAndRouter
from landlab.components import FlowAccumulator
from landlab.grid.create_network import AtMostNodes
from landlab.grid.create_network import network_grid_from_raster


def check_vals(vals, vals_e):
    np.testing.assert_allclose(vals, vals_e, atol=1e-4)


@pytest.fixture
def grid():
    """small watershed used for tests, drains to the south"""
    mg = RasterModelGrid((8, 7), 10)
    mg.at_node["topographic__elevation"] = np.array(
        [
            [5, 3, 2, 1, 2, 3, 5],
            [7, 4, 1.5, 3, 3, 4, 7],
            [7.5, 5, 4, 2, 4, 5, 7.5],
            [7.5, 5, 3, 5, 6, 6.5, 9],
            [9, 5.5, 4, 4.1, 5.5, 8, 10],
            [10.5, 6, 8, 7, 7, 7, 11],
            [12, 7, 9, 8, 9, 10, 12],
            [14, 12, 11.5, 11, 11.5, 12, 14],
        ]
    ).reshape(-1)
    mg.at_node["id"] = mg.nodes
    # add flow direction and drainage area and correct for any pits
    fr = FlowAccumulator(mg, "topographic__elevation", flow_director="D8")  # add
    fr.run_one_step()
    df_4 = DepressionFinderAndRouter(mg)
    df_4.map_depressions()
    return mg


@pytest.fixture
def grid_up():
    """small watershed used for tests, flipped upside down, drains to the north"""
    mg_up = RasterModelGrid((8, 7), 10)
    mg_up.at_node["topographic__elevation"] = np.flipud(
        np.array(
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
    )
    mg_up.at_node["id"] = mg_up.nodes
    # add flow direction and drainage area and correct for any pits
    fr = FlowAccumulator(mg_up, "topographic__elevation", flow_director="D8")  # add
    fr.run_one_step()
    df_4 = DepressionFinderAndRouter(mg_up)
    df_4.map_depressions()
    return mg_up


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
def nmgrid_f_up(grid_up):
    """fine-scale network model grid representation of watershed, flipped
    upside down"""
    nmg_f_up = network_grid_from_raster(
        grid_up,
        reducer=AtMostNodes(count=4),
        minimum_channel_threshold=200,  # upstream drainage area to truncate network, in m^2
        include=["drainage_area", "topographic__elevation"],
    )
    nmg_f_up.at_link["drainage_area"] = nmg_f_up.map_mean_of_link_nodes_to_link(
        "drainage_area"
    )
    return nmg_f_up


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


@pytest.fixture
def nmgrid_c_up(grid_up):
    """coarse-scale network model grid representation of watershed, flipped
    upside down"""
    nmg_c_up = network_grid_from_raster(
        grid_up,
        reducer=AtMostNodes(count=2),
        minimum_channel_threshold=300,
        include=["drainage_area", "topographic__elevation"],
    )
    nmg_c_up.at_link["drainage_area"] = nmg_c_up.map_mean_of_link_nodes_to_link(
        "drainage_area"
    )
    return nmg_c_up


@pytest.fixture
def nmgrid_s(grid):
    """single-channel-reach network model grid representation of watershed"""
    nmg_s = network_grid_from_raster(
        grid,
        reducer=AtMostNodes(count=2),
        minimum_channel_threshold=1000,
        include=["drainage_area", "topographic__elevation"],
    )
    nmg_s.at_link["drainage_area"] = nmg_s.map_mean_of_link_nodes_to_link(
        "drainage_area"
    )
    return nmg_s


class TestMapNmg1LinksToNmg2Links:

    def test_map_nmg1_links_to_nmg2_links_1(self, nmgrid_c, nmgrid_f):
        """map fine grid to coarse grid"""
        link_mapper = gt.map_nmg1_links_to_nmg2_links(nmgrid_c, nmgrid_f)
        v = pd.DataFrame.from_dict(link_mapper, orient="index").values
        v_e = np.array([[0], [3], [4]])
        check_vals(v, v_e)

    def test_map_nmg1_links_to_nmg2_links_2(self, nmgrid_c, nmgrid_f):
        """map coarse grid to fine grid"""
        link_mapper = gt.map_nmg1_links_to_nmg2_links(nmgrid_f, nmgrid_c)
        v = pd.DataFrame.from_dict(link_mapper, orient="index").values
        v_e = np.array([[0], [0], [0], [1], [2], [2], [1]])
        check_vals(v, v_e)

    def test_map_nmg1_links_to_nmg2_links_3(self, nmgrid_c_up, nmgrid_f_up):
        """map upside down fine grid to the upside-down coarse grids"""
        link_mapper_up = gt.map_nmg1_links_to_nmg2_links(nmgrid_c_up, nmgrid_f_up)
        v = pd.DataFrame.from_dict(link_mapper_up, orient="index").values
        v_e = np.array([[2], [3], [6]])
        check_vals(v, v_e)

    def test_map_nmg1_links_to_nmg2_links_4(self, nmgrid_c_up, nmgrid_f_up):
        """map upside down coarse grid to the upside-down fine grids"""
        link_mapper_up = gt.map_nmg1_links_to_nmg2_links(nmgrid_f_up, nmgrid_c_up)
        v = pd.DataFrame.from_dict(link_mapper_up, orient="index").values
        v_e = np.array([[1], [0], [0], [1], [2], [2], [2]])
        check_vals(v, v_e)

    def test_map_nmg1_links_to_nmg2_links_5(self, nmgrid_c, nmgrid_f):
        """map coarse grid to fine grid, but set number_of_points low to purposely
        cause function to fail test. When only a few points are used to estimate
        the mean distance between links, numerous links appear to be the same
        distance from each other"""
        link_mapper = gt.map_nmg1_links_to_nmg2_links(
            nmgrid_f, nmgrid_c, number_of_points=3
        )
        v = pd.DataFrame.from_dict(link_mapper, orient="index").values
        v_e = np.array([[0], [0], [0], [1], [2], [2], [1]])
        check_vals(v, v_e)


class TestCreateDfOfLinkPoints:

    def test_create_df_of_link_points_1(self, nmgrid_s):
        """for single-channel grid"""
        number_of_points = 3
        link_nodes_s = np.array([[0, 1]])
        df = gt.create_df_of_link_points(nmgrid_s, link_nodes_s, number_of_points)
        v_e = np.array([[0.0, 30.0, 0.0], [0.0, 25.0, 15.0], [0.0, 20.0, 30.0]])

        check_vals(df.values, v_e)

    def test_create_df_of_link_points_2(self, nmgrid_c):
        """for coarse grid"""
        number_of_points = 3
        link_nodes_c = np.array([[0, 1], [1, 2], [1, 3]])
        df = gt.create_df_of_link_points(nmgrid_c, link_nodes_c, number_of_points)
        v_e = np.array(
            [
                [0.0, 30.0, 0.0],
                [0.0, 25.0, 15.0],
                [0.0, 20.0, 30.0],
                [1.0, 20.0, 30.0],
                [1.0, 20.0, 35.0],
                [1.0, 20.0, 40.0],
                [2.0, 20.0, 30.0],
                [2.0, 30.0, 35.0],
                [2.0, 40.0, 40.0],
            ]
        )

        check_vals(df.values, v_e)

    def test_create_df_of_link_points_3(self, nmgrid_c_up):
        """for upside-down coarse grid"""
        number_of_points = 3
        link_nodes_c_up = np.array([[2, 0], [2, 1], [3, 2]])
        df = gt.create_df_of_link_points(nmgrid_c_up, link_nodes_c_up, number_of_points)
        v_e = np.array(
            [
                [0.0, 40.0, 40.0],
                [0.0, 30.0, 35.0],
                [0.0, 20.0, 30.0],
                [1.0, 40.0, 40.0],
                [1.0, 40.0, 35.0],
                [1.0, 40.0, 30.0],
                [2.0, 30.0, 70.0],
                [2.0, 35.0, 55.0],
                [2.0, 40.0, 40.0],
            ]
        )

        check_vals(df.values, v_e)
