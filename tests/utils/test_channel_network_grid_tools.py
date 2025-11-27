import numpy as np
import pandas as pd
import pytest

import landlab.utils.channel_network_grid_tools as gt
from landlab import RasterModelGrid
from landlab.components import DepressionFinderAndRouter
from landlab.components import FlowAccumulator


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
def nmg_link_to_rmg_coincident_nodes_mapper():
    """example link to rmg coincident nodes mapper output by map_nmg_links_to_rmg_coincident_nodes"""
    nmg_link_to_rmg_coincident_nodes_mapper = pd.DataFrame(
        {
            "linkID": {0: 0, 1: 0, 2: 0, 3: 0, 4: 1, 5: 2, 6: 2, 7: 2},
            "coincident_node": {0: 3, 1: 10, 2: 16, 3: 23, 4: 30, 5: 24, 6: 31, 7: 32},
            "x": {
                0: 30.0,
                1: 30.0,
                2: 20.0,
                3: 20.0,
                4: 20.0,
                5: 30.0,
                6: 30.0,
                7: 40.0,
            },
            "y": {
                0: 0.0,
                1: 10.0,
                2: 20.0,
                3: 30.0,
                4: 40.0,
                5: 30.0,
                6: 40.0,
                7: 40.0,
            },
            "dist": {
                0: 31.622776601683793,
                1: 26.336486619220135,
                2: 15.795561085325538,
                3: 5.25463555143094,
                4: 4.994994994994997,
                5: 16.764914065538967,
                6: 11.169148356080033,
                7: 5.5733826466210985,
            },
            "drainage_area": {
                0: 2300.0,
                1: 2300.0,
                2: 2300.0,
                3: 2300.0,
                4: 1100.0,
                5: 1050.0,
                6: 1050.0,
                7: 1050.0,
            },
        }
    )

    return nmg_link_to_rmg_coincident_nodes_mapper


class TestMapChannelNodesToNmgLinks:

    def test_map_channel_nodes_to_nmg_links_1(
        self, grid, nmg_link_to_rmg_coincident_nodes_mapper
    ):
        """normal grid"""

        # example channel node ids, output by extract_channel_nodes
        acn = np.array([3, 9, 17, 23, 30, 31, 32, 36])

        cn_to_nmg_link_mapper = gt.map_rmg_nodes_to_nmg_links(
            grid, nmg_link_to_rmg_coincident_nodes_mapper, acn
        )

        # expected values
        l0_e = [np.int64(3), np.int64(9), np.int64(17), np.int64(23)]
        l1_e = [np.int64(30), np.int64(36)]
        l2_e = [np.int64(31), np.int64(32)]

        check_vals(
            cn_to_nmg_link_mapper["node"][cn_to_nmg_link_mapper["linkID"] == 0].values,
            l0_e,
        )
        check_vals(
            cn_to_nmg_link_mapper["node"][cn_to_nmg_link_mapper["linkID"] == 1].values,
            l1_e,
        )
        check_vals(
            cn_to_nmg_link_mapper["node"][cn_to_nmg_link_mapper["linkID"] == 2].values,
            l2_e,
        )

        # check that each channel node has only one equivalent coincident node.
        assert len(np.unique(cn_to_nmg_link_mapper["node"].values)) == len(
            cn_to_nmg_link_mapper["node"].values
        )

    def test_map_channel_nodes_to_nmg_links_2(self, grid_up):
        """upside down grid"""

        # example link to rmg coincident nodes mapper output by map_nmg_links_to_rmg_coincident_nodes
        nmg_link_to_rmg_coincident_nodes_mapper_up = pd.DataFrame(
            {
                "linkID": {0: 0, 1: 0, 2: 0, 3: 1, 4: 2, 5: 2, 6: 2, 7: 2},
                "coincident_node": {
                    0: 31,
                    1: 24,
                    2: 23,
                    3: 25,
                    4: 52,
                    5: 45,
                    6: 39,
                    7: 32,
                },
                "x": {
                    0: 30.0,
                    1: 30.0,
                    2: 20.0,
                    3: 40.0,
                    4: 30.0,
                    5: 30.0,
                    6: 40.0,
                    7: 40.0,
                },
                "y": {
                    0: 40.0,
                    1: 30.0,
                    2: 30.0,
                    3: 30.0,
                    4: 70.0,
                    5: 60.0,
                    6: 50.0,
                    7: 40.0,
                },
                "dist": {
                    0: 16.764914065538967,
                    1: 11.169148356080033,
                    2: 5.5733826466210985,
                    3: 4.994994994994997,
                    4: 31.622776601683793,
                    5: 26.33648661922014,
                    6: 15.795561085325536,
                    7: 5.254635551430944,
                },
                "drainage_area": {
                    0: 1050.0,
                    1: 1050.0,
                    2: 1050.0,
                    3: 1100.0,
                    4: 2300.0,
                    5: 2300.0,
                    6: 2300.0,
                    7: 2300.0,
                },
            }
        )
        # example channel node ids, output by extract_channel_nodes
        acn_up = np.array([19, 23, 24, 25, 32, 38, 46, 52])

        cn_to_nmg_link_mapper_up = gt.map_rmg_nodes_to_nmg_links(
            grid_up, nmg_link_to_rmg_coincident_nodes_mapper_up, acn_up
        )

        # expected values
        l0_e = [np.int64(23), np.int64(24)]
        l1_e = [np.int64(19), np.int64(25)]
        l2_e = [np.int64(32), np.int64(38), np.int64(46), np.int64(52)]

        check_vals(
            cn_to_nmg_link_mapper_up["node"][
                cn_to_nmg_link_mapper_up["linkID"] == 0
            ].values,
            l0_e,
        )
        check_vals(
            cn_to_nmg_link_mapper_up["node"][
                cn_to_nmg_link_mapper_up["linkID"] == 1
            ].values,
            l1_e,
        )
        check_vals(
            cn_to_nmg_link_mapper_up["node"][
                cn_to_nmg_link_mapper_up["linkID"] == 2
            ].values,
            l2_e,
        )

        # check that each channel and terrace node has only one equivalent nmg rmg node.
        assert len(np.unique(cn_to_nmg_link_mapper_up["node"].values)) == len(
            cn_to_nmg_link_mapper_up["node"].values
        )


class TestRemoveSmallTribs:

    def test_remove_small_tribs_1(self, nmg_link_to_rmg_coincident_nodes_mapper):
        """small trib at node 17, mapped to coincident node 16"""
        rmg_nodes_to_nmg_links_mapper = pd.DataFrame(
            {
                "node": {0: 3, 1: 9, 2: 17, 3: 23},
                "linkID": {0: 0, 1: 0, 2: 0, 3: 0},
                "dist": {
                    0: 31.622776601683793,
                    1: 26.336486619220135,
                    2: 26.336486619220135,
                    3: 5.25463555143094,
                },
                "coincident_node": {0: 3, 1: 10, 2: 16, 3: 23},
                "drainage_area": {0: 2300.0, 1: 2300.0, 2: 2300.0, 3: 2300.0},
                "node_drainage_area": {0: 2900.0, 1: 2300.0, 2: 100.0, 3: 1600.0},
            }
        )

        rmg_nodes_to_nmg_links_mapper = gt._remove_small_tribs(
            rmg_nodes_to_nmg_links_mapper, nmg_link_to_rmg_coincident_nodes_mapper
        )  # cn_to_nmg_link_mapper)

        rmg_nodes_to_nmg_links_mapper_e = np.array(
            [
                [3, 0.0, 31.6227766, 3.0, 2300.0, 2900.0],
                [9.0, 0.0, 26.33648662, 10.0, 2300.0, 2300.0],
                [23.0, 0.0, 5.25463555, 23.0, 2300.0, 1600.0],
            ]
        )
        check_vals(rmg_nodes_to_nmg_links_mapper, rmg_nodes_to_nmg_links_mapper_e)

    def test_remove_small_tribs_2(self, nmg_link_to_rmg_coincident_nodes_mapper):
        """small trib linked to outlet node"""
        rmg_nodes_to_nmg_links_mapper = pd.DataFrame(
            {
                "node": {0: 3, 1: 4, 2: 9, 3: 23},
                "linkID": {0: 0, 1: 0, 2: 0, 3: 0},
                "dist": {
                    0: 31.622776601683793,
                    1: 31.622776601683793,
                    2: 26.336486619220135,
                    3: 5.25463555143094,
                },
                "coincident_node": {0: 3, 1: 3, 2: 10, 3: 23},
                "drainage_area": {0: 2300.0, 1: 2300, 2: 2300.0, 3: 2300.0},
                "node_drainage_area": {0: 2900.0, 1: 100, 2: 2300, 3: 1600.0},
            }
        )

        rmg_nodes_to_nmg_links_mapper = gt._remove_small_tribs(
            rmg_nodes_to_nmg_links_mapper, nmg_link_to_rmg_coincident_nodes_mapper
        )  # cn_to_nmg_link_mapper)

        rmg_nodes_to_nmg_links_mapper_e = np.array(
            [
                [3.0, 0.0, 31.6227766, 3.0, 2300.0, 2900.0],
                [9.0, 0.0, 26.33648662, 10.0, 2300.0, 2300.0],
                [23.0, 0.0, 5.25463555, 23.0, 2300.0, 1600.0],
            ]
        )

        check_vals(rmg_nodes_to_nmg_links_mapper, rmg_nodes_to_nmg_links_mapper_e)

    def test_remove_small_tribs_3(self, nmg_link_to_rmg_coincident_nodes_mapper):
        """small trib linked to the inlet node (node 24 is mapped to coincident
        node 23,"""
        rmg_nodes_to_nmg_links_mapper = pd.DataFrame(
            {
                "node": {0: 3, 1: 9, 2: 24, 3: 23},
                "linkID": {0: 0, 1: 0, 2: 0, 3: 0},
                "dist": {
                    0: 31.622776601683793,
                    1: 26.336486619220135,
                    2: 5.25463555143094,
                    3: 5.25463555143094,
                },
                "coincident_node": {0: 3, 1: 10, 2: 23, 3: 23},
                "drainage_area": {0: 2300.0, 1: 2300.0, 2: 2300.0, 3: 2300.0},
                "node_drainage_area": {0: 2900.0, 1: 2300.0, 2: 100.0, 3: 1600.0},
            }
        )

        rmg_nodes_to_nmg_links_mapper_e = np.array(
            [
                [27, 0, 31.6227766, 3, 2300, 2900],
                [25, 0, 26.33648662, 10, 2300, 2300],
                [31, 0, 5.25463555, 23, 2300, 1600],
            ]
        )

        rmg_nodes_to_nmg_links_mapper = gt._remove_small_tribs(
            rmg_nodes_to_nmg_links_mapper, nmg_link_to_rmg_coincident_nodes_mapper
        )  # cn_to_nmg_link_mapper)
        # in this case, remove_small_tribs doesn't remove node 24 from the mapper
        # and the test fails.
        # This result is expected and the test is included to highlight that there
        # are situations where the remove_small_tribs function (as presently coded)
        # may not work as expected
        check_vals(rmg_nodes_to_nmg_links_mapper, rmg_nodes_to_nmg_links_mapper_e)
