import numpy as np
import pytest

import landlab.utils.channel_network_grid_tools as gt


def check_vals(vals, vals_e):
    np.testing.assert_allclose(vals, vals_e, atol=1e-4)


def test_create_lol():
    a_list = [1, 2, 3, 3, 2, 4, 5, 6]
    linkID = [0, 0, 0, 1, 1, 1, 2, 2]
    lol_e = [[1, 2, 3], [3, 2, 4], [5, 6]]
    assert lol_e == gt._create_lol(linkID, a_list)


def test_flaten_lol():
    lol = [[1, 2, 3], [3, 2, 4], [5, 6]]
    a_list_e = [1, 2, 3, 3, 2, 4, 5, 6]
    assert a_list_e == gt._flatten_lol(lol)


def test_dist_func():
    check_vals(gt._dist_func(0, 1, 0, 1), 1.414213)


class TestGetLinknodes:

    def test_get_link_nodes_coarse(self, nmgrid_c):
        """test using simple nmgrid"""
        link_nodes_c = gt.get_link_nodes(nmgrid_c)
        link_nodes_c_e = np.array([[0, 1], [1, 2], [1, 3]])
        check_vals(link_nodes_c, link_nodes_c_e)

    def test_get_link_nodes_fine(self, nmgrid_f):
        """test using a more complex nmgrid"""
        link_nodes_f = gt.get_link_nodes(nmgrid_f)
        link_nodes_f_e = np.array(
            [[0, 1], [1, 2], [2, 3], [3, 4], [3, 5], [5, 6], [4, 7]]
        )
        check_vals(link_nodes_f, link_nodes_f_e)

    def test_get_link_nodes_slope_is_zero(self, nmgrid_c):
        """what happens if no slope?"""
        nmgrid_c.at_node["topographic__elevation"][
            1
        ] = 1  # set elevation of second node equal to first, so no slope
        link_nodes_c = gt.get_link_nodes(nmgrid_c)
        link_nodes_c_e = np.array(
            [
                [
                    -1,
                    -1,
                ],  # there are no downstream or upstream nodes, both listed as -1
                [1, 2],
                [1, 3],
            ]
        )
        check_vals(link_nodes_c, link_nodes_c_e)


class TestLinkToPointsAndDist:
    def test_link_to_points_and_dist_1(self):
        """from 0,0 to 1,1"""
        x0 = 0
        x1 = 1
        y0 = 0
        y1 = 1
        x_, y_, dist_ = gt._link_to_points_and_dist(
            (x0, y0), (x1, y1), number_of_points=5
        )

        v_e = (
            np.array([0.0, 0.25, 0.5, 0.75, 1.0]),
            np.array([0.0, 0.25, 0.5, 0.75, 1.0]),
            np.array([0.0, 0.35355339, 0.70710678, 1.06066017, 1.41421356]),
        )
        check_vals(x_, v_e[0])
        check_vals(y_, v_e[1])
        check_vals(dist_, v_e[2])

    def test_link_to_points_and_dist_2(self):
        """from 1,1 to 0,0"""
        x0 = 1
        x1 = 0
        y0 = 1
        y1 = 0
        x_, y_, dist_ = gt._link_to_points_and_dist(
            (x0, y0), (x1, y1), number_of_points=5
        )
        v_e = (
            np.array([1.0, 0.75, 0.5, 0.25, 0.0]),
            np.array([1.0, 0.75, 0.5, 0.25, 0.0]),
            np.array([0.0, 0.35355339, 0.70710678, 1.06066017, 1.41421356]),
        )
        check_vals(x_, v_e[0])
        check_vals(y_, v_e[1])
        check_vals(dist_, v_e[2])

    def test_link_to_points_and_dist_3(self):
        """from 0,0 to 1,-1"""
        x0 = 0
        x1 = 1
        y0 = 0
        y1 = -1
        x_, y_, dist_ = gt._link_to_points_and_dist(
            (x0, y0), (x1, y1), number_of_points=5
        )
        v_e = (
            np.array([0.0, 0.25, 0.5, 0.75, 1.0]),
            np.array([0.0, -0.25, -0.5, -0.75, -1.0]),
            np.array([0.0, 0.35355339, 0.70710678, 1.06066017, 1.41421356]),
        )
        check_vals(x_, v_e[0])
        check_vals(y_, v_e[1])
        check_vals(dist_, v_e[2])

    def test_link_to_points_and_dist_4(self):
        """from 1,-1 to 0,0"""
        x0 = 1
        x1 = 0
        y0 = -1
        y1 = 0
        x_, y_, dist_ = gt._link_to_points_and_dist(
            (x0, y0), (x1, y1), number_of_points=5
        )
        v_e = (
            np.array([1.0, 0.75, 0.5, 0.25, 0.0]),
            np.array([-1.0, -0.75, -0.5, -0.25, 0.0]),
            np.array([0.0, 0.35355339, 0.70710678, 1.06066017, 1.41421356]),
        )
        check_vals(x_, v_e[0])
        check_vals(y_, v_e[1])
        check_vals(dist_, v_e[2])

    def test_link_to_points_and_dist_5(self):
        """from 0,0 to -1,-1"""
        x0 = 0
        x1 = -1
        y0 = 0
        y1 = -1
        x_, y_, dist_ = gt._link_to_points_and_dist(
            (x0, y0), (x1, y1), number_of_points=5
        )
        v_e = (
            np.array([0.0, -0.25, -0.5, -0.75, -1.0]),
            np.array([0.0, -0.25, -0.5, -0.75, -1.0]),
            np.array([0.0, 0.35355339, 0.70710678, 1.06066017, 1.41421356]),
        )
        check_vals(x_, v_e[0])
        check_vals(y_, v_e[1])
        check_vals(dist_, v_e[2])

    def test_link_to_points_and_dist_6(self):
        """from -1,-1 to 0,0"""
        x0 = -1
        x1 = 0
        y0 = -1
        y1 = 0
        x_, y_, dist_ = gt._link_to_points_and_dist(
            (x0, y0), (x1, y1), number_of_points=5
        )
        v_e = (
            np.array([-1.0, -0.75, -0.5, -0.25, 0.0]),
            np.array([-1.0, -0.75, -0.5, -0.25, 0.0]),
            np.array([0.0, 0.35355339, 0.70710678, 1.06066017, 1.41421356]),
        )
        check_vals(x_, v_e[0])
        check_vals(y_, v_e[1])
        check_vals(dist_, v_e[2])

    def test_link_to_points_and_dist_7(self):
        """from 0,0 to -1,1"""
        x0 = 0
        x1 = -1
        y0 = 0
        y1 = 1
        x_, y_, dist_ = gt._link_to_points_and_dist(
            (x0, y0), (x1, y1), number_of_points=5
        )
        v_e = (
            np.array([0.0, -0.25, -0.5, -0.75, -1.0]),
            np.array([0.0, 0.25, 0.5, 0.75, 1.0]),
            np.array([0.0, 0.35355339, 0.70710678, 1.06066017, 1.41421356]),
        )
        check_vals(x_, v_e[0])
        check_vals(y_, v_e[1])
        check_vals(dist_, v_e[2])

    def test_link_to_points_and_dist_8(self):
        """from -1,1 to 0,0"""
        x0 = -1
        x1 = 0
        y0 = 1
        y1 = 0
        x_, y_, dist_ = gt._link_to_points_and_dist(
            (x0, y0), (x1, y1), number_of_points=5
        )
        v_e = (
            np.array([-1.0, -0.75, -0.5, -0.25, 0.0]),
            np.array([1.0, 0.75, 0.5, 0.25, 0.0]),
            np.array([0.0, 0.35355339, 0.70710678, 1.06066017, 1.41421356]),
        )
        check_vals(x_, v_e[0])
        check_vals(y_, v_e[1])
        check_vals(dist_, v_e[2])

    def test_link_to_points_and_dist_9(self):
        """from 0,0 to 1,0"""
        x0 = 0
        x1 = 1
        y0 = 0
        y1 = 0
        x_, y_, dist_ = gt._link_to_points_and_dist(
            (x0, y0), (x1, y1), number_of_points=5
        )
        v_e = (
            np.array([0.0, 0.25, 0.5, 0.75, 1.0]),
            np.array([0.0, 0.0, 0.0, 0.0, 0.0]),
            np.array([0.0, 0.25, 0.5, 0.75, 1.0]),
        )
        check_vals(x_, v_e[0])
        check_vals(y_, v_e[1])
        check_vals(dist_, v_e[2])

    def test_link_to_points_and_dist_10(self):
        """from 1,0 to 0,0"""
        x0 = 1
        x1 = 0
        y0 = 0
        y1 = 0
        x_, y_, dist_ = gt._link_to_points_and_dist(
            (x0, y0), (x1, y1), number_of_points=5
        )
        v_e = (
            np.array([1.0, 0.75, 0.5, 0.25, 0.0]),
            np.array([0.0, 0.0, 0.0, 0.0, 0.0]),
            np.array([0.0, 0.25, 0.5, 0.75, 1.0]),
        )
        check_vals(x_, v_e[0])
        check_vals(y_, v_e[1])
        check_vals(dist_, v_e[2])

    def test_link_to_points_and_dist_11(self):
        """from 0,0 to 0,-1"""
        x0 = 0
        x1 = 0
        y0 = 0
        y1 = -1
        x_, y_, dist_ = gt._link_to_points_and_dist(
            (x0, y0), (x1, y1), number_of_points=5
        )
        v_e = (
            np.array([0.0, 0.0, 0.0, 0.0, 0.0]),
            np.array([0.0, -0.25, -0.5, -0.75, -1.0]),
            np.array([0.0, 0.25, 0.5, 0.75, 1.0]),
        )
        check_vals(x_, v_e[0])
        check_vals(y_, v_e[1])
        check_vals(dist_, v_e[2])

    def test_link_to_points_and_dist_12(self):
        """from 0,-1 to 0,0"""
        x0 = 0
        x1 = 0
        y0 = -1
        y1 = 0
        x_, y_, dist_ = gt._link_to_points_and_dist(
            (x0, y0), (x1, y1), number_of_points=5
        )
        v_e = (
            np.array([0.0, 0.0, 0.0, 0.0, 0.0]),
            np.array([-1.0, -0.75, -0.5, -0.25, 0.0]),
            np.array([0.0, 0.25, 0.5, 0.75, 1.0]),
        )
        check_vals(x_, v_e[0])
        check_vals(y_, v_e[1])
        check_vals(dist_, v_e[2])

    def test_link_to_points_and_dist_13(self):
        """from 0,0 to -1,0"""
        x0 = 0
        x1 = -1
        y0 = 0
        y1 = 0
        x_, y_, dist_ = gt._link_to_points_and_dist(
            (x0, y0), (x1, y1), number_of_points=5
        )
        v_e = (
            np.array([0.0, -0.25, -0.5, -0.75, -1.0]),
            np.array([0.0, 0.0, 0.0, 0.0, 0.0]),
            np.array([0.0, 0.25, 0.5, 0.75, 1.0]),
        )
        check_vals(x_, v_e[0])
        check_vals(y_, v_e[1])
        check_vals(dist_, v_e[2])

    def test_link_to_points_and_dist_14(self):
        """from -1,0 to 0,0"""
        x0 = -1
        x1 = 0
        y0 = 0
        y1 = 0
        x_, y_, dist_ = gt._link_to_points_and_dist(
            (x0, y0), (x1, y1), number_of_points=5
        )
        v_e = (
            np.array([-1.0, -0.75, -0.5, -0.25, 0.0]),
            np.array([0.0, 0.0, 0.0, 0.0, 0.0]),
            np.array([0.0, 0.25, 0.5, 0.75, 1.0]),
        )
        check_vals(x_, v_e[0])
        check_vals(y_, v_e[1])
        check_vals(dist_, v_e[2])

    def test_link_to_points_and_dist_15(self):
        """from 0,0 to 0,1"""
        x0 = 0
        x1 = 0
        y0 = 0
        y1 = 1
        x_, y_, dist_ = gt._link_to_points_and_dist(
            (x0, y0), (x1, y1), number_of_points=5
        )
        v_e = (
            np.array([0.0, 0.0, 0.0, 0.0, 0.0]),
            np.array([0.0, 0.25, 0.5, 0.75, 1.0]),
            np.array([0.0, 0.25, 0.5, 0.75, 1.0]),
        )
        check_vals(x_, v_e[0])
        check_vals(y_, v_e[1])
        check_vals(dist_, v_e[2])

    def test_link_to_points_and_dist_16(self):
        """from 0,1 to 0,0"""
        x0 = 0
        x1 = 0
        y0 = 1
        y1 = 0
        x_, y_, dist_ = gt._link_to_points_and_dist(
            (x0, y0), (x1, y1), number_of_points=5
        )
        v_e = (
            np.array([0.0, 0.0, 0.0, 0.0, 0.0]),
            np.array([1.0, 0.75, 0.5, 0.25, 0.0]),
            np.array([0.0, 0.25, 0.5, 0.75, 1.0]),
        )
        check_vals(x_, v_e[0])
        check_vals(y_, v_e[1])
        check_vals(dist_, v_e[2])


class TestExtractChannelNodes:

    def test_extract_channel_nodes_1(self, grid):
        """large contributing area"""
        fcn = gt.extract_channel_nodes(grid, 1200)
        fcn_e = np.array([3, 9, 17, 23])
        check_vals(fcn, fcn_e)

    def test_extract_channel_nodes_2(self, grid):
        """small contributing area"""
        acn = gt.extract_channel_nodes(grid, 300)
        cn_e = np.array([3, 9, 17, 23, 30, 31, 32, 36])
        check_vals(acn, cn_e)

    def test_extract_channel_nodes_3(self, grid):
        """contributing area is zero"""
        acn = gt.extract_channel_nodes(grid, 0)
        cn_e = grid.nodes.flatten()  # should return all nodes
        check_vals(acn, cn_e)


class TestExtractTerraceNodes:

    def test_extract_terrace_nodes_1(self, grid):
        """terrace width = 1"""
        terrace_width = 1
        acn = np.array([3, 9, 17, 23, 30, 31, 32, 36])
        fcn = np.array([3, 9, 17, 23])
        tn = gt.extract_terrace_nodes(grid, terrace_width, acn, fcn)
        tn_e = np.array([1, 2, 4, 8, 10, 11, 15, 16, 18, 22, 24, 25, 29])
        check_vals(tn, tn_e)

    def test_extract_terrace_nodes_2(self, grid):
        """terrace width = 2"""
        terrace_width = 2
        acn = np.array([3, 9, 17, 23, 30, 31, 32, 36])
        fcn = np.array([3, 9, 17, 23])
        tn = gt.extract_terrace_nodes(grid, terrace_width, acn, fcn)
        tn_e = np.array(
            [
                0,
                1,
                2,
                4,
                5,
                7,
                8,
                10,
                11,
                12,
                14,
                15,
                16,
                18,
                19,
                21,
                22,
                24,
                25,
                26,
                28,
                29,
                33,
                35,
                37,
            ]
        )
        check_vals(tn, tn_e)

    def test_extract_terrace_nodes_3(self, grid):
        """terrace width = 0"""
        terrace_width = 0
        acn = np.array([3, 9, 17, 23, 30, 31, 32, 36])
        fcn = np.array([3, 9, 17, 23])
        with pytest.raises(ValueError) as exc_info:
            gt.extract_terrace_nodes(grid, terrace_width, acn, fcn)
        assert exc_info.match("terrace width must be 1 or greater")

    def test_extract_terrace_nodes_4(self, grid):
        """terrace width = 1.5 - should round to 2"""
        terrace_width = 1.5
        acn = np.array([3, 9, 17, 23, 30, 31, 32, 36])
        fcn = np.array([3, 9, 17, 23])
        tn = gt.extract_terrace_nodes(grid, terrace_width, acn, fcn)
        tn_e = np.array(
            [
                0,
                1,
                2,
                4,
                5,
                7,
                8,
                10,
                11,
                12,
                14,
                15,
                16,
                18,
                19,
                21,
                22,
                24,
                25,
                26,
                28,
                29,
                33,
                35,
                37,
            ]
        )
        check_vals(tn, tn_e)


class TestMinDistToNetwork:

    def test_min_dist_to_network_1(self, grid):
        """a normal test, only one closest node"""
        node_id = 47
        acn = np.array([3, 9, 17, 23, 30, 31, 32, 36])
        dist, mdn = gt.min_distance_to_network(grid, acn, node_id)
        dist_e = 22.36068
        mdn_e = 32
        check_vals(dist, dist_e)
        check_vals(mdn, mdn_e)

    def test_min_dist_to_network_2(self, grid):
        """two closest nodes, will pick first"""
        node_id = 2
        acn = np.array([3, 9, 17, 23, 30, 31, 32, 36])
        dist, mdn = gt.min_distance_to_network(grid, acn, node_id)
        dist_e = 10
        mdn_e = 3
        check_vals(dist, dist_e)
        check_vals(mdn, mdn_e)

    def test_min_dist_to_network_3(self, grid):
        """node is coincides with a channel node"""
        node_id = 3
        acn = np.array([3, 9, 17, 23, 30, 31, 32, 36])
        dist, mdn = gt.min_distance_to_network(grid, acn, node_id)
        dist_e = 0
        mdn_e = 3
        check_vals(dist, dist_e)
        check_vals(mdn, mdn_e)


class TestMapNMGLinksToRMGCoincidentNodes:

    def test_map_nmg_links_to_rmg_coincident_nodes_1(self, grid, nmgrid_c):
        """remove duplicates"""
        link_nodes_c = gt.get_link_nodes(nmgrid_c)
        nmg_link_to_rmg_coincident_nodes_mapper = (
            gt.map_nmg_links_to_rmg_coincident_nodes(
                grid, nmgrid_c, link_nodes_c, remove_duplicates=True
            )
        )

        nmg_link_to_rmg_coincident_nodes_mapper_e = np.array(
            [
                [0, 3, 30.0, 0.0, 31.622777, 2300.0],
                [0, 10, 30.0, 10.0, 26.336487, 2300.0],
                [0, 16, 20.0, 20.0, 15.795561, 2300.0],
                [0, 23, 20.0, 30.0, 5.254636, 2300.0],
                [1, 30, 20.0, 40.0, 4.994995, 1100.0],
            ]
        )

        check_vals(
            nmg_link_to_rmg_coincident_nodes_mapper.iloc[0:5].to_numpy(),
            nmg_link_to_rmg_coincident_nodes_mapper_e,
        )

        # check no duplicate nodes
        assert len(
            np.unique(nmg_link_to_rmg_coincident_nodes_mapper["coincident_node"].values)
        ) == len(nmg_link_to_rmg_coincident_nodes_mapper["coincident_node"].values)

    def test_map_nmg_links_to_rmg_coincident_nodes_2(self, grid, nmgrid_c):
        """don't remove duplicate nodes"""
        link_nodes_c = gt.get_link_nodes(nmgrid_c)
        nmg_link_to_rmg_coincident_nodes_mapper = (
            gt.map_nmg_links_to_rmg_coincident_nodes(
                grid, nmgrid_c, link_nodes_c, remove_duplicates=False
            )
        )

        nmg_link_to_rmg_coincident_nodes_mapper_e = np.array(
            [
                [0, 3, 30.0, 0.0, 31.622777, 2300.0],
                [0, 10, 30.0, 10.0, 26.336487, 2300.0],
                [0, 16, 20.0, 20.0, 15.795561, 2300.0],
                [0, 23, 20.0, 30.0, 5.254636, 2300.0],
                [1, 23, 20.0, 30.0, 10.00000, 1100.0],
            ]
        )

        check_vals(
            nmg_link_to_rmg_coincident_nodes_mapper.iloc[0:5].to_numpy(),
            nmg_link_to_rmg_coincident_nodes_mapper_e,
        )

        # check that duplicate nodes are included
        assert len(
            np.unique(nmg_link_to_rmg_coincident_nodes_mapper["coincident_node"].values)
        ) != len(nmg_link_to_rmg_coincident_nodes_mapper["coincident_node"].values)

    def test_map_nmg_links_to_rmg_coincident_nodes_3(self, grid, nmgrid_f):
        """use the fine-scale network model grid, remove duplicates"""
        link_nodes_f = gt.get_link_nodes(nmgrid_f)
        nmg_link_to_rmg_coincident_nodes_mapper = (
            gt.map_nmg_links_to_rmg_coincident_nodes(
                grid, nmgrid_f, link_nodes_f, remove_duplicates=True
            )
        )

        nmg_link_to_rmg_coincident_nodes_mapper_e = np.array(
            [
                [0, 3, 30.0, 0.0, 14.142136, 2750.0],
                [0, 9, 20.0, 10.0, 7.063990, 2750.0],
                [1, 17, 30.0, 20.0, 7.063990, 2400.0],
                [2, 23, 20.0, 30.0, 7.063990, 1950.0],
                [3, 30, 20.0, 40.0, 4.994995, 1100.0],
            ]
        )

        check_vals(
            nmg_link_to_rmg_coincident_nodes_mapper.iloc[0:5].to_numpy(),
            nmg_link_to_rmg_coincident_nodes_mapper_e,
        )

        # check no duplicate nodes
        assert len(
            np.unique(nmg_link_to_rmg_coincident_nodes_mapper["coincident_node"].values)
        ) == len(nmg_link_to_rmg_coincident_nodes_mapper["coincident_node"].values)
