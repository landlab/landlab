import os

import numpy as np
import pytest

from landlab import RasterModelGrid

XX = RasterModelGrid.BAD_INDEX


@pytest.fixture
def dans_grid1():
    """
    Create a 5x5 test grid.
    This is a sheet flow test.
    """
    mg = RasterModelGrid((5, 5), xy_spacing=(10.0, 10.0))

    this_dir = os.path.abspath(os.path.dirname(__file__))
    infile = os.path.join(this_dir, "test_fr_input.txt")

    z = mg.node_x.copy()

    A_target = (
        np.array(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [3.0, 3.0, 2.0, 1.0, 0.0],
                [3.0, 3.0, 2.0, 1.0, 0.0],
                [3.0, 3.0, 2.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0],
            ]
        ).flatten()
        * 100.0
    )

    frcvr_target = np.array(
        [
            [0, 1, 2, 3, 4],
            [5, 5, 6, 7, 9],
            [10, 10, 11, 12, 14],
            [15, 15, 16, 17, 19],
            [20, 21, 22, 23, 24],
        ]
    ).flatten()

    upids_target = np.array(
        [
            [0, 1, 2, 3, 4],
            [5, 6, 7, 8, 9],
            [10, 11, 12, 13, 14],
            [15, 16, 17, 18, 19],
            [20, 21, 22, 23, 24],
        ]
    ).flatten()

    links2rcvr_target = np.full(25, XX)
    links2rcvr_target[mg.core_nodes] = np.array([9, 10, 11, 18, 19, 20, 27, 28, 29])

    Q_target = A_target * 2.0  # only once Q_in is used

    steepest_target = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 1.0, 1.0, 0.0],
            [0.0, 1.0, 1.0, 1.0, 0.0],
            [0.0, 1.0, 1.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    ).flatten()

    mg.add_field("topographic__elevation", z, at="node", units="-")

    class DansGrid(object):
        pass

    dans_grid = DansGrid()
    dans_grid.mg = mg
    dans_grid.z = z
    dans_grid.infile = infile
    dans_grid.A_target = A_target
    dans_grid.frcvr_target = frcvr_target
    dans_grid.upids_target = upids_target
    dans_grid.Q_target = Q_target
    dans_grid.steepest_target = steepest_target
    dans_grid.links2rcvr_target = links2rcvr_target

    return dans_grid


@pytest.fixture
def internal_closed():
    """
    Create a 6x5 test grid, but with two internal nodes closed.
    This is a sheet flow test.
    """
    mg = RasterModelGrid((6, 5), xy_spacing=(10.0, 10.0))

    mg.set_closed_boundaries_at_grid_edges(True, True, False, True)
    mg.status_at_node[7] = mg.BC_NODE_IS_CLOSED
    mg.status_at_node[16] = mg.BC_NODE_IS_CLOSED

    z = mg.node_x.copy()

    Q_in = np.full(25, 2.0)

    A_target = (
        np.array(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [1.0, 1.0, 0.0, 1.0, 0.0],
                [6.0, 6.0, 3.0, 1.0, 0.0],
                [0.0, 0.0, 2.0, 1.0, 0.0],
                [3.0, 3.0, 2.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0],
            ]
        ).flatten()
        * 100.0
    )

    frcvr_target = np.array(
        [
            [0, 1, 2, 3, 4],
            [5, 5, 7, 12, 9],
            [10, 10, 11, 12, 14],
            [15, 16, 11, 17, 19],
            [20, 20, 21, 22, 24],
            [25, 26, 27, 28, 29],
        ]
    ).flatten()

    links2rcvr_target = np.full(mg.number_of_nodes, XX)
    links2rcvr_target[mg.core_nodes] = np.array([9, 62, 18, 19, 20, 67, 29, 36, 37, 38])

    steepest_target = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 1.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 1.0, 1.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    ).flatten()

    steepest_target[np.array([8, 17])] = 1.0 / np.sqrt(2.0)

    mg.add_field("topographic__elevation", z, at="node", units="-")

    class DansGrid(object):
        pass

    dans_grid = DansGrid()
    dans_grid.mg = mg
    dans_grid.z = z
    dans_grid.Q_in = Q_in
    dans_grid.A_target = A_target
    dans_grid.frcvr_target = frcvr_target
    dans_grid.steepest_target = steepest_target
    dans_grid.links2rcvr_target = links2rcvr_target

    return dans_grid


@pytest.fixture
def dans_grid2():
    """
    Create a 5x5 test grid.
    This tests more complex routing, with diffs between D4 & D8.
    """
    mg = RasterModelGrid((5, 5), xy_spacing=(10.0, 10.0))

    this_dir = os.path.abspath(os.path.dirname(__file__))
    infile = os.path.join(this_dir, "test_fr_input.txt")

    z = np.array(
        [
            [7.0, 7.0, 7.0, 7.0, 7.0],
            [7.0, 5.0, 3.2, 6.0, 7.0],
            [7.0, 2.0, 3.0, 5.0, 7.0],
            [7.0, 1.0, 1.9, 4.0, 7.0],
            [7.0, 0.0, 7.0, 7.0, 7.0],
        ]
    ).flatten()

    A_target_D8 = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 100.0, 200.0, 100.0, 0.0],
            [0.0, 400.0, 100.0, 100.0, 0.0],
            [0.0, 600.0, 300.0, 100.0, 0.0],
            [0.0, 900.0, 0.0, 0.0, 0.0],
        ]
    ).flatten()

    A_target_D4 = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 100.0, 200.0, 100.0, 0.0],
            [0.0, 200.0, 400.0, 100.0, 0.0],
            [0.0, 900.0, 600.0, 100.0, 0.0],
            [0.0, 900.0, 0.0, 0.0, 0.0],
        ]
    ).flatten()

    frcvr_target_D8 = np.array(
        [
            [0, 1, 2, 3, 4],
            [5, 11, 11, 7, 9],
            [10, 16, 16, 17, 14],
            [15, 21, 21, 17, 19],
            [20, 21, 22, 23, 24],
        ]
    ).flatten()

    frcvr_target_D4 = np.array(
        [
            [0, 1, 2, 3, 4],
            [5, 11, 12, 7, 9],
            [10, 16, 17, 12, 14],
            [15, 21, 16, 17, 19],
            [20, 21, 22, 23, 24],
        ]
    ).flatten()

    upids_target_D8 = np.array(
        [
            [0, 1, 2, 3, 4],
            [5, 9, 10, 14, 15],
            [19, 20, 21, 16, 11],
            [6, 7, 8, 12, 17],
            [13, 18, 22, 23, 24],
        ]
    ).flatten()

    upids_target_D4 = np.array(
        [
            [0, 1, 2, 3, 4],
            [5, 9, 10, 14, 15],
            [19, 20, 21, 16, 11],
            [6, 17, 12, 7, 8],
            [13, 18, 22, 23, 24],
        ]
    ).flatten()

    links2rcvr_target_D8 = np.full(25, XX)
    links2rcvr_target_D8[mg.core_nodes] = np.array([14, 51, 11, 23, 59, 61, 32, 67, 29])

    links2rcvr_target_D4 = np.full(25, XX)
    links2rcvr_target_D4[mg.core_nodes] = np.array([14, 15, 11, 23, 24, 20, 32, 28, 29])

    steepest_target_D8 = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.3, 0.08485281, 0.28, 0.0],
            [0.0, 0.1, 0.14142136, 0.21920310, 0.0],
            [0.0, 0.1, 0.13435029, 0.21, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    ).flatten()

    steepest_target_D4 = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.3, 0.02, 0.28, 0.0],
            [0.0, 0.1, 0.11, 0.2, 0.0],
            [0.0, 0.1, 0.09, 0.21, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    ).flatten()

    mg.add_field("topographic__elevation", z, at="node", units="-")

    class DansGrid(object):
        pass

    dans_grid = DansGrid()
    dans_grid.mg = mg
    dans_grid.z = z
    dans_grid.infile = infile
    dans_grid.A_target_D8 = A_target_D8
    dans_grid.A_target_D4 = A_target_D4
    dans_grid.frcvr_target_D8 = frcvr_target_D8
    dans_grid.frcvr_target_D4 = frcvr_target_D4
    dans_grid.upids_target_D8 = upids_target_D8
    dans_grid.upids_target_D4 = upids_target_D4
    dans_grid.steepest_target_D8 = steepest_target_D8
    dans_grid.steepest_target_D4 = steepest_target_D4
    dans_grid.links2rcvr_target_D8 = links2rcvr_target_D8
    dans_grid.links2rcvr_target_D4 = links2rcvr_target_D4

    return dans_grid
