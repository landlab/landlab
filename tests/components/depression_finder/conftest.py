import numpy as np
import pytest

from landlab import RasterModelGrid
from landlab.components import DepressionFinderAndRouter
from landlab.components import FlowAccumulator

XX = RasterModelGrid.BAD_INDEX


@pytest.fixture
def dans_grid3():
    """
    Create a 7x7 test grid with a well defined hole in it.
    """
    mg = RasterModelGrid((7, 7))

    z = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0],
            [0.0, 2.0, 1.6, 1.5, 1.6, 2.0, 0.0],
            [0.0, 2.0, 1.7, 1.6, 1.7, 2.0, 0.0],
            [0.0, 2.0, 1.8, 2.0, 2.0, 2.0, 0.0],
            [0.0, 1.0, 0.6, 1.0, 1.0, 1.0, 0.0],
            [0.0, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0],
        ]
    ).flatten()

    r_old = np.array(
        [
            [0, 1, 2, 3, 4, 5, 6],
            [7, 1, 2, 3, 4, 5, 13],
            [14, 14, 17, 17, 17, 20, 20],
            [21, 21, 17, 17, 17, 27, 27],
            [28, 28, 37, 38, 39, 34, 34],
            [35, 44, 44, 44, 46, 41, 41],
            [42, 43, 44, 45, 46, 47, 48],
        ]
    ).flatten()

    r_new = np.array(
        [
            [0, 1, 2, 3, 4, 5, 6],
            [7, 1, 2, 3, 4, 5, 13],
            [14, 14, 23, 24, 24, 20, 20],
            [21, 21, 30, 30, 24, 27, 27],
            [28, 28, 37, 38, 39, 34, 34],
            [35, 44, 44, 44, 46, 41, 41],
            [42, 43, 44, 45, 46, 47, 48],
        ]
    ).flatten()

    A_old = np.array(
        [
            [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0],
            [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0, 6.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 1.0],
            [0.0, 0.0, 5.0, 0.0, 2.0, 0.0, 0.0],
        ]
    ).flatten()

    A_new = np.array(
        [
            [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0],
            [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 2.0, 4.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 7.0, 1.0, 1.0, 1.0, 1.0],
            [0.0, 1.0, 8.0, 2.0, 2.0, 1.0, 1.0],
            [0.0, 0.0, 11.0, 0.0, 2.0, 0.0, 0.0],
        ]
    ).flatten()

    s_new = np.array(
        [
            [0, 1, 8, 2, 9, 3, 10],
            [4, 11, 5, 12, 6, 7, 13],
            [14, 15, 20, 19, 21, 22, 27],
            [26, 28, 29, 34, 33, 35, 41],
            [40, 42, 43, 44, 36, 37, 30],
            [23, 16, 24, 17, 18, 25, 38],
            [31, 45, 46, 39, 32, 47, 48],
        ]
    ).flatten()

    links_old = np.array(
        [
            [-1, -1, -1, -1, -1, -1, -1],
            [-1, 7, 8, 9, 10, 11, -1],
            [-1, 26, 28, -1, 29, 31, -1],
            [-1, 39, 113, 35, 114, 44, -1],
            [-1, 52, 60, 61, 62, 57, -1],
            [-1, 146, 73, 149, 75, 70, -1],
            [-1, -1, -1, -1, -1, -1, -1],
        ]
    ).flatten()

    links_new = np.array(
        [
            [-1, -1, -1, -1, -1, -1, -1],
            [-1, 7, 8, 9, 10, 11, -1],
            [-1, 26, 34, 35, 115, 31, -1],
            [-1, 39, 47, 125, 42, 44, -1],
            [-1, 52, 60, 61, 62, 57, -1],
            [-1, 146, 73, 149, 75, 70, -1],
            [-1, -1, -1, -1, -1, -1, -1],
        ]
    ).flatten()

    depr_outlet_target = np.array(
        [
            [XX, XX, XX, XX, XX, XX, XX],
            [XX, XX, XX, XX, XX, XX, XX],
            [XX, XX, 30, 30, 30, XX, XX],
            [XX, XX, 30, 30, 30, XX, XX],
            [XX, XX, XX, XX, XX, XX, XX],
            [XX, XX, XX, XX, XX, XX, XX],
            [XX, XX, XX, XX, XX, XX, XX],
        ]
    ).flatten()

    mg.add_field("topographic__elevation", z, at="node", units="-")

    fr = FlowAccumulator(mg, flow_director="D8")
    lf = DepressionFinderAndRouter(mg)

    class DansGrid:
        pass

    dans_grid = DansGrid()
    dans_grid.mg = mg
    dans_grid.fr = fr
    dans_grid.lf = lf
    dans_grid.z = z
    dans_grid.r_new = r_new
    dans_grid.r_old = r_old
    dans_grid.A_new = A_new
    dans_grid.A_old = A_old
    dans_grid.s_new = s_new
    dans_grid.depr_outlet_target = depr_outlet_target
    dans_grid.links_old = links_old
    dans_grid.links_new = links_new

    return dans_grid


@pytest.fixture
def d4_grid():
    """Test functionality of routing when D4 is specified.

    The elevation field in this test looks like::

    1   2   3   4   5   6   7

    1   2   3   0   5   0   7

    1   2   3   4   0   0   7

    1   2   3   0   5   6   7

    1   2   0   0   0   6   7

    1   2   3   0   5   6   7

    1   2   3   4   5   6   7
    """
    mg1 = RasterModelGrid((7, 7))
    mg2 = RasterModelGrid((7, 7))
    z = mg1.node_x.copy() + 1.0
    lake_nodes = np.array([10, 16, 17, 18, 24, 32, 33, 38, 40])
    z[lake_nodes] = 0.0
    mg1.add_field("topographic__elevation", z, at="node", units="-")
    mg2.add_field("topographic__elevation", z, at="node", units="-")

    frD8 = FlowAccumulator(mg1, flow_director="D8")
    frD4 = FlowAccumulator(mg2, flow_director="D4")
    lfD8 = DepressionFinderAndRouter(mg1, routing="D8")
    lfD4 = DepressionFinderAndRouter(mg2, routing="D4")

    class DansGrid:
        pass

    d4_grid = DansGrid()
    d4_grid.mg1 = mg1
    d4_grid.mg2 = mg2
    d4_grid.z = z
    d4_grid.lake_nodes = lake_nodes
    d4_grid.frD8 = frD8
    d4_grid.frD4 = frD4
    d4_grid.lfD8 = lfD8
    d4_grid.lfD4 = lfD4

    return d4_grid
