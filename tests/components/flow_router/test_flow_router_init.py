import numpy as np
from numpy.testing import (
    assert_array_equal,
    assert_array_almost_equal,
    assert_almost_equal,
)
from landlab.components import FlowRouter

from landlab import RasterModelGrid, HexModelGrid
from landlab import NetworkModelGrid


def test_init_raster():
    spacing = 10
    g = RasterModelGrid((4, 4), (spacing, spacing))
    g.status_at_node[g.perimeter_nodes] = g.BC_NODE_IS_CLOSED
    nodes_n = g.number_of_nodes
    self = FlowRouter(g)
    g.at_node["topographic__elevation"] = np.array(
        [10., 20., 10., 10.] +
        [10., 0., 5., 10.] +
        [20., 20., 5., 10.] +
        [10., 20., 25., 15.])

    assert_array_equal(self._nodes, np.arange(16))

    # 1. Initialization of the input parameters
    ###########################################
    # 1.1. Surface where the flow is directed
    assert self._surface == "topographic__elevation"

    # 1.2. Water influx external to grid
    assert self._uniform_water_external_influx
    assert_almost_equal(
        g.at_node["water__unit_flux_in"],
        np.array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.] +
            [1.,  1.,  1.]),
    )

    # 1.3. Options
    assert self._diagonals

    # 2. Boundary conditions
    ########################
    # 2.1. Boundary settings: guarantee of one base-level node (at least)
    assert_array_equal(
        g.closed_boundary_nodes,
        np.array([ 0,  1,  2,  4,  7,  8, 11, 12, 13, 14, 15]),
    )
    assert_array_equal(g.open_boundary_nodes, np.array([3]))
    assert_array_equal(
        self._closed_nodes,
        np.arrayarray([ 0,  1,  2,  4,  7,  8, 11, 12, 13, 14, 15]),
    )
    assert_array_equal(self._base_level_nodes, np.array([3]))
    assert_array_equal(
        self._base_level_and_closed_nodes,
        np.array([ 3,  0,  1,  2,  4,  7,  8, 11, 12, 13, 14, 15]),
    )

    # 2.2 Cell area at boundary nodes = 0, NetworkModelGrid cell area = 1
    # used in run_flow_accumulations()
    assert_array_almost_equal(
        g.at_node["cell_area_at_node"],
        np.array([   0.,    0.,    0.,    0.,    0.,  100.,  100.,    0.,    0.] +
            [100.,  100.,    0.,    0.,    0.,    0.,    0.]),
    )

    # 2.3. Max number of nodes (for sort head/tails/links)
    assert self._max_number_of_nodes == 1e9

    # 3. Determination of stable input grid data necessary to calculate
    # run_flow_directions
    ###################################################################

    assert self._neighbors_max_number == 8

    # Link infos(tail, head, link id, gradient) sorted by head id.
    # These infos are voluntarily duplicate
    assert_array_equal(
        self._dupli_links,
        np.array(
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            + [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36]
            + [37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53]
            + [54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70]
            + [71, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
            + [19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
            + [36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52]
            + [53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69]
            + [70, 71]
        ),
    )
    assert_array_equal(
        self._pseudo_head_nodes,
        np.array(
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 6, 7, 8, 9, 10, 11, 12, 13, 14, 11, 12]
            + [13, 14, 15, 16, 17, 18, 19, 16, 17, 18, 19, 20, 21, 22, 23, 24, 21]
            + [22, 23, 24, 6, 5, 7, 6, 8, 7, 9, 8, 11, 10, 12, 11, 13, 12, 14, 13]
            + [16, 15, 17, 16, 18, 17, 19, 18, 21, 20, 22, 21, 23, 22, 24, 23, 0]
            + [1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7, 8, 5, 6, 7, 8, 9, 10, 11, 12, 13]
            + [10, 11, 12, 13, 14, 15, 16, 17, 18, 15, 16, 17, 18, 19, 20, 21, 22]
            + [23, 0, 1, 1, 2, 2, 3, 3, 4, 5, 6, 6, 7, 7, 8, 8, 9, 10, 11, 11, 12]
            + [12, 13, 13, 14, 15, 16, 16, 17, 17, 18, 18, 19]
        ),
    )
    assert_array_equal(
        self._pseudo_tail_nodes,
        np.array(
            [0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7, 8, 5, 6, 7, 8, 9, 10, 11, 12, 13]
            + [10, 11, 12, 13, 14, 15, 16, 17, 18, 15, 16, 17, 18, 19, 20, 21, 22]
            + [23, 0, 1, 1, 2, 2, 3, 3, 4, 5, 6, 6, 7, 7, 8, 8, 9, 10, 11, 11, 12]
            + [12, 13, 13, 14, 15, 16, 16, 17, 17, 18, 18, 19, 1, 2, 3, 4, 5, 6, 7]
            + [8, 9, 6, 7, 8, 9, 10, 11, 12, 13, 14, 11, 12, 13, 14, 15, 16, 17, 18]
            + [19, 16, 17, 18, 19, 20, 21, 22, 23, 24, 21, 22, 23, 24, 6, 5, 7, 6]
            + [8, 7, 9, 8, 11, 10, 12, 11, 13, 12, 14, 13, 16, 15, 17, 16, 18, 17]
            + [19, 18, 21, 20, 22, 21, 23, 22, 24, 23]
        ),
    )

    assert_array_equal(
        self._sorted_pseudo_heads,
        np.array(
            [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5]
            + [5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8]
            + [8, 8, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11]
            + [11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13]
            + [13, 13, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16]
            + [16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18]
            + [18, 18, 19, 19, 19, 19, 19, 20, 20, 20, 21, 21, 21, 21, 21, 22, 22]
            + [22, 22, 22, 23, 23, 23, 23, 23, 24, 24, 24]
        ),
    )
    assert_array_equal(
        self._sorted_pseudo_tails,
        np.array(
            [6, 5, 1, 0, 7, 5, 6, 2, 1, 6, 7, 3, 8, 9, 8, 4, 7, 2, 8, 9, 3, 6, 0]
            + [1, 10, 11, 0, 12, 5, 2, 1, 7, 10, 11, 1, 12, 8, 11, 6, 3, 2, 13, 12]
            + [7, 3, 13, 14, 9, 2, 4, 13, 3, 4, 8, 14, 6, 5, 11, 15, 16, 10, 16, 6]
            + [15, 5, 17, 7, 12, 11, 17, 8, 18, 13, 16, 7, 6, 9, 14, 19, 17, 18, 7]
            + [12, 8, 19, 9, 13, 8, 18, 10, 11, 21, 20, 16, 20, 15, 17, 21, 22, 10]
            + [12, 11, 12, 16, 18, 11, 23, 21, 13, 22, 23, 17, 19, 24, 22, 13, 12]
            + [14, 23, 24, 14, 18, 13, 21, 16, 15, 22, 17, 20, 15, 16, 23, 21, 16]
            + [18, 17, 19, 22, 24, 17, 18, 19, 18, 23]
        ),
    )
    assert_array_equal(
        self._sorted_dupli_links,
        np.array(
            [40, 4, 0, 0, 42, 41, 5, 1, 1, 43, 6, 2, 44, 46, 7, 3, 45, 2, 47, 8]
            + [3, 9, 4, 41, 13, 48, 40, 50, 9, 43, 5, 10, 49, 14, 42, 15, 11, 51]
            + [10, 45, 6, 52, 53, 11, 7, 16, 54, 12, 44, 47, 55, 46, 8, 12, 17, 49]
            + [13, 18, 22, 56, 18, 23, 14, 57, 48, 58, 51, 19, 19, 24, 53, 60, 20]
            + [59, 15, 50, 55, 21, 62, 61, 25, 52, 20, 16, 26, 17, 21, 54, 63, 22]
            + [57, 64, 31, 27, 65, 27, 28, 32, 66, 56, 59, 23, 24, 28, 29, 58, 68]
            + [67, 61, 33, 34, 29, 30, 70, 69, 25, 60, 63, 71, 35, 26, 30, 62, 36]
            + [65, 31, 37, 67, 36, 64, 32, 38, 37, 66, 69, 33, 71, 38, 39, 68, 34]
            + [35, 70, 39]
        ),
    )

    assert_array_equal(
        np.asarray(self._head_start_end_indexes),
        np.array(
            [
                [0, 3, 8, 13, 18, 21, 26, 34, 42, 50, 55, 60, 68, 76, 84, 89, 94]
                + [102, 110, 118, 123, 126, 131, 136, 141],
                [2, 7, 12, 17, 20, 25]
                + [33, 41, 49, 54, 59, 67, 75, 83, 88, 93, 101, 109, 117, 122, 125]
                + [130, 135, 140, 143],
            ]
        ),
    )
    assert_array_equal(
        self._link_idx_sorted_by_heads,
        np.array(
            [112, 76, 72, 0, 114, 113, 77, 73, 1, 115, 78, 74, 116, 118, 79, 75]
            + [117, 2, 119, 80, 3, 81, 4, 41, 85, 120, 40, 122, 9, 43, 5, 82, 121]
            + [86, 42, 87, 83, 123, 10, 45, 6, 124, 125, 11, 7, 88, 126, 84, 44, 47]
            + [127, 46, 8, 12, 89, 49, 13, 90, 94, 128, 18, 95, 14, 129, 48, 130]
            + [51, 91, 19, 96, 53, 132, 92, 131, 15, 50, 55, 93, 134, 133, 97, 52]
            + [20, 16, 98, 17, 21, 54, 135, 22, 57, 136, 103, 99, 137, 27, 100, 104]
            + [138, 56, 59, 23, 24, 28, 101, 58, 140, 139, 61, 105, 106, 29, 102]
            + [142, 141, 25, 60, 63, 143, 107, 26, 30, 62, 108, 65, 31, 109, 67, 36]
            + [64, 32, 110, 37, 66, 69, 33, 71, 38, 111, 68, 34, 35, 70, 39]
        ),
    )


def test_init_hex():
    spacing = 10

    g = HexModelGrid((5, 3), spacing, node_layout="hex")
    g.status_at_node[g.perimeter_nodes] = g.BC_NODE_IS_FIXED_VALUE
    g.status_at_node[0] = g.BC_NODE_IS_CLOSED
    nodes_n = g.number_of_nodes

    self = FlowRouter(g, surface="soil__elevation", diagonals=True, runoff_rate=2.0)
    g.at_node["topographic__elevation"] = np.array(
        [10., 20., 10.] +
        [10., 0., 5., 10.] +
        [20., 10., 5., 10., 20.] +
        [10., 20., 25., 15.] +
        [5., 0., 5.])

    assert_array_equal(self._nodes, np.arange(19))

    # 1. Initialization of the input parameters
    ###########################################
    # 1.1. Surface where the flow is directed
    assert self._surface == "soil__elevation"

    # 1.2. Water influx external to grid
    assert self._uniform_water_external_influx
    assert_almost_equal(g.at_node["water__unit_flux_in"][0], 2.0)

    # 1.3. Options
    assert not self._diagonals

    # 2. Boundary conditions
    ########################
    # 2.1. Boundary settings: cell area = 0 and guarantee of one base-level node (at
    # least)
    assert_array_equal(g.closed_boundary_nodes, np.array([0]))
    assert_array_equal(
        g.open_boundary_nodes,
        np.array([ 1,  2,  3,  6,  7, 11, 12, 15, 16, 17, 18]),
    )
    assert_array_equal(self._closed_nodes, np.array([0]))
    assert_array_equal(
        self._base_level_nodes,
        np.array([ 1,  2,  3,  6,  7, 11, 12, 15, 16, 17, 18]),
    )
    assert_array_equal(
        self._base_level_and_closed_nodes,
        np.array([ 1,  2,  3,  6,  7, 11, 12, 15, 16, 17, 18,  0]),
    )

    # 2.2 Cell area at boundary nodes = 0, NetworkModelGrid cell area = 1
    # used in run_flow_accumulations()
    assert_array_almost_equal(
        g.at_node["cell_area_at_node"],
        np.array([0., 0., 0., 0., 86.60254, 86.60254] +
            [0., 0., 86.60254, 86.60254, 86.60254, 0.] +
            [0., 86.60254, 86.60254, 0., 0., 0., 0.])
    )

    # 2.3. Max number of nodes (for sort head/tails/links)
    assert self._max_number_of_nodes == 1e9

    # 3. Determination of stable input grid data necessary to calculate
    # run_flow_directions
    ###################################################################################

    assert self._neighbors_max_number == 6

    # Link infos(tail, head, link id, gradient) sorted by head id.
    # These infos are voluntarily duplicate
    assert_array_equal(
        self._dupli_links,
        np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16] +
            [17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33] +
            [34, 35, 36, 37, 38, 39, 40, 41,  0,  1,  2,  3,  4,  5,  6,  7,  8] +
            [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25] +
            [26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41]
        ),
    )
    assert_array_equal(
        self._pseudo_head_nodes,
        np.array(
            [ 1,  2,  3,  4,  4,  5,  5,  6,  4,  5,  6,  7,  8,  8,  9,  9, 10] +
            [10, 11,  8,  9, 10, 11, 12, 12, 13, 13, 14, 14, 15, 15, 13, 14, 15] +
            [16, 16, 17, 17, 18, 18, 17, 18,  0,  1,  0,  0,  1,  1,  2,  2,  3] +
            [4,  5,  3,  3,  4,  4,  5,  5,  6,  6,  7,  8,  9, 10,  7,  8,  8] +
            [9,  9, 10, 10, 11, 12, 13, 14, 12, 13, 13, 14, 14, 15, 16, 17]
        ),
    )
    assert_array_equal(
        self._pseudo_tail_nodes,
        np.array([ 0,  1,  0,  0,  1,  1,  2,  2,  3,  4,  5,  3,  3,  4,  4,  5,  5] +
            [6,  6,  7,  8,  9, 10,  7,  8,  8,  9,  9, 10, 10, 11, 12, 13, 14] +
            [12, 13, 13, 14, 14, 15, 16, 17,  1,  2,  3,  4,  4,  5,  5,  6,  4] +
            [5,  6,  7,  8,  8,  9,  9, 10, 10, 11,  8,  9, 10, 11, 12, 12, 13] +
            [13, 14, 14, 15, 15, 13, 14, 15, 16, 16, 17, 17, 18, 18, 17, 18]),
    )

    assert_array_equal(
        self._sorted_pseudo_heads,
        np.array([ 0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  3,  3,  3,  3,  4,  4,  4] +
            [4,  4,  4,  5,  5,  5,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  8] +
            [8,  8,  8,  8,  8,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 10] +
            [11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14] +
            [14, 14, 15, 15, 15, 15, 16, 16, 16, 17, 17, 17, 17, 18, 18, 18]),
    )
    assert_array_equal(
        self._sorted_pseudo_tails,
        np.array([ 4,  3,  1,  0,  5,  4,  2,  1,  6,  5,  4,  7,  8,  0,  5,  0,  1] +
            [3,  9,  8, 10,  6,  1,  2,  9,  4, 11,  5, 10,  2, 12,  8,  3, 13] +
            [12,  3,  9,  7,  4, 10, 13,  8,  5,  4, 14, 14, 15,  9,  6,  5, 11] +
            [15, 10,  6, 13, 16,  7,  8, 16, 17, 14,  8,  9, 12,  9, 10, 17, 13] +
            [18, 15, 10, 14, 18, 11, 12, 13, 17, 18, 14, 13, 16, 15, 14, 17]),
    )
    assert_array_equal(
        self._sorted_dupli_links,
        np.array([ 3,  2,  0,  0,  5,  4,  1,  1,  7,  6,  8, 11, 12,  2,  9,  3,  4] +
                [8, 14, 13, 16, 10,  5,  6, 15,  9, 18, 10, 17,  7, 23, 19, 11, 25] +
                [24, 12, 20, 19, 13, 21, 26, 20, 15, 14, 27, 28, 29, 21, 17, 16, 22] +
                [30, 22, 18, 31, 34, 23, 24, 35, 36, 32, 25, 26, 31, 27, 28, 37, 32] +
                [38, 33, 29, 33, 39, 30, 34, 35, 40, 41, 37, 36, 40, 39, 38, 41]
        ),
    )

    assert_array_equal(
        np.asarray(self._head_start_end_indexes),
        np.array([
            [0,  3,  7, 10, 14, 20, 26, 30, 33, 39, 45, 51, 54, 58, 64, 70, 74] +
            [77, 81],
            [2,  6,  9, 13, 19, 25, 29, 32, 38, 44, 50, 53, 57, 63, 69, 73, 76] +
            [80, 83]]),
    )
    assert_array_equal(
        self._link_idx_sorted_by_heads,
        np.array(
            [45, 44, 42,  0, 47, 46, 43,  1, 49, 48, 50, 53, 54,  2, 51,  3,  4] +
            [8, 56, 55, 58, 52,  5,  6, 57,  9, 60, 10, 59,  7, 65, 61, 11, 67] +
            [66, 12, 62, 19, 13, 63, 68, 20, 15, 14, 69, 70, 71, 21, 17, 16, 64] +
            [72, 22, 18, 73, 76, 23, 24, 77, 78, 74, 25, 26, 31, 27, 28, 79, 32] +
            [80, 75, 29, 33, 81, 30, 34, 35, 82, 83, 37, 36, 40, 39, 38, 41]),
    )


def test_init_network():
    params = {
        "yx_of_node": (
            (0, 100, 200, 200, 300, 400, 400, 125),
            (0, 0, 100, -50, -100, 50, -150, -100),
        ),
        "links": ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6)),
    }
    g = NetworkModelGrid(**params)
    self = FlowRouter(g, runoff_rate="water__unit_flux_in")
    g.at_node["topographic__elevation"] = [0.0, 0.08, 0.25, 0.15, 0.25, 0.4, 0.8, 0.8]
    g.at_node["water__unit_flux_in"] = [1.0, 2.0, 3.0, 4.0, 3.0, 3.2, 4.5, 1.0]

    g.status_at_node[7] = g.BC_NODE_IS_CLOSED

    assert_array_equal(self._nodes, np.arange(8))

    # 1. Initialization of the input parameters
    ###########################################
    # 1.1. Surface where the flow is directed
    assert self._surface == "topographic__elevation"

    # 1.2. Water influx external to grid
    assert not self._uniform_water_external_influx
    assert_almost_equal(
        g.at_node["water__unit_flux_in"], [1.0, 2.0, 3.0, 4.0, 3.0, 3.2, 4.5, 1.0]
    )

    # 1.3. Options
    assert not self._diagonals

    # 2. Boundary conditions
    ########################
    # 2.1. Boundary settings: cell area = 0 and guarantee of one base-level
    # node (at least)
    assert_array_equal(g.status_at_node[g.BC_NODE_IS_CLOSED], np.array([]))
    assert_array_equal(g.status_at_node[g.BC_NODE_IS_FIXED_VALUE], np.array([0]))
    assert_array_equal(self._closed_nodes, np.array([]))
    assert_array_equal(self._base_level_nodes, np.array([0]))
    assert_array_equal(self._base_level_and_closed_nodes, np.array([0]))

    # 2.2 Cell area at boundary nodes = 0, NetworkModelGrid cell area = 1
    # used in run_flow_accumulations()
    assert_array_equal(
        g.at_node["cell_area_at_node"],
        np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]),
    )

    # 2.3. Max number of nodes (for sort head/tails/links)
    assert self._max_number_of_nodes == 1e9

    # 3. Determination of stable input grid data necessary to calculate
    # run_flow_directions
    ###################################################################

    # Link infos(tail, head, link id, gradient) sorted by head id.
    # These infos are voluntarily duplicate
    assert_array_equal(
        self._dupli_links, np.array([0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6])
    )
    assert_array_equal(
        self._pseudo_head_nodes, np.array([1, 1, 3, 4, 5, 6, 7, 0, 2, 1, 1, 3, 5, 5])
    )
    assert_array_equal(
        self._pseudo_tail_nodes, np.array([0, 2, 1, 1, 3, 5, 5, 1, 1, 3, 4, 5, 6, 7])
    )

    assert_array_equal(
        self._sorted_pseudo_heads, np.array([0, 1, 1, 1, 1, 2, 3, 3, 4, 5, 5, 5, 6, 7])
    )
    assert_array_equal(
        self._sorted_pseudo_tails, np.array([1, 0, 2, 3, 4, 1, 1, 5, 1, 3, 6, 7, 5, 5])
    )
    assert_array_equal(
        self._sorted_dupli_links, np.array([0, 0, 1, 2, 3, 1, 2, 4, 3, 4, 5, 6, 5, 6])
    )

    assert_array_equal(
        np.asarray(self._head_start_end_indexes),
        np.array([[0, 1, 5, 6, 8, 9, 12, 13], [0, 4, 5, 7, 8, 11, 12, 13]]),
    )
    assert_array_equal(
        self._link_idx_sorted_by_heads,
        np.array([7, 0, 1, 9, 10, 8, 2, 11, 3, 4, 12, 13, 5, 6]),
    )

    assert g.status_at_node[0] == g.BC_NODE_IS_FIXED_VALUE
    assert self._neighbors_max_number == 4
