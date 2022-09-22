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
    self = FlowRouter(g)
    g.at_node["topographic__elevation"] = np.array(
        [10.0, 20.0, 10.0, 10.0]
        + [10.0, 0.0, 5.0, 10.0]
        + [20.0, 20.0, 5.0, 10.0]
        + [10.0, 20.0, 25.0, 15.0]
    )

    assert_array_equal(self._nodes, np.arange(16))

    # 1. Initialization of the input parameters
    ###########################################
    # 1.1. Surface where the flow is directed
    assert self._surface == "topographic__elevation"

    # 1.2. Water influx external to grid
    assert self._uniform_water_external_influx
    assert_almost_equal(
        g.at_node["water__unit_flux_in"],
        np.array(
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
            + [1.0, 1.0, 1.0]
        ),
    )

    # 1.3. Options
    assert self._diagonals

    # 2. Boundary conditions
    ########################
    # 2.1. Boundary settings: guarantee of one base-level node (at least)
    assert_array_equal(
        g.closed_boundary_nodes,
        np.array([0, 1, 2, 4, 7, 8, 11, 12, 13, 14, 15]),
    )
    assert_array_equal(g.open_boundary_nodes, np.array([3]))
    assert_array_equal(
        self._closed_nodes,
        np.array([0, 1, 2, 4, 7, 8, 11, 12, 13, 14, 15]),
    )
    assert_array_equal(self._base_level_nodes, np.array([3]))
    assert_array_equal(
        self._base_level_and_closed_nodes,
        np.array([3, 0, 1, 2, 4, 7, 8, 11, 12, 13, 14, 15]),
    )

    # 2.2 Cell area at boundary nodes = 0, NetworkModelGrid cell area = 1
    # used in run_flow_accumulations()
    assert_array_almost_equal(
        g.at_node["cell_area_at_node"],
        np.array(
            [0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 100.0, 0.0, 0.0]
            + [100.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ),
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
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
            + [17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33]
            + [34, 35, 36, 37, 38, 39, 40, 41, 0, 1, 2, 3, 4, 5, 6, 7, 8]
            + [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
            + [26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41]
        ),
    )
    assert_array_equal(
        self._pseudo_head_nodes,
        np.array(
            [1, 2, 3, 4, 5, 6, 7, 5, 6, 7, 8, 9, 10, 11, 9, 10, 11]
            + [12, 13, 14, 15, 13, 14, 15, 5, 4, 6, 5, 7, 6, 9, 8, 10, 9]
            + [11, 10, 13, 12, 14, 13, 15, 14, 0, 1, 2, 0, 1, 2, 3, 4, 5]
            + [6, 4, 5, 6, 7, 8, 9, 10, 8, 9, 10, 11, 12, 13, 14, 0, 1]
            + [1, 2, 2, 3, 4, 5, 5, 6, 6, 7, 8, 9, 9, 10, 10, 11]
        ),
    )
    assert_array_equal(
        self._pseudo_tail_nodes,
        np.array(
            [0, 1, 2, 0, 1, 2, 3, 4, 5, 6, 4, 5, 6, 7, 8, 9, 10]
            + [8, 9, 10, 11, 12, 13, 14, 0, 1, 1, 2, 2, 3, 4, 5, 5, 6]
            + [6, 7, 8, 9, 9, 10, 10, 11, 1, 2, 3, 4, 5, 6, 7, 5, 6]
            + [7, 8, 9, 10, 11, 9, 10, 11, 12, 13, 14, 15, 13, 14, 15, 5, 4]
            + [6, 5, 7, 6, 9, 8, 10, 9, 11, 10, 13, 12, 14, 13, 15, 14]
        ),
    )

    assert_array_equal(
        self._sorted_pseudo_heads,
        np.array(
            [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4]
            + [4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6]
            + [6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9]
            + [9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11]
            + [12, 12, 12, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 15, 15, 15]
        ),
    )
    assert_array_equal(
        self._sorted_pseudo_tails,
        np.array(
            [5, 4, 1, 5, 2, 4, 6, 0, 3, 5, 7, 1, 6, 6, 2, 7, 1]
            + [9, 0, 8, 5, 4, 10, 8, 0, 1, 2, 9, 6, 7, 2, 10, 5, 1]
            + [9, 11, 3, 2, 11, 3, 6, 10, 9, 12, 5, 4, 13, 6, 13, 10, 5]
            + [8, 12, 4, 14, 9, 15, 13, 6, 11, 14, 5, 7, 14, 7, 15, 10, 6]
            + [8, 9, 13, 14, 10, 8, 9, 12, 10, 11, 13, 9, 15, 10, 11, 14]
        ),
    )
    assert_array_equal(
        self._sorted_dupli_links,
        np.array(
            [24, 3, 0, 4, 1, 25, 26, 0, 2, 27, 28, 1, 5, 29, 2, 6, 25]
            + [30, 3, 10, 7, 7, 32, 31, 24, 4, 27, 11, 8, 9, 5, 12, 8, 26]
            + [33, 34, 29, 28, 13, 6, 9, 35, 14, 17, 31, 10, 36, 33, 18, 15, 11]
            + [14, 37, 30, 38, 15, 40, 39, 12, 16, 19, 32, 35, 41, 13, 20, 16, 34]
            + [17, 37, 21, 22, 39, 36, 18, 21, 19, 41, 22, 38, 23, 40, 20, 23]
        ),
    )

    assert_array_equal(
        np.asarray(self._head_start_end_indexes),
        np.array(
            [
                [0, 3, 8, 13, 16, 21, 29, 37, 42, 47, 55, 63, 68, 71, 76, 81],
                [2, 7, 12, 15, 20, 28, 36, 41, 46, 54, 62, 67, 70, 75, 80, 83],
            ]
        ),
    )
    assert_array_equal(
        self._link_idx_sorted_by_heads,
        np.array(
            [66, 45, 42, 46, 43, 67, 68, 0, 44, 69, 70, 1, 47, 71, 2, 48, 25]
            + [72, 3, 52, 49, 7, 74, 73, 24, 4, 27, 53, 50, 51, 5, 54, 8, 26]
            + [75, 76, 29, 28, 55, 6, 9, 77, 56, 59, 31, 10, 78, 33, 60, 57, 11]
            + [14, 79, 30, 80, 15, 82, 81, 12, 58, 61, 32, 35, 83, 13, 62, 16, 34]
            + [17, 37, 63, 64, 39, 36, 18, 21, 19, 41, 22, 38, 65, 40, 20, 23]
        ),
    )


def test_init_hex():
    spacing = 10

    g = HexModelGrid((5, 3), spacing, node_layout="hex")
    g.status_at_node[g.perimeter_nodes] = g.BC_NODE_IS_FIXED_VALUE
    g.status_at_node[0] = g.BC_NODE_IS_CLOSED

    self = FlowRouter(g, surface="soil__elevation", diagonals=True, runoff_rate=2.0)
    g.at_node["soil__elevation"] = np.array(
        [10.0, 20.0, 10.0]
        + [10.0, 0.0, 5.0, 10.0]
        + [20.0, 10.0, 5.0, 10.0, 20.0]
        + [10.0, 20.0, 25.0, 15.0]
        + [5.0, 0.0, 5.0]
    )

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
        np.array([1, 2, 3, 6, 7, 11, 12, 15, 16, 17, 18]),
    )
    assert_array_equal(self._closed_nodes, np.array([0]))
    assert_array_equal(
        self._base_level_nodes,
        np.array([1, 2, 3, 6, 7, 11, 12, 15, 16, 17, 18]),
    )
    assert_array_equal(
        self._base_level_and_closed_nodes,
        np.array([1, 2, 3, 6, 7, 11, 12, 15, 16, 17, 18, 0]),
    )

    # 2.2 Cell area at boundary nodes = 0, NetworkModelGrid cell area = 1
    # used in run_flow_accumulations()
    assert_array_almost_equal(
        g.at_node["cell_area_at_node"],
        np.array(
            [0.0, 0.0, 0.0, 0.0, 86.60254, 86.60254]
            + [0.0, 0.0, 86.60254, 86.60254, 86.60254, 0.0]
            + [0.0, 86.60254, 86.60254, 0.0, 0.0, 0.0, 0.0]
        ),
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
        np.array(
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
            + [17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33]
            + [34, 35, 36, 37, 38, 39, 40, 41, 0, 1, 2, 3, 4, 5, 6, 7, 8]
            + [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
            + [26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41]
        ),
    )
    assert_array_equal(
        self._pseudo_head_nodes,
        np.array(
            [1, 2, 3, 4, 4, 5, 5, 6, 4, 5, 6, 7, 8, 8, 9, 9, 10]
            + [10, 11, 8, 9, 10, 11, 12, 12, 13, 13, 14, 14, 15, 15, 13, 14, 15]
            + [16, 16, 17, 17, 18, 18, 17, 18, 0, 1, 0, 0, 1, 1, 2, 2, 3]
            + [4, 5, 3, 3, 4, 4, 5, 5, 6, 6, 7, 8, 9, 10, 7, 8, 8]
            + [9, 9, 10, 10, 11, 12, 13, 14, 12, 13, 13, 14, 14, 15, 16, 17]
        ),
    )
    assert_array_equal(
        self._pseudo_tail_nodes,
        np.array(
            [0, 1, 0, 0, 1, 1, 2, 2, 3, 4, 5, 3, 3, 4, 4, 5, 5]
            + [6, 6, 7, 8, 9, 10, 7, 8, 8, 9, 9, 10, 10, 11, 12, 13, 14]
            + [12, 13, 13, 14, 14, 15, 16, 17, 1, 2, 3, 4, 4, 5, 5, 6, 4]
            + [5, 6, 7, 8, 8, 9, 9, 10, 10, 11, 8, 9, 10, 11, 12, 12, 13]
            + [13, 14, 14, 15, 15, 13, 14, 15, 16, 16, 17, 17, 18, 18, 17, 18]
        ),
    )

    assert_array_equal(
        self._sorted_pseudo_heads,
        np.array(
            [0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4]
            + [4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8]
            + [8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10]
            + [11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14]
            + [14, 14, 15, 15, 15, 15, 16, 16, 16, 17, 17, 17, 17, 18, 18, 18]
        ),
    )
    assert_array_equal(
        self._sorted_pseudo_tails,
        np.array(
            [4, 3, 1, 0, 5, 4, 2, 1, 6, 5, 4, 7, 8, 0, 5, 0, 1]
            + [3, 9, 8, 10, 6, 1, 2, 9, 4, 11, 5, 10, 2, 12, 8, 3, 13]
            + [12, 3, 9, 7, 4, 10, 13, 8, 5, 4, 14, 14, 15, 9, 6, 5, 11]
            + [15, 10, 6, 13, 16, 7, 8, 16, 17, 14, 8, 9, 12, 9, 10, 17, 13]
            + [18, 15, 10, 14, 18, 11, 12, 13, 17, 18, 14, 13, 16, 15, 14, 17]
        ),
    )
    assert_array_equal(
        self._sorted_dupli_links,
        np.array(
            [3, 2, 0, 0, 5, 4, 1, 1, 7, 6, 8, 11, 12, 2, 9, 3, 4]
            + [8, 14, 13, 16, 10, 5, 6, 15, 9, 18, 10, 17, 7, 23, 19, 11, 25]
            + [24, 12, 20, 19, 13, 21, 26, 20, 15, 14, 27, 28, 29, 21, 17, 16, 22]
            + [30, 22, 18, 31, 34, 23, 24, 35, 36, 32, 25, 26, 31, 27, 28, 37, 32]
            + [38, 33, 29, 33, 39, 30, 34, 35, 40, 41, 37, 36, 40, 39, 38, 41]
        ),
    )

    assert_array_equal(
        np.asarray(self._head_start_end_indexes),
        np.array(
            [
                [0, 3, 7, 10, 14, 20, 26, 30, 33, 39, 45, 51, 54, 58, 64, 70, 74]
                + [77, 81],
                [2, 6, 9, 13, 19, 25, 29, 32, 38, 44, 50, 53, 57, 63, 69, 73, 76]
                + [80, 83],
            ]
        ),
    )
    assert_array_equal(
        self._link_idx_sorted_by_heads,
        np.array(
            [45, 44, 42, 0, 47, 46, 43, 1, 49, 48, 50, 53, 54, 2, 51, 3, 4]
            + [8, 56, 55, 58, 52, 5, 6, 57, 9, 60, 10, 59, 7, 65, 61, 11, 67]
            + [66, 12, 62, 19, 13, 63, 68, 20, 15, 14, 69, 70, 71, 21, 17, 16, 64]
            + [72, 22, 18, 73, 76, 23, 24, 77, 78, 74, 25, 26, 31, 27, 28, 79, 32]
            + [80, 75, 29, 33, 81, 30, 34, 35, 82, 83, 37, 36, 40, 39, 38, 41]
        ),
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
