import numpy as np
from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

from landlab import HexModelGrid
from landlab import NetworkModelGrid
from landlab import RasterModelGrid
from landlab.components import FlowRouter


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
        self._dupli_links[np.r_[0:2, 40:44, 82:84]],
        np.array([0, 1, 40, 41, 0, 1, 40, 41]),
    )

    assert_array_equal(
        self._pseudo_head_nodes[0:4],
        np.array([1, 2, 3, 4]),
    )
    assert_array_equal(
        self._pseudo_tail_nodes[0:4],
        np.array([0, 1, 2, 0]),
    )

    assert_array_equal(
        self._sorted_pseudo_heads[np.r_[0:9, 75:84]],
        np.array([0, 0, 0, 1, 1, 1, 1, 1, 2] + [13, 14, 14, 14, 14, 14, 15, 15, 15]),
    )
    assert_array_equal(
        self._sorted_pseudo_tails[np.r_[0:9, 75:84]],
        np.array([5, 4, 1, 5, 2, 4, 6, 0, 3] + [12, 10, 11, 13, 9, 15, 10, 11, 14]),
    )
    assert_array_equal(
        self._sorted_dupli_links[np.r_[0:3, 82:84]],
        np.array([24, 3, 0, 20, 23]),
    )

    assert_array_equal(
        np.asarray(self._head_start_end_indexes)[0, 0:4],
        np.array([0, 3, 8, 13]),
    )
    assert_array_equal(
        self._link_idx_sorted_by_heads[0:4],
        np.array([66, 45, 42, 46]),
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
        self._dupli_links[np.r_[0:2, 40:44, 82:84]],
        np.array([0, 1, 40, 41, 0, 1, 40, 41]),
    )
    assert_array_equal(
        self._pseudo_head_nodes[0:4],
        np.array([1, 2, 3, 4]),
    )
    assert_array_equal(
        self._pseudo_tail_nodes[0:4],
        np.array([0, 1, 0, 0]),
    )

    assert_array_equal(
        self._sorted_pseudo_heads[np.r_[0:7, 77:84]],
        np.array([0, 0, 0, 1, 1, 1, 1] + [17, 17, 17, 17, 18, 18, 18]),
    )
    assert_array_equal(
        self._sorted_pseudo_tails[np.r_[0:7, 77:84]],
        np.array([4, 3, 1, 0, 5, 4, 2] + [18, 14, 13, 16, 15, 14, 17]),
    )
    assert_array_equal(
        self._sorted_dupli_links[np.r_[0:3, 82:84]],
        np.array([3, 2, 0] + [38, 41]),
    )

    assert_array_equal(
        np.asarray(self._head_start_end_indexes)[0, 0:3],
        np.array([0, 3, 7]),
    )
    assert_array_equal(
        self._link_idx_sorted_by_heads[0:4],
        np.array([45, 44, 42, 0]),
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
    g.at_node["topographic__elevation"] = [0.0, 0.08, 0.25, 0.15, 0.25, 0.4, 0.8, 0.8]
    g.at_node["water__unit_flux_in"] = [1.0, 2.0, 3.0, 4.0, 3.0, 3.2, 4.5, 1.0]

    self = FlowRouter(g, runoff_rate="water__unit_flux_in")

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
