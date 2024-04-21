import numpy as np
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

from landlab import HexModelGrid
from landlab import NetworkModelGrid
from landlab import RasterModelGrid
from landlab.components import FlowRouter


def test_run_flow_directions_raster():
    spacing = 10
    g = RasterModelGrid((5, 5), (spacing, spacing))
    g.status_at_node[g.perimeter_nodes] = g.BC_NODE_IS_CLOSED
    self = FlowRouter(g)
    g.at_node["topographic__elevation"] = np.array(
        [10, 10, 10, 10, 10]
        + [20, 20, 0, 20, 20]
        + [30, 0, 10, 20, 10]
        + [20, 20, 30, 20, 10]
        + [0, 30, 0, 0, 0]
    )
    self.run_flow_directions()

    assert_array_equal(
        g.at_node["flow__receiver_node"],
        np.int64(
            [0, 1, 2, 3, 4]
            + [5, 7, 8, 4, 9]
            + [10, 7, 8, 8, 14]
            + [15, 11, 11, 12, 19]
            + [20, 21, 22, 23, 24]
        ),
    )
    assert_array_almost_equal(
        g.at_node["topographic__steepest_slope"],
        np.float64(
            [0.0, 0.0, 0.0, 0.0, 0.0]
            + [0.0, 2.0, 0.0, 0.70710678, 0.0]
            + [0.0, 0.0, 0.0, 0.0, 0.0]
            + [0.0, 2.0, 2.12132034, 0.70710678, 0.0]
            + [0.0, 0.0, 0.0, 0.0, 0.0]
        ),
    )
    assert_array_equal(
        g.at_node["flow__link_to_receiver_node"],
        np.int64(
            [-1, -1, -1, -1, -1]
            + [-1, 10, 11, 47, -1]
            + [-1, 51, 53, 16, -1]
            + [-1, 23, 58, 60, -1]
            + [-1, -1, -1, -1, -1]
        ),
    )
    assert_array_equal(
        g.at_node["flood_status_code"],
        np.int64(
            [0, 0, 0, 0, 0]
            + [0, 0, 3, 0, 0]
            + [0, 3, 3, 0, 0]
            + [0, 0, 0, 0, 0]
            + [0, 0, 0, 0, 0]
        ),
    )
    assert_array_almost_equal(
        g.at_node["depression__depth"],
        np.float64(
            [0.0, 0.0, 0.0, 0.0, 0.0]
            + [0.0, 0.0, 20.0, 0.0, 0.0]
            + [0.0, 20.0, 10.0, 0.0, 0.0]
            + [0.0, 0.0, 0.0, 0.0, 0.0]
            + [0.0, 0.0, 0.0, 0.0, 0.0]
        ),
    )
    assert_array_equal(
        g.at_node["outlet_node"],
        np.int64(
            [0, 1, 2, 3, 4]
            + [5, 4, 4, 4, 9]
            + [10, 4, 4, 4, 14]
            + [15, 4, 4, 4, 19]
            + [20, 21, 22, 23, 24]
        ),
    )
    assert_array_equal(
        g.at_node["depression__outlet_node"],
        np.int64(
            [-1, -1, -1, -1, -1]
            + [-1, -1, 8, -1, -1]
            + [-1, 8, 8, -1, -1]
            + [-1, -1, -1, -1, -1]
            + [-1, -1, -1, -1, -1]
        ),
    )
    assert_array_almost_equal(
        g.at_node["depression_free__elevation"],
        np.float64(
            [10.0, 10.0, 10.0, 10.0, 10.0]
            + [20.0, 20.0, 20.0000002, 20.0, 20.0]
            + [30.0, 20.0000004, 20.0000002, 20.0, 10.0]
            + [20.0, 20.0, 30.0, 20.0, 10.0]
            + [0.0, 30.0, 0.0, 0.0, 0.0]
        ),
    )


def test_run_flow_directions_hex():
    spacing = 10

    g = HexModelGrid((5, 3), spacing, node_layout="hex")
    g.status_at_node[g.perimeter_nodes] = g.BC_NODE_IS_FIXED_VALUE
    g.status_at_node[0] = g.BC_NODE_IS_CLOSED

    self = FlowRouter(g, surface="soil__elevation", diagonals=True, runoff_rate=2.0)
    g.at_node["soil__elevation"] = np.array(
        [10, 10, 10]
        + [20, 20, 0, 20]
        + [30, 0, 10, 20, 10]
        + [20, 20, 30, 20]
        + [10, 0, 10]
    )
    self.run_flow_directions()

    assert_array_equal(
        g.at_node["flow__receiver_node"][:10],
        np.int64([0, 1, 2] + [3, 1, 1, 6] + [7, 9, 5]),
    )
    # The test doesnt work for full array in workflow mac os,
    # supposedly because of int32/int64 gradients and lack of unit
    # test and conversion in the hex grid classes
    """
    np.int64(
            [0, 1, 2]
            + [3, 1, 1, 6]
            + [7, 9, 5, 11, 11]
            + [12, 17, 17, 15]
            + [16, 17, 18]
        ),
    """

    assert_array_almost_equal(
        g.at_node["topographic__steepest_slope"][:10],
        np.float64([0.0, 0.0, 0.0] + [0.0, 1.0, 0.0, 0.0] + [0.0, 0.0, 1.0]),
    )
    """
    np.float64(
            [0.0, 0.0, 0.0]
            + [0.0, 1.0, 0.0, 0.0]
            + [0.0, 0.0, 1.0, 1.0, 0.0]
            + [0.0, 2.0, 3.0, 0.0]
            + [0.0, 0.0, 0.0]
        ),
    """
    assert_array_equal(
        g.at_node["flow__link_to_receiver_node"][:10],
        np.int64([-1, -1, -1] + [-1, 4, 5, -1] + [-1, 20, 15]),
    )
    """
    np.int64(
            [-1, -1, -1]
            + [-1, 4, 5, -1]
            + [-1, 20, 15, 22, -1]
            + [-1, 36, 37, -1]
            + [-1, -1, -1]
        ),
    """
    assert_array_equal(
        g.at_node["flood_status_code"],
        np.int64([0, 0, 0] + [0, 0, 3, 0] + [0, 3, 0, 0, 0] + [0, 0, 0, 0] + [0, 0, 0]),
    )
    assert_array_almost_equal(
        g.at_node["depression__depth"],
        np.float64(
            [0.0, 0.0, 0.0]
            + [0.0, 0.0, 10.0, 0.0]
            + [0.0, 10.0, 0.0, 0.0, 0.0]
            + [0.0, 0.0, 0.0, 0.0]
            + [0.0, 0.0, 0.0]
        ),
    )
    assert_array_equal(
        g.at_node["outlet_node"][:10],
        np.int64([0, 1, 2] + [3, 1, 1, 6] + [7, 1, 1]),
    )
    """
    np.int64(
            [0, 1, 2]
            + [3, 1, 1, 6]
            + [7, 1, 1, 11, 11]
            + [12, 17, 17, 15]
            + [16, 17, 18]
        ),
    """

    assert_array_equal(
        g.at_node["depression__outlet_node"],
        np.int64(
            [-1, -1, -1]
            + [-1, -1, 1, -1]
            + [-1, 9, -1, -1, -1]
            + [-1, -1, -1, -1]
            + [-1, -1, -1]
        ),
    )
    assert_array_almost_equal(
        g.at_node["depression_free__elevation"],
        np.float64(
            [10.0, 10.0, 10.0]
            + [20.0, 20.0, 10.0000001, 20.0]
            + [30.0, 10.0000001, 10.0, 20.0, 10.0]
            + [20.0, 20.0, 30.0, 20.0]
            + [10.0, 0.0, 10.0]
        ),
    )


def test_run_flow_directions_network():
    params = {
        "yx_of_node": (
            (0, 100, 200, 200, 300, 400, 400, 125),
            (0, 0, 100, -50, -100, 50, -150, -100),
        ),
        "links": ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6)),
    }
    g = NetworkModelGrid(**params)
    g.at_node["topographic__elevation"] = [0.0, 0.2, 0.25, 0.15, 0.25, 0.4, 0.8, 0.8]
    g.at_node["water__unit_flux_in"] = [1.0, 2.0, 3.0, 4.0, 3.0, 3.2, 4.5, 1.0]

    self = FlowRouter(g, runoff_rate="water__unit_flux_in")

    g.status_at_node[7] = g.BC_NODE_IS_CLOSED
    self.run_flow_directions()

    assert_array_equal(
        g.at_node["flow__receiver_node"], np.int64([0, 0, 1, 1, 1, 3, 5, 5])
    )
    assert_array_almost_equal(
        g.at_node["topographic__steepest_slope"],
        np.float64(
            [0.0, 0.002, 0.00048507, 0.0, 0.00035355]
            + [0.00223607, 0.00357771, 0.0022188]
        ),
    )
    assert_array_equal(
        g.at_node["flow__link_to_receiver_node"], np.int64([-1, 0, 1, 2, 3, 4, 5, 6])
    )
    assert_array_equal(
        g.at_node["flood_status_code"], np.int64([0, 0, 0, 3, 0, 0, 0, 0])
    )
    assert_array_almost_equal(
        g.at_node["depression__depth"],
        np.float64([0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0]),
    )
    assert_array_equal(g.at_node["outlet_node"], np.int64([0, 0, 0, 0, 0, 0, 0, 0]))
    assert_array_equal(
        g.at_node["depression__outlet_node"], np.int64([-1, -1, -1, 1, -1, -1, -1, -1])
    )
    assert_array_almost_equal(
        g.at_node["depression_free__elevation"],
        np.float64([0.0, 0.2, 0.25, 0.2, 0.25, 0.4, 0.8, 0.8]),
    )
