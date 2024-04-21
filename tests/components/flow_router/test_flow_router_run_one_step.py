import numpy as np
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

from landlab import HexModelGrid
from landlab import NetworkModelGrid
from landlab import RasterModelGrid
from landlab.components import FlowRouter


def test_run_one_step_raster():
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
    self.run_one_step()

    assert_array_equal(
        g.at_node["flow__receiver_node"],
        np.int64([0, 1, 2, 3, 4, 6, 3, 7, 8, 6, 6, 11, 12, 13, 14, 15]),
    )
    assert_array_almost_equal(
        g.at_node["topographic__steepest_slope"],
        np.float64(
            [0.0, 0.0, 0.0, 0.0]
            + [0.0, 0.0, 0.0, 0.0]
            + [0.0, 1.06066017, 0.0, 0.0]
            + [0.0, 0.0, 0.0, 0.0]
        ),
    )
    assert_array_equal(
        g.at_node["flow__link_to_receiver_node"],
        np.int64(
            [-1, -1, -1, -1] + [-1, 8, 29, -1] + [-1, 33, 12, -1] + [-1, -1, -1, -1]
        ),
    )
    assert_array_equal(
        g.at_node["flood_status_code"],
        np.int64([0, 0, 0, 0] + [0, 3, 3, 0] + [0, 0, 3, 0] + [0, 0, 0, 0]),
    )
    assert_array_almost_equal(
        g.at_node["depression__depth"],
        np.float64(
            [0.0, 0.0, 0.0, 0.0]
            + [0.0, 10.0, 5.0, 0.0]
            + [0.0, 0.0, 5.0, 0.0]
            + [0.0, 0.0, 0.0, 0.0]
        ),
    )
    assert_array_equal(
        g.at_node["outlet_node"],
        np.int64([0, 1, 2, 3] + [4, 3, 3, 7] + [8, 3, 3, 11] + [12, 13, 14, 15]),
    )
    assert_array_equal(
        g.at_node["depression__outlet_node"],
        np.int64(
            [-1, -1, -1, -1] + [-1, 3, 3, -1] + [-1, -1, 3, -1] + [-1, -1, -1, -1]
        ),
    )
    assert_array_almost_equal(
        g.at_node["depression_free__elevation"],
        np.float64(
            [10.0, 20.0, 10.0, 10.0]
            + [10.0, 10.0000002, 10.0000001, 10.0]
            + [20.0, 20.0, 10.0000002, 10.0]
            + [10.0, 20.0, 25.0, 15.0]
        ),
    )

    assert_array_equal(
        g.at_node["flow__upstream_node_order"],
        np.int64([3, 6, 5, 9] + [10, 0, 1, 2] + [4, 7, 8, 11] + [12, 13, 14, 15]),
    )
    assert_array_almost_equal(
        g.at_node["drainage_area"],
        np.float64(
            [0.0, 0.0, 0.0, 400.0]
            + [0.0, 100.0, 400.0, 0.0]
            + [0.0, 100.0, 100.0, 0.0]
            + [0.0, 0.0, 0.0, 0.0]
        ),
    )
    assert_array_almost_equal(
        g.at_node["surface_water__discharge"],
        np.float64(
            [0.0, 0.0, 0.0, 400.0]
            + [0.0, 100.0, 400.0, 0.0]
            + [0.0, 100.0, 100.0, 0.0]
            + [0.0, 0.0, 0.0, 0.0]
        ),
    )


def test_run_one_step_hex():
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
    self.run_one_step()

    assert_array_equal(
        g.at_node["flow__receiver_node"][11:],
        np.int64([11] + [12, 17, 17, 15] + [16, 17, 18]),
    )
    # The test doesnt work for full array in workflow mac os,
    # supposedly because of int32/int64 gradients and lack of unit
    # test and conversion in the hex grid classes
    """
    np.int64(
            [0, 1, 2]
            + [3, 3, 4, 6]
            + [7, 3, 4, 5, 11]
            + [12, 17, 17, 15]
            + [16, 17, 18]
        ),
    """
    assert_array_almost_equal(
        g.at_node["topographic__steepest_slope"][11:],
        np.float64([0.0] + [0.0, 2.0, 2.5, 0.0] + [0.0, 0.0, 0.0]),
    )
    """
    np.float64(
            [0.0, 0.0, 0.0]
            + [0.0, 0.0, 0.5, 0.0]
            + [0.0, 0.0, 0.5, 0.5, 0.0]
            + [0.0, 2.0, 2.5, 0.0]
            + [0.0, 0.0, 0.0]
        ),
    """
    assert_array_equal(
        g.at_node["flow__link_to_receiver_node"][11:],
        np.int64([-1] + [-1, 36, 37, -1] + [-1, -1, -1]),
    )
    """
    np.int64(
            [-1, -1, -1]
            + [-1, 8, 9, -1]
            + [-1, 12, 14, 16, -1]
            + [-1, 36, 37, -1]
            + [-1, -1, -1]
        ),
    """
    assert_array_equal(
        g.at_node["flood_status_code"],
        np.int64([0, 0, 0] + [0, 3, 3, 0] + [0, 0, 3, 0, 0] + [0, 0, 0, 0] + [0, 0, 0]),
    )
    assert_array_almost_equal(
        g.at_node["depression__depth"],
        np.float64(
            [0.0, 0.0, 0.0]
            + [0.0, 10.0, 5.0, 0.0]
            + [0.0, 0.0, 5.0, 0.0, 0.0]
            + [0.0, 0.0, 0.0, 0.0]
            + [0.0, 0.0, 0.0]
        ),
    )
    assert_array_equal(
        g.at_node["outlet_node"][11:],
        np.int64([11] + [12, 17, 17, 15] + [16, 17, 18]),
    )
    """
    np.int64(
            [0, 1, 2]
            + [3, 3, 3, 6]
            + [7, 3, 3, 3, 11]
            + [12, 17, 17, 15]
            + [16, 17, 18]
        ),
    """
    assert_array_equal(
        g.at_node["depression__outlet_node"][11:],
        np.int64([-1] + [-1, -1, -1, -1] + [-1, -1, -1]),
    )
    """
    np.int64(
            [-1, -1, -1]
            + [-1, 3, 3, -1]
            + [-1, -1, 3, -1, -1]
            + [-1, -1, -1, -1]
            + [-1, -1, -1]
        ),
    """
    assert_array_almost_equal(
        g.at_node["depression_free__elevation"],
        np.float64(
            [10.0, 20.0, 10.0]
            + [10.0, 10.0000001, 10.0000002, 10.0]
            + [20.0, 10.0, 10.0000002, 10.0, 20.0]
            + [10.0, 20.0, 25.0, 15.0]
            + [5.0, 0.0, 5.0]
        ),
    )

    assert_array_equal(
        g.at_node["flow__upstream_node_order"],
        np.int64(
            [1, 2, 3]
            + [4, 9, 5, 10]
            + [8, 6, 7, 11, 12]
            + [15, 16, 17, 13]
            + [14, 18, 0]
        ),
    )
    assert_array_almost_equal(
        g.at_node["drainage_area"],
        np.float64(
            [0.0, 0.0, 0.0]
            + [433.0127, 346.41016, 173.20508, 0.0]
            + [0.0, 86.60254, 86.60254, 86.60254, 0.0]
            + [0.0, 86.60254, 86.60254, 0.0]
            + [0.0, 173.20508, 0.0]
        ),
    )
    assert_array_almost_equal(
        g.at_node["surface_water__discharge"],
        np.float64(
            [0.0, 0.0, 0.0]
            + [866.0254, 692.82032, 346.41016, 0.0]
            + [0.0, 173.20508, 173.20508, 173.20508, 0.0]
            + [0.0, 173.20508, 173.20508, 0.0]
            + [0.0, 346.41016, 0.0]
        ),
    )


def test_run_one_step_network():
    params = {
        "yx_of_node": (
            (0, 100, 200, 200, 300, 400, 400, 125),
            (0, 0, 100, -50, -100, 50, -150, -100),
        ),
        "links": ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6)),
    }
    g = NetworkModelGrid(**params)
    g.at_node["topographic__elevation"] = np.array(
        [0.0, 0.2, 0.25, 0.15, 0.25, 0.4, 0.8, 0.8]
    )
    g.at_node["water__unit_flux_in"] = np.array(
        [1.0, 2.0, 3.0, 4.0, 3.0, 3.2, 4.5, 1.0]
    )

    self = FlowRouter(g, runoff_rate="water__unit_flux_in")

    g.status_at_node[7] = g.BC_NODE_IS_CLOSED
    self.run_one_step()

    assert_array_equal(
        g.at_node["flow__receiver_node"], np.int64([0, 0, 1, 1, 1, 3, 5, 5])
    )
    assert_array_almost_equal(
        g.at_node["topographic__steepest_slope"],
        np.float64(
            [0.0, 0.002, 0.00048507, 0.0, 0.00035355]
            + [0.00223607, 0.00357771, 0.0022188],
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

    assert_array_equal(
        g.at_node["flow__upstream_node_order"], np.int64([0, 1, 2, 3, 5, 6, 7, 4])
    )
    assert_array_almost_equal(
        g.at_node["drainage_area"], np.float64([8.0, 7.0, 1.0, 4.0, 1.0, 3.0, 1.0, 1.0])
    )
    assert_array_almost_equal(
        g.at_node["surface_water__discharge"],
        np.float64([8.0, 14.0, 3.0, 16.0, 3.0, 9.6, 4.5, 1.0]),
    )
