import numpy as np
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

from landlab import HexModelGrid
from landlab import NetworkModelGrid
from landlab import RasterModelGrid
from landlab.components import FlowRouter


def test_run_flow_accumulations_raster():
    spacing = 100
    g = RasterModelGrid((5, 5), (spacing, spacing))
    g.status_at_node[g.perimeter_nodes] = g.BC_NODE_IS_CLOSED
    self = FlowRouter(g)
    g.at_node["topographic__elevation"] = np.float64(
        [10, 10, 10, 10, 10]
        + [20, 20, 0, 20, 20]
        + [30, 0, 10, 20, 10]
        + [20, 20, 30, 20, 10]
        + [0, 30, 0, 0, 0]
    )
    self.run_flow_directions()
    self.run_flow_accumulations()

    assert_array_equal(
        g.at_node["flow__upstream_node_order"],
        np.int64(
            [4, 8, 13, 12, 18, 7, 6, 11, 16, 17, 0, 1, 2, 3, 5, 9, 10]
            + [14, 15, 19, 20, 21, 22, 23, 24]
        ),
    )
    assert_array_almost_equal(
        g.at_node["drainage_area"],
        np.float64(
            [0.0, 0.0, 0.0, 0.0, 90000.0]
            + [0.0, 10000.0, 50000.0, 90000.0, 0.0]
            + [0.0, 30000.0, 20000.0, 10000.0, 0.0]
            + [0.0, 10000.0, 10000.0, 10000.0, 0.0]
            + [0.0, 0.0, 0.0, 0.0, 0.0]
        ),
    )
    assert_array_almost_equal(
        g.at_node["surface_water__discharge"],
        np.float64(
            [0.0, 0.0, 0.0, 0.0, 90000.0]
            + [0.0, 10000.0, 50000.0, 90000.0, 0.0]
            + [0.0, 30000.0, 20000.0, 10000.0, 0.0]
            + [0.0, 10000.0, 10000.0, 10000.0, 0.0]
            + [0.0, 0.0, 0.0, 0.0, 0.0]
        ),
    )


def test_run_flow_accumulations_hex():
    spacing = 10

    g = HexModelGrid((5, 3), spacing, node_layout="hex")
    g.status_at_node[g.perimeter_nodes] = g.BC_NODE_IS_FIXED_VALUE
    g.status_at_node[0] = g.BC_NODE_IS_CLOSED

    self = FlowRouter(g, surface="soil__elevation", diagonals=True, runoff_rate=2.0)
    g.at_node["soil__elevation"] = np.float64(
        [10.0, 20.0, 10.0]
        + [10.0, 0.0, 5.0, 10.0]
        + [20.0, 10.0, 5.0, 10.0, 20.0]
        + [10.0, 20.0, 25.0, 15.0]
        + [5.0, 0.0, 5.0]
    )
    self.run_flow_directions()
    self.run_flow_accumulations()

    assert_array_equal(
        g.at_node["flow__upstream_node_order"],
        np.int64([1, 2, 3, 4, 9, 5, 10, 8, 6, 7, 11, 12, 15, 16, 17, 13, 14] + [18, 0]),
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


def test_run_flow_accumulations_network():
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
    self.run_flow_accumulations()

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
