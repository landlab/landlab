import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from landlab.components import FlowRouter

from landlab import RasterModelGrid, HexModelGrid
from landlab import NetworkModelGrid


def test_run_flow_accumulations_raster():
    spacing = 10
    g = RasterModelGrid((5, 5), (spacing, spacing))
    g.status_at_node[g.perimeter_nodes] = g.BC_NODE_IS_CLOSED
    nodes_n = g.number_of_nodes
    self = FlowRouter(g)
    random_generator = np.random.Generator(np.random.PCG64(seed=500))
    g.at_node["topographic__elevation"] = 10 * random_generator.random(nodes_n)
    self.run_flow_directions()
    self.run_flow_accumulations()

    assert_array_equal(
        g.at_node["flow__upstream_node_order"],
        np.array(
            [4, 8, 13, 17, 18, 12, 11, 6, 16, 7, 0, 1, 2, 3, 5, 9, 10]
            + [14, 15, 19, 20, 21, 22, 23, 24]
        ),
    )
    assert_array_almost_equal(
        g.at_node["drainage_area"],
        np.array(
            [0.0, 0.0, 0.0, 0.0, 900.0, 0.0, 100.0, 100.0, 900.0]
            + [0.0, 0.0, 100.0, 400.0, 300.0, 0.0, 0.0, 100.0, 100.0]
            + [100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ),
    )
    assert_array_almost_equal(
        g.at_node["surface_water__discharge"],
        np.array(
            [0.0, 0.0, 0.0, 0.0, 900.0, 0.0, 100.0, 100.0, 900.0]
            + [0.0, 0.0, 100.0, 400.0, 300.0, 0.0, 0.0, 100.0, 100.0]
            + [100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ),
    )


def test_run_flow_accumulations_hex():
    spacing = 10

    g = HexModelGrid((5, 3), spacing, node_layout="hex")
    g.status_at_node[g.perimeter_nodes] = g.BC_NODE_IS_FIXED_VALUE
    g.status_at_node[0] = g.BC_NODE_IS_CLOSED
    nodes_n = g.number_of_nodes

    self = FlowRouter(g, surface="soil__elevation", diagonals=True, runoff_rate=2.0)
    random_generator = np.random.Generator(np.random.PCG64(seed=500))
    g.at_node["soil__elevation"] = 10 * random_generator.random(nodes_n)
    self.run_flow_directions()
    self.run_flow_accumulations()

    assert_array_equal(
        g.at_node["flow__upstream_node_order"],
        np.array([1, 2, 3, 6, 7, 11, 12, 15, 16, 13, 17, 18, 14, 9, 4, 5, 8] + [10, 0]),
    )
    assert_array_almost_equal(
        g.at_node["drainage_area"],
        np.array(
            [0.0, 0.0, 0.0, 0.0, 86.60254]
            + [86.60254, 0.0, 0.0, 86.60254, 346.41016]
            + [86.60254, 0.0, 0.0, 86.60254, 519.61524]
            + [0.0, 86.60254, 0.0, 519.61524]
        ),
    )
    assert_array_almost_equal(
        g.at_node["surface_water__discharge"],
        np.array(
            [0.0, 0.0, 0.0, 0.0, 173.20508]
            + [173.20508, 0.0, 0.0, 173.20508, 692.82032]
            + [173.20508, 0.0, 0.0, 173.20508, 1039.23048]
            + [0.0, 173.20508, 0.0, 1039.23048]
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
    self = FlowRouter(g, runoff_rate="water__unit_flux_in")
    g.at_node["topographic__elevation"] = [0.0, 0.2, 0.25, 0.15, 0.25, 0.4, 0.8, 0.8]
    g.at_node["water__unit_flux_in"] = [1.0, 2.0, 3.0, 4.0, 3.0, 3.2, 4.5, 1.0]

    g.status_at_node[7] = g.BC_NODE_IS_CLOSED
    self.run_flow_directions()
    self.run_flow_accumulations()

    assert_array_equal(
        g.at_node["flow__upstream_node_order"], np.array([0, 1, 2, 3, 5, 6, 7, 4])
    )
    assert_array_almost_equal(
        g.at_node["drainage_area"], np.array([8.0, 7.0, 1.0, 4.0, 1.0, 3.0, 1.0, 1.0])
    )
    assert_array_almost_equal(
        g.at_node["surface_water__discharge"],
        np.array([8.0, 14.0, 3.0, 16.0, 3.0, 9.6, 4.5, 1.0]),
    )
