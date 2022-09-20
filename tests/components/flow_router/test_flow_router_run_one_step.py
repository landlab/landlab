import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from landlab.components import FlowRouter

from landlab import RasterModelGrid, HexModelGrid
from landlab import NetworkModelGrid


def test_run_one_step_raster():
    spacing = 10
    g = RasterModelGrid((5, 5), (spacing, spacing))
    g.status_at_node[g.perimeter_nodes] = g.BC_NODE_IS_CLOSED
    nodes_n = g.number_of_nodes
    self = FlowRouter(g)
    random_generator = np.random.Generator(np.random.PCG64(seed=500))
    g.at_node["topographic__elevation"] = 10 * random_generator.random(nodes_n)
    self.run_one_step()

    assert_array_equal(
        g.at_node["flow__receiver_node"],
        np.array(
            [0, 1, 2, 3, 4, 5, 12, 8, 4, 9, 10, 12, 8, 8, 14, 15, 12]
            + [13, 13, 19, 20, 21, 22, 23, 24]
        ),
    )
    assert_array_almost_equal(
        g.at_node["topographic__steepest_slope"],
        np.array(
            [0.0, 0.0, 0.0, 0.0, 0.0]
            + [0.0, 0.0, 0.35067063, 0.0, 0.0]
            + [0.0, 0.33992738, 0.08589512, 0.09309407, 0.0]
            + [0.0, 0.0, 0.19606383, 0.0, 0.0]
            + [0.0, 0.0, 0.0, 0.0, 0.0],
        ),
    )
    assert_array_equal(
        g.at_node["flow__link_to_receiver_node"],
        np.array(
            [-1, -1, -1, -1, -1, -1, 50, 11, 47, -1, -1, 19, 53, 16, -1, -1, 59]
            + [61, 25, -1, -1, -1, -1, -1, -1],
        ),
    )
    assert_array_equal(
        g.at_node["flood_status_code"],
        np.array(
            [0, 0, 0, 0, 0, 0, 3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 3, 0, 3, 0, 0, 0, 0]
            + [0, 0],
        ),
    )
    assert_array_almost_equal(
        g.at_node["depression__depth"],
        np.array(
            [0.0, 0.0, 0.0, 0.0, 0.0]
            + [0.0, 2.44501033, 0.0, 0.32481788, 0.0]
            + [0.0, 0.0, 0.0, 0.0, 0.0]
            + [0.0, 2.97824959, 0.0, 4.07782519, 0.0]
            + [0.0, 0.0, 0.0, 0.0, 0.0],
        ),
    )
    assert_array_equal(
        g.at_node["outlet_node"],
        np.array(
            [0, 1, 2, 3, 4, 5, 4, 4, 4, 9, 10, 4, 4, 4, 14, 15, 4]
            + [4, 4, 19, 20, 21, 22, 23, 24],
        ),
    )
    assert_array_equal(
        g.at_node["depression__outlet_node"],
        np.array(
            [-1, -1, -1, -1, -1, -1, 12, -1, 4, -1, -1, -1, -1, -1, -1, -1, 12]
            + [-1, 13, -1, -1, -1, -1, -1, -1],
        ),
    )
    assert_array_almost_equal(
        g.at_node["depression_free__elevation"],
        np.array(
            [5.66743143, 8.53977988, 6.45357199, 4.11156813, 4.68031945]
            + [8.21361221, 5.57024205, 7.8622079, 4.68031945, 2.47630211]
            + [4.42145537, 8.96951584, 5.57024205, 5.28644224, 2.75982042]
            + [8.74862188, 5.57024205, 8.05920356, 5.28644224, 4.31320515]
            + [3.49732191, 5.16016994, 3.69042619, 4.72153783, 1.93599482],
        ),
    )

    assert_array_equal(
        g.at_node["flow__upstream_node_order"],
        np.array(
            [4, 8, 13, 17, 18, 12, 11, 6, 16, 7, 0, 1, 2, 3, 5, 9, 10]
            + [14, 15, 19, 20, 21, 22, 23, 24],
        ),
    )
    assert_array_almost_equal(
        g.at_node["drainage_area"],
        np.array(
            [0.0, 0.0, 0.0, 0.0, 900.0, 0.0, 100.0, 100.0, 900.0]
            + [0.0, 0.0, 100.0, 400.0, 300.0, 0.0, 0.0, 100.0, 100.0]
            + [100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ),
    )
    assert_array_almost_equal(
        g.at_node["surface_water__discharge"],
        np.array(
            [0.0, 0.0, 0.0, 0.0, 900.0, 0.0, 100.0, 100.0, 900.0]
            + [0.0, 0.0, 100.0, 400.0, 300.0, 0.0, 0.0, 100.0, 100.0]
            + [100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ),
    )


def test_run_one_step_hex():
    spacing = 10

    g = HexModelGrid((7, 4), spacing, node_layout="hex")
    g.status_at_node[g.perimeter_nodes] = g.BC_NODE_IS_FIXED_VALUE
    g.status_at_node[0] = g.BC_NODE_IS_CLOSED
    nodes_n = g.number_of_nodes

    self = FlowRouter(g, surface="soil__elevation", diagonals=True, runoff_rate=2.0)
    random_generator = np.random.Generator(np.random.PCG64(seed=500))
    g.at_node["soil__elevation"] = 10 * random_generator.random(nodes_n)
    self.run_one_step()

    assert_array_equal(
        g.at_node["flow__receiver_node"],
        np.array(
            [0, 1, 2, 3, 4, 10, 12, 3, 8, 9, 9, 18, 19, 14, 14, 15, 9]
            + [16, 19, 20, 27, 21, 22, 16, 18, 19, 27, 27, 28, 24, 24, 36, 32, 33]
            + [34, 35, 36],
        ),
    )
    assert_array_almost_equal(
        g.at_node["topographic__steepest_slope"],
        np.array(
            [0.0, 0.0, 0.0, 0.0, 0.0]
            + [0.37921568, 0.0, 0.37506398, 0.0, 0.0]
            + [0.19451533, 0.77608988, 0.12570369, 0.25266218, 0.0]
            + [0.0, 0.01156903, 0.54672111, 0.0, 0.08158832]
            + [0.07843021, 0.0, 0.0, 0.21295454, 0.07273778]
            + [0.0, 0.22328804, 0.0, 0.0, 0.44582296]
            + [0.74436412, 0.78271885, 0.0, 0.0, 0.0]
            + [0.0, 0.0],
        ),
    )
    assert_array_equal(
        g.at_node["flow__link_to_receiver_node"],
        np.array(
            [-1, -1, -1, -1, -1, 17, 20, 9, -1, -1, 25, 35, 37, 29, -1, -1, 31]
            + [43, 45, 46, 58, -1, -1, 50, 53, 55, 64, -1, -1, 68, 69, 85, -1, -1]
            + [-1, -1, -1],
        ),
    )
    assert_array_equal(
        g.at_node["flood_status_code"],
        np.array(
            [0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0]
            + [0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ),
    )
    assert_array_almost_equal(
        g.at_node["depression__depth"],
        np.array(
            [0.0, 0.0, 0.0, 0.0, 0.0]
            + [0.0, 2.44501033, 0.0, 0.0, 0.0]
            + [0.0, 0.0, 0.0, 0.0, 0.0]
            + [0.0, 0.0, 0.0, 3.1045881, 0.0]
            + [0.0, 0.0, 0.0, 0.0, 2.37721033]
            + [1.8170932, 0.0, 0.0, 0.0, 0.0]
            + [0.0, 0.0, 0.0, 0.0, 0.0]
            + [0.0, 0.0],
        ),
    )
    assert_array_equal(
        g.at_node["outlet_node"],
        np.array(
            [0, 1, 2, 3, 4, 9, 27, 3, 8, 9, 9, 27, 27, 14, 14, 15, 9]
            + [9, 27, 27, 27, 21, 22, 9, 27, 27, 27, 27, 28, 27, 27, 36, 32, 33]
            + [34, 35, 36],
        ),
    )
    assert_array_equal(
        g.at_node["depression__outlet_node"],
        np.array(
            [-1, -1, -1, -1, -1, -1, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
            + [-1, 19, -1, -1, -1, -1, -1, 19, 19, -1, -1, -1, -1, -1, -1, -1, -1]
            + [-1, -1, -1],
        ),
    )
    assert_array_almost_equal(
        g.at_node["depression_free__elevation"],
        np.array(
            [5.66743143, 8.53977988, 6.45357199, 4.11156813, 4.68031945]
            + [8.21361221, 5.57024205, 7.8622079, 4.35550157, 2.47630211]
            + [4.42145537, 8.96951584, 5.57024205, 5.28644224, 2.75982042]
            + [8.74862188, 2.59199246, 8.05920356, 4.31320515, 4.31320515]
            + [3.49732191, 5.16016994, 3.69042619, 4.72153783, 4.31320515]
            + [4.31320515, 4.9459002, 2.71301982, 4.8779288, 6.39422443]
            + [9.37963598, 8.29564356, 6.35885287, 6.11800132, 9.41483585]
            + [8.8676027, 0.46845509],
        ),
    )

    assert_array_equal(
        g.at_node["flow__upstream_node_order"],
        np.array(
            [1, 2, 3, 7, 4, 8, 9, 16, 17, 23, 10, 5, 14, 13, 15, 21, 22]
            + [27, 26, 20, 19, 25, 18, 11, 24, 29, 30, 12, 6, 28, 32, 33, 34, 35]
            + [36, 31, 0],
        ),
    )
    assert_array_almost_equal(
        g.at_node["drainage_area"],
        np.array(
            [0.0, 0.0, 0.0, 86.60254, 0.0]
            + [86.60254, 86.60254, 86.60254, 0.0, 433.0127]
            + [173.20508, 86.60254, 173.20508, 86.60254, 86.60254]
            + [0.0, 259.80762, 86.60254, 433.01271, 779.42287]
            + [866.02541, 0.0, 0.0, 86.60254, 259.80763]
            + [86.60254, 86.60254, 952.62795, 0.0, 86.602545]
            + [86.602545, 86.602545, 0.0, 0.0, 0.0]
            + [0.0, 86.602545],
        ),
    )
    assert_array_almost_equal(
        g.at_node["surface_water__discharge"],
        np.array(
            [0.0, 0.0, 0.0, 173.20508, 0.0]
            + [173.20508, 173.20508, 173.20508, 0.0, 866.0254]
            + [346.41016, 173.20508, 346.41016, 173.20508, 173.20508]
            + [0.0, 519.61524, 173.20508, 866.02542, 1558.84574]
            + [1732.05082, 0.0, 0.0, 173.20508, 519.61526]
            + [173.20508, 173.20508, 1905.2559, 0.0, 173.20509]
            + [173.20509, 173.20509, 0.0, 0.0, 0.0]
            + [0.0, 173.20509],
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
    self = FlowRouter(g, runoff_rate="water__unit_flux_in")
    g.at_node["topographic__elevation"] = [0.0, 0.2, 0.25, 0.15, 0.25, 0.4, 0.8, 0.8]
    g.at_node["water__unit_flux_in"] = [1.0, 2.0, 3.0, 4.0, 3.0, 3.2, 4.5, 1.0]

    g.status_at_node[7] = g.BC_NODE_IS_CLOSED
    self.run_one_step()

    assert_array_equal(
        g.at_node["flow__receiver_node"], np.array([0, 0, 1, 1, 1, 3, 5, 5])
    )
    assert_array_almost_equal(
        g.at_node["topographic__steepest_slope"],
        np.array(
            [0.0, 0.002, 0.00048507, 0.0, 0.00035355]
            + [0.00223607, 0.00357771, 0.0022188],
        ),
    )
    assert_array_equal(
        g.at_node["flow__link_to_receiver_node"], np.array([-1, 0, 1, 2, 3, 4, 5, 6])
    )
    assert_array_equal(
        g.at_node["flood_status_code"], np.array([0, 0, 0, 3, 0, 0, 0, 0])
    )
    assert_array_almost_equal(
        g.at_node["depression__depth"],
        np.array([0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0]),
    )
    assert_array_equal(g.at_node["outlet_node"], np.array([0, 0, 0, 0, 0, 0, 0, 0]))
    assert_array_equal(
        g.at_node["depression__outlet_node"], np.array([-1, -1, -1, 1, -1, -1, -1, -1])
    )
    assert_array_almost_equal(
        g.at_node["depression_free__elevation"],
        np.array([0.0, 0.2, 0.25, 0.2, 0.25, 0.4, 0.8, 0.8]),
    )

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
