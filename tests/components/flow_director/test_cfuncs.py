import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from landlab import RasterModelGrid, HexModelGrid
import landlab.components.flow_director.cfuncs as _cfuncs


def test_raster_plane_dipping_north():
    grid = RasterModelGrid((3, 4))
    z = np.array(
        [
            [2.0, 2.0, 2.0, 2.0],
            [1.0, 1.0, 1.0, 1.0],
            [0.0, 0.0, 0.0, 0.0],
        ]
    ).flatten()

    steepest_slope = np.zeros(grid.number_of_nodes)
    receiver = np.arange(grid.number_of_nodes)
    receiver_link = np.full(grid.number_of_nodes, -1, dtype=int)
    active_links = np.arange(grid.number_of_links, dtype=int)

    _cfuncs.adjust_flow_receivers(
        grid.node_at_link_head,
        grid.node_at_link_tail,
        z,
        grid.calc_grad_at_link(z),
        active_links,
        receiver,
        receiver_link,
        steepest_slope,
    )

    assert_array_equal(
        receiver.reshape(grid.shape), [[4, 5, 6, 7], [8, 9, 10, 11], [8, 9, 10, 11]]
    )
    assert_array_equal(
        receiver_link.reshape(grid.shape),
        [[3, 4, 5, 6], [10, 11, 12, 13], [-1, -1, -1, -1]],
    )
    assert_array_almost_equal(
        steepest_slope.reshape(grid.shape), [[1, 1, 1, 1], [1, 1, 1, 1], [0, 0, 0, 0]]
    )


def test_raster_plane_dipping_south():
    grid = RasterModelGrid((3, 4))
    z = np.array(
        [
            [0.0, 0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0, 1.0],
            [2.0, 2.0, 2.0, 2.0],
        ]
    ).flatten()

    steepest_slope = np.zeros(grid.number_of_nodes)
    receiver = np.arange(grid.number_of_nodes)
    receiver_link = np.full(grid.number_of_nodes, -1, dtype=int)
    active_links = np.arange(grid.number_of_links, dtype=int)

    _cfuncs.adjust_flow_receivers(
        grid.node_at_link_head,
        grid.node_at_link_tail,
        z,
        grid.calc_grad_at_link(z),
        active_links,
        receiver,
        receiver_link,
        steepest_slope,
    )

    assert_array_equal(
        receiver.reshape(grid.shape), [[0, 1, 2, 3], [0, 1, 2, 3], [4, 5, 6, 7]]
    )
    assert_array_equal(
        receiver_link.reshape(grid.shape),
        [[-1, -1, -1, -1], [3, 4, 5, 6], [10, 11, 12, 13]],
    )
    assert_array_almost_equal(
        steepest_slope.reshape(grid.shape), [[0, 0, 0, 0], [1, 1, 1, 1], [1, 1, 1, 1]]
    )


def test_adjust_flow_receivers():
    params = {"shape": (7, 4), "spacing": 10, "node_layout": "hex"}
    g = HexModelGrid(**params)
    nodes_n = g.number_of_nodes
    random_generator = np.random.Generator(np.random.PCG64(seed=500))
    z = g.add_field("topographic__elevation", 10 * random_generator.random(nodes_n))
    receiver = -1 * np.ones(nodes_n, dtype=int)
    receiver_link = -1 * np.ones(nodes_n, dtype=int)
    steepest_slope = np.zeros(nodes_n, dtype=float)
    src_nodes = g.node_at_link_tail[g.active_links]
    dst_nodes = g.node_at_link_head[g.active_links]
    link_slope = -g.calc_grad_at_link(z)[g.active_links]

    _cfuncs.adjust_flow_receivers(
        src_nodes,
        dst_nodes,
        z,
        link_slope,
        g.active_links,
        receiver,
        receiver_link,
        steepest_slope,
    )
    assert_array_equal(src_nodes, g.node_at_link_tail[g.active_links])
    assert_array_equal(dst_nodes, g.node_at_link_head[g.active_links])
    assert_array_equal(
        receiver,
        np.array(
            [-1, 6, 6, -1, 10, 6, -1, 6, -1, -1, 9, 18, 18, 14, -1, 16, 9]
            + [18, -1, 18, 27, 20, 16, 24, 18, 18, 25, -1, 23, 24, 24, 36]
            + [26, -1, 29, 31, -1]
        ),
    )
    assert_array_equal(
        receiver_link,
        np.array(
            [-1, 6, 7, -1, 16, 12, -1, 13, -1, -1, 25, 35, 36, 29, -1, 42]
            + [31, 44, -1, 45, 58, 47, 49, 61, 53, 54, 63, -1, 66, 68, 69]
            + [85, 73, -1, 81, 84, -1]
        ),
    )
    assert_array_almost_equal(
        steepest_slope,
        np.array(
            [0.0, 0.54145482, 0.33283403, 0.0, 0.02588641, 0.50883805, 0.0]
            + [0.47369762, 0.0, 0.0, 0.19451533, 0.77608988, 0.4361625, 0.25266218]
            + [0.0, 0.61566294, 0.01156903, 0.68505865, 0.0, 0.31045881, 0.07843021]
            + [0.1662848, 0.10984337, 0.2785543, 0.07273778, 0.12874949, 0.24497883]
            + [0.0, 0.0156391, 0.44582296, 0.74436412, 0.78271885, 0.14129527, 0.0]
            + [0.30206114, 0.05719591, 0.0]
        ),
    )
