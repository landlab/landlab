import numpy as np
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

import landlab.components.flow_director.cfuncs as _cfuncs
from landlab import HexModelGrid
from landlab import RasterModelGrid


def test_hex():
    grid = HexModelGrid((3, 3))
    z = grid.x_of_node

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

    assert_array_equal(receiver, [3, 0, 1, 3, 3, 4, 5, 3, 7, 8])
    assert_array_equal(
        receiver_link,
        [2, 0, 1, -1, 8, 9, 10, 11, 17, 18],
    )
    assert_array_almost_equal(steepest_slope, [0.5, 1, 1, 0, 1, 1, 1, 0.5, 1, 1])


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
    # Grid of 5 nodes with these coordinates:
    # x = [0., 0., 3., 2., 3.]; y = [0., 4., 0., 2., 4.]
    # All nodes are linked (7 links)
    nodes_n = 5
    z = np.array([0.0, 1.0, 0.0, 4.0, 5.0])  # elevation
    receiver = -1 * np.ones(nodes_n, dtype=int)
    receiver_link = -1 * np.ones(nodes_n, dtype=int)
    steepest_slope = np.zeros(nodes_n, dtype=float)
    src_nodes = np.array([0, 0, 1, 0, 1, 3, 2, 3])  # tails
    dst_nodes = np.array([1, 2, 2, 3, 4, 2, 4, 4])  # heads
    link_slope = np.array(
        [-0.33333333, -0.0, 0.4472136, -1.0, -1.0]
        + [1.41421356, -2.23606798, -0.33333333]
    )

    _cfuncs.adjust_flow_receivers(
        src_nodes,
        dst_nodes,
        z,
        link_slope,
        np.arange(7),
        receiver,
        receiver_link,
        steepest_slope,
    )
    assert_array_equal(src_nodes, np.array([0, 0, 1, 0, 1, 3, 2, 3]))
    assert_array_equal(dst_nodes, np.array([1, 2, 2, 3, 4, 2, 4, 4]))
    assert_array_equal(receiver, np.array([-1, 2, -1, 2, 2]))
    assert_array_equal(receiver_link, np.array([-1, 2, -1, 5, 6]))
    assert_array_almost_equal(
        steepest_slope, np.array([0.0, 0.4472136, 0.0, 1.41421356, 2.23606798])
    )
