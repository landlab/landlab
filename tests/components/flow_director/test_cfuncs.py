import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
import landlab.components.flow_director.cfuncs as _cfuncs


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
