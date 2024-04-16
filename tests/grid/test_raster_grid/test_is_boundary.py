import numpy as np
from numpy.testing import assert_array_equal

from landlab import NodeStatus
from landlab import RasterModelGrid


def test_id_as_int():
    rmg = RasterModelGrid((4, 5))
    assert rmg.node_is_boundary(0)


def test_id_as_small_list():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(rmg.node_is_boundary([0]), np.array([True]))


def test_id_as_array():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(
        rmg.node_is_boundary(np.arange(20)),
        np.array(
            [
                True,
                True,
                True,
                True,
                True,
                True,
                False,
                False,
                False,
                True,
                True,
                False,
                False,
                False,
                True,
                True,
                True,
                True,
                True,
                True,
            ],
            dtype=bool,
        ),
    )


def test_id_as_list():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(rmg.node_is_boundary([8, 9]), np.array([False, True]))


def test_boundary_flag():
    rmg = RasterModelGrid((4, 5))
    rmg.status_at_node[0] = NodeStatus.FIXED_GRADIENT
    assert_array_equal(
        rmg.node_is_boundary(np.arange(20)),
        np.array(
            [
                True,
                True,
                True,
                True,
                True,
                True,
                False,
                False,
                False,
                True,
                True,
                False,
                False,
                False,
                True,
                True,
                True,
                True,
                True,
                True,
            ],
            dtype=bool,
        ),
    )

    assert_array_equal(
        rmg.node_is_boundary(np.arange(20), boundary_flag=NodeStatus.FIXED_GRADIENT),
        np.array(
            [
                True,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
            ],
            dtype=bool,
        ),
    )
