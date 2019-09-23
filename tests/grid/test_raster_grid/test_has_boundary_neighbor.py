import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid


def test_boundary_node():
    rmg = RasterModelGrid((5, 6))
    assert rmg.node_has_boundary_neighbor(0)
    assert not rmg.node_has_boundary_neighbor(14)


def test_last_index():
    rmg = RasterModelGrid((4, 5))
    assert rmg.node_has_boundary_neighbor(-1)


def test_id_as_list():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(rmg.node_has_boundary_neighbor([-1, 0]), np.array([True, True]))


def test_id_as_array():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(
        rmg.node_has_boundary_neighbor(np.arange(20)),
        np.array(
            [
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
            ]
        ),
    )


def test_id_as_array_with_one_interior():
    rmg = RasterModelGrid((5, 5))
    assert_array_equal(
        rmg.node_has_boundary_neighbor(np.arange(25)),
        np.array(
            [
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                False,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
            ]
        ),
    )
