import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid


def test_inactive_boundaries():
    rmg = RasterModelGrid(3, 4)
    assert_array_equal(rmg.active_links, np.array([4, 5, 7, 8, 9, 11, 12]))
    with pytest.deprecated_call():
        rmg.set_inactive_boundaries(True, True, True, True)
    assert_array_equal(rmg.active_links, np.array([8]))
    assert_array_equal(
        rmg._active_links_at_node(),
        np.array(
            [
                [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
                [-1, -1, -1, -1, -1, -1, 8, -1, -1, -1, -1, -1],
                [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
                [-1, -1, -1, -1, -1, 8, -1, -1, -1, -1, -1, -1],
            ]
        ),
    )


def test_inactive_interiors():
    rmg = RasterModelGrid((4, 5))
    rmg.set_closed_nodes([6, 12])
    assert_array_equal(
        rmg._active_links_at_node(),
        np.array(
            [
                [
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    6,
                    7,
                    -1,
                    -1,
                    -1,
                    -1,
                    16,
                    -1,
                    -1,
                    23,
                    -1,
                    25,
                    -1,
                ],
                [
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    11,
                    12,
                    -1,
                    18,
                    -1,
                    -1,
                    21,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                ],
                [
                    -1,
                    -1,
                    6,
                    7,
                    -1,
                    -1,
                    -1,
                    -1,
                    16,
                    -1,
                    -1,
                    23,
                    -1,
                    25,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                ],
                [
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    11,
                    12,
                    -1,
                    18,
                    -1,
                    -1,
                    21,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                ],
            ]
        ),
    )
