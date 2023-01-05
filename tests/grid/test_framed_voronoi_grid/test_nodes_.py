import numpy as np
from numpy.testing import assert_array_equal

from landlab import FramedVoronoiGrid


def test_rect_patches_at_node():
    grid = FramedVoronoiGrid((3, 3))
    assert_array_equal(
        grid.patches_at_node,
        np.array(
            [
                [0, -1, -1, -1, -1, -1],
                [1, 2, 0, -1, -1, -1],
                [3, 1, -1, -1, -1, -1],
                [4, 0, 2, -1, -1, -1],
                [5, 6, 4, 2, 1, 3],
                [7, 5, 3, -1, -1, -1],
                [4, 6, -1, -1, -1, -1],
                [6, 5, 7, -1, -1, -1],
                [7, -1, -1, -1, -1, -1],
            ]
        ),
    )
