from numpy.testing import assert_array_equal

from landlab import HexModelGrid


def test_patches_at_node():
    grid = HexModelGrid((3, 3))
    assert_array_equal(
        grid.patches_at_node,
        [
            [0, 2, -1, -1, -1, -1],
            [1, 3, 0, -1, -1, -1],
            [4, 1, -1, -1, -1, -1],
            [5, 2, -1, -1, -1, -1],
            [6, 8, 5, 2, 0, 3],
            [7, 9, 6, 3, 1, 4],
            [7, 4, -1, -1, -1, -1],
            [5, 8, -1, -1, -1, -1],
            [8, 6, 9, -1, -1, -1],
            [9, 7, -1, -1, -1, -1],
        ],
    )
