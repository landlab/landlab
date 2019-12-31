import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid


def test_id_as_int():
    grid = RasterModelGrid((4, 5))
    assert_array_equal(grid.faces_at_cell[0], np.array([4, 7, 3, 0]))


def test_length_of_face():
    grid = RasterModelGrid((3, 3))
    assert_array_equal(grid.length_of_face, [1.0, 1.0, 1.0, 1.0])

    grid = RasterModelGrid((3, 3), xy_spacing=2.0)
    assert_array_equal(grid.length_of_face, [2.0, 2.0, 2.0, 2.0])

    grid = RasterModelGrid((3, 3), xy_spacing=(3.0, 2.0))
    assert_array_equal(grid.length_of_face, [3.0, 2.0, 2.0, 3.0])


def test_id_as_array():
    grid = RasterModelGrid((4, 5))
    assert_array_equal(
        grid.faces_at_cell[[0, 1]], np.array([[4, 7, 3, 0], [5, 8, 4, 1]])
    )
