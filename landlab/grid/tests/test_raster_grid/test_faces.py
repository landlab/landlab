import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid


def test_horizontally_adjacent_cells():
    grid = RasterModelGrid((4, 5))
    assert_array_equal(grid.face_connecting_cell_pair(0, 1), np.array([4]))


def test_vertically_adjacent_cells():
    grid = RasterModelGrid((4, 5))
    assert_array_equal(grid.face_connecting_cell_pair(0, 3), np.array([7]))


def test_diagonally_adjacent_cells():
    grid = RasterModelGrid((4, 5))
    assert_array_equal(grid.face_connecting_cell_pair(1, 5), np.array([]))


def test_non_adjacent_cells():
    grid = RasterModelGrid((4, 5))
    assert_array_equal(grid.face_connecting_cell_pair(0, 2), np.array([]))


def test_id_as_int():
    grid = RasterModelGrid((4, 5))
    assert_array_equal(grid.faces_at_cell[0], np.array([4, 7, 3, 0]))


def test_width_of_face():
    grid = RasterModelGrid((3, 3))
    assert_array_equal(grid.width_of_face, [1., 1., 1., 1.])

    grid = RasterModelGrid((3, 3), spacing=2.)
    assert_array_equal(grid.width_of_face, [2., 2., 2., 2.])

    grid = RasterModelGrid((3, 3), spacing=(2., 3.))
    assert_array_equal(grid.width_of_face, [3., 2., 2., 3.])


def test_id_as_array():
    grid = RasterModelGrid((4, 5))
    assert_array_equal(
        grid.faces_at_cell[[0, 1]], np.array([[4, 7, 3, 0], [5, 8, 4, 1]])
    )
