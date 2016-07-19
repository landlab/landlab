import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import (assert_equal, assert_raises, raises, assert_true,
                        assert_false)
try:
    from nose.tools import assert_is
except ImportError:
    from landlab.testing.tools import assert_is

from landlab import RasterModelGrid


class TestRasterModelGridConnectingFaces():

    def setup(self):
        self.rmg = RasterModelGrid(4, 5)

    def test_horizontally_adjacent_cells(self):
        assert_array_equal(self.rmg.face_connecting_cell_pair(0, 1),
                           np.array([4]))

    def test_vertically_adjacent_cells(self):
        assert_array_equal(self.rmg.face_connecting_cell_pair(0, 3),
                           np.array([7]))

    def test_diagonally_adjacent_cells(self):
        assert_array_equal(self.rmg.face_connecting_cell_pair(1, 5),
                           np.array([]))

    def test_non_adjacent_cells(self):
        assert_array_equal(self.rmg.face_connecting_cell_pair(0, 2),
                           np.array([]))


class TestRasterModelGridCellFaces():

    def setup(self):
        self.rmg = RasterModelGrid(4, 5)

    def test_id_as_int(self):
        assert_array_equal(self.rmg.faces_at_cell[0], np.array([4, 7, 3, 0]))

    def test_id_as_array(self):
        assert_array_equal(self.rmg.faces_at_cell[[0, 1]],
                           np.array([[4, 7, 3, 0], [5, 8, 4, 1]]))
