#! /usr/bin/env python
"""
Unit tests for landlab.model_grid
"""

import unittest
import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.model_grid import BAD_INDEX_VALUE


_SPACING = 1.
(_NUM_ROWS, _NUM_COLS) = (4, 5)

class TestRasterModelGrid(unittest.TestCase):
    def setUp(self):
        """
        These tests use a grid that 4x5 nodes::

         15------16------17------18------19
          |       |       |       |       |
          |       |       |       |       |
          |       |       |       |       |
         10------11------12------13------14
          |       |       |       |       |
          |       |       |       |       |   
          |       |       |       |       |
          5-------6-------7-------8-------9
          |       |       |       |       |
          |       |       |       |       |
          |       |       |       |       |
          0-------1-------2-------3-------4
        """
        self.grid = RasterModelGrid(num_rows=_NUM_ROWS, num_cols=_NUM_COLS,
                                    dx=_SPACING)

    def test_create_grid(self):
        """
        Create a raster grid.
        """
        grid = RasterModelGrid(num_rows=4, num_cols=5)

    def test_grid_dimensions(self):
        """
        Check the width of the grid. The x-dimension is the columns, and the
        y-dimension is the rows.
        """
        SPACING = 2.
        (NUM_ROWS, NUM_COLS) = (4, 5)

        grid = RasterModelGrid(num_rows=NUM_ROWS, num_cols=NUM_COLS,
                               dx=SPACING)

        self.assertEqual(grid.get_grid_ydimension(), NUM_ROWS * SPACING)
        self.assertEqual(grid.get_grid_xdimension(), NUM_COLS * SPACING)

    def test_grid_dimensions_default_spacing(self):
        """
        Use the default spacing of 1.
        """
        (NUM_ROWS, NUM_COLS) = (4, 5)

        grid = RasterModelGrid(num_rows=NUM_ROWS, num_cols=NUM_COLS)

        self.assertEqual(grid.get_grid_ydimension(), NUM_ROWS * 1.)
        self.assertEqual(grid.get_grid_xdimension(), NUM_COLS * 1.)

    def test_nodes_around_point(self):
        surrounding_ids = self.grid.get_nodes_around_point(2.1, 1.1)
        self.assertListEqual(
            list(surrounding_ids), [7, 8, 12, 13],
            'incorrect surrounding ids in test_nodes_around_point')

        surrounding_ids = self.grid.get_nodes_around_point(2.1, .9)
        self.assertListEqual(
            list(surrounding_ids), [2, 3, 7, 8],
            'incorrect surrounding ids in test_nodes_around_point')

    def test_neighbor_list(self):
        neighbors = self.grid.get_neighbor_list(id=6)
        self.assertEqual(list(neighbors), [7, 11, 5, 1],
                         'incorrect neighbors in test_neighbor_list')

    def test_neighbor_list_boundary(self):
        """
        All of the neighbor IDs for a boundary cell are -1.
        """
        neighbors = self.grid.get_neighbor_list(id=0)

        self.assertEqual(list(neighbors), [-1, -1, -1, -1])

    def test_cell_x(self):
        expected_x = [0., 1., 2., 3., 4.,
                      0., 1., 2., 3., 4.,
                      0., 1., 2., 3., 4.,
                      0., 1., 2., 3., 4.]

        for (cell_id, expected) in zip(xrange(12), expected_x):
            cell_x = self.grid.get_node_x(cell_id)
            self.assertEqual(cell_x, expected)

            cell_x = self.grid.node_x[cell_id]
            self.assertEqual(cell_x, expected)

    def test_cell_y(self):
        expected_y = [0., 0., 0., 0., 0.,
                      1., 1., 1., 1., 1.,
                      2., 2., 2., 2., 2.,
                      3., 3., 3., 3., 3.]

        for (cell_id, expected) in zip(xrange(12), expected_y):
            cell_y = self.grid.get_node_y(cell_id)
            self.assertEqual(cell_y, expected)

            cell_y = self.grid.node_y[cell_id]
            self.assertEqual(cell_y, expected)

    def test_node_x_coordinates(self):
        x_coords = self.grid.get_node_x_coords()

        self.assertEqual(list(x_coords), [0., 1., 2., 3., 4.,
                                          0., 1., 2., 3., 4.,
                                          0., 1., 2., 3., 4.,
                                          0., 1., 2., 3., 4.])

    def test_cell_y_coordinates(self):
        y_coords = self.grid.get_node_y_coords()

        self.assertEqual(list(y_coords), [0., 0., 0., 0., 0.,
                                          1., 1., 1., 1., 1.,
                                          2., 2., 2., 2., 2.,
                                          3., 3., 3., 3., 3.])

    def test_diagonal_list(self):
        diagonals = self.grid.get_diagonal_list(id=6)
        self.assertEqual(list(diagonals), [12, 10, 0, 2])

    def test_diagonal_list_boundary(self):
        diagonals = self.grid.get_diagonal_list(id=0)
        self.assertEqual(list(diagonals), [-1, -1, -1, -1])

    def test_is_interior(self):
        for cell_id in [0, 1, 2, 3, 4, 5, 9, 10, 14, 15, 16, 17, 18, 19]:
            self.assertFalse(self.grid.is_interior(cell_id))

        for cell_id in [6, 7, 8, 11, 12, 13]:
            self.assertTrue(self.grid.is_interior(cell_id))

    def test_get_interior_cells(self):
        interiors = self.grid.get_active_cell_node_ids()
        self.assertEqual(list(interiors), [6, 7, 8, 11, 12, 13])

    def test_active_links(self):
        self.assertTrue(isinstance(self.grid.active_links, np.ndarray))
        self.assertEqual(self.grid.num_active_links, 17)

        assert_array_equal(self.grid.active_links,
                           np.array([1, 2, 3, 6, 7, 8, 11, 12, 13,
                                     19, 20, 21, 22, 23, 24, 25, 26]))

    def test_active_link_fromnode(self):
        assert_array_equal(self.grid.activelink_fromnode,
                           np.array([1, 2, 3, 6, 7, 8, 11, 12, 13,
                                     5, 6, 7, 8, 10, 11, 12, 13]))

    def test_active_link_tonode(self):
        assert_array_equal(self.grid.activelink_tonode,
                           np.array([6, 7, 8, 11, 12, 13, 16, 17, 18,
                                     6, 7, 8, 9, 11, 12, 13, 14]))

    def test_active_link_num_inlink(self):
        assert_array_equal(self.grid.node_numactiveinlink,
                           np.array([0, 0, 0, 0, 0,
                                     0, 2, 2, 2, 1,
                                     0, 2, 2, 2, 1,
                                     0, 1, 1, 1, 0]))

    def test_active_link_num_outlink(self):
        assert_array_equal(self.grid.node_numactiveoutlink,
                           np.array([0, 1, 1, 1, 0,
                                     1, 2, 2, 2, 0,
                                     1, 2, 2, 2, 0,
                                     0, 0, 0, 0, 0]))

    def test_active_inlink_matrix(self):
        self.assertTupleEqual(self.grid.node_active_inlink_matrix.shape,
                              (2, 20))

        assert_array_equal(self.grid.node_active_inlink_matrix[0],
                           np.array([-1, -1, -1, -1, -1,
                                     -1,  0,  1,  2, 12,
                                     -1,  3,  4,  5, 16,
                                     -1,  6,  7,  8, -1]))
        assert_array_equal(self.grid.node_active_inlink_matrix[1],
                           np.array([-1, -1, -1, -1, -1,
                                     -1,  9, 10, 11, -1,
                                     -1, 13, 14, 15, -1,
                                     -1, -1, -1, -1, -1]))

    def test_active_outlink_matrix(self):
        self.assertTupleEqual(self.grid.node_active_outlink_matrix.shape,
                              (2, 20))

        assert_array_equal(self.grid.node_active_outlink_matrix[0],
                           np.array([-1,  0,  1,  2, -1,
                                      9,  3,  4,  5, -1,
                                     13,  6,  7,  8, -1,
                                     -1, -1, -1, -1, -1]))
        assert_array_equal(self.grid.node_active_outlink_matrix[1],
                           np.array([-1, -1, -1, -1, -1,
                                     -1, 10, 11, 12, -1,
                                     -1, 14, 15, 16, -1,
                                     -1, -1, -1, -1, -1]))

    def test_link_fromnode(self):
        assert_array_equal(self.grid.link_fromnode,
                           np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                     10, 11, 12, 13, 14,  0,  1,  2,  3,  5,
                                      6,  7,  8, 10, 11, 12, 13, 15, 16, 17,
                                     18]))

    def test_link_tonode(self):
        assert_array_equal(self.grid.link_tonode,
                           np.array([ 5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
                                     15, 16, 17, 18, 19,  1,  2,  3,  4,  6,
                                      7,  8,  9, 11, 12, 13, 14, 16, 17, 18,
                                     19]))

    def test_link_num_inlink(self):
        assert_array_equal(self.grid.node_numinlink,
                           np.array([0, 1, 1, 1, 1,
                                     1, 2, 2, 2, 2,
                                     1, 2, 2, 2, 2,
                                     1, 2, 2, 2, 2]))

    def test_link_num_outlink(self):
        assert_array_equal(self.grid.node_numoutlink,
                           np.array([2, 2, 2, 2, 1,
                                     2, 2, 2, 2, 1,
                                     2, 2, 2, 2, 1,
                                     1, 1, 1, 1, 0]))

    def test_node_inlink_matrix(self):
        self.assertTupleEqual(self.grid.node_inlink_matrix.shape,
                              (2, 20))
        assert_array_equal(self.grid.node_inlink_matrix[0],
                           np.array([-1, 15, 16, 17, 18,
                                      0,  1,  2,  3,  4,
                                      5,  6,  7,  8,  9,
                                     10, 11, 12, 13, 14]))
        assert_array_equal(self.grid.node_inlink_matrix[1],
                           np.array([-1, -1, -1, -1, -1,
                                     -1, 19, 20, 21, 22,
                                     -1, 23, 24, 25, 26,
                                     -1, 27, 28, 29, 30]))

    def test_node_outlink_matrix(self):
        self.assertTupleEqual(self.grid.node_outlink_matrix.shape,
                              (2, 20))
        assert_array_equal(self.grid.node_outlink_matrix[0],
                           np.array([ 0,  1,  2,  3,  4,
                                     5,  6,  7,  8,  9,
                                     10, 11, 12, 13, 14,
                                     27, 28, 29, 30, -1]))
        assert_array_equal(self.grid.node_outlink_matrix[1],
                           np.array([15, 16, 17, 18, -1,
                                     19, 20, 21, 22, -1,
                                     23, 24, 25, 26, -1,
                                     -1, -1, -1, -1, -1]))

    def test_link_face(self):
        BAD = BAD_INDEX_VALUE
        assert_array_equal(self.grid.link_face,
                           np.array([BAD, 0, 1, 2, BAD,
                                     BAD, 3, 4, 5, BAD,
                                     BAD, 6, 7, 8, BAD,
                                     BAD, BAD, BAD, BAD,
                                     9, 10, 11, 12,
                                     13, 14, 15, 16,
                                     BAD, BAD, BAD, BAD]))


if __name__ == '__main__':
    unittest.main()
