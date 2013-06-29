#! /usr/bin/env python
"""
Unit tests for landlab.model_grid
"""

import unittest

from landlab import RasterModelGrid

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
        self.assertEqual(surrounding_ids, [7, 8, 12, 13],
                         'incorrect surrounding ids in test_nodes_around_point')

        surrounding_ids = self.grid.get_nodes_around_point(2.1, .9)
        self.assertEqual(surrounding_ids, [2, 3, 7, 8],
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
            cell_x = self.grid.x(cell_id)
            self.assertEqual(cell_x, expected)

    def test_cell_y(self):
        expected_y = [0., 0., 0., 0., 0.,
                      1., 1., 1., 1., 1.,
                      2., 2., 2., 2., 2.,
                      3., 3., 3., 3., 3.]

        for (cell_id, expected) in zip(xrange(12), expected_y):
            cell_y = self.grid.y(cell_id)
            self.assertEqual(cell_y, expected)

    def test_node_x_coordinates(self):
        x_coords = self.grid.get_node_x_coords()

        self.assertEqual(list(x_coords), [0., 1., 2., 3., 4.,
                                          0., 1., 2., 3., 4.,
                                          0., 1., 2., 3., 4.,
                                          0., 1., 2., 3., 4.])

    def test_cell_y_coordinates(self):
        y_coords = self.grid.get_cell_y_coords()

        self.assertEqual(list(y_coords), [0., 0., 0., 0., 0.,
                                          1., 1., 1., 1., 1.,
                                          2., 2., 2., 2., 2.,
                                          3., 3., 3., 3., 3.])

    def test_diagonal_list(self):
        diagonals = self.grid.get_diagonal_list(id=6)
        self.assertEqual(list(diagonals), [12, 10, 0, 2])

    def test_diagonal_list_boundary(self):
        diagonals = self.grid.get_diagonal_list(id=0)
        self.assertEqual(list(diagonals), [1, -1, -1, -1])

    def test_is_interior(self):
        for cell_id in [0, 1, 2, 3, 4, 5, 9, 10, 14, 15, 16, 17, 18, 19]:
            self.assertFalse(self.grid.is_interior(cell_id))

        for cell_id in [6, 7, 8, 11, 12, 13]:
            self.assertTrue(self.grid.is_interior(cell_id))

    def test_get_interior_cells(self):
        interiors = self.grid.get_interior_cells()
        self.assertEqual(list(interiors), [6, 7, 8, 11, 12, 13])


if __name__ == '__main__':
    unittest.main()
