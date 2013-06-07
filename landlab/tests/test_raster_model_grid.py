#! /usr/bin/env python
"""
Unit tests for landlab.model_grid
"""

import unittest

from landlab import RasterModelGrid

_SPACING = 1.
(_NUM_ROWS, _NUM_COLS) = (3, 4)

class TestRasterModelGrid(unittest.TestCase):
    def setUp(self):
        """
        These tests use a grid that 3x4 cells::

        |-------|-------|-------|-------|
        |       |       |       |       |
        |   8   |   9   |  10   |  11   |
        |       |       |       |       |
        |-------|---5---|---6---|-------|
        |       |       |       |       |
        |   4   0   5   1   6   2   7   |
        |       |       |       |       |
        |-------|---3---|---4---|-------|
        |       |       |       |       |
        |   0   |   1   |   2   |   3   |
        |       |       |       |       |
        |-------|-------|-------|-------|
        """
        self.grid = RasterModelGrid(num_rows=_NUM_ROWS, num_cols=_NUM_COLS,
                                    dx=_SPACING)

    def test_create_grid(self):
        """
        Create a raster grid.
        """
        grid = RasterModelGrid(num_rows=3, num_cols=4)

    def test_grid_dimensions(self):
        """
        Check the width of the grid. The x-dimension is the columns, and the
        y-dimension is the rows.
        """
        SPACING = 2.
        (NUM_ROWS, NUM_COLS) = (3, 4)

        grid = RasterModelGrid(num_rows=NUM_ROWS, num_cols=NUM_COLS,
                               dx=SPACING)

        self.assertEqual(grid.get_grid_ydimension(), NUM_ROWS * SPACING)
        self.assertEqual(grid.get_grid_xdimension(), NUM_COLS * SPACING)

    def test_grid_dimensions_default_spacing(self):
        """
        Use the default spacing of 1.
        """
        (NUM_ROWS, NUM_COLS) = (3, 4)

        grid = RasterModelGrid(num_rows=NUM_ROWS, num_cols=NUM_COLS)

        self.assertEqual(grid.get_grid_ydimension(), NUM_ROWS * 1.)
        self.assertEqual(grid.get_grid_xdimension(), NUM_COLS * 1.)

    def test_nodes_around_point(self):
        surrounding_ids = self.grid.get_nodes_around_point(2.1, 1.1)
        self.assertEqual(surrounding_ids, [6, 7, 10, 11])

        surrounding_ids = self.grid.get_nodes_around_point(2.1, .9)
        self.assertEqual(surrounding_ids, [2, 3, 6, 7])

    def test_neighbor_list(self):
        neighbors = self.grid.get_neighbor_list(id=6)
        self.assertEqual(list(neighbors), [7, 10, 5, 2])

    def test_neighbor_list_boundary(self):
        """
        All of the neighbor IDs for a boundary cell are -1.
        """
        neighbors = self.grid.get_neighbor_list(id=0)

        self.assertEqual(list(neighbors), [-1, -1, -1, -1])

    def test_cell_x(self):
        expected_x = [0., 1., 2., 3.,
                      0., 1., 2., 3.,
                      0., 1., 2., 3.]

        for (cell_id, expected) in zip(xrange(12), expected_x):
            cell_x = self.grid.x(cell_id)
            self.assertEqual(cell_x, expected)

    def test_cell_y(self):
        expected_y = [0., 0., 0., 0.,
                      1., 1., 1., 1.,
                      2., 2., 2., 2.]

        for (cell_id, expected) in zip(xrange(12), expected_y):
            cell_y = self.grid.y(cell_id)
            self.assertEqual(cell_y, expected)

    def test_cell_x_coordinates(self):
        x_coords = self.grid.get_cell_x_coords()

        self.assertEqual(list(x_coords), [0., 1., 2., 3.,
                                          0., 1., 2., 3.,
                                          0., 1., 2., 3.])

    def test_cell_y_coordinates(self):
        y_coords = self.grid.get_cell_y_coords()

        self.assertEqual(list(y_coords), [0., 0., 0., 0.,
                                          1., 1., 1., 1.,
                                          2., 2., 2., 2.])

    def test_diagonal_list(self):
        diagonals = self.grid.get_diagonal_list(id=5)
        self.assertEqual(list(diagonals), [10, 8, 0, 2])

    def test_diagonal_list_boundary(self):
        diagonals = self.grid.get_diagonal_list(id=0)
        self.assertEqual(list(diagonals), [-1, -1, -1, -1])

    def test_is_interior(self):
        for cell_id in [0, 1, 2, 3, 4, 7, 8, 9, 10, 11]:
            self.assertFalse(self.grid.is_interior(cell_id))

        for cell_id in [5, 6]:
            self.assertTrue(self.grid.is_interior(cell_id))

    def test_get_interior_cells(self):
        interiors = self.grid.get_interior_cells()
        self.assertEqual(list(interiors), [5, 6])


if __name__ == '__main__':
    unittest.main()
