#! /usr/bin/env python
"""
Unit tests for landlab.model_grid
"""

import unittest
import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid, BAD_INDEX_VALUE, FIXED_GRADIENT_BOUNDARY


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
            list(surrounding_ids), [7, 12, 13, 8],
            #list(surrounding_ids), [7, 8, 12, 13],
            'incorrect surrounding ids in test_nodes_around_point')

        surrounding_ids = self.grid.get_nodes_around_point(2.1, .9)
        self.assertListEqual(
            list(surrounding_ids), [2, 7, 8, 3],
            #list(surrounding_ids), [2, 3, 7, 8],
            'incorrect surrounding ids in test_nodes_around_point')

    def test_neighbor_list(self):
        neighbors = self.grid.get_neighbor_list(6)
        assert_array_equal(neighbors, np.array([7, 11, 5, 1]))

        neighbors = self.grid.get_neighbor_list(-1)
        assert_array_equal(neighbors, np.array([-1, -1, -1, -1]))

        neighbors = self.grid.get_neighbor_list([6, -1])
        assert_array_equal(neighbors,
                           np.array([[7, 11, 5, 1], [-1, -1, -1, -1]]))

        neighbors = self.grid.get_neighbor_list()
        assert_array_equal(
            neighbors,
            np.array([[-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1],
                      [-1, -1, -1, -1], [ 7, 11,  5,  1], [ 8, 12,  6,  2], [ 9, 13,  7,  3], [-1, -1, -1, -1],
                      [-1, -1, -1, -1], [12, 16, 10,  6], [13, 17, 11,  7], [14, 18, 12,  8], [-1, -1, -1, -1],
                      [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1]]))

    def test_neighbor_list_boundary(self):
        """
        All of the neighbor IDs for a boundary cell are -1.
        """
        neighbors = self.grid.get_neighbor_list(0)

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
        diagonals = self.grid.get_diagonal_list(6)
        assert_array_equal(diagonals, np.array([12, 10, 0, 2]))

        diagonals = self.grid.get_diagonal_list(-1)
        assert_array_equal(diagonals, np.array([-1, -1, -1, -1]))

        diagonals = self.grid.get_diagonal_list([6, -1])
        assert_array_equal(diagonals,
                           np.array([[12, 10, 0, 2], [-1, -1, -1, -1]]))

        diagonals = self.grid.get_diagonal_list()
        assert_array_equal(
            diagonals,
            np.array([[-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1],
                      [-1, -1, -1, -1], [12, 10,  0,  2], [13, 11,  1,  3], [14, 12,  2,  4], [-1, -1, -1, -1],
                      [-1, -1, -1, -1], [17, 15,  5,  7], [18, 16,  6,  8], [19, 17,  7,  9], [-1, -1, -1, -1],
                      [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1]]))

    def test_diagonal_list_boundary(self):
        diagonals = self.grid.get_diagonal_list(0)
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

    def test_grid_coords_to_node_id(self):
        self.assertEqual(self.grid.grid_coords_to_node_id(3, 4), 19)
        assert_array_equal(self.grid.grid_coords_to_node_id((3, 2), (4, 1)),
                           np.array([19, 11]))

        with self.assertRaises(ValueError):
            self.grid.grid_coords_to_node_id(5, 0)

    def test_create_diagonal_list(self):
        self.grid.create_diagonal_list()

        assert_array_equal(
            self.grid.get_diagonal_list(),
            np.array([[-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1],
                      [-1, -1, -1, -1], [-1, -1, -1, -1], [12, 10,  0,  2], [13, 11,  1,  3],
                      [14, 12,  2,  4], [-1, -1, -1, -1], [-1, -1, -1, -1], [17, 15,  5,  7],
                      [18, 16,  6,  8], [19, 17,  7,  9], [-1, -1, -1, -1], [-1, -1, -1, -1],
                      [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1]]))


class TestZerosArray(unittest.TestCase):
    def test_default(self):
        mg = RasterModelGrid(4, 5)
        v = mg.zeros()
        assert_array_equal(v, np.zeros(20., dtype=np.float))

    def test_int_array(self):
        mg = RasterModelGrid(4, 5)
        v = mg.zeros(dtype=np.int)
        assert_array_equal(v, np.zeros(20., dtype=np.int))

    def test_at_cells(self):
        mg = RasterModelGrid(4, 5)
        v = mg.zeros(centering='cell')
        assert_array_equal(v, np.zeros(6, dtype=np.float))

    def test_at_links(self):
        mg = RasterModelGrid(4, 5)
        v = mg.zeros(centering='link')
        assert_array_equal(v, np.zeros(17, dtype=np.float))

    def test_at_faces(self):
        mg = RasterModelGrid(4, 5)
        v = mg.zeros(centering='face')
        assert_array_equal(v, np.zeros(17, dtype=np.float))


class TestRasterModelGridZeroArrays(unittest.TestCase):
    def test_default(self):
        mg = RasterModelGrid(4, 5)
        assert_array_equal(mg.zeros(), np.zeros((20, ), dtype=np.float))

    def test_at_nodes(self):
        mg = RasterModelGrid(4, 5)
        assert_array_equal(mg.zeros(centering='node'),
                           np.zeros((20, ), dtype=np.float))

    def test_at_cells(self):
        mg = RasterModelGrid(4, 5)
        assert_array_equal(mg.zeros(centering='cell'),
                           np.zeros((6, ), dtype=np.float))

    def test_at_links(self):
        mg = RasterModelGrid(4, 5)
        assert_array_equal(mg.zeros(centering='link'),
                           np.zeros((17, ), dtype=np.float))

    def test_at_faces(self):
        mg = RasterModelGrid(4, 5)
        assert_array_equal(mg.zeros(centering='face'),
                           np.zeros((17, ), dtype=np.float))

    def test_bad_centering(self):
        mg = RasterModelGrid(4, 5)
        with self.assertRaises(TypeError):
            mg.zeros(centering='bad_centering_string')

    def test_int(self):
        mg = RasterModelGrid(4, 5)
        assert_array_equal(mg.zeros(dtype=np.int),
                           np.zeros((20, ), dtype=np.int))

    def test_int8(self):
        mg = RasterModelGrid(4, 5)
        assert_array_equal(mg.zeros(dtype=np.int8),
                           np.zeros((20, ), dtype=np.int8))

class TestRasterModelGridIsBoundary(unittest.TestCase):
    def test_id_as_int(self):
        mg = RasterModelGrid(4, 5)
        self.assertTrue(mg.is_boundary(0))

    def test_id_as_small_list(self):
        mg = RasterModelGrid(4, 5)
        assert_array_equal(
            mg.is_boundary([0]),
            np.array([True]))

    def test_id_as_array(self):
        mg = RasterModelGrid(4, 5)
        assert_array_equal(
            mg.is_boundary(np.arange(20)),
            np.array([True,  True,  True,  True,  True,
                      True, False, False, False,  True,  
                      True, False, False, False,  True,  
                      True,  True,  True,  True,  True], dtype=bool))

    def test_id_as_list(self):
        mg = RasterModelGrid(4, 5)
        assert_array_equal(
            mg.is_boundary([8, 9]),
            np.array([False, True]))

    def test_boundary_flag(self):
        mg = RasterModelGrid(4, 5)
        mg.node_status[0] = FIXED_GRADIENT_BOUNDARY
        assert_array_equal(
            mg.is_boundary(np.arange(20)),
            np.array([True,  True,  True,  True,  True,
                      True, False, False, False,  True,  
                      True, False, False, False,  True,  
                      True,  True,  True,  True,  True], dtype=bool))

        assert_array_equal(
            mg.is_boundary(np.arange(20), boundary_flag=FIXED_GRADIENT_BOUNDARY),
            np.array([ True, False, False, False, False,
                      False, False, False, False, False,
                      False, False, False, False, False,
                      False, False, False, False, False], dtype=bool))


class TestRasterModelGridHasBoundaryNeighbor(unittest.TestCase):
    def test_boundary_node(self):
        mg = RasterModelGrid(5, 6)
        self.assertTrue(mg.has_boundary_neighbor(0))
        self.assertFalse(mg.has_boundary_neighbor(14))

    def test_last_index(self):
        mg = RasterModelGrid(4, 5)
        self.assertTrue(mg.has_boundary_neighbor(-1))

    def test_id_as_list(self):
        mg = RasterModelGrid(4, 5)
        assert_array_equal(
            mg.has_boundary_neighbor([-1, 0]),
            np.array([True, True])
        )

    def test_id_as_array(self):
        mg = RasterModelGrid(4, 5)
        assert_array_equal(
            mg.has_boundary_neighbor(np.arange(20)),
            np.array(
                [True, True, True, True, True,
                 True, True, True, True, True,
                 True, True, True, True, True,
                 True, True, True, True, True]))

        mg = RasterModelGrid(5, 5)
        assert_array_equal(
            mg.has_boundary_neighbor(np.arange(25)),
            np.array(
                [True, True,  True, True, True,
                 True, True,  True, True, True,
                 True, True, False, True, True,
                 True, True,  True, True, True,
                 True, True,  True, True, True]))


class TestRasterModelGridConnectingFaces(unittest.TestCase):
    def test_horizontally_adjacent_cells(self):
        mg = RasterModelGrid(4, 5)
        assert_array_equal(
            mg.get_face_connecting_cell_pair(0, 1),
            np.array([10]))

    def test_vertically_adjacent_cells(self):
        mg = RasterModelGrid(4, 5)
        assert_array_equal(
            mg.get_face_connecting_cell_pair(0, 3),
            np.array([3]))

    def test_diagonally_adjacent_cells(self):
        mg = RasterModelGrid(4, 5)
        assert_array_equal(
            mg.get_face_connecting_cell_pair(1, 5),
            np.array([]))

    def test_non_adjacent_cells(self):
        mg = RasterModelGrid(4, 5)
        assert_array_equal(
            mg.get_face_connecting_cell_pair(0, 2),
            np.array([]))


class TestRasterModelGridCellFaces(unittest.TestCase):
    def test_id_as_int(self):
        mg = RasterModelGrid(4, 5)
        assert_array_equal(
            mg.cell_faces(0),
            np.array([0, 9, 3, 10]))

    def test_id_as_array(self):
        mg = RasterModelGrid(4, 5)
        assert_array_equal(
            mg.cell_faces(np.array([0, 1])),
            np.array([[0, 9, 3, 10],
                      [1, 10, 4, 11]]))


class TestRasterModelGridNodesAroundPoint(unittest.TestCase):
    def test_lower_left_cell(self):
        mg = RasterModelGrid(3, 3)
        assert_array_equal(
            mg.get_nodes_around_point(.1, .1),
            np.array([0, 3, 4, 1]))
        assert_array_equal(
            mg.get_nodes_around_point(.1, .9),
            np.array([0, 3, 4, 1]))
        assert_array_equal(
            mg.get_nodes_around_point(.9, .9),
            np.array([0, 3, 4, 1]))
        assert_array_equal(
            mg.get_nodes_around_point(.9, .1),
            np.array([0, 3, 4, 1]))


class TestRasterModelGridCellAreas(unittest.TestCase):
    def test_unit_grid_all_cells(self):
        mg = RasterModelGrid(4, 4)
        assert_array_equal(mg.get_cell_areas(), np.ones(4.))

    def test_unit_grid_one_cell(self):
        mg = RasterModelGrid(4, 4)
        self.assertEqual(mg.get_cell_areas(0), 1.)

    def test_unit_grid_last_cell(self):
        mg = RasterModelGrid(4, 4)
        self.assertEqual(mg.get_cell_areas(-1), 1.)

    def test_out_of_range(self):
        mg = RasterModelGrid(4, 4)
        with self.assertRaises(IndexError):
            mg.get_cell_areas(5)

    def test_all_cells_with_spacing(self):
        mg = RasterModelGrid(4, 4, 10.)
        assert_array_equal(mg.get_cell_areas(), 100. * np.ones(4.))


if __name__ == '__main__':
    unittest.main()
