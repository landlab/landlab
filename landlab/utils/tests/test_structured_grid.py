#! /usr/bin/env python

import unittest
import numpy as np
from numpy.testing import assert_array_equal

import landlab.utils.structured_grid as sgrid
from landlab.utils.structured_grid import BAD_INDEX_VALUE


class TestGetNodeCoords(unittest.TestCase):
    def test_node_x_2d(self):
        (x, _) = sgrid.node_coords((3,2))

        assert_array_equal(x, np.array([0., 1.,
                                        0., 1.,
                                        0., 1.,]))

    def test_node_x_2d_with_spacing(self):
        (x, _) = sgrid.node_coords((3,2), (2., 10.))

        assert_array_equal(x, np.array([0., 10.,
                                        0., 10.,
                                        0., 10.,]))

    def test_node_x_2d_with_origin(self):
        (x, _) = sgrid.node_coords((3,2), (2., 10.), (-1., 1.))

        assert_array_equal(x, np.array([1., 11.,
                                        1., 11.,
                                        1., 11.,]))

    def test_node_y_2d(self):
        (_, y) = sgrid.node_coords((3,2))

        assert_array_equal(y, np.array([0., 0.,
                                        1., 1.,
                                        2., 2.,]))

    def test_node_y_2d_with_spacing(self):
        (_, y) = sgrid.node_coords((3,2), (2., 10.))

        assert_array_equal(y, np.array([0., 0.,
                                        2., 2.,
                                        4., 4.,]))

    def test_node_y_2d_with_origin(self):
        (_, y) = sgrid.node_coords((3,2), (2., 10.), (-1., 1.))

        assert_array_equal(y, np.array([-1., -1.,
                                         1.,  1.,
                                         3.,  3.,]))


class TestGetCellNode(unittest.TestCase):
    def test_2d_shape_2_by_3(self):
        cell_nodes = sgrid.cell_node_index((2, 3))

        assert_array_equal(cell_nodes, np.array([]))

    def test_2d_shape_3_by_3(self):
        cell_nodes = sgrid.cell_node_index((3, 3))

        assert_array_equal(cell_nodes, np.array([4]))

    def test_shape_4_by_5(self):
        cell_nodes = sgrid.cell_node_index((4, 5))

        assert_array_equal(cell_nodes, np.array([ 6,  7,  8,
                                                 11, 12, 13]))


class TestGetNodeLinks(unittest.TestCase):
    def test_2d_3_by_2_from_links(self):
        (from_indices, _) = sgrid.node_link_index((3, 2))

        assert_array_equal(from_indices,
                           np.array([0, 1, 2, 3,
                                     0, 2, 4]))

    def test_2d_3_by_2_to_links(self):
        (_, to_indices) = sgrid.node_link_index((3, 2))

        assert_array_equal(to_indices,
                           np.array([2, 3, 4, 5,
                                     1, 3, 5]))


class TestNodeActiveCell(unittest.TestCase):
    def test_one_active_cell(self):
        active_cells = sgrid.node_active_cell((3, 3))

        assert_array_equal(
            active_cells,
            np.array([BAD_INDEX_VALUE, BAD_INDEX_VALUE, BAD_INDEX_VALUE,
                      BAD_INDEX_VALUE,               0, BAD_INDEX_VALUE,
                      BAD_INDEX_VALUE, BAD_INDEX_VALUE, BAD_INDEX_VALUE])
        )

    def test_no_active_cells(self):
        active_cells = sgrid.node_active_cell((3, 2))

        assert_array_equal(
            active_cells,
            np.array([BAD_INDEX_VALUE, BAD_INDEX_VALUE,
                      BAD_INDEX_VALUE, BAD_INDEX_VALUE,
                      BAD_INDEX_VALUE, BAD_INDEX_VALUE])
        )
                                     
                                     
class TestActiveCells(unittest.TestCase):
    def test_one_active_cell(self):
        active_cells = sgrid.active_cells((3, 3))

        assert_array_equal(active_cells, np.array([0]))

    def test_no_active_cells(self):
        active_cells = sgrid.active_cells((3, 2))

        assert_array_equal(active_cells,
            np.array([]))


class TestCellCount(unittest.TestCase):
    def test_one_cell(self):
        n_cells = sgrid.cell_count((3, 3))
        self.assertEqual(n_cells, 1)

    def test_no_cells(self):
        n_cells = sgrid.cell_count((2, 3))
        self.assertEqual(n_cells, 0)


class TestInteriorCellCount(unittest.TestCase):
    def test_one_cell(self):
        n_cells = sgrid.interior_cell_count((3, 3))
        self.assertEqual(n_cells, 1)

    def test_no_cells(self):
        n_cells = sgrid.interior_cell_count((2, 3))
        self.assertEqual(n_cells, 0)


class TestActiveCellCount(unittest.TestCase):
    def test_one_cell(self):
        n_cells = sgrid.active_cell_count((3, 3))
        self.assertEqual(n_cells, 1)

    def test_no_cells(self):
        n_cells = sgrid.active_cell_count((2, 3))
        self.assertEqual(n_cells, 0)


class TestInteriorNodes(unittest.TestCase):
    def test_4_by_5(self):
        interiors = sgrid.interior_nodes((4, 5))
        assert_array_equal(interiors, np.array([6, 7, 8, 11, 12, 13]))

    def test_no_interiors(self):
        interiors = sgrid.interior_nodes((2, 3))
        assert_array_equal(interiors, np.array([]))


class TestNodeStatus(unittest.TestCase):
    def test_4_by_5(self):
        status = sgrid.node_status((4, 5))
        self.assertEqual(status.dtype, np.int8)
        assert_array_equal(status,
                           np.array([1, 1, 1, 1, 1,
                                     1, 0, 0, 0, 1,
                                     1, 0, 0, 0, 1,
                                     1, 1, 1, 1, 1, ]))

    def test_no_interiors(self):
        status = sgrid.node_status((2, 3))
        self.assertEqual(status.dtype, np.int8)
        assert_array_equal(status,
                           np.array([1, 1, 1,
                                     1, 1, 1,]))


class TestActiveLinks(unittest.TestCase):
    """
    *--27-->*--28-->*--29-->*--30-->*
    ^       ^       ^       ^       ^
    10      11      12      13      14
    |       |       |       |       |
    *--23-->*--24-->*--25-->*--26-->*
    ^       ^       ^       ^       ^
    5       6       7       8       9   
    |       |       |       |       |
    *--19-->*--20-->*--21-->*--22-->*
    ^       ^       ^       ^       ^
    0       1       2       3       4
    |       |       |       |       |
    *--15-->*--16-->*--17-->*--18-->*
    """
    def test_4_by_5(self):
        active_links = sgrid.active_links((4, 5))
        assert_array_equal(active_links,
                           np.array([1, 2, 3, 6, 7, 8, 11, 12, 13,
                                     19, 20, 21, 22, 23, 24, 25, 26]))
        self.assertEqual(len(active_links), sgrid.active_link_count((4, 5)))

    def test_with_node_status(self):
        status = sgrid.node_status((4, 5))
        status[6] = sgrid.INACTIVE_BOUNDARY
        active_links = sgrid.active_links((4, 5), node_status_array=status)

        assert_array_equal(active_links,
                           np.array([2, 3, 7, 8, 11, 12, 13,
                                     21, 22, 23, 24, 25, 26]))

    def test_with_link_nodes(self):
        link_nodes = sgrid.node_link_index((4, 5))
        active_links = sgrid.active_links((4, 5), link_nodes=link_nodes)

        assert_array_equal(active_links,
                           np.array([1, 2, 3, 6, 7, 8, 11, 12, 13,
                                     19, 20, 21, 22, 23, 24, 25, 26]))
        self.assertEqual(len(active_links), sgrid.active_link_count((4, 5)))


class TestLinkFaces(unittest.TestCase):
    def test_4_by_5(self):
        link_faces = sgrid.link_faces((4, 5))

        BAD = sgrid.BAD_INDEX_VALUE

        assert_array_equal(link_faces, np.array([BAD, 0, 1, 2, BAD,
                                                 BAD, 3, 4, 5, BAD,
                                                 BAD, 6, 7, 8, BAD,
                                                 BAD, BAD, BAD, BAD,
                                                 9, 10, 11, 12,
                                                 13, 14, 15, 16,
                                                 BAD, BAD, BAD, BAD]))

    def test_with_active_links(self):
        active_links = sgrid.active_links((4, 5))
        active_links = active_links[:-1]
        link_faces = sgrid.link_faces((4, 5), actives=active_links)

        BAD = sgrid.BAD_INDEX_VALUE

        assert_array_equal(link_faces, np.array([BAD, 0, 1, 2, BAD,
                                                 BAD, 3, 4, 5, BAD,
                                                 BAD, 6, 7, 8, BAD,
                                                 BAD, BAD, BAD, BAD,
                                                 9, 10, 11, 12,
                                                 13, 14, 15, BAD,
                                                 BAD, BAD, BAD, BAD]))


class TestReshapeArray(unittest.TestCase):
    def test_default(self):
        x = np.arange(12.)
        y = sgrid.reshape_array((3, 4), x)

        self.assertEqual(y.shape, (3, 4))
        assert_array_equal(x, y.flat)
        self.assertTrue(y.flags['C_CONTIGUOUS'])
        self.assertIs(y.base, x)

    def test_copy(self):
        x = np.arange(12.)
        y = sgrid.reshape_array((3, 4), x, copy=True)

        self.assertEqual(y.shape, (3, 4))
        assert_array_equal(x, y.flat)
        self.assertTrue(y.flags['C_CONTIGUOUS'])
        self.assertIsNone(y.base)

        y[0][0] = 0.
        assert_array_equal(x, np.array([ 0.,  1.,  2., 3.,
                                         4.,  5.,  6., 7.,
                                         8.,  9., 10., 11.]))

    def test_flip(self):
        x = np.arange(12.)
        y = sgrid.reshape_array((3, 4), x, flip_vertically=True)

        self.assertEqual(y.shape, (3, 4))
        assert_array_equal(y, np.array([[ 8.,  9., 10., 11.],
                                        [ 4.,  5.,  6., 7.],
                                        [ 0.,  1.,  2., 3.]]))
        self.assertFalse(y.flags['C_CONTIGUOUS'])
        self.assertIsNotNone(y.base)

        y[0][0] = 0.
        assert_array_equal(x, np.array([ 0.,  1.,  2., 3.,
                                         4.,  5.,  6., 7.,
                                         0.,  9., 10., 11.]))

    def test_flip_copy(self):
        x = np.arange(12.)
        y = sgrid.reshape_array((3, 4), x, flip_vertically=True, copy=True)

        self.assertEqual(y.shape, (3, 4))
        assert_array_equal(y, np.array([[ 8.,  9., 10., 11.],
                                        [ 4.,  5.,  6., 7.],
                                        [ 0.,  1.,  2., 3.]]))
        self.assertTrue(y.flags['C_CONTIGUOUS'])
        self.assertIsNot(y.base, x)


class TestDiagonalArray(unittest.TestCase):
    def test_default(self):
        diags = sgrid.diagonal_node_array((2, 3), out_of_bounds=-1)
        assert_array_equal(diags,
                           np.array([[ 4, -1, -1, -1],
                                     [ 5,  3, -1, -1],
                                     [-1,  4, -1, -1],
                                     [-1, -1, -1,  1],
                                     [-1, -1,  0,  2],
                                     [-1, -1,  1, -1]]))

        self.assertTrue(diags.base is None)
        self.assertTrue(diags.flags['C_CONTIGUOUS'])

    def test_non_contiguous(self):
        diags = sgrid.diagonal_node_array((2, 3), out_of_bounds=-1,
                                          contiguous=False)
        assert_array_equal(diags,
                           np.array([[ 4, -1, -1, -1],
                                     [ 5,  3, -1, -1],
                                     [-1,  4, -1, -1],
                                     [-1, -1, -1,  1],
                                     [-1, -1,  0,  2],
                                     [-1, -1,  1, -1]]))

        self.assertTrue(isinstance(diags.base, np.ndarray))
        self.assertFalse(diags.flags['C_CONTIGUOUS'])

    def test_boundary_node_mask_no_actives(self):
        diags = sgrid.diagonal_node_array((2, 3), out_of_bounds=-1,
                                          boundary_node_mask=-2)
        assert_array_equal(diags, - 2 * np.ones((6, 4)))

    def test_boundary_node_mask(self):
        diags = sgrid.diagonal_node_array((3, 3), out_of_bounds=-1,
                                          boundary_node_mask=-2)
        assert_array_equal(diags, 
                           np.array([[-2, -2, -2, -2], [-2, -2, -2, -2], [-2, -2, -2, -2],
                                     [-2, -2, -2, -2], [ 8,  6,  0,  2], [-2, -2, -2, -2],
                                     [-2, -2, -2, -2], [-2, -2, -2, -2], [-2, -2, -2, -2]]))


class TestNeighborArray(unittest.TestCase):
    def test_default(self):
        neighbors = sgrid.neighbor_node_array((2, 3))

        BAD = sgrid.BAD_INDEX_VALUE

        assert_array_equal(
            neighbors,
            np.array([[  1,   3, BAD, BAD],
                      [  2,   4,   0, BAD],
                      [BAD,   5,   1, BAD],
                      [  4, BAD, BAD,   0],
                      [  5, BAD,   3,   1],
                      [BAD, BAD,   4,   2]]))

        self.assertTrue(neighbors.flags['C_CONTIGUOUS'])
        self.assertTrue(neighbors.base is None)

    def test_set_out_of_bounds(self):
        neighbors = sgrid.neighbor_node_array((2, 3), out_of_bounds=-1)
        assert_array_equal(neighbors,
                           np.array([[ 1,  3, -1, -1],
                                     [ 2,  4,  0, -1],
                                     [-1,  5,  1, -1],
                                     [ 4, -1, -1,  0],
                                     [ 5, -1,  3,  1],
                                     [-1, -1,  4,  2]]))

    def test_as_view(self):
        neighbors = sgrid.neighbor_node_array((2, 3), out_of_bounds=-1,
                                              contiguous=False)
        assert_array_equal(neighbors,
                           np.array([[ 1,  3, -1, -1],
                                     [ 2,  4,  0, -1],
                                     [-1,  5,  1, -1],
                                     [ 4, -1, -1,  0],
                                     [ 5, -1,  3,  1],
                                     [-1, -1,  4,  2]]))

        self.assertFalse(neighbors.flags['C_CONTIGUOUS'])
        self.assertTrue(isinstance(neighbors.base, np.ndarray))

    def test_boundary_node_mask_no_actives(self):
        neighbors = sgrid.neighbor_node_array((2, 3), out_of_bounds=-1,
                                          boundary_node_mask=-2)
        assert_array_equal(neighbors, - 2 * np.ones((6, 4)))

    def test_boundary_node_mask(self):
        neighbors = sgrid.neighbor_node_array((3, 3), out_of_bounds=-1,
                                              boundary_node_mask=-2)
        assert_array_equal(neighbors, 
                           np.array([[-2, -2, -2, -2], [-2, -2, -2, -2], [-2, -2, -2, -2],
                                     [-2, -2, -2, -2], [ 5,  7,  3,  1], [-2, -2, -2, -2],
                                     [-2, -2, -2, -2], [-2, -2, -2, -2], [-2, -2, -2, -2]]))


if __name__ == '__main__':
    unittest.main()
