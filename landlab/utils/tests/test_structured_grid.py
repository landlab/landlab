#! /usr/bin/env python

import unittest
import numpy as np
from numpy.testing import assert_array_equal

import landlab.utils.structured_grid as sgrid
from landlab.utils.structured_grid import BAD_INDEX_VALUE


class TestGetNodeCoords(unittest.TestCase):
    def test_node_x_2d(self):
        (x, _, _) = sgrid.node_xyz((3,2))

        assert_array_equal(x, np.array([0., 1.,
                                        0., 1.,
                                        0., 1.,]))

    def test_node_x_2d_with_spacing(self):
        (x, _, _) = sgrid.node_xyz((3,2), (2., 10.))

        assert_array_equal(x, np.array([0., 10.,
                                        0., 10.,
                                        0., 10.,]))

    def test_node_x_2d_with_origin(self):
        (x, _, _) = sgrid.node_xyz((3,2), (2., 10.), (-1., 1.))

        assert_array_equal(x, np.array([1., 11.,
                                        1., 11.,
                                        1., 11.,]))

    def test_node_y_2d(self):
        (_, y, _) = sgrid.node_xyz((3,2))

        assert_array_equal(y, np.array([0., 0.,
                                        1., 1.,
                                        2., 2.,]))

    def test_node_y_2d_with_spacing(self):
        (_, y, _) = sgrid.node_xyz((3,2), (2., 10.))

        assert_array_equal(y, np.array([0., 0.,
                                        2., 2.,
                                        4., 4.,]))

    def test_node_y_2d_with_origin(self):
        (_, y, _) = sgrid.node_xyz((3,2), (2., 10.), (-1., 1.))

        assert_array_equal(y, np.array([-1., -1.,
                                         1.,  1.,
                                         3.,  3.,]))

    def test_node_z_2d(self):
        (_, _, z) = sgrid.node_xyz((3,2))

        assert_array_equal(z, np.zeros(3 * 2))


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
                      BAD_INDEX_VALUE,               4, BAD_INDEX_VALUE,
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
        status = sgrid.node_boundary_status((4, 5))
        assert_array_equal(status,
                           np.array([1, 1, 1, 1, 1,
                                     1, 0, 0, 0, 1,
                                     1, 0, 0, 0, 1,
                                     1, 1, 1, 1, 1, ]))

    def test_no_interiors(self):
        status = sgrid.node_boundary_status((2, 3))
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

    def test_with_node_status(self):
        status = sgrid.node_boundary_status((4, 5))
        status[6] = sgrid.INACTIVE_BOUNDARY
        active_links = sgrid.active_links((4, 5), node_status=status)

        assert_array_equal(active_links,
                           np.array([2, 3, 7, 8, 11, 12, 13,
                                     21, 22, 23, 24, 25, 26]))

    def test_with_link_nodes(self):
        link_nodes = sgrid.node_link_index((4, 5))
        active_links = sgrid.active_links((4, 5), link_nodes=link_nodes)

        assert_array_equal(active_links,
                           np.array([1, 2, 3, 6, 7, 8, 11, 12, 13,
                                     19, 20, 21, 22, 23, 24, 25, 26]))


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


if __name__ == '__main__':
    unittest.main()
