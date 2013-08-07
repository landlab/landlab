#! /usr/bin/env python

import unittest
import numpy as np
from numpy.testing import assert_array_equal

import landlab.utils.structured_grid as sgrid


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


if __name__ == '__main__':
    unittest.main()
