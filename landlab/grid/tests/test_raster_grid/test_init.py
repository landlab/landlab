import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import (assert_is, assert_equal, assert_raises, raises,
                        assert_true, assert_false)

from landlab import RasterModelGrid
from landlab import BAD_INDEX_VALUE


class TestRasterModelGrid(object):
    def setup(self):
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
        self.num_rows = 4
        self.num_cols = 5
        self.spacing = 1.
        self.rmg = RasterModelGrid(self.num_rows, self.num_cols,
                                   dx=self.spacing)

    def test_create_grid(self):
        """
        Create a raster grid.
        """
        grid = RasterModelGrid(num_rows=4, num_cols=5)

    def test_grid_dimensions(self):
        """
        Use the default spacing of 1.
        """
        assert_equal(self.rmg.get_grid_ydimension(), self.num_rows)
        assert_equal(self.rmg.get_grid_xdimension(), self.num_cols)

    def test_grid_dimensions_non_unit_spacing(self):
        """
        Check the width of the grid. The x-dimension is the columns, and the
        y-dimension is the rows.
        """
        rmg = RasterModelGrid(4, 5, dx=2.)
        assert_equal(rmg.get_grid_ydimension(), 8.)
        assert_equal(rmg.get_grid_xdimension(), 10.)

    def test_nodes_around_point(self):
        surrounding_ids = self.rmg.get_nodes_around_point(2.1, 1.1)
        assert_array_equal(surrounding_ids, np.array([7, 12, 13, 8]))

        surrounding_ids = self.rmg.get_nodes_around_point(2.1, .9)
        assert_array_equal(surrounding_ids, np.array([2, 7, 8, 3]))

    def test_neighbor_list_with_scalar_arg(self):
        assert_array_equal(self.rmg.get_neighbor_list(6),
                           np.array([7, 11, 5, 1]))
        assert_array_equal(self.rmg.get_neighbor_list(-1),
                           np.array([-1, -1, -1, -1]))

    def test_neighbor_list_with_array_arg(self):
        assert_array_equal(self.rmg.get_neighbor_list([6, -1]),
                           np.array([[7, 11, 5, 1], [-1, -1, -1, -1]]))

    def test_neighbor_list_with_no_args(self):
        assert_array_equal(
            self.rmg.get_neighbor_list(),
            np.array([[-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1],
                      [-1, -1, -1, -1], [ 7, 11,  5,  1], [ 8, 12,  6,  2], [ 9, 13,  7,  3], [-1, -1, -1, -1],
                      [-1, -1, -1, -1], [12, 16, 10,  6], [13, 17, 11,  7], [14, 18, 12,  8], [-1, -1, -1, -1],
                      [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1]]))

    def test_neighbor_list_boundary(self):
        """
        All of the neighbor IDs for a boundary cell are -1.
        """
        import landlab.utils.structured_grid as sgrid
        for node_id in sgrid.boundary_iter(self.rmg.shape):
            assert_array_equal(self.rmg.get_neighbor_list(node_id),
                               np.array([-1, -1, -1, -1]))

    def test_node_x(self):
        assert_array_equal(self.rmg.node_x,
                           np.array([0., 1., 2., 3., 4.,
                                     0., 1., 2., 3., 4.,
                                     0., 1., 2., 3., 4.,
                                     0., 1., 2., 3., 4.]))

    def test_node_y(self):
        assert_array_equal(self.rmg.node_y,
                              np.array([0., 0., 0., 0., 0.,
                                        1., 1., 1., 1., 1.,
                                        2., 2., 2., 2., 2.,
                                        3., 3., 3., 3., 3.]))

    @raises(ValueError)
    def test_node_x_is_immutable(self):
        self.rmg.node_x[0] = 0

    @raises(ValueError)
    def test_node_y_is_immutable(self):
        self.rmg.node_y[0] = 0

    def test_node_axis_coordinates(self):
        assert_is(self.rmg.node_axis_coordinates(axis=0).base,
                  self.rmg.node_y.base)
        assert_is(self.rmg.node_axis_coordinates(axis=1).base,
                  self.rmg.node_x.base)
        assert_is(self.rmg.node_axis_coordinates(axis=-1).base,
                  self.rmg.node_x.base)
        assert_is(self.rmg.node_axis_coordinates(axis=-2).base,
                  self.rmg.node_y.base)

    def test_diagonal_list(self):
        assert_array_equal(self.rmg.get_diagonal_list(6),
                           np.array([12, 10, 0, 2]))
        assert_array_equal(self.rmg.get_diagonal_list(-1),
                           np.array([-1, -1, -1, -1]))
        assert_array_equal(self.rmg.get_diagonal_list([6, -1]),
                           np.array([[12, 10, 0, 2], [-1, -1, -1, -1]]))
        assert_array_equal(
            self.rmg.get_diagonal_list(),
            np.array([[-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1],
                        [-1, -1, -1, -1], [-1, -1, -1, -1],
                      [-1, -1, -1, -1], [12, 10,  0,  2], [13, 11,  1,  3],
                        [14, 12,  2,  4], [-1, -1, -1, -1],
                      [-1, -1, -1, -1], [17, 15,  5,  7], [18, 16,  6,  8],
                        [19, 17,  7,  9], [-1, -1, -1, -1],
                      [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1],
                        [-1, -1, -1, -1], [-1, -1, -1, -1]]))

    def test_diagonal_list_boundary(self):
        assert_array_equal(self.rmg.get_diagonal_list(0),
                           np.array([-1, -1, -1, -1]))

    def test_is_interior(self):
        for cell_id in [0, 1, 2, 3, 4, 5, 9, 10, 14, 15, 16, 17, 18, 19]:
            assert_false(self.rmg.is_interior(cell_id))
        for cell_id in [6, 7, 8, 11, 12, 13]:
            assert_true(self.rmg.is_interior(cell_id))

    def test_get_interior_cells(self):
        assert_array_equal(self.rmg.get_active_cell_node_ids(),
                           np.array([6, 7, 8, 11, 12, 13]))

    def test_active_links(self):
        assert_equal(self.rmg.number_of_active_links, 17)
        assert_array_equal(self.rmg.active_links,
                           np.array([1, 2, 3, 6, 7, 8, 11, 12, 13,
                                     19, 20, 21, 22, 23, 24, 25, 26]))

    def test_active_link_fromnode(self):
        assert_array_equal(self.rmg.activelink_fromnode,
                           np.array([1, 2, 3, 6, 7, 8, 11, 12, 13,
                                     5, 6, 7, 8, 10, 11, 12, 13]))

    def test_active_link_tonode(self):
        assert_array_equal(self.rmg.activelink_tonode,
                           np.array([6, 7, 8, 11, 12, 13, 16, 17, 18,
                                     6, 7, 8, 9, 11, 12, 13, 14]))

    def test_active_link_num_inlink(self):
        assert_array_equal(self.rmg.node_numactiveinlink,
                           np.array([0, 0, 0, 0, 0,
                                     0, 2, 2, 2, 1,
                                     0, 2, 2, 2, 1,
                                     0, 1, 1, 1, 0]))

    def test_active_link_num_outlink(self):
        assert_array_equal(self.rmg.node_numactiveoutlink,
                           np.array([0, 1, 1, 1, 0,
                                     1, 2, 2, 2, 0,
                                     1, 2, 2, 2, 0,
                                     0, 0, 0, 0, 0]))

    def test_active_inlink_matrix(self):
        assert_array_equal(
            self.rmg.node_active_inlink_matrix,
            np.array([[-1, -1, -1, -1, -1,
                       -1,  0,  1,  2, -1,
                       -1,  3,  4,  5, -1,
                       -1,  6, 7,  8, -1],
                      [-1, -1, -1, -1, -1,
                       -1,  9, 10, 11, 12,
                       -1, 13, 14, 15, 16,
                       -1, -1, -1, -1, -1]]))

    def test_active_outlink_matrix(self):
        assert_array_equal(
            self.rmg.node_active_outlink_matrix,
            np.array([[-1,  0,  1,  2, -1,
                       -1,  3,  4,  5, -1,
                       -1,  6,  7,  8, -1,
                       -1, -1, -1, -1, -1],
                      [-1, -1, -1, -1, -1,
                        9, 10, 11, 12, -1,
                       13, 14, 15, 16, -1,
                       -1, -1, -1, -1, -1]]))

    def test_active_node_links_scalar_interior(self):
        assert_array_equal(
            self.rmg.active_node_links([6]),
            np.array([[0, 9, 3, 10]]).T)

    def test_active_node_links_scalar_boundary(self):
        assert_array_equal(
            self.rmg.active_node_links([1]),
            np.array([[-1, -1, 0, -1]]).T)

    def test_active_node_with_array_arg(self):
        assert_array_equal(
            self.rmg.active_node_links([6, 7]),
            np.array([[0,  9, 3, 10],
                      [1, 10, 4, 11]]).T)

    def test_active_node_links_with_no_args(self):
        assert_array_equal(
            self.rmg.active_node_links(),
            np.array([[-1, -1, -1, -1, -1, -1,  0,  1,  2, -1,
                       -1,  3,  4,  5, -1, -1,  6,  7,  8, -1],
                      [-1, -1, -1, -1, -1, -1,  9, 10, 11, 12,
                       -1, 13, 14, 15, 16, -1, -1, -1, -1, -1],
                      [-1,  0,  1,  2, -1, -1,  3,  4,  5, -1,
                       -1,  6,  7,  8, -1, -1, -1, -1, -1, -1],
                      [-1, -1, -1, -1, -1,  9, 10, 11, 12, -1,
                       13, 14, 15, 16, -1, -1, -1, -1, -1, -1]]))

    def test_link_fromnode(self):
        assert_array_equal(
            self.rmg.link_fromnode,
            np.array(
                [0, 1, 2, 3, 4, 5, 6, 7,  8,  9, 10, 11, 12, 13, 14,
                 0, 1, 2, 3, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 18]))

    def test_link_tonode(self):
        assert_array_equal(
            self.rmg.link_tonode,
            np.array(
                [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                 1, 2, 3, 4, 6,  7,  8,  9, 11, 12, 13, 14, 16, 17, 18, 19]))

    def test_link_num_inlink(self):
        assert_array_equal(self.rmg.node_numinlink,
                           np.array([0, 1, 1, 1, 1,
                                     1, 2, 2, 2, 2,
                                     1, 2, 2, 2, 2,
                                     1, 2, 2, 2, 2]))

    def test_link_num_outlink(self):
        assert_array_equal(self.rmg.node_numoutlink,
                           np.array([2, 2, 2, 2, 1,
                                     2, 2, 2, 2, 1,
                                     2, 2, 2, 2, 1,
                                     1, 1, 1, 1, 0]))

    def test_node_inlink_matrix(self):
        assert_array_equal(self.rmg.node_inlink_matrix,
                           np.array([[-1, -1, -1, -1, -1,
                                       0,  1,  2,  3,  4,
                                       5,  6,  7,  8,  9,
                                      10, 11, 12, 13, 14],
                                     [-1, 15, 16, 17, 18,
                                      -1, 19, 20, 21, 22,
                                      -1, 23, 24, 25, 26,
                                      -1, 27, 28, 29, 30]]))

    def test_node_outlink_matrix(self):
        assert_array_equal(self.rmg.node_outlink_matrix,
                           np.array([[ 0,  1,  2,  3,  4,
                                       5,  6,  7,  8,  9,
                                      10, 11, 12, 13, 14,
                                      -1, -1, -1, -1, -1],
                                     [15, 16, 17, 18, -1,
                                      19, 20, 21, 22, -1,
                                      23, 24, 25, 26, -1,
                                      27, 28, 29, 30, -1]]))

    def test_node_links_with_scalar_interior(self):
        assert_array_equal(self.rmg.node_links([6]),
                           np.array([[1, 19, 6, 20]]).T)

    def test_node_links_with_scalar_boundary(self):
        assert_array_equal(self.rmg.node_links([1]),
                           np.array([[-1, 15, 1, 16]]).T)

    def test_node_links_with_array_arg(self):
        assert_array_equal(self.rmg.node_links([6, 7]),
                           np.array([[1, 19, 6, 20], [2, 20, 7, 21]]).T)

    def test_node_links_with_no_args(self):
        assert_array_equal(
            self.rmg.node_links(),
            np.array([[-1, -1, -1, -1, -1,  0,  1,  2,  3,  4,
                        5,  6,  7,  8,  9, 10, 11, 12, 13, 14],
                      [-1, 15, 16, 17, 18, -1, 19, 20, 21, 22,
                       -1, 23, 24, 25, 26, -1, 27, 28, 29, 30],
                      [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                       10, 11, 12, 13, 14, -1, -1, -1, -1, -1],
                      [15, 16, 17, 18, -1, 19, 20, 21, 22, -1,
                       23, 24, 25, 26, -1, 27, 28, 29, 30, -1]]))

    def test_link_face(self):
        BAD = BAD_INDEX_VALUE
        assert_array_equal(self.rmg.link_face,
                           np.array([BAD, 0, 1, 2, BAD,
                                     BAD, 3, 4, 5, BAD,
                                     BAD, 6, 7, 8, BAD,
                                     BAD, BAD, BAD, BAD,
                                     9, 10, 11, 12,
                                     13, 14, 15, 16,
                                     BAD, BAD, BAD, BAD]))

    def test_grid_coords_to_node_id_with_scalar(self):
        assert_equal(self.rmg.grid_coords_to_node_id(3, 4), 19)

    def test_grid_coords_to_node_id_with_array(self):
        assert_array_equal(self.rmg.grid_coords_to_node_id((3, 2), (4, 1)),
                           np.array([19, 11]))

    def test_grid_coords_to_node_id_outside_of_grid(self):
        assert_raises(ValueError, self.rmg.grid_coords_to_node_id, 5, 0)

    def test_create_diagonal_list(self):
        self.rmg.create_diagonal_list()

        assert_array_equal(
            self.rmg.get_diagonal_list(),
            np.array([[-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1],
                      [-1, -1, -1, -1], [-1, -1, -1, -1], [12, 10,  0,  2], [13, 11,  1,  3],
                      [14, 12,  2,  4], [-1, -1, -1, -1], [-1, -1, -1, -1], [17, 15,  5,  7],
                      [18, 16,  6,  8], [19, 17,  7,  9], [-1, -1, -1, -1], [-1, -1, -1, -1],
                      [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1]]))


