#! /usr/bin/env python

from nose.tools import assert_equal, assert_true, assert_false
try:
    from nose.tools import (assert_tuple_equal, assert_is_not, assert_is_none,
                            assert_is_not_none, assert_is)
except ImportError:
    from landlab.testing.tools import (assert_tuple_equal, assert_is_not,
                                       assert_is_none, assert_is_not_none,
                                       assert_is)
from numpy.testing import assert_array_equal
import numpy as np

from landlab.testing import NumpyArrayTestingMixIn

import landlab.utils.structured_grid as sgrid
from landlab.grid.base import BAD_INDEX_VALUE, CLOSED_BOUNDARY


def test_node_x_2d():
    (x, _) = sgrid.node_coords((3, 2))

    assert_array_equal(x, np.array([0., 1.,
                                    0., 1.,
                                    0., 1., ]))


def test_node_x_2d_with_spacing():
    (x, _) = sgrid.node_coords((3, 2), (2., 10.))

    assert_array_equal(x, np.array([0., 10.,
                                    0., 10.,
                                    0., 10., ]))


def test_node_x_2d_with_origin():
    (x, _) = sgrid.node_coords((3, 2), (2., 10.), (-1., 1.))

    assert_array_equal(x, np.array([1., 11.,
                                    1., 11.,
                                    1., 11., ]))


def test_node_y_2d():
    (_, y) = sgrid.node_coords((3, 2))

    assert_array_equal(y, np.array([0., 0.,
                                    1., 1.,
                                    2., 2., ]))


def test_node_y_2d_with_spacing():
    (_, y) = sgrid.node_coords((3, 2), (2., 10.))

    assert_array_equal(y, np.array([0., 0.,
                                    2., 2.,
                                    4., 4., ]))


def test_node_y_2d_with_origin():
    (_, y) = sgrid.node_coords((3, 2), (2., 10.), (-1., 1.))

    assert_array_equal(y, np.array([-1., -1.,
                                    1.,  1.,
                                    3.,  3., ]))


def test_round_off_error():
    (x, y) = sgrid.node_coords((135, 127),
                               (5.4563957090392, 5.4563957090392),
                               (0., 0.))

    assert_tuple_equal(x.shape, (135 * 127, ))
    assert_tuple_equal(y.shape, (135 * 127, ))


def test_2d_shape_2_by_3():
    cell_nodes = sgrid.node_at_cell((2, 3))

    assert_array_equal(cell_nodes, np.array([]))


def test_2d_shape_3_by_3():
    cell_nodes = sgrid.node_at_cell((3, 3))

    assert_array_equal(cell_nodes, np.array([4]))


def test_shape_4_by_5():
    cell_nodes = sgrid.node_at_cell((4, 5))

    assert_array_equal(cell_nodes, np.array([6,  7,  8,
                                             11, 12, 13]))


# class TestGetNodeLinks(unittest.TestCase, NumpyArrayTestingMixIn):
def test_2d_3_by_2_from_links():
    (from_indices, _) = sgrid.node_index_at_link_ends((3, 2))

    assert_array_equal(from_indices, np.array([0, 1, 2, 3, 0, 2, 4]))


def test_2d_3_by_2_to_links():
    (_, to_indices) = sgrid.node_index_at_link_ends((3, 2))

    assert_array_equal(to_indices, np.array([2, 3, 4, 5, 1, 3, 5]))


def test_west_links():
    links = sgrid.west_links((3, 4))
    assert_array_equal(links, np.array([[-1,  8,  9, 10],
                                        [-1, 11, 12, 13],
                                        [-1, 14, 15, 16]]))

    links = sgrid.west_links((1, 4))
    assert_array_equal(links, np.array([[-1, 0, 1, 2]]))

    links = sgrid.west_links((4, 1))
    assert_array_equal(links, np.array([[-1], [-1], [-1], [-1]]))


def test_east_links():
    links = sgrid.east_links((3, 4))
    assert_array_equal(links, np.array([[8,  9, 10, -1],
                                        [11, 12, 13, -1],
                                        [14, 15, 16, -1]]))

    links = sgrid.east_links((1, 4))
    assert_array_equal(links, np.array([[0, 1, 2, -1]]))

    links = sgrid.east_links((4, 1))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]).T)


def test_north_links():
    links = sgrid.north_links((3, 4))
    assert_array_equal(links, np.array([[0,  1,  2,  3],
                                        [4,  5,  6,  7],
                                        [-1, -1, -1, -1]]))

    links = sgrid.north_links((1, 4))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]))

    links = sgrid.north_links((4, 1))
    assert_array_equal(links, np.array([[0, 1, 2, -1]]).T)


def test_south_links():
    links = sgrid.south_links((3, 4))
    assert_array_equal(links, np.array([[-1, -1, -1, -1],
                                        [0,  1,  2,  3],
                                        [4,  5,  6,  7]]))

    links = sgrid.south_links((1, 4))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]))

    links = sgrid.south_links((4, 1))
    assert_array_equal(links, np.array([[-1, 0, 1, 2]]).T)


def test_inlinks():
    links = sgrid.inlinks((3, 4))
    assert_array_equal(
        np.array([[-1, -1, -1, -1,  0,  1,  2,  3,  4,  5,  6,  7],
                  [-1,  8,  9, 10, -1, 11, 12, 13, -1, 14, 15, 16]]),
        links)


def test_outlinks():
    links = sgrid.outlinks((3, 4))
    assert_array_equal(
        np.array([[0,  1,  2,  3,  4,  5,  6,  7, -1, -1, -1, -1],
                  [8,  9, 10, -1, 11, 12, 13, -1, 14, 15, 16, -1]]),
        links)


# class TestNodeActiveCell(unittest.TestCase, NumpyArrayTestingMixIn):
def test_one_active_cell():
    active_cells = sgrid.active_cell_index_at_nodes((3, 3))

    assert_array_equal(
        active_cells,
        np.array([BAD_INDEX_VALUE, BAD_INDEX_VALUE, BAD_INDEX_VALUE,
                  BAD_INDEX_VALUE,               0, BAD_INDEX_VALUE,
                  BAD_INDEX_VALUE, BAD_INDEX_VALUE, BAD_INDEX_VALUE])
    )


def test_no_active_cells():
    active_cells = sgrid.active_cell_index_at_nodes((3, 2))

    assert_array_equal(
        active_cells,
        np.array([BAD_INDEX_VALUE, BAD_INDEX_VALUE,
                  BAD_INDEX_VALUE, BAD_INDEX_VALUE,
                  BAD_INDEX_VALUE, BAD_INDEX_VALUE])
    )


# class TestActiveCells(unittest.TestCase, NumpyArrayTestingMixIn):
def test_one_active_cell():
    active_cells = sgrid.active_cell_index((3, 3))

    assert_array_equal(active_cells, np.array([0]))


def test_no_active_cells():
    active_cells = sgrid.active_cell_index((3, 2))

    assert_array_equal(active_cells,
                       np.array([]))


# class TestCellCount(unittest.TestCase, NumpyArrayTestingMixIn):
def test_one_cell():
    n_cells = sgrid.cell_count((3, 3))
    assert_equal(n_cells, 1)


def test_no_cells():
    n_cells = sgrid.cell_count((2, 3))
    assert_equal(n_cells, 0)


# class TestInteriorCellCount(unittest.TestCase, NumpyArrayTestingMixIn):
def test_one_cell():
    n_cells = sgrid.interior_cell_count((3, 3))
    assert_equal(n_cells, 1)


def test_no_cells():
    n_cells = sgrid.interior_cell_count((2, 3))
    assert_equal(n_cells, 0)


# class TestActiveCellCount(unittest.TestCase, NumpyArrayTestingMixIn):
def test_one_cell():
    n_cells = sgrid.active_cell_count((3, 3))
    assert_equal(n_cells, 1)


def test_no_cells():
    n_cells = sgrid.active_cell_count((2, 3))
    assert_equal(n_cells, 0)


# class TestInteriorNodes(unittest.TestCase, NumpyArrayTestingMixIn):
def test_4_by_5():
    interiors = sgrid.interior_nodes((4, 5))
    assert_array_equal(interiors, np.array([6, 7, 8, 11, 12, 13]))


def test_no_interiors():
    interiors = sgrid.interior_nodes((2, 3))
    assert_array_equal(interiors, np.array([]))


# class TestNodeStatus(unittest.TestCase, NumpyArrayTestingMixIn):
def test_4_by_5():
    status = sgrid.stat_at_node((4, 5))
    assert_equal(status.dtype, np.int8)
    assert_array_equal(status,
                       np.array([1, 1, 1, 1, 1,
                                 1, 0, 0, 0, 1,
                                 1, 0, 0, 0, 1,
                                 1, 1, 1, 1, 1, ]))


def test_no_interiors():
    status = sgrid.status_at_node((2, 3))
    assert_equal(status.dtype, np.int8)
    assert_array_equal(status,
                       np.array([1, 1, 1,
                                 1, 1, 1, ]))


# class TestActiveLinks(unittest.TestCase, NumpyArrayTestingMixIn):
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


def test_4_by_5():
    active_links = sgrid.active_links((4, 5))
    assert_array_equal(active_links,
                       np.array([1, 2, 3, 6, 7, 8, 11, 12, 13,
                                 19, 20, 21, 22, 23, 24, 25, 26]))
    assert_equal(len(active_links), sgrid.active_link_count((4, 5)))


def test_with_status_at_node():
    status = sgrid.status_at_node((4, 5))
    status[6] = sgrid.CLOSED_BOUNDARY
    active_links = sgrid.active_links((4, 5), node_status_array=status)

    assert_array_equal(active_links,
                       np.array([2, 3, 7, 8, 11, 12, 13,
                                 21, 22, 23, 24, 25, 26]))


def test_with_link_nodes():
    link_nodes = sgrid.node_index_at_link_ends((4, 5))
    active_links = sgrid.active_links((4, 5), link_nodes=link_nodes)

    assert_array_equal(active_links,
                       np.array([1, 2, 3, 6, 7, 8, 11, 12, 13,
                                 19, 20, 21, 22, 23, 24, 25, 26]))
    assert_equal(len(active_links), sgrid.active_link_count((4, 5)))


def test_vertical_active_link_count():
    link_count = sgrid.vertical_active_link_count((3, 4))
    assert_equal(4, link_count)

    link_count = sgrid.vertical_active_link_count((3, 2))
    assert_equal(0, link_count)

    node_status = np.ones((4, 5), dtype=np.int)
    link_count = sgrid.vertical_active_link_count((4, 5),
                                                  node_status=node_status)
    assert_equal(9, link_count)

    link_count = sgrid.vertical_active_link_count((4, 5),
                                                  node_status=node_status)
    node_status[0, 1] = 0
    link_count = sgrid.vertical_active_link_count((4, 5),
                                                  node_status=node_status)
    assert_equal(8, link_count)

    node_status[2, 1] = 0
    link_count = sgrid.vertical_active_link_count((4, 5),
                                                  node_status=node_status)
    assert_equal(6, link_count)

    node_status[2, 2] = 0
    link_count = sgrid.vertical_active_link_count((4, 5),
                                                  node_status=node_status)
    assert_equal(4, link_count)

    node_status[1, 1] = 0
    link_count = sgrid.vertical_active_link_count((4, 5),
                                                  node_status=node_status)
    assert_equal(4, link_count)


def test_horizontal_active_link_count():
    link_count = sgrid.horizontal_active_link_count((3, 4))
    assert_equal(3, link_count)

    link_count = sgrid.horizontal_active_link_count((2, 3))
    assert_equal(0, link_count)

    node_status = np.ones((4, 5), dtype=np.int)
    link_count = sgrid.horizontal_active_link_count(
        (4, 5), node_status=node_status)
    assert_equal(8, link_count)

    link_count = sgrid.horizontal_active_link_count(
        (4, 5), node_status=node_status)
    node_status[0, 1] = 0
    link_count = sgrid.horizontal_active_link_count((4, 5),
                                                    node_status=node_status)
    assert_equal(8, link_count)

    node_status[2, 1] = 0
    link_count = sgrid.horizontal_active_link_count((4, 5),
                                                    node_status=node_status)
    assert_equal(6, link_count)

    node_status[2, 2] = 0
    link_count = sgrid.horizontal_active_link_count((4, 5),
                                                    node_status=node_status)
    assert_equal(5, link_count)

    node_status[1, 1] = 0
    link_count = sgrid.horizontal_active_link_count((4, 5),
                                                    node_status=node_status)
    assert_equal(3, link_count)


def test_horizontal_active_link_ids():
    links = sgrid.horizontal_active_link_ids((3, 4))
    assert_array_equal(links, np.array([[4, 5, 6]]))

    links = sgrid.horizontal_active_link_ids((1, 4))
    expected = np.array([], ndmin=2, dtype=np.int64)
    expected.shape = (0, 3)
    assert_array_equal(expected, links)

    links = sgrid.horizontal_active_link_ids((4, 1))
    expected.shape = (2, 0)
    assert_array_equal(expected, links)

    node_status = np.ones((4, 5), dtype=int)
    links = sgrid.horizontal_active_link_ids((4, 5),
                                             node_status=node_status)
    assert_array_equal(links, np.array([[9, 10, 11, 12],
                                        [13, 14, 15, 16]]))

    node_status = np.ones((4, 5), dtype=int)
    node_status[1, 1] = 0
    links = sgrid.horizontal_active_link_ids((4, 5),
                                             node_status=node_status)
    assert_array_equal(links, np.array([[-1, -1,  7,  8],
                                        [9, 10, 11, 12]]))

    node_status[2, 1] = 0
    links = sgrid.horizontal_active_link_ids((4, 5),
                                             node_status=node_status)
    assert_array_equal(links, np.array([[-1, -1, 6, 7],
                                        [-1, -1, 8, 9]]))

    node_status[0, 0] = 0
    links = sgrid.horizontal_active_link_ids((4, 5),
                                             node_status=node_status)
    assert_array_equal(links, np.array([[-1, -1, 6, 7],
                                        [-1, -1, 8, 9]]))


def test_vertical_active_link_ids():
    links = sgrid.vertical_active_link_ids((3, 4))
    assert_array_equal(links, np.array([[0, 1], [2, 3]]))

    links = sgrid.vertical_active_link_ids((1, 4))
    expected = np.array([], ndmin=2, dtype=np.int64)
    expected.shape = (0, 2)
    assert_array_equal(expected, links)

    links = sgrid.vertical_active_link_ids((4, 1))
    expected.shape = (3, 0)
    assert_array_equal(expected, links)

    node_status = np.ones((4, 5), dtype=int)
    links = sgrid.vertical_active_link_ids((4, 5), node_status=node_status)
    assert_array_equal(links, np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]]))

    node_status = np.ones((4, 5), dtype=int)
    node_status[1, 1] = 0
    links = sgrid.vertical_active_link_ids((4, 5), node_status=node_status)
    assert_array_equal(links,
                       np.array([[-1, 0, 1], [-1, 2, 3], [4, 5, 6]]))

    node_status[2, 1] = 0
    links = sgrid.vertical_active_link_ids((4, 5), node_status=node_status)
    assert_array_equal(links,
                       np.array([[-1, 0, 1], [-1, 2, 3], [-1, 4, 5]]))

    node_status[0, 0] = 0
    links = sgrid.vertical_active_link_ids((4, 5), node_status=node_status)
    assert_array_equal(links,
                       np.array([[-1, 0, 1], [-1, 2, 3], [-1, 4, 5]]))


def test_west_links():
    links = sgrid.active_west_links((3, 4))
    assert_array_equal(links, np.array([[-1, -1, -1, -1],
                                        [-1,  4,  5,  6],
                                        [-1, -1, -1, -1]]))

    links = sgrid.active_west_links((1, 4))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]))

    links = sgrid.active_west_links((4, 1))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]).T)


def test_east_links():
    links = sgrid.active_east_links((3, 4))
    assert_array_equal(links, np.array([[-1, -1, -1, -1],
                                        [4,  5,  6, -1],
                                        [-1, -1, -1, -1]]))

    links = sgrid.active_east_links((1, 4))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]))

    links = sgrid.active_east_links((4, 1))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]).T)

    links = sgrid.horizontal_active_link_ids((4, 5))
    assert_array_equal(np.array([[9, 10, 11, 12],
                                 [13, 14, 15, 16]]),
                       links)

    links = sgrid.active_east_links((4, 5))
    assert_array_equal(np.array([[-1, -1, -1, -1, -1],
                                 [9, 10, 11, 12, -1],
                                 [13, 14, 15, 16, -1],
                                 [-1, -1, -1, -1, -1]]),
                       links)


def test_north_links():
    links = sgrid.active_north_links((3, 4))
    assert_array_equal(links, np.array([[-1,  0,  1, -1],
                                        [-1,  2,  3, -1],
                                        [-1, -1, -1, -1]]))

    links = sgrid.active_north_links((1, 4))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]))

    links = sgrid.active_north_links((4, 1))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]).T)


def test_south_links():
    links = sgrid.active_south_links((3, 4))
    assert_array_equal(links, np.array([[-1, -1, -1, -1],
                                        [-1,  0,  1, -1],
                                        [-1,  2,  3, -1]]))

    links = sgrid.active_south_links((1, 4))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]))

    links = sgrid.active_south_links((4, 1))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]).T)


def test_inlinks():
    links = sgrid.active_inlinks((3, 4))
    assert_array_equal(
        np.array([[-1, -1, -1, -1, -1,  0,  1, -1, -1,  2,  3, -1],
                  [-1, -1, -1, -1, -1,  4,  5,  6, -1, -1, -1, -1]]),
        links)


def test_outlinks():
    links = sgrid.active_outlinks((3, 4))
    assert_array_equal(
        np.array([[-1,  0,  1, -1, -1,  2,  3, -1, -1, -1, -1, -1],
                  [-1, -1, -1, -1,  4,  5,  6, -1, -1, -1, -1, -1]]),
        links)


def test_outlinks_4x5():
    links = sgrid.active_outlinks((4, 5))

    assert_array_equal(np.array([[-1,  0,  1,  2, -1,
                                  -1,  3,  4,  5, -1,
                                  -1,  6,  7,  8, -1,
                                  -1, -1, -1, -1, -1],
                                 [-1, -1, -1, -1, -1,
                                  9, 10, 11, 12, -1,
                                  13, 14, 15, 16, -1,
                                  -1, -1, -1, -1, -1]]),
                       links)


def test_inlinks_4x5():
    links = sgrid.active_inlinks((4, 5))

    assert_array_equal(np.array([[-1, -1, -1, -1, -1,
                                  -1,  0,  1,  2, -1,
                                  -1,  3,  4,  5, -1,
                                  -1,  6, 7,  8, -1],
                                 [-1, -1, -1, -1, -1,
                                  -1,  9, 10, 11, 12,
                                  -1, 13, 14, 15, 16,
                                  -1, -1, -1, -1, -1]]),
                       links)


# class TestFaces(unittest.TestCase, NumpyArrayTestingMixIn):
def test_face_count():
    assert_equal(17, sgrid.face_count((4, 5)))
    assert_equal(4, sgrid.face_count((3, 3)))
    assert_equal(0, sgrid.face_count((2, 100)))
    assert_equal(0, sgrid.face_count((100, 2)))
    assert_equal(0, sgrid.face_count((100, 1)))


def test_active_face_count():
    assert_equal(17, sgrid.active_face_count((4, 5)))
    assert_equal(4, sgrid.active_face_count((3, 3)))
    assert_equal(0, sgrid.active_face_count((2, 100)))
    assert_equal(0, sgrid.active_face_count((100, 2)))
    assert_equal(0, sgrid.active_face_count((100, 1)))


def test_active_faces():
    active_faces = sgrid.active_face_index((4, 5))
    assert_array_equal(np.array([0,  1,  2,
                                 3,  4,  5,
                                 6,  7,  8,
                                 9, 10, 11, 12,
                                 13, 14, 15, 16]),
                       active_faces)


# class TestLinkFaces(unittest.TestCase, NumpyArrayTestingMixIn):
def test_4_by_5():
    link_faces = sgrid.face_at_link((4, 5))

    BAD = sgrid.BAD_INDEX_VALUE

    assert_array_equal(link_faces, np.array([BAD, 0, 1, 2, BAD,
                                             BAD, 3, 4, 5, BAD,
                                             BAD, 6, 7, 8, BAD,
                                             BAD, BAD, BAD, BAD,
                                             9, 10, 11, 12,
                                             13, 14, 15, 16,
                                             BAD, BAD, BAD, BAD]))


def test_with_active_links():
    active_links = sgrid.active_links((4, 5))
    active_links = active_links[:-1]
    link_faces = sgrid.face_at_link((4, 5), actives=active_links)

    BAD = sgrid.BAD_INDEX_VALUE

    assert_array_equal(link_faces, np.array([BAD, 0, 1, 2, BAD,
                                             BAD, 3, 4, 5, BAD,
                                             BAD, 6, 7, 8, BAD,
                                             BAD, BAD, BAD, BAD,
                                             9, 10, 11, 12,
                                             13, 14, 15, BAD,
                                             BAD, BAD, BAD, BAD]))


# class TestReshapeArray(unittest.TestCase, NumpyArrayTestingMixIn):
def test_default():
    x = np.arange(12.)
    y = sgrid.reshape_array((3, 4), x)

    assert_equal(y.shape, (3, 4))
    assert_array_equal(x, y.flat)
    assert_true(y.flags['C_CONTIGUOUS'])
    assert_is(y.base, x)


def test_copy():
    x = np.arange(12.)
    y = sgrid.reshape_array((3, 4), x, copy=True)

    assert_equal(y.shape, (3, 4))
    assert_array_equal(x, y.flat)
    assert_true(y.flags['C_CONTIGUOUS'])
    assert_is_none(y.base)

    y[0][0] = 0.
    assert_array_equal(x, np.array([0.,  1.,  2., 3.,
                                    4.,  5.,  6., 7.,
                                    8.,  9., 10., 11.]))


def test_flip():
    x = np.arange(12.)
    y = sgrid.reshape_array((3, 4), x, flip_vertically=True)

    assert_equal(y.shape, (3, 4))
    assert_array_equal(y, np.array([[8.,  9., 10., 11.],
                                    [4.,  5.,  6., 7.],
                                    [0.,  1.,  2., 3.]]))
    assert_false(y.flags['C_CONTIGUOUS'])
    assert_is_not_none(y.base)

    y[0][0] = 0.
    assert_array_equal(x, np.array([0.,  1.,  2., 3.,
                                    4.,  5.,  6., 7.,
                                    0.,  9., 10., 11.]))


def test_flip_copy():
    x = np.arange(12.)
    y = sgrid.reshape_array((3, 4), x, flip_vertically=True, copy=True)

    assert_equal(y.shape, (3, 4))
    assert_array_equal(y, np.array([[8.,  9., 10., 11.],
                                    [4.,  5.,  6., 7.],
                                    [0.,  1.,  2., 3.]]))
    assert_true(y.flags['C_CONTIGUOUS'])
    assert_is_not(y.base, x)


# class TestDiagonalArray(unittest.TestCase, NumpyArrayTestingMixIn):
def test_default():
    diags = sgrid.diagonal_node_array((2, 3), out_of_bounds=-1)
    assert_array_equal(diags,
                       np.array([[4, -1, -1, -1],
                                 [5,  3, -1, -1],
                                 [-1,  4, -1, -1],
                                 [-1, -1, -1,  1],
                                 [-1, -1,  0,  2],
                                 [-1, -1,  1, -1]]))

    assert_true(diags.base is None)
    assert_true(diags.flags['C_CONTIGUOUS'])


def test_non_contiguous():
    diags = sgrid.diagonal_node_array((2, 3), out_of_bounds=-1,
                                      contiguous=False)
    assert_array_equal(diags,
                       np.array([[4, -1, -1, -1],
                                 [5,  3, -1, -1],
                                 [-1,  4, -1, -1],
                                 [-1, -1, -1,  1],
                                 [-1, -1,  0,  2],
                                 [-1, -1,  1, -1]]))

    assert_true(isinstance(diags.base, np.ndarray))
    assert_false(diags.flags['C_CONTIGUOUS'])


def test_boundary_node_mask_no_actives():
    diags = sgrid.diagonal_node_array((2, 3), out_of_bounds=-1,
                                      boundary_node_mask=-2)
    assert_array_equal(diags, - 2 * np.ones((6, 4)))


def test_boundary_node_mask():
    diags = sgrid.diagonal_node_array((3, 3), out_of_bounds=-1,
                                      boundary_node_mask=-2)
    assert_array_equal(diags,
                       np.array([[-2, -2, -2, -2], [-2, -2, -2, -2], [-2, -2, -2, -2],
                                 [-2, -2, -2, -2], [8,  6,
                                                    0,  2], [-2, -2, -2, -2],
                                 [-2, -2, -2, -2], [-2, -2, -2, -2], [-2, -2, -2, -2]]))


# class TestNeighborArray(unittest.TestCase, NumpyArrayTestingMixIn):
def test_default():
    neighbors = sgrid.neighbor_node_array((2, 3))

    BAD = sgrid.BAD_INDEX_VALUE

    assert_array_equal(
        neighbors,
        np.array([[1,   3, BAD, BAD],
                  [2,   4,   0, BAD],
                  [BAD,   5,   1, BAD],
                  [4, BAD, BAD,   0],
                  [5, BAD,   3,   1],
                  [BAD, BAD,   4,   2]]).T)

    assert_true(neighbors.flags['C_CONTIGUOUS'])
    assert_true(neighbors.base is None)


def test_set_out_of_bounds():
    neighbors = sgrid.neighbor_node_array((2, 3), inactive=-1)
    assert_array_equal(neighbors,
                       np.array([[1,  3, -1, -1],
                                 [2,  4,  0, -1],
                                 [-1,  5,  1, -1],
                                 [4, -1, -1,  0],
                                 [5, -1,  3,  1],
                                 [-1, -1,  4,  2]]).T)


# class TestInlinkMatrix(unittest.TestCase, NumpyArrayTestingMixIn):
def test_no_inactive():
    inlinks = sgrid.setup_active_inlink_matrix((4, 5), return_count=False)
    assert_array_equal(inlinks,
                       np.array([[-1, -1, -1, -1, -1,
                                  -1,  0,  1,  2, -1,
                                  -1,  3,  4,  5, -1,
                                  -1,  6,  7,  8, -1],
                                 [-1, -1, -1, -1, -1,
                                  -1,  9, 10, 11, 12,
                                  -1, 13, 14, 15, 16,
                                  -1, -1, -1, -1, -1]]))


def test_inactive():
    status = np.ones((4, 5))
    status[1, 1] = 0
    inlinks = sgrid.setup_active_inlink_matrix((4, 5), return_count=False,
                                               node_status=status)
    assert_array_equal(inlinks,
                       np.array([[-1, -1, -1, -1, -1,
                                  -1, -1,  0,  1, -1,
                                  -1, -1,  2,  3, -1,
                                  -1,  4,  5,  6, -1],
                                 [-1, -1, -1, -1, -1,
                                  -1, -1, -1,  7,  8,
                                  -1,  9, 10, 11, 12,
                                  -1, -1, -1, -1, -1]]))


def test_out_link_ids_at_nodes():
    links_ids = sgrid.outlink_index_at_node((4, 5))
    assert_array_equal(
        np.array([[ 4,  5,  6,  7,  8,
                   13, 14, 15, 16, 17,
                   22, 23, 24, 25, 26,
                   -1, -1, -1, -1, -1],
                  [ 0,  1,  2,  3, -1,
                    9, 10, 11, 12, -1,
                   18, 19, 20, 21, -1,
                   27, 28, 29, 30, -1]]),
        links_ids)


def test_in_link_ids_at_nodes():
    links_ids = sgrid.inlink_index_at_node((4, 5))
    assert_array_equal(
        np.array([[-1, -1, -1, -1, -1,  
                    4,  5,  6,  7,  8,
                   13, 14, 15, 16, 17, 
                   22, 23, 24, 25, 26],
                  [-1,  0,  1,  2,  3, 
                   -1,  9, 10, 11, 12,
                   -1, 18, 19, 20, 21,
                   -1, 27, 28, 29, 30]]),
        links_ids)
