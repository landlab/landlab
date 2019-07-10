#! /usr/bin/env python

import numpy as np
from numpy.testing import assert_array_equal

import landlab.utils.structured_grid as sgrid


def test_node_x_2d():
    (x, _) = sgrid.node_coords((3, 2))

    assert_array_equal(x, np.array([0.0, 1.0, 0.0, 1.0, 0.0, 1.0]))


def test_node_x_2d_with_spacing():
    (x, _) = sgrid.node_coords((3, 2), (2.0, 10.0))

    assert_array_equal(x, np.array([0.0, 10.0, 0.0, 10.0, 0.0, 10.0]))


def test_node_x_2d_with_origin():
    (x, _) = sgrid.node_coords((3, 2), (2.0, 10.0), (-1.0, 1.0))

    assert_array_equal(x, np.array([1.0, 11.0, 1.0, 11.0, 1.0, 11.0]))


def test_node_y_2d():
    (_, y) = sgrid.node_coords((3, 2))

    assert_array_equal(y, np.array([0.0, 0.0, 1.0, 1.0, 2.0, 2.0]))


def test_node_y_2d_with_spacing():
    (_, y) = sgrid.node_coords((3, 2), (2.0, 10.0))

    assert_array_equal(y, np.array([0.0, 0.0, 2.0, 2.0, 4.0, 4.0]))


def test_node_y_2d_with_origin():
    (_, y) = sgrid.node_coords((3, 2), (2.0, 10.0), (-1.0, 1.0))

    assert_array_equal(y, np.array([-1.0, -1.0, 1.0, 1.0, 3.0, 3.0]))


def test_round_off_error():
    (x, y) = sgrid.node_coords(
        (135, 127), (5.4563957090392, 5.4563957090392), (0.0, 0.0)
    )

    assert x.shape == (135 * 127,)
    assert y.shape == (135 * 127,)


def test_2d_shape_2_by_3():
    cell_nodes = sgrid.node_at_cell((2, 3))

    assert_array_equal(cell_nodes, np.array([]))


def test_2d_shape_3_by_3():
    cell_nodes = sgrid.node_at_cell((3, 3))

    assert_array_equal(cell_nodes, np.array([4]))


def test_shape_4_by_5():
    cell_nodes = sgrid.node_at_cell((4, 5))

    assert_array_equal(cell_nodes, np.array([6, 7, 8, 11, 12, 13]))


def test_2d_3_by_2_from_links():
    (from_indices, _) = sgrid.node_index_at_link_ends((3, 2))

    assert_array_equal(from_indices, np.array([0, 1, 2, 3, 0, 2, 4]))


def test_2d_3_by_2_to_links():
    (_, to_indices) = sgrid.node_index_at_link_ends((3, 2))

    assert_array_equal(to_indices, np.array([2, 3, 4, 5, 1, 3, 5]))


def test_west_links():
    links = sgrid.west_links((3, 4))
    assert_array_equal(
        links, np.array([[-1, 0, 1, 2], [-1, 7, 8, 9], [-1, 14, 15, 16]])
    )

    links = sgrid.west_links((1, 4))
    assert_array_equal(links, np.array([[-1, 0, 1, 2]]))

    links = sgrid.west_links((4, 1))
    assert_array_equal(links, np.array([[-1], [-1], [-1], [-1]]))


def test_east_links():
    links = sgrid.east_links((3, 4))
    assert_array_equal(
        links, np.array([[0, 1, 2, -1], [7, 8, 9, -1], [14, 15, 16, -1]])
    )

    links = sgrid.east_links((1, 4))
    assert_array_equal(links, np.array([[0, 1, 2, -1]]))

    links = sgrid.east_links((4, 1))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]).T)


def test_north_links():
    links = sgrid.north_links((3, 4))
    assert_array_equal(
        links, np.array([[3, 4, 5, 6], [10, 11, 12, 13], [-1, -1, -1, -1]])
    )

    links = sgrid.north_links((1, 4))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]))

    links = sgrid.north_links((4, 1))
    assert_array_equal(links, np.array([[0, 1, 2, -1]]).T)


def test_south_links():
    links = sgrid.south_links((3, 4))
    assert_array_equal(
        links, np.array([[-1, -1, -1, -1], [3, 4, 5, 6], [10, 11, 12, 13]])
    )

    links = sgrid.south_links((1, 4))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]))

    links = sgrid.south_links((4, 1))
    assert_array_equal(links, np.array([[-1, 0, 1, 2]]).T)


def test_inlinks():
    links = sgrid.inlinks((3, 4))
    assert_array_equal(
        np.array(
            [
                [-1, -1, -1, -1, 3, 4, 5, 6, 10, 11, 12, 13],
                [-1, 0, 1, 2, -1, 7, 8, 9, -1, 14, 15, 16],
            ]
        ),
        links,
    )


def test_outlinks():
    links = sgrid.outlinks((3, 4))
    assert_array_equal(
        np.array(
            [
                [3, 4, 5, 6, 10, 11, 12, 13, -1, -1, -1, -1],
                [0, 1, 2, -1, 7, 8, 9, -1, 14, 15, 16, -1],
            ]
        ),
        links,
    )


def test_cell_count_one_cell():
    n_cells = sgrid.cell_count((3, 3))
    assert n_cells == 1


def test_no_cells():
    n_cells = sgrid.cell_count((2, 3))
    assert n_cells == 0


def test_interior_cell_count_one_cell():
    n_cells = sgrid.interior_cell_count((3, 3))
    assert n_cells == 1


def test_interior_cell_count_no_cells():
    n_cells = sgrid.interior_cell_count((2, 3))
    assert n_cells == 0


def test_active_cell_count_one_cell():
    n_cells = sgrid.active_cell_count((3, 3))
    assert n_cells == 1


def test_active_cell_count_no_cells():
    n_cells = sgrid.active_cell_count((2, 3))
    assert n_cells == 0


def test_interior_nodes_4_by_5():
    interiors = sgrid.interior_nodes((4, 5))
    assert_array_equal(interiors, np.array([6, 7, 8, 11, 12, 13]))


def test_no_interiors():
    interiors = sgrid.interior_nodes((2, 3))
    assert_array_equal(interiors, np.array([]))


def test_node_status_4_by_5():
    status = sgrid.status_at_node((4, 5))
    assert status.dtype == np.int8
    assert_array_equal(
        status,
        np.array(
            [[1, 1, 1, 1, 1], [1, 0, 0, 0, 1], [1, 0, 0, 0, 1], [1, 1, 1, 1, 1]]
        ).flatten(),
    )


def test_node_status_no_interiors():
    status = sgrid.status_at_node((2, 3))
    assert status.dtype == np.int8
    assert_array_equal(status, np.array([1, 1, 1, 1, 1, 1]))


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
    assert_array_equal(
        active_links,
        np.array([1, 2, 3, 6, 7, 8, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26]),
    )
    assert len(active_links) == sgrid.active_link_count((4, 5))


def test_with_status_at_node():
    status = sgrid.status_at_node((4, 5))
    status[6] = sgrid.CLOSED_BOUNDARY
    active_links = sgrid.active_links((4, 5), node_status_array=status)

    assert_array_equal(
        active_links, np.array([2, 3, 7, 8, 11, 12, 13, 21, 22, 23, 24, 25, 26])
    )


def test_with_link_nodes():
    link_nodes = sgrid.node_index_at_link_ends((4, 5))
    active_links = sgrid.active_links((4, 5), link_nodes=link_nodes)

    assert_array_equal(
        active_links,
        np.array([1, 2, 3, 6, 7, 8, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26]),
    )
    assert len(active_links) == sgrid.active_link_count((4, 5))


def test_vertical_active_link_count():
    link_count = sgrid.vertical_active_link_count((3, 4))
    assert 4 == link_count

    link_count = sgrid.vertical_active_link_count((3, 2))
    assert 0 == link_count

    node_status = np.ones((4, 5), dtype=np.int)
    link_count = sgrid.vertical_active_link_count((4, 5), node_status=node_status)
    assert 9 == link_count

    link_count = sgrid.vertical_active_link_count((4, 5), node_status=node_status)
    node_status[0, 1] = 0
    link_count = sgrid.vertical_active_link_count((4, 5), node_status=node_status)
    assert 8 == link_count

    node_status[2, 1] = 0
    link_count = sgrid.vertical_active_link_count((4, 5), node_status=node_status)
    assert 6 == link_count

    node_status[2, 2] = 0
    link_count = sgrid.vertical_active_link_count((4, 5), node_status=node_status)
    assert 4 == link_count

    node_status[1, 1] = 0
    link_count = sgrid.vertical_active_link_count((4, 5), node_status=node_status)
    assert 4 == link_count


def test_horizontal_active_link_count():
    link_count = sgrid.horizontal_active_link_count((3, 4))
    assert 3 == link_count

    link_count = sgrid.horizontal_active_link_count((2, 3))
    assert 0 == link_count

    node_status = np.ones((4, 5), dtype=np.int)
    link_count = sgrid.horizontal_active_link_count((4, 5), node_status=node_status)
    assert 8 == link_count

    link_count = sgrid.horizontal_active_link_count((4, 5), node_status=node_status)
    node_status[0, 1] = 0
    link_count = sgrid.horizontal_active_link_count((4, 5), node_status=node_status)
    assert 8 == link_count

    node_status[2, 1] = 0
    link_count = sgrid.horizontal_active_link_count((4, 5), node_status=node_status)
    assert 6 == link_count

    node_status[2, 2] = 0
    link_count = sgrid.horizontal_active_link_count((4, 5), node_status=node_status)
    assert 5 == link_count

    node_status[1, 1] = 0
    link_count = sgrid.horizontal_active_link_count((4, 5), node_status=node_status)
    assert 3 == link_count


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
    links = sgrid.horizontal_active_link_ids((4, 5), node_status=node_status)
    assert_array_equal(links, np.array([[9, 10, 11, 12], [13, 14, 15, 16]]))

    node_status = np.ones((4, 5), dtype=int)
    node_status[1, 1] = 0
    links = sgrid.horizontal_active_link_ids((4, 5), node_status=node_status)
    assert_array_equal(links, np.array([[-1, -1, 7, 8], [9, 10, 11, 12]]))

    node_status[2, 1] = 0
    links = sgrid.horizontal_active_link_ids((4, 5), node_status=node_status)
    assert_array_equal(links, np.array([[-1, -1, 6, 7], [-1, -1, 8, 9]]))

    node_status[0, 0] = 0
    links = sgrid.horizontal_active_link_ids((4, 5), node_status=node_status)
    assert_array_equal(links, np.array([[-1, -1, 6, 7], [-1, -1, 8, 9]]))


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
    assert_array_equal(links, np.array([[-1, 0, 1], [-1, 2, 3], [4, 5, 6]]))

    node_status[2, 1] = 0
    links = sgrid.vertical_active_link_ids((4, 5), node_status=node_status)
    assert_array_equal(links, np.array([[-1, 0, 1], [-1, 2, 3], [-1, 4, 5]]))

    node_status[0, 0] = 0
    links = sgrid.vertical_active_link_ids((4, 5), node_status=node_status)
    assert_array_equal(links, np.array([[-1, 0, 1], [-1, 2, 3], [-1, 4, 5]]))


def test_active_west_links():
    links = sgrid.active_west_links((3, 4))
    assert_array_equal(
        links, np.array([[-1, -1, -1, -1], [-1, 4, 5, 6], [-1, -1, -1, -1]])
    )

    links = sgrid.active_west_links((1, 4))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]))

    links = sgrid.active_west_links((4, 1))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]).T)


def test_active_east_links():
    links = sgrid.active_east_links((3, 4))
    assert_array_equal(
        links, np.array([[-1, -1, -1, -1], [4, 5, 6, -1], [-1, -1, -1, -1]])
    )

    links = sgrid.active_east_links((1, 4))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]))

    links = sgrid.active_east_links((4, 1))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]).T)

    links = sgrid.horizontal_active_link_ids((4, 5))
    assert_array_equal(np.array([[9, 10, 11, 12], [13, 14, 15, 16]]), links)

    links = sgrid.active_east_links((4, 5))
    assert_array_equal(
        np.array(
            [
                [-1, -1, -1, -1, -1],
                [9, 10, 11, 12, -1],
                [13, 14, 15, 16, -1],
                [-1, -1, -1, -1, -1],
            ]
        ),
        links,
    )


def test_active_north_links():
    links = sgrid.active_north_links((3, 4))
    assert_array_equal(
        links, np.array([[-1, 0, 1, -1], [-1, 2, 3, -1], [-1, -1, -1, -1]])
    )

    links = sgrid.active_north_links((1, 4))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]))

    links = sgrid.active_north_links((4, 1))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]).T)


def test_active_south_links():
    links = sgrid.active_south_links((3, 4))
    assert_array_equal(
        links, np.array([[-1, -1, -1, -1], [-1, 0, 1, -1], [-1, 2, 3, -1]])
    )

    links = sgrid.active_south_links((1, 4))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]))

    links = sgrid.active_south_links((4, 1))
    assert_array_equal(links, np.array([[-1, -1, -1, -1]]).T)


def test_active_inlinks():
    links = sgrid.active_inlinks((3, 4))
    assert_array_equal(
        np.array(
            [
                [-1, -1, -1, -1, -1, 0, 1, -1, -1, 2, 3, -1],
                [-1, -1, -1, -1, -1, 4, 5, 6, -1, -1, -1, -1],
            ]
        ),
        links,
    )


def test_active_outlinks():
    links = sgrid.active_outlinks((3, 4))
    assert_array_equal(
        np.array(
            [
                [-1, 0, 1, -1, -1, 2, 3, -1, -1, -1, -1, -1],
                [-1, -1, -1, -1, 4, 5, 6, -1, -1, -1, -1, -1],
            ]
        ),
        links,
    )


def test_active_outlinks_4x5():
    links = sgrid.active_outlinks((4, 5))

    assert_array_equal(
        np.array(
            [
                [-1, 0, 1, 2, -1, -1, 3, 4, 5, -1, -1, 6, 7, 8, -1, -1, -1, -1, -1, -1],
                [
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    9,
                    10,
                    11,
                    12,
                    -1,
                    13,
                    14,
                    15,
                    16,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                ],
            ]
        ),
        links,
    )


def test_active_inlinks_4x5():
    links = sgrid.active_inlinks((4, 5))

    assert_array_equal(
        np.array(
            [
                [-1, -1, -1, -1, -1, -1, 0, 1, 2, -1, -1, 3, 4, 5, -1, -1, 6, 7, 8, -1],
                [
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    9,
                    10,
                    11,
                    12,
                    -1,
                    13,
                    14,
                    15,
                    16,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                ],
            ]
        ),
        links,
    )


def test_face_count():
    assert 17 == sgrid.face_count((4, 5))
    assert 4 == sgrid.face_count((3, 3))
    assert 0 == sgrid.face_count((2, 100))
    assert 0 == sgrid.face_count((100, 2))
    assert 0 == sgrid.face_count((100, 1))


def test_active_face_count():
    assert 17 == sgrid.active_face_count((4, 5))
    assert 4 == sgrid.active_face_count((3, 3))
    assert 0 == sgrid.active_face_count((2, 100))
    assert 0 == sgrid.active_face_count((100, 2))
    assert 0 == sgrid.active_face_count((100, 1))


def test_active_faces():
    active_faces = sgrid.active_face_index((4, 5))
    assert_array_equal(
        np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]),
        active_faces,
    )


def test_link_faces_4_by_5():
    link_faces = sgrid.face_at_link((4, 5))

    BAD = sgrid.BAD_INDEX_VALUE

    assert_array_equal(
        link_faces,
        np.array(
            [
                BAD,
                0,
                1,
                2,
                BAD,
                BAD,
                3,
                4,
                5,
                BAD,
                BAD,
                6,
                7,
                8,
                BAD,
                BAD,
                BAD,
                BAD,
                BAD,
                9,
                10,
                11,
                12,
                13,
                14,
                15,
                16,
                BAD,
                BAD,
                BAD,
                BAD,
            ]
        ),
    )


def test_with_active_links():
    active_links = sgrid.active_links((4, 5))
    active_links = active_links[:-1]
    link_faces = sgrid.face_at_link((4, 5), actives=active_links)

    BAD = sgrid.BAD_INDEX_VALUE

    assert_array_equal(
        link_faces,
        np.array(
            [
                BAD,
                0,
                1,
                2,
                BAD,
                BAD,
                3,
                4,
                5,
                BAD,
                BAD,
                6,
                7,
                8,
                BAD,
                BAD,
                BAD,
                BAD,
                BAD,
                9,
                10,
                11,
                12,
                13,
                14,
                15,
                BAD,
                BAD,
                BAD,
                BAD,
                BAD,
            ]
        ),
    )


def test_reshape_array_default():
    x = np.arange(12.0)
    y = sgrid.reshape_array((3, 4), x)

    assert y.shape == (3, 4)
    assert_array_equal(x, y.flat)
    assert y.flags["C_CONTIGUOUS"]
    assert y.base is x


def test_copy():
    x = np.arange(12.0)
    y = sgrid.reshape_array((3, 4), x, copy=True)

    assert y.shape == (3, 4)
    assert_array_equal(x, y.flat)
    assert y.flags["C_CONTIGUOUS"]
    assert y.base is None

    y[0][0] = 0.0
    assert_array_equal(
        x, np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0])
    )


def test_flip():
    x = np.arange(12.0)
    y = sgrid.reshape_array((3, 4), x, flip_vertically=True)

    assert y.shape == (3, 4)
    assert_array_equal(
        y,
        np.array([[8.0, 9.0, 10.0, 11.0], [4.0, 5.0, 6.0, 7.0], [0.0, 1.0, 2.0, 3.0]]),
    )
    assert not y.flags["C_CONTIGUOUS"]
    assert y.base is not None

    y[0][0] = 0.0
    assert_array_equal(
        x, np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 0.0, 9.0, 10.0, 11.0])
    )


def test_flip_copy():
    x = np.arange(12.0)
    y = sgrid.reshape_array((3, 4), x, flip_vertically=True, copy=True)

    assert y.shape == (3, 4)
    assert_array_equal(
        y,
        np.array([[8.0, 9.0, 10.0, 11.0], [4.0, 5.0, 6.0, 7.0], [0.0, 1.0, 2.0, 3.0]]),
    )
    assert y.flags["C_CONTIGUOUS"]
    assert y.base is not x


def test_diagonal_array_default():
    diags = sgrid.diagonal_node_array((2, 3), out_of_bounds=-1)
    assert_array_equal(
        diags,
        np.array(
            [
                [4, -1, -1, -1],
                [5, 3, -1, -1],
                [-1, 4, -1, -1],
                [-1, -1, -1, 1],
                [-1, -1, 0, 2],
                [-1, -1, 1, -1],
            ]
        ),
    )

    assert diags.base is None
    assert diags.flags["C_CONTIGUOUS"]


def test_non_contiguous():
    diags = sgrid.diagonal_node_array((2, 3), out_of_bounds=-1, contiguous=False)
    assert_array_equal(
        diags,
        np.array(
            [
                [4, -1, -1, -1],
                [5, 3, -1, -1],
                [-1, 4, -1, -1],
                [-1, -1, -1, 1],
                [-1, -1, 0, 2],
                [-1, -1, 1, -1],
            ]
        ),
    )

    assert isinstance(diags.base, np.ndarray)
    assert not diags.flags["C_CONTIGUOUS"]


def test_boundary_node_mask_no_actives():
    diags = sgrid.diagonal_node_array((2, 3), out_of_bounds=-1, boundary_node_mask=-2)
    assert_array_equal(diags, -2 * np.ones((6, 4)))


def test_boundary_node_mask():
    diags = sgrid.diagonal_node_array((3, 3), out_of_bounds=-1, boundary_node_mask=-2)
    assert_array_equal(
        diags,
        np.array(
            [
                [-2, -2, -2, -2],
                [-2, -2, -2, -2],
                [-2, -2, -2, -2],
                [-2, -2, -2, -2],
                [8, 6, 0, 2],
                [-2, -2, -2, -2],
                [-2, -2, -2, -2],
                [-2, -2, -2, -2],
                [-2, -2, -2, -2],
            ]
        ),
    )


def test_neighbor_array_default():
    neighbors = sgrid.neighbor_node_array((2, 3))

    BAD = sgrid.BAD_INDEX_VALUE

    assert_array_equal(
        neighbors,
        np.array(
            [
                [1, 3, BAD, BAD],
                [2, 4, 0, BAD],
                [BAD, 5, 1, BAD],
                [4, BAD, BAD, 0],
                [5, BAD, 3, 1],
                [BAD, BAD, 4, 2],
            ]
        ).T,
    )

    assert neighbors.flags["C_CONTIGUOUS"]
    assert neighbors.base is None


def test_set_out_of_bounds():
    neighbors = sgrid.neighbor_node_array((2, 3), inactive=-1)
    assert_array_equal(
        neighbors,
        np.array(
            [
                [1, 3, -1, -1],
                [2, 4, 0, -1],
                [-1, 5, 1, -1],
                [4, -1, -1, 0],
                [5, -1, 3, 1],
                [-1, -1, 4, 2],
            ]
        ).T,
    )


def test_no_inactive():
    inlinks = sgrid.setup_active_inlink_matrix((4, 5), return_count=False)
    assert_array_equal(
        inlinks,
        np.array(
            [
                [-1, -1, -1, -1, -1, -1, 0, 1, 2, -1, -1, 3, 4, 5, -1, -1, 6, 7, 8, -1],
                [
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    9,
                    10,
                    11,
                    12,
                    -1,
                    13,
                    14,
                    15,
                    16,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                ],
            ]
        ),
    )


def test_inactive():
    status = np.ones((4, 5))
    status[1, 1] = 0
    inlinks = sgrid.setup_active_inlink_matrix(
        (4, 5), return_count=False, node_status=status
    )
    assert_array_equal(
        inlinks,
        np.array(
            [
                [
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    0,
                    1,
                    -1,
                    -1,
                    -1,
                    2,
                    3,
                    -1,
                    -1,
                    4,
                    5,
                    6,
                    -1,
                ],
                [
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    7,
                    8,
                    -1,
                    9,
                    10,
                    11,
                    12,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                ],
            ]
        ),
    )


def test_out_link_ids_at_nodes():
    links_ids = sgrid.outlink_index_at_node((4, 5))
    assert_array_equal(
        np.array(
            [
                [
                    4,
                    5,
                    6,
                    7,
                    8,
                    13,
                    14,
                    15,
                    16,
                    17,
                    22,
                    23,
                    24,
                    25,
                    26,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                ],
                [
                    0,
                    1,
                    2,
                    3,
                    -1,
                    9,
                    10,
                    11,
                    12,
                    -1,
                    18,
                    19,
                    20,
                    21,
                    -1,
                    27,
                    28,
                    29,
                    30,
                    -1,
                ],
            ]
        ),
        links_ids,
    )


def test_in_link_ids_at_nodes():
    links_ids = sgrid.inlink_index_at_node((4, 5))
    assert_array_equal(
        np.array(
            [
                [
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    4,
                    5,
                    6,
                    7,
                    8,
                    13,
                    14,
                    15,
                    16,
                    17,
                    22,
                    23,
                    24,
                    25,
                    26,
                ],
                [
                    -1,
                    0,
                    1,
                    2,
                    3,
                    -1,
                    9,
                    10,
                    11,
                    12,
                    -1,
                    18,
                    19,
                    20,
                    21,
                    -1,
                    27,
                    28,
                    29,
                    30,
                ],
            ]
        ),
        links_ids,
    )
