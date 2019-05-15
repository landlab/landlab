import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import BAD_INDEX_VALUE as X, RasterModelGrid


@pytest.mark.parametrize("params", [[(3, 4), {}], [(), {"num_rows": 3, "num_cols": 4}]])
def test_parse_parameters_deprecated_shape(params):
    with pytest.deprecated_call():
        args, kwds = RasterModelGrid._parse_parameters(*params)
    assert args == ((3, 4),)
    assert kwds == {}


@pytest.mark.parametrize(
    "params",
    [[((3, 4),), {"dx": 20.0}], [((3, 4),), {"spacing": 20.0}], [(3, 4, 20.0), {}]],
)
def test_parse_parameters_deprecated_spacing(params):
    with pytest.deprecated_call():
        args, kwds = RasterModelGrid._parse_parameters(*params)
    assert args == ((3, 4),)
    assert kwds == {"xy_spacing": 20.0}


def test_parse_parameters_deprecated_xy_of_lower_left():
    with pytest.deprecated_call():
        args, kwds = RasterModelGrid._parse_parameters([(3, 4)], {"origin": (1.0, 2.0)})
    assert args == ((3, 4),)
    assert kwds == {"xy_of_lower_left": (1.0, 2.0)}


@pytest.mark.parametrize(
    "kwds",
    [
        dict(xy_of_lower_left=(3.0, 4.0)),
        dict(xy_spacing=(30.0, 40.0)),
        dict(xy_spacing=(30.0, 40.0), xy_of_lower_left=(3.0, 4.0)),
        dict(xy_spacing=(30.0, 40.0), xy_of_lower_left=(3.0, 4.0), another_kwd=True),
    ],
)
def test_parse_parameters_new_style(kwds):
    args, new_kwds = RasterModelGrid._parse_parameters([(3, 4)], kwds)
    assert args == ((3, 4),)
    assert new_kwds == kwds


def test_init_with_kwds_classic():

    with pytest.deprecated_call():
        grid = RasterModelGrid(num_rows=4, num_cols=5, xy_spacing=1.0)

    assert grid.number_of_node_rows == 4
    assert grid.number_of_node_columns == 5
    assert grid.dy == 1
    assert grid.dx == 1

    with pytest.deprecated_call():
        grid = RasterModelGrid(3, 7, 2)

    assert grid.number_of_node_rows == 3
    assert grid.number_of_node_columns == 7
    assert grid.dy == 2.0
    assert grid.dx == 2.0


def test_init_new_style():
    grid = RasterModelGrid((4, 5), xy_spacing=2)

    assert grid.number_of_node_rows == 4
    assert grid.number_of_node_columns == 5
    assert grid.dy == 2.0
    assert grid.dx == 2.0

    grid = RasterModelGrid((4, 5))

    assert grid.number_of_node_rows == 4
    assert grid.number_of_node_columns == 5
    assert grid.dy == 1.0
    assert grid.dx == 1.0


def test_spacing_is_float():
    grid = RasterModelGrid((4, 5))
    assert grid.dy == 1.0
    assert isinstance(grid.dy, float)
    assert grid.dx == 1.0
    assert isinstance(grid.dx, float)

    grid = RasterModelGrid((4, 5), xy_spacing=2)
    assert grid.dy == 2.0
    assert isinstance(grid.dy, float)
    assert grid.dx == 2.0
    assert isinstance(grid.dx, float)


def test_grid_dimensions():
    """Test extent of grid with unit spacing."""
    rmg = RasterModelGrid((4, 5))
    assert rmg.extent[0] == rmg.number_of_node_rows - 1
    assert rmg.extent[1] == rmg.number_of_node_columns - 1


def test_grid_dimensions_non_unit_spacing():
    """Test extent of grid with non-unit spacing."""
    rmg = RasterModelGrid((4, 5), xy_spacing=2.0)
    assert rmg.extent[0] == 6.0
    assert rmg.extent[1] == 8.0


def test_nodes_around_point():
    rmg = RasterModelGrid((4, 5))
    surrounding_ids = rmg.nodes_around_point(2.1, 1.1)
    assert_array_equal(surrounding_ids, np.array([7, 12, 13, 8]))

    surrounding_ids = rmg.nodes_around_point(2.1, 0.9)
    assert_array_equal(surrounding_ids, np.array([2, 7, 8, 3]))


def test_neighbor_list_with_scalar_arg():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(rmg.active_adjacent_nodes_at_node[6], np.array([7, 11, 5, 1]))
    assert_array_equal(rmg.active_adjacent_nodes_at_node[-1], np.array([X, X, X, X]))
    assert_array_equal(rmg.active_adjacent_nodes_at_node[-2], np.array([X, X, X, 13]))


def test_neighbor_list_with_array_arg():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(
        rmg.active_adjacent_nodes_at_node[[6, -1]],
        np.array([[7, 11, 5, 1], [X, X, X, X]]),
    )


def test_neighbor_list_with_no_args():
    rmg = RasterModelGrid((4, 5))
    expected = np.array(
        [
            [X, X, X, X],
            [X, 6, X, X],
            [X, 7, X, X],
            [X, 8, X, X],
            [X, X, X, X],
            [6, X, X, X],
            [7, 11, 5, 1],
            [8, 12, 6, 2],
            [9, 13, 7, 3],
            [X, X, 8, X],
            [11, X, X, X],
            [12, 16, 10, 6],
            [13, 17, 11, 7],
            [14, 18, 12, 8],
            [X, X, 13, X],
            [X, X, X, X],
            [X, X, X, 11],
            [X, X, X, 12],
            [X, X, X, 13],
            [X, X, X, X],
        ]
    )

    assert_array_equal(rmg.active_adjacent_nodes_at_node, expected)


def test_node_x():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(
        rmg.node_x,
        np.array(
            [
                0.0,
                1.0,
                2.0,
                3.0,
                4.0,
                0.0,
                1.0,
                2.0,
                3.0,
                4.0,
                0.0,
                1.0,
                2.0,
                3.0,
                4.0,
                0.0,
                1.0,
                2.0,
                3.0,
                4.0,
            ]
        ),
    )


def test_node_y():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(
        rmg.node_y,
        np.array(
            [
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                2.0,
                2.0,
                2.0,
                2.0,
                2.0,
                3.0,
                3.0,
                3.0,
                3.0,
                3.0,
            ]
        ),
    )


def test_node_x_is_immutable():
    rmg = RasterModelGrid((4, 5))
    with pytest.raises(ValueError):
        rmg.node_x[0] = 0


def test_node_y_is_immutable():
    rmg = RasterModelGrid((4, 5))
    with pytest.raises(ValueError):
        rmg.node_y[0] = 0


def test_node_axis_coordinates():
    rmg = RasterModelGrid((4, 5))
    assert rmg.node_axis_coordinates(axis=0).base is rmg.node_y.base
    assert rmg.node_axis_coordinates(axis=1).base is rmg.node_x.base
    assert rmg.node_axis_coordinates(axis=-1).base is rmg.node_x.base
    assert rmg.node_axis_coordinates(axis=-2).base is rmg.node_y.base


def test_diagonal_list():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(rmg.diagonal_adjacent_nodes_at_node[6], np.array([12, 10, 0, 2]))
    assert_array_equal(rmg.diagonal_adjacent_nodes_at_node[-1], np.array([X, X, 13, X]))
    assert_array_equal(
        rmg.diagonal_adjacent_nodes_at_node[[6, -1]],
        np.array([[12, 10, 0, 2], [X, X, 13, X]]),
    )
    assert_array_equal(
        rmg.diagonal_adjacent_nodes_at_node,
        np.array(
            [
                [6, X, X, X],
                [7, 5, X, X],
                [8, 6, X, X],
                [9, 7, X, X],
                [X, 8, X, X],
                [11, X, X, 1],
                [12, 10, 0, 2],
                [13, 11, 1, 3],
                [14, 12, 2, 4],
                [X, 13, 3, X],
                [16, X, X, 6],
                [17, 15, 5, 7],
                [18, 16, 6, 8],
                [19, 17, 7, 9],
                [X, 18, 8, X],
                [X, X, X, 11],
                [X, X, 10, 12],
                [X, X, 11, 13],
                [X, X, 12, 14],
                [X, X, 13, X],
            ]
        ),
    )


def test_diagonal_list_boundary():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(rmg.diagonal_adjacent_nodes_at_node[0], np.array([6, X, X, X]))


def test_node_is_core():
    rmg = RasterModelGrid((4, 5))
    for cell_id in [0, 1, 2, 3, 4, 5, 9, 10, 14, 15, 16, 17, 18, 19]:
        assert not rmg.node_is_core(cell_id)
    for cell_id in [6, 7, 8, 11, 12, 13]:
        assert rmg.node_is_core(cell_id)


def test_get_interior_cells():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(rmg.node_at_core_cell, np.array([6, 7, 8, 11, 12, 13]))


def test_active_links():
    rmg = RasterModelGrid((4, 5))
    assert rmg.number_of_active_links == 17
    assert_array_equal(
        rmg.active_links,
        np.array([5, 6, 7, 9, 10, 11, 12, 14, 15, 16, 18, 19, 20, 21, 23, 24, 25]),
    )


# def test_active_link_fromnode():
#    rmg = RasterModelGrid((4, 5))
#    assert_array_equal(rmg._activelink_fromnode,
#                       np.array([1, 2, 3, 6, 7, 8, 11, 12, 13,
#                                 5, 6, 7, 8, 10, 11, 12, 13]))
#
#
# def test_active_link_tonode():
#    rmg = RasterModelGrid((4, 5))
#    assert_array_equal(rmg._activelink_tonode,
#                       np.array([6, 7, 8, 11, 12, 13, 16, 17, 18,
#                                 6, 7, 8, 9, 11, 12, 13, 14]))


def test__active_links_at_node_scalar_interior():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(rmg._active_links_at_node([6]), np.array([[5, 9, 14, 10]]).T)


def test__active_links_at_node_scalar_boundary():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(rmg._active_links_at_node([1]), np.array([[-1, -1, 5, -1]]).T)


def test_active_node_with_array_arg():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(
        rmg._active_links_at_node([6, 7]), np.array([[5, 9, 14, 10], [6, 10, 15, 11]]).T
    )


def test__active_links_at_node_with_no_args():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(
        rmg._active_links_at_node(),
        np.array(
            [
                [
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    5,
                    6,
                    7,
                    -1,
                    -1,
                    14,
                    15,
                    16,
                    -1,
                    -1,
                    23,
                    24,
                    25,
                    -1,
                ],
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
                    18,
                    19,
                    20,
                    21,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                ],
                [
                    -1,
                    5,
                    6,
                    7,
                    -1,
                    -1,
                    14,
                    15,
                    16,
                    -1,
                    -1,
                    23,
                    24,
                    25,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                ],
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
                    18,
                    19,
                    20,
                    21,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                ],
            ]
        ),
    )


def test_nodes_at_link():
    """Test nodes_at_link shares data with tail and head."""
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(rmg.nodes_at_link[:, 0], rmg.node_at_link_tail)
    assert_array_equal(rmg.nodes_at_link[:, 1], rmg.node_at_link_head)

    assert np.may_share_memory(rmg.nodes_at_link, rmg.node_at_link_tail)
    assert np.may_share_memory(rmg.nodes_at_link, rmg.node_at_link_head)


def test_node_at_link_tail():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(
        rmg.node_at_link_tail,
        np.array(
            [
                0,
                1,
                2,
                3,
                0,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                13,
                10,
                11,
                12,
                13,
                14,
                15,
                16,
                17,
                18,
            ]
        ),
    )


def test_node_at_link_head():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(
        rmg.node_at_link_head,
        np.array(
            [
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                13,
                14,
                11,
                12,
                13,
                14,
                15,
                16,
                17,
                18,
                19,
                16,
                17,
                18,
                19,
            ]
        ),
    )


def test_links_at_node_with_scalar_interior():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(rmg.links_at_node[6], np.array([10, 14, 9, 5]))


def test_links_at_node_with_scalar_boundary():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(rmg.links_at_node[1], np.array([1, 5, 0, -1]))


def test_links_at_node_with_array_arg():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(
        rmg.links_at_node[6:8], np.array([[10, 14, 9, 5], [11, 15, 10, 6]])
    )


def test_links_at_node_with_no_args():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(
        rmg.links_at_node,
        np.array(
            [
                [0, 4, -1, -1],
                [1, 5, 0, -1],
                [2, 6, 1, -1],
                [3, 7, 2, -1],
                [-1, 8, 3, -1],
                [9, 13, -1, 4],
                [10, 14, 9, 5],
                [11, 15, 10, 6],
                [12, 16, 11, 7],
                [-1, 17, 12, 8],
                [18, 22, -1, 13],
                [19, 23, 18, 14],
                [20, 24, 19, 15],
                [21, 25, 20, 16],
                [-1, 26, 21, 17],
                [27, -1, -1, 22],
                [28, -1, 27, 23],
                [29, -1, 28, 24],
                [30, -1, 29, 25],
                [-1, -1, 30, 26],
            ]
        ),
    )


def test_face_at_link():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(
        rmg.face_at_link,
        np.array(
            [
                X,
                X,
                X,
                X,
                X,
                0,
                1,
                2,
                X,
                3,
                4,
                5,
                6,
                X,
                7,
                8,
                9,
                X,
                10,
                11,
                12,
                13,
                X,
                14,
                15,
                16,
                X,
                X,
                X,
                X,
                X,
            ]
        ),
    )


def test_grid_coords_to_node_id_with_scalar():
    rmg = RasterModelGrid((4, 5))
    assert rmg.grid_coords_to_node_id(3, 4) == 19


def test_grid_coords_to_node_id_with_array():
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(rmg.grid_coords_to_node_id((3, 2), (4, 1)), np.array([19, 11]))


def test_grid_coords_to_node_id_outside_of_grid():
    rmg = RasterModelGrid((4, 5))
    with pytest.raises(ValueError):
        rmg.grid_coords_to_node_id(5, 0)


def test_diagonal_adjacent_nodes_at_node():
    """Test diagonally adjacent nodes."""
    rmg = RasterModelGrid((4, 5))
    assert_array_equal(
        rmg.diagonal_adjacent_nodes_at_node,
        np.array(
            [
                [6, X, X, X],
                [7, 5, X, X],
                [8, 6, X, X],
                [9, 7, X, X],
                [X, 8, X, X],
                [11, X, X, 1],
                [12, 10, 0, 2],
                [13, 11, 1, 3],
                [14, 12, 2, 4],
                [X, 13, 3, X],
                [16, X, X, 6],
                [17, 15, 5, 7],
                [18, 16, 6, 8],
                [19, 17, 7, 9],
                [X, 18, 8, X],
                [X, X, X, 11],
                [X, X, 10, 12],
                [X, X, 11, 13],
                [X, X, 12, 14],
                [X, X, 13, X],
            ]
        ),
    )
