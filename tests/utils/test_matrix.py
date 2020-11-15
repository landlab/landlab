import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import HexModelGrid, RasterModelGrid
from landlab.utils import get_core_node_matrix


@pytest.mark.parametrize("diff", (-1, 1))
def test_bad_value_at_node(diff):
    grid = RasterModelGrid((5, 6))
    value_at_node = np.full(grid.number_of_nodes + diff, 1.0)

    with pytest.raises(ValueError):
        get_core_node_matrix(grid, value_at_node)


@pytest.mark.parametrize("diff", (-1, 1))
def test_bad_coef_at_link(diff):
    grid = RasterModelGrid((5, 6))
    value_at_node = np.full(grid.number_of_nodes, 1.0)
    coef_at_link = np.full(grid.number_of_links + diff, 1.0)

    with pytest.raises(ValueError):
        get_core_node_matrix(grid, value_at_node, coef_at_link=coef_at_link)


def test_default_bc():
    grid = RasterModelGrid((5, 6))
    value_at_node = np.full(grid.number_of_nodes, 1.0, dtype=float)

    mat, rhs = get_core_node_matrix(grid, value_at_node)

    assert_array_equal(
        mat.toarray(),
        [
            [-4.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [1.0, -4.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, -4.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, -4.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, -4.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 1.0, -4.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, 1.0, -4.0, 1.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, -4.0, 0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -4.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, -4.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, -4.0, 1.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, -4.0],
        ],
    )
    assert_array_equal(
        rhs,
        [
            [-2.0],
            [-1.0],
            [-1.0],
            [-2.0],
            [-1.0],
            [0.0],
            [0.0],
            [-1.0],
            [-2.0],
            [-1.0],
            [-1.0],
            [-2.0],
        ],
    )


def test_small_grids():
    grid = RasterModelGrid((3, 3))
    value_at_node = np.full(grid.number_of_nodes, 1.0)
    mat, rhs = get_core_node_matrix(grid, value_at_node)
    assert_array_equal(mat.toarray(), [[-4.0]])
    assert_array_equal(rhs, [[-4.0]])

    grid = RasterModelGrid((4, 4))
    value_at_node = np.full(grid.number_of_nodes, 1.0)
    mat, rhs = get_core_node_matrix(grid, value_at_node)

    assert_array_equal(
        mat.toarray(),
        [
            [-4.0, 1.0, 1.0, 0.0],
            [1.0, -4.0, 0.0, 1.0],
            [1.0, 0.0, -4.0, 1.0],
            [0.0, 1.0, 1.0, -4.0],
        ],
    )
    assert_array_equal(rhs, [[-2.0], [-2.0], [-2.0], [-2.0]])


def test_with_fixed_value_bc():
    grid = RasterModelGrid((4, 5))
    grid.status_at_node[13] = grid.BC_NODE_IS_FIXED_VALUE
    grid.status_at_node[2] = grid.BC_NODE_IS_CLOSED

    value_at_node = np.arange(grid.number_of_nodes, dtype=float)

    mat, rhs = get_core_node_matrix(grid, value_at_node)
    assert_array_equal(
        mat.toarray(),
        [
            [-4.0, 1.0, 0.0, 1.0, 0.0],
            [1.0, -3.0, 1.0, 0.0, 1.0],
            [0.0, 1.0, -4.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, -4.0, 1.0],
            [0.0, 1.0, 0.0, 1.0, -4.0],
        ],
    )
    assert_array_equal(
        rhs, [[-6.0], [0.0], [-25.0], [-26.0], [-30.0]],
    )


def test_with_coef():
    grid = RasterModelGrid((4, 5))
    grid.status_at_node[13] = grid.BC_NODE_IS_FIXED_VALUE
    grid.status_at_node[2] = grid.BC_NODE_IS_CLOSED

    value_at_node = np.arange(grid.number_of_nodes, dtype=float)

    coef_at_link = np.arange(grid.number_of_links, dtype=np.double)
    mat, rhs = get_core_node_matrix(grid, value_at_node, coef_at_link=coef_at_link)
    assert_array_equal(
        mat.toarray(),
        [
            [-38.0, 10.0, 0.0, 14.0, 0.0],
            [10.0, -36.0, 11.0, 0.0, 15.0],
            [0.0, 11.0, -46.0, 0.0, 0.0],
            [14.0, 0.0, 0.0, -74.0, 19.0],
            [0.0, 15.0, 0.0, 19.0, -78.0],
        ],
    )
    assert_array_equal(rhs, [[-6.0], [0.0], [-25.0], [-26.0], [-30.0]])


@pytest.mark.parametrize("Grid", (RasterModelGrid, HexModelGrid))
def test_scalar_value_at_node(Grid):
    grid = Grid((4, 5))

    _, expected = get_core_node_matrix(grid, np.full(grid.number_of_nodes, 11.0))
    _, actual = get_core_node_matrix(grid, 11.0)

    assert_array_equal(expected, actual)


@pytest.mark.parametrize("Grid", (RasterModelGrid, HexModelGrid))
def test_scalar_coef_at_link(Grid):
    grid = Grid((4, 5))

    expected, _ = get_core_node_matrix(
        grid, 1.0, coef_at_link=np.full(grid.number_of_links, 11.0)
    )
    actual, _ = get_core_node_matrix(grid, 1.0, coef_at_link=11.0)

    assert_array_equal(expected.toarray(), actual.toarray())


@pytest.mark.parametrize("Grid", (RasterModelGrid, HexModelGrid))
def test_default_coef_is_one(Grid):
    grid = Grid((4, 5))
    grid.status_at_node[13] = grid.BC_NODE_IS_FIXED_VALUE
    grid.status_at_node[2] = grid.BC_NODE_IS_CLOSED

    value_at_node = np.arange(grid.number_of_nodes, dtype=float)

    coef_at_link = np.full(grid.number_of_links, 1.0, dtype=float)

    expected, _ = get_core_node_matrix(grid, value_at_node, coef_at_link=coef_at_link)
    actual, _ = get_core_node_matrix(grid, value_at_node)

    assert_array_equal(expected.toarray(), actual.toarray())


@pytest.mark.parametrize("Grid", (RasterModelGrid, HexModelGrid))
def test_value_at_node(Grid):
    grid = Grid((4, 5))

    value_at_node = np.full(grid.number_of_nodes, 1.0)
    coef_at_link = np.full(grid.number_of_links, 1.0)

    mat, rhs = get_core_node_matrix(grid, value_at_node, coef_at_link=coef_at_link)
    actual, actual_rhs = get_core_node_matrix(grid, value_at_node * 2, coef_at_link)

    assert_array_equal(mat.toarray(), actual.toarray())
    assert_array_equal(rhs, actual_rhs / 2.0)


@pytest.mark.parametrize("Grid", (RasterModelGrid, HexModelGrid))
def test_coef_at_link(Grid):
    grid = Grid((4, 5))

    value_at_node = np.full(grid.number_of_nodes, 1.0)
    coef_at_link = np.full(grid.number_of_links, 1.0)

    mat, rhs = get_core_node_matrix(grid, value_at_node, coef_at_link=coef_at_link)
    actual, actual_rhs = get_core_node_matrix(grid, value_at_node, coef_at_link * 2.0)

    assert_array_equal(mat.toarray(), actual.toarray() / 2.0)
    assert_array_equal(rhs, actual_rhs)


def test_hex_grid():
    grid = HexModelGrid((4, 3))
    value_at_node = np.full(grid.number_of_nodes, 1.0)
    coef_at_link = np.full(grid.number_of_links, 1.0)

    mat, rhs = get_core_node_matrix(grid, value_at_node, coef_at_link=coef_at_link)

    assert_array_equal(
        mat.toarray(),
        [
            [-6.0, 1.0, 1.0, 1.0, 0.0],
            [1.0, -6.0, 0.0, 1.0, 1.0],
            [1.0, 0.0, -6.0, 1.0, 0.0],
            [1.0, 1.0, 1.0, -6.0, 1.0],
            [0.0, 1.0, 0.0, 1.0, -6.0],
        ],
    )
    assert_array_equal(rhs, [[-3.0], [-3.0], [-4.0], [-2.0], [-4.0]])
