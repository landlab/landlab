import numpy as np
import pytest
from numpy.testing import assert_equal

from landlab.graph._links_at_node import _get_links_at_node


@pytest.mark.parametrize("dtype", (np.int32, np.int64))
def test_links_at_node_2d(dtype):
    # 3 - 4 - 5
    # |   |   |
    # 0 - 1 - 2
    nodes_at_link = [
        [0, 1],
        [1, 2],
        [0, 3],
        [1, 4],
        [2, 5],
        [3, 4],
        [4, 5],
    ]
    links_at_node = np.full((6, 4), -2, dtype=dtype)
    link_dirs_at_node = np.full((6, 4), -2, dtype=np.int8)

    _get_links_at_node(
        np.asarray(nodes_at_link, dtype=dtype), links_at_node, link_dirs_at_node
    )

    sorted_by_link = np.argsort(links_at_node, axis=1)

    assert_equal(
        np.take_along_axis(links_at_node, sorted_by_link),
        [
            [-2, -2, 0, 2],
            [-2, 0, 1, 3],
            [-2, -2, 1, 4],
            [-2, -2, 2, 5],
            [-2, 3, 5, 6],
            [-2, -2, 4, 6],
        ],
    )
    assert_equal(
        np.take_along_axis(link_dirs_at_node, sorted_by_link),
        [
            [-2, -2, -1, -1],
            [-2, 1, -1, -1],
            [-2, -2, 1, -1],
            [-2, -2, 1, -1],
            [-2, 1, 1, -1],
            [-2, -2, 1, 1],
        ],
    )


@pytest.mark.parametrize("dtype", (np.int32, np.int64))
def test_links_at_node_1d(dtype):
    # 0 - 1 - 2 - 3
    nodes_at_link = [[0, 1], [1, 2], [2, 3]]

    links_at_node = np.full((4, 2), -2, dtype=dtype)
    link_dirs_at_node = np.full((4, 2), -2, dtype=np.int8)

    _get_links_at_node(
        np.asarray(nodes_at_link, dtype=dtype), links_at_node, link_dirs_at_node
    )

    sorted_by_link = np.argsort(links_at_node, axis=1)
    assert_equal(
        np.take_along_axis(links_at_node, sorted_by_link),
        [[-2, 0], [0, 1], [1, 2], [-2, 2]],
    )
    assert_equal(
        np.take_along_axis(link_dirs_at_node, sorted_by_link),
        [[-2, -1], [1, -1], [1, -1], [-2, 1]],
    )


@pytest.mark.parametrize("dtype", (np.int32, np.int64))
def test_links_at_node_missing_node(dtype):
    # 0 - 1 - 3 - 4
    nodes_at_link = [[0, 1], [1, 3], [3, 4]]

    links_at_node = np.full((5, 2), -2, dtype=dtype)
    link_dirs_at_node = np.full((5, 2), -2, dtype=np.int8)

    _get_links_at_node(
        np.asarray(nodes_at_link, dtype=dtype), links_at_node, link_dirs_at_node
    )

    sorted_by_link = np.argsort(links_at_node, axis=1)
    assert_equal(
        np.take_along_axis(links_at_node, sorted_by_link),
        [[-2, 0], [0, 1], [-2, -2], [1, 2], [-2, 2]],
    )
    assert_equal(
        np.take_along_axis(link_dirs_at_node, sorted_by_link),
        [[-2, -1], [1, -1], [-2, -2], [1, -1], [-2, 1]],
    )


@pytest.mark.parametrize("dtype", (np.int32, np.int64))
def test_links_at_node_disconnected_node(dtype):
    # 0 ... 1 - 2 - 3 - 4
    nodes_at_link = [[1, 2], [2, 3], [3, 4]]

    links_at_node = np.full((5, 2), -2, dtype=dtype)
    link_dirs_at_node = np.full((5, 2), -2, dtype=np.int8)

    _get_links_at_node(
        np.asarray(nodes_at_link, dtype=dtype), links_at_node, link_dirs_at_node
    )

    sorted_by_link = np.argsort(links_at_node, axis=1)
    assert_equal(
        np.take_along_axis(links_at_node, sorted_by_link),
        [[-2, -2], [-2, 0], [0, 1], [1, 2], [-2, 2]],
    )
    assert_equal(
        np.take_along_axis(link_dirs_at_node, sorted_by_link),
        [[-2, -2], [-2, -1], [1, -1], [1, -1], [-2, 1]],
    )
