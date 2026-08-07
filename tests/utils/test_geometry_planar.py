import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.utils.geometry.planar import NearestNodeFinder
from landlab.utils.geometry.planar import calc_cumulative_path_length
from landlab.utils.geometry.planar import calc_path_length
from landlab.utils.geometry.planar import calc_path_segment_lengths
from landlab.utils.geometry.planar import find_nearest_node


def test_nearest_node_with_bad_shape():
    with pytest.raises(ValueError):
        NearestNodeFinder([0.0, 1.0, 2.0])

    with pytest.raises(ValueError):
        NearestNodeFinder(np.zeros((2, 2, 2)))


def test_nearest_node_with_one_target():
    find_nearest = NearestNodeFinder([[0.0, 0.0], [2.0, 2.0], [5.0, 5.0]])
    node = find_nearest([1.9, 2.1])
    assert isinstance(node, (int, np.integer))
    assert node == 1


def test_nearest_node_with_multiple_targets():
    find_nodes = NearestNodeFinder([[0.0, 0.0], [1.0, 1.0], [5.0, 5.0]])
    nodes = find_nodes([[0.2, 0.2], [3.9, 4.1], [0.8, 0.7]])
    assert isinstance(nodes, np.ndarray)
    assert nodes.shape == (3,)
    assert np.array_equal(nodes, [0, 2, 1])


def test_nearest_node_agrees_with_bruteforce_argmin():
    rng = np.random.default_rng(42)
    coords_of_node = rng.random((1000, 2))
    coords_of_target = rng.random((20, 2))
    find_nearest = NearestNodeFinder(coords_of_node)

    actual = find_nearest(coords_of_target)

    expected = np.argmin(
        np.sum(
            (coords_of_node[None, :, :] - coords_of_target[:, None, :]) ** 2, axis=2
        ),
        axis=1,
    )

    assert np.array_equal(actual, expected)


def test_nearest_node_with_a_tie():
    coords_of_node = np.array([[0.0, 0.0], [0.0, 1.0], [2.0, 2.0]])
    coords_of_target = np.array([0.0, 0.5])

    tie_nodes = {0, 1}

    find_nearest = NearestNodeFinder(coords_of_node)
    node = find_nearest(coords_of_target)

    assert node in tie_nodes

    distance_to_node = np.sum((coords_of_node - coords_of_target) ** 2, axis=1)
    assert np.isclose(distance_to_node[node], distance_to_node.min())


def test_function_matches_class():
    coords_of_node = [[0.0, 0.0], [1.0, 1.0], [5.0, 5.0]]
    coords_of_target = [[0.2, 0.2], [3.9, 4.1], [0.8, 0.7]]
    expected = NearestNodeFinder(coords_of_node)(coords_of_target)
    actual = find_nearest_node(coords_of_node, coords_of_target)

    assert np.array_equal(actual, expected)


@pytest.mark.parametrize(
    "points",
    (
        [[0.0, 0.0], [1.0, 0.0]],
        [[0.0, 0.0], [0.0, -1.0]],
        [[1.0, 2.0], [1.5, 2.0], [1.5, 2.5]],
        [[1.0, 0.0, 2.0], [1.5, 0.0, 2.0], [1.5, 0.0, 2.5]],
    ),
)
def test_path_length(points):
    assert calc_path_length(points) == 1.0
    assert_array_equal(calc_cumulative_path_length(points)[[0, -1]], (0.0, 1.0))
    assert len(calc_cumulative_path_length(points)) == len(points)


def test_path_with_one_point():
    points = [[1.0, 2.0]]
    assert_array_equal(calc_path_segment_lengths(points), [])
    assert_array_equal(calc_cumulative_path_length(points), [0.0])
    assert calc_path_length(points) == 0.0


@pytest.mark.parametrize("points", (np.empty((0, 2)), np.empty((0, 3)), []))
def test_path_with_no_points(points):
    assert_array_equal(calc_path_segment_lengths(points), [])
    assert_array_equal(calc_cumulative_path_length(points), [])
    assert calc_path_length(points) == 0.0


def test_cumulative_path_length():
    assert_array_equal(
        calc_cumulative_path_length([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [2.0, 1.0]]),
        [0.0, 1.0, 2.0, 3.0],
    )


def test_path_segment_lengths():
    assert_array_equal(
        calc_path_segment_lengths([[0.0, 0.0], [1.0, 0.0], [1.0, 2.0], [4.0, 2.0]]),
        [1.0, 2.0, 3.0],
    )


@pytest.mark.parametrize("ndim", (1, 2, 3, 4))
def test_wrong_dimensionality(ndim):
    points = np.zeros((1,) * ndim)
    if ndim == 2:
        assert_array_equal(calc_path_segment_lengths(points), [])
    else:
        with pytest.raises(ValueError) as excinfo:
            calc_path_segment_lengths(points)
        assert str(excinfo.value).startswith("Expected a 2D array")
