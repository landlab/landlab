import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.grid.raster_funcs import neighbor_to_arrow


@pytest.mark.parametrize("ndim", (1, 2))
@pytest.mark.parametrize("diagonals_ok", (True, False))
def test_vector(ndim, diagonals_ok):
    array = np.array([0, 0, 1], ndmin=ndim)

    actual = neighbor_to_arrow(array, diagonals_ok=diagonals_ok)

    assert actual.ndim == ndim
    assert_array_equal(actual.squeeze(), ["○", "←", "←"])


@pytest.mark.parametrize("bad_neighbor", (2, -1))
def test_bad_neighbor(bad_neighbor):
    assert_array_equal(neighbor_to_arrow([bad_neighbor, 0, 1]), ["?", "←", "←"])


def test_too_many_dimensions():
    with pytest.raises(ValueError):
        neighbor_to_arrow([[[0, 1, 2]]])


@pytest.mark.parametrize("diagonals_ok", (True, False))
@pytest.mark.parametrize("direction", ("left", "right", "up", "down", "center"))
def test_d4_directions(diagonals_ok, direction):
    array = np.arange(12).reshape((3, 4))
    expected = np.full(array.shape, "○")

    if direction == "left":
        array[:, 1:] = array[:, :-1]
        expected[:, 1:] = "←"
    elif direction == "right":
        array[:, :-1] = array[:, 1:]
        expected[:, :-1] = "→"
    elif direction == "up":
        array[1:, :] = array[:-1, :]
        expected[1:, :] = "↑"
    elif direction == "down":
        array[:-1, :] = array[1:, :]
        expected[:-1, :] = "↓"

    assert_array_equal(neighbor_to_arrow(array, diagonals_ok=diagonals_ok), expected)


@pytest.mark.parametrize(
    "direction", ("up-right", "up-left", "down-left", "down-right", "center")
)
def test_d8_directions(direction):
    array = np.arange(12).reshape((3, 4))
    expected = np.full(array.shape, "○")

    if direction == "up-left":
        array[1:, 1:] = array[:-1, :-1]
        expected[1:, 1:] = "↖"
    elif direction == "up-right":
        array[1:, :-1] = array[:-1, 1:]
        expected[1:, :-1] = "↗"
    elif direction == "down-left":
        array[:-1, 1:] = array[1:, :-1]
        expected[:-1, 1:] = "↙"
    elif direction == "down-right":
        array[:-1, :-1] = array[1:, 1:]
        expected[:-1, :-1] = "↘"

    assert_array_equal(neighbor_to_arrow(array, diagonals_ok=True), expected)
