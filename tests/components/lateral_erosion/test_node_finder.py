import numpy as np
import pytest

from landlab import RasterModelGrid
from landlab.components.lateral_erosion.node_finder import angle_finder


@pytest.mark.parametrize(
    "node_1,node_2",
    [(12, 8), (13, 3), (8, 2), (3, 1), (6, 2), (11, 1), (12, 6), (13, 11)],
)
def test_angle_finder_90(node_1, node_2):
    grid = RasterModelGrid((4, 5))
    assert angle_finder(grid, node_1, 7, node_2) == pytest.approx(np.pi * 0.5)
    assert angle_finder(grid, node_2, 7, node_1) == pytest.approx(np.pi * 0.5)


@pytest.mark.parametrize(
    "node_1,node_2",
    [(10, 6), (6, 2), (2, 1), (1, 0), (0, 4), (4, 8), (8, 9), (9, 10)],
)
def test_angle_finder_45(node_1, node_2):
    grid = RasterModelGrid((3, 4))
    assert angle_finder(grid, node_1, 5, node_2) == pytest.approx(np.pi * 0.25)
    assert angle_finder(grid, node_2, 5, node_1) == pytest.approx(np.pi * 0.25)


@pytest.mark.parametrize("node", [6, 10, 9, 8, 4, 0, 1, 2])
def test_angle_finder_0(node):
    grid = RasterModelGrid((3, 4))
    assert angle_finder(grid, node, 5, node) == pytest.approx(0.0)


def test_angle_finder_array():
    grid = RasterModelGrid((3, 4))
    assert angle_finder(grid, (10, 6, 2), 5, (6, 2, 1)) == pytest.approx(np.pi * 0.25)


@pytest.mark.parametrize("dx", [0.5, 1.0, 2.0])
@pytest.mark.parametrize("dy", [0.5, 1.0, 2.0])
def test_unequal_spacing(dx, dy):
    grid = RasterModelGrid((3, 4), xy_spacing=(dx, dy))
    assert angle_finder(grid, (6, 9, 4, 1), 5, (1, 6, 9, 4)) == pytest.approx(
        np.pi * 0.5
    )
    assert angle_finder(grid, (10, 8, 0, 2), 5, (6, 4, 4, 6)) == pytest.approx(
        np.arctan(dy / dx)
    )
    assert angle_finder(grid, (9, 1, 1, 9), 5, (8, 0, 2, 10)) == pytest.approx(
        np.arctan(dx / dy)
    )
