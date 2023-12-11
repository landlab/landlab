"""Tests for the TriangleModelGrid class."""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from landlab.graph.triangle.mesh import TriangleMesh
from landlab.grid.triangle import TriangleModelGrid

if not TriangleMesh.validate_triangle():
    pytestmark = pytest.mark.skip(reason="triangle is not installed")


@pytest.fixture(scope="session")
def square_grid():
    return TriangleModelGrid(
        ([-1.0, -1.0, 11.0, 11.0], [0.0, 10.0, 10.0, 0.0]), triangle_opts="pqa1Devjz"
    )


@pytest.mark.parametrize("point", ("corner", "node"))
def test_all_points_in_box(square_grid, point):
    x, y = getattr(square_grid, f"x_of_{point}"), getattr(square_grid, f"y_of_{point}")

    assert np.all(x >= 0.0) and np.all(x <= 10.0)
    assert np.all(y >= -1.0) and np.all(y <= 11.0)


@pytest.mark.parametrize("edge", ("face", "link"))
def test_no_zero_length_edges(square_grid, edge):
    assert np.all(getattr(square_grid, f"length_of_{edge}") >= 0.0)


@pytest.mark.parametrize("polygon", ("cell", "patch"))
def test_no_zero_area_polygons(square_grid, polygon):
    assert np.all(getattr(square_grid, f"area_of_{polygon}") >= 0.0)


def test_boundary_nodes_on_boundary(square_grid):
    x, y = square_grid.x_of_node, square_grid.y_of_node
    assert np.all(
        (x[square_grid.boundary_nodes] == 0.0)
        | (x[square_grid.boundary_nodes] == 10.0)
        | (y[square_grid.boundary_nodes] == -1.0)
        | (y[square_grid.boundary_nodes] == 11.0)
    )
    assert square_grid.number_of_cells == square_grid.number_of_nodes - len(
        square_grid.boundary_nodes
    )


def test_grid_init():
    grid = TriangleModelGrid(
        ([-1.0, -1.0, 11.0, 11.0], [0.0, 10.0, 10.0, 0.0]), triangle_opts="pqa1Devjz"
    )
    assert grid.number_of_corners == grid.number_of_patches


def test_plot_nodes_and_links(square_grid):
    plot = square_grid.plot_nodes_and_links()

    assert len(plot.axes) == 1
    assert plot.axes[0].has_data()

    plt.close()
