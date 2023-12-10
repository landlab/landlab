"""Tests for the TriangleModelGrid class."""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from landlab.graph.triangle.mesh import TriangleMesh
from landlab.grid.triangle import TriangleModelGrid

if not TriangleMesh.validate_triangle():
    pytestmark = pytest.mark.skip(reason="triangle is not installed")


ys = np.array([0, 0, 10, 10])
xs = np.array([0, 10, 10, 0])


def test_grid_init():
    grid = TriangleModelGrid((ys, xs), triangle_opts="pqa1Devjz")

    assert grid.number_of_nodes == 89


def test_plot_nodes_and_links():
    grid = TriangleModelGrid((ys, xs), triangle_opts="pqa1Devjz")
    plot = grid.plot_nodes_and_links()

    assert len(plot.axes) == 1
    assert plot.axes[0].has_data()

    plt.close()


def test_circular_polygon(datadir):
    grid = TriangleModelGrid.from_shapefile(
        datadir / "polygon_circular.geojson", triangle_opts="pqDjevz"
    )

    assert grid.number_of_nodes == 104
    assert len(grid.boundary_nodes) == 64
    assert grid.number_of_cells == grid.number_of_nodes - len(grid.boundary_nodes)
    assert len(grid.x_of_node) == grid.number_of_nodes
    assert len(grid.y_of_node) == grid.number_of_nodes
    assert grid.number_of_links == 245
    assert grid.number_of_patches == 142
    assert grid.number_of_corners == grid.number_of_patches
    assert grid.number_of_faces == 181
    assert grid.number_of_cells == 40


def test_concave_polygon(datadir):
    grid = TriangleModelGrid.from_shapefile(
        datadir / "polygon_concave.geojson", triangle_opts="pqa10Djevz"
    )

    assert grid.number_of_nodes == 25
    assert len(grid.boundary_nodes) == 24
    assert grid.number_of_cells == grid.number_of_nodes - len(grid.boundary_nodes)
    assert len(grid.x_of_node) == grid.number_of_nodes
    assert len(grid.y_of_node) == grid.number_of_nodes
    assert grid.number_of_links == 51
    assert grid.number_of_patches == 26
    assert grid.number_of_corners == grid.number_of_patches
    assert grid.number_of_faces == 27
    assert grid.number_of_cells == 1


def test_multiple_interior_rings(datadir):
    grid = TriangleModelGrid.from_shapefile(
        datadir / "polygon_two_interior_rings.geojson", triangle_opts="pqa10Djevz"
    )

    assert grid.number_of_nodes == 26
    assert len(grid.boundary_nodes) == 20
    assert grid.number_of_cells == grid.number_of_nodes - len(grid.boundary_nodes)
    assert len(grid.x_of_node) == grid.number_of_nodes
    assert len(grid.y_of_node) == grid.number_of_nodes
    assert grid.number_of_links == 61
    assert grid.number_of_patches == 34
    assert grid.number_of_corners == grid.number_of_patches
    assert grid.number_of_faces == 41
    assert grid.number_of_cells == 6
