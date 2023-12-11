"""Tests for the TriangleModelGrid class."""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from landlab.graph.triangle.mesh import TriangleMesh
from landlab.grid.triangle import TriangleModelGrid

if not TriangleMesh.validate_triangle():
    pytestmark = pytest.mark.skip(reason="triangle is not installed")


def test_grid_init():
    grid = TriangleModelGrid(
        ([-1.0, -1.0, 11.0, 11.0], [0.0, 10.0, 10.0, 0.0]), triangle_opts="pqa1Devjz"
    )

    assert np.all(grid.x_of_node >= 0.0) and np.all(grid.x_of_node <= 10.0)
    assert np.all(grid.y_of_node >= -1.0) and np.all(grid.y_of_node <= 11.0)
    assert np.all(grid.length_of_link >= 0.0)
    assert np.all(grid.length_of_face >= 0.0)
    assert np.all(grid.area_of_cell >= 0.0)
    assert np.all(grid.area_of_patch >= 0.0)


def test_plot_nodes_and_links():
    grid = TriangleModelGrid(
        ([0.0, 0.0, 11.0, 11.0], [0.0, 10.0, 10.0, 0.0]), triangle_opts="pqa1Devjz"
    )
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
        datadir / "polygon_two_interior_rings.geojson", triangle_opts="pq10a10Djevz"
    )

    assert np.all(grid.x_of_node >= 0.0) and np.all(grid.x_of_node <= 10.0)
    assert np.all(grid.y_of_node >= -1.0) and np.all(grid.y_of_node <= 11.0)
    assert np.all(grid.length_of_link >= 0.0)
    assert np.all(grid.length_of_face >= 0.0)
    assert np.all(grid.area_of_cell >= 0.0)
    assert np.all(grid.area_of_patch >= 0.0)

    assert grid.number_of_cells == grid.number_of_nodes - len(grid.boundary_nodes)
    assert grid.number_of_corners == grid.number_of_patches
