"""Tests for the TriangleMeshGrid class."""

import numpy as np
import matplotlib.pyplot as plt
import pytest
from numpy.testing import assert_array_equal

from landlab.graph.triangle import TriangleMesh
from landlab.grid.triangle import TriangleMeshGrid
from landlab.plot.graph import plot_graph

if not TriangleMesh.validate_triangle():
    pytestmark = pytest.mark.skip(reason="triangle is not installed")


ys = [0, 0, 10, 10]
xs = [0, 10, 10, 0]
    

def test_grid_init():
    grid = TriangleMeshGrid((ys, xs), triangle_opts="pqa1Devjz")


def test_grid_from_shapefile(datadir):
    grid = TriangleMeshGrid.from_shapefile(
        datadir / "example_geojson_for_grid_no_holes.geojson", triangle_opts="pq10Djevz"
    )

    assert grid.number_of_nodes == 1068


def test_grid_from_file_with_holes(datadir):
    grid = TriangleMeshGrid.from_shapefile(
        datadir / "example_geojson_for_grid.geojson", triangle_opts="pqDjevz"
    )
    

def test_plot_nodes_and_links():
    grid = TriangleMeshGrid((ys, xs), triangle_opts="pqa1Devjz")
    plot = grid.plot_nodes_and_links()

    assert len(plot.axes) == 1
    assert plot.axes[0].has_data()

    plt.close()
