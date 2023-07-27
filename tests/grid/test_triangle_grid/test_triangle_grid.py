"""Tests for the TriangleMeshGrid class."""

import numpy as np
import matplotlib.pyplot as plt
import pytest
from numpy.testing import assert_array_equal

from landlab.graph.triangle import TriangleGraph, TriangleMesh
from landlab.grid.triangle import TriangleMeshGrid
from landlab.plot.graph import plot_graph

if not TriangleMesh.validate_triangle():
    pytestmark = pytest.mark.skip(reason="triangle is not installed")


ys = [0, 0, 10, 10]
xs = [0, 10, 10, 0]


def test_grid_init():
    grid = TriangleMeshGrid(xs, ys, triangle_opts='pqa1Devjz')

def test_plot_nodes_and_links():
    grid = TriangleMeshGrid(xs, ys, triangle_opts='pqa1Devjz')
    plot = grid.plot_nodes_and_links(
        nodes_args={'s': 10},
        links_args={'linewidth': 1, 'linestyle': ':'}
    )
    plt.close()

def test_plot_graph():
    grid = TriangleMeshGrid(xs, ys, triangle_opts='pqa1Devjz')

    _ = plot_graph(grid, at = 'node')
    plt.show()

    _ = plot_graph(grid, at = 'link')
    plt.close()

    _ = plot_graph(grid, at = 'patch')
    plt.close()

