import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.spatial import SphericalVoronoi

from landlab.graph.quasi_spherical.dual_icosphere import DualIcosphereGraph


@pytest.mark.parametrize("level", [0, 1, 2, 3])
@pytest.mark.parametrize("radius", [1.0, 6.371e6])
def test_area_of_cell_sums_to_sphere(level, radius):
    ico = DualIcosphereGraph(radius=radius, mesh_densification_level=level)
    assert_allclose(ico.area_of_cell.sum(), 4.0 * np.pi * radius**2, rtol=1e-12)
    assert_allclose(ico.area_of_patch.sum(), 4.0 * np.pi * radius**2, rtol=1e-12)


@pytest.mark.parametrize("level", [1, 2, 3])
def test_area_of_cell_matches_spherical_voronoi(level):
    # https://github.com/landlab/landlab/issues/2509
    ico = DualIcosphereGraph(radius=1.0, mesh_densification_level=level)
    points = np.column_stack([ico.x_of_node, ico.y_of_node, ico.z_of_node])
    expected = SphericalVoronoi(points, radius=1.0).calculate_areas()
    assert_allclose(ico.area_of_cell, expected[ico.node_at_cell], rtol=1e-10)


def test_area_of_cell_base_icosahedron():
    ico = DualIcosphereGraph()
    assert ico.number_of_cells == 12
    assert_allclose(ico.area_of_cell, 4.0 * np.pi / 12.0)
