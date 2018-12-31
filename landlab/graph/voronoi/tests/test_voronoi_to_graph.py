import inspect

import numpy as np
from pytest import approx
from scipy.spatial import Delaunay, Voronoi

from landlab.graph.voronoi.voronoi_to_graph import VoronoiDelaunay, VoronoiDelaunayToGraph


def pytest_generate_tests(metafunc):
    if "at_property" in metafunc.fixturenames:
        props = dict(
            inspect.getmembers(VoronoiDelaunayToGraph, lambda o: isinstance(o, property))
        )
        metafunc.parametrize(
            "at_property", [name for name in props.keys() if "_at_" in name]
        )
    elif "of_property" in metafunc.fixturenames:
        props = dict(
            inspect.getmembers(VoronoiDelaunayToGraph, lambda o: isinstance(o, property))
        )
        metafunc.parametrize(
            "of_property", [name for name in props.keys() if "_of_" in name and not name.startswith("number")]
        )


def test_voronoi_name_mapping():
    """Test scipy Voronoi names are mapped to landlab-style names."""
    xy_of_node = [
        [0.0, 0.0],
        [1.0, 0.0],
        [2.0, 0.0],
        [0.5, 1.0],
        [1.5, 1.0],
        [2.5, 1.0],
        [0.0, 2.0],
        [1.0, 2.0],
        [2.0, 2.0],
    ]
    voronoi = Voronoi(xy_of_node)
    delaunay = Delaunay(xy_of_node)
    graph = VoronoiDelaunay(xy_of_node)

    assert np.all(graph.x_of_node == approx(voronoi.points[:, 0]))
    assert np.all(graph.y_of_node == approx(voronoi.points[:, 1]))

    assert np.all(graph.x_of_corner == approx(voronoi.vertices[:, 0]))
    assert np.all(graph.y_of_corner == approx(voronoi.vertices[:, 1]))

    assert np.all(graph.nodes_at_link == voronoi.ridge_points)

    assert tuple(graph.n_corners_at_cell) == tuple(
        len(region) for region in voronoi.regions
    )
    for cell, corners in enumerate(graph.corners_at_cell):
        assert np.all(corners[: graph.n_corners_at_cell[cell]] == voronoi.regions[cell])
        assert np.all(corners[graph.n_corners_at_cell[cell] :] == -1)
    assert np.all(graph.corners_at_face == voronoi.ridge_vertices)
    assert np.all(graph.nodes_at_face == voronoi.ridge_points)
    assert np.all(graph.cell_at_node == voronoi.point_region)

    assert np.all(graph.nodes_at_patch == delaunay.simplices)


def test_at_array_is_int(at_property):
    """Test that _at_ properties are arrays of int."""
    xy_of_node = [
        [0.0, 0.0],
        [1.0, 0.0],
        [2.0, 0.0],
        [0.5, 1.0],
        [1.5, 1.0],
        [2.5, 1.0],
        [0.0, 2.0],
        [1.0, 2.0],
        [2.0, 2.0],
    ]
    graph = VoronoiDelaunayToGraph(xy_of_node)
    assert getattr(graph, at_property).dtype == int


def test_of_array_is_float(of_property):
    """Test that _of_ properties are arrays of float."""
    xy_of_node = [
        [0, 0],
        [1, 0],
        [2, 0],
        [1, 1],
        [2, 1],
        [3, 1],
        [0, 2],
        [1, 2],
        [2, 2],
    ]
    graph = VoronoiDelaunayToGraph(xy_of_node)
    assert getattr(graph, of_property).dtype == float


def test_without_perimeter_nodes():
    xy_of_node = [
        [0.0, 0.0],
        [1.0, 0.0],
        [2.0, 0.0],
        [0.5, 1.0],
        [1.5, 1.0],
        [2.5, 1.0],
        [0.0, 2.0],
        [1.0, 2.0],
        [2.0, 2.0],
    ]
    perimeter_links = None
    graph = VoronoiDelaunayToGraph(xy_of_node, perimeter_links=perimeter_links)
    assert graph.number_of_nodes == 9
    assert graph.number_of_links == 17
    assert graph.number_of_patches == 9
    assert graph.number_of_corners == 9
    assert graph.number_of_faces == 10
    assert graph.number_of_cells == 2


def test_with_perimeter_nodes():
    xy_of_node = [
        [0.0, 0.0],
        [1.0, 0.0],
        [2.0, 0.0],
        [0.5, 1.0],
        [1.5, 1.0],
        [2.5, 1.0],
        [0.0, 2.0],
        [1.0, 2.0],
        [2.0, 2.0],
    ]
    perimeter_links = [[0, 1], [1, 2], [2, 5], [5, 8], [8, 7], [7, 6], [6, 3], [3, 0]]
    graph = VoronoiDelaunayToGraph(xy_of_node, perimeter_links=perimeter_links)
    assert graph.number_of_nodes == 9
    assert graph.number_of_links == 16
    assert graph.number_of_patches == 8
    assert graph.number_of_corners == 8
    assert graph.number_of_faces == 8
    assert graph.number_of_cells == 1
