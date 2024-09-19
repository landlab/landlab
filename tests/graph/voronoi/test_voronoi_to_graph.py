import inspect

import numpy as np
import pytest
from pytest import approx
from scipy.spatial import Delaunay
from scipy.spatial import Voronoi

from landlab.graph.voronoi.voronoi_to_graph import VoronoiDelaunay
from landlab.graph.voronoi.voronoi_to_graph import VoronoiDelaunayToGraph

XY_OF_NODE = {
    "rect-horizontal-3-3": [
        [0.0, 0.0],
        [1.0, 0.0],
        [2.0, 0.0],
        [0.5, 1.0],
        [1.5, 1.0],
        [2.5, 1.0],
        [0.0, 2.0],
        [1.0, 2.0],
        [2.0, 2.0],
    ],
    "rect-horizontal-4-3": [
        [0.0, 0.0],
        [1.0, 0.0],
        [2.0, 0.0],
        [0.5, 1.0],
        [1.5, 1.0],
        [2.5, 1.0],
        [0.0, 2.0],
        [1.0, 2.0],
        [2.0, 2.0],
        [0.0, 3.0],
        [1.0, 3.0],
        [2.0, 3.0],
    ],
    "rect-vertical-3-3": [
        [0.0, 0.0],
        [2.0, 0.0],
        [1.0, 0.5],
        [0.0, 1.0],
        [2.0, 1.0],
        [1.0, 1.5],
        [0.0, 2.0],
        [2.0, 2.0],
        [1.0, 2.5],
    ],
    "rect-vertical-3-4": [
        [0.0, 0.0],
        [2.0, 0.0],
        [1.0, 0.5],
        [3.0, 0.5],
        [0.0, 1.0],
        [2.0, 1.0],
        [1.0, 1.5],
        [3.0, 1.5],
        [0.0, 2.0],
        [2.0, 2.0],
        [1.0, 2.5],
        [3.0, 2.5],
    ],
}


@pytest.fixture
def hex_graph():
    return VoronoiDelaunayToGraph(XY_OF_NODE["rect-horizontal-3-3"])


@pytest.fixture
def xy_of_hex():
    return XY_OF_NODE["rect-horizontal-3-3"]


def pytest_generate_tests(metafunc):
    if "at_property" in metafunc.fixturenames:
        props = dict(
            inspect.getmembers(
                VoronoiDelaunayToGraph, lambda o: isinstance(o, property)
            )
        )
        metafunc.parametrize(
            "at_property", [name for name in props.keys() if "_at_" in name]
        )
    if "of_property" in metafunc.fixturenames:
        props = dict(
            inspect.getmembers(
                VoronoiDelaunayToGraph, lambda o: isinstance(o, property)
            )
        )
        metafunc.parametrize(
            "of_property",
            [
                name
                for name in props.keys()
                if "_of_" in name and not name.startswith("number")
            ],
        )


def test_voronoi_name_mapping(xy_of_hex):
    """Test scipy Voronoi names are mapped to landlab-style names."""
    voronoi = Voronoi(xy_of_hex)
    delaunay = Delaunay(xy_of_hex)
    graph = VoronoiDelaunay(xy_of_hex)

    voronoi.regions, voronoi.point_region = VoronoiDelaunay._remove_empty_regions(
        voronoi.regions, voronoi.point_region
    )
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


def test_at_array_is_int(hex_graph, at_property):
    """Test that _at_ properties are arrays of int."""
    assert getattr(hex_graph, at_property).dtype == int


def test_degenerate_case():
    xy_of_node = np.array(
        [[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1], [0, 2], [1, 2], [2, 2]],
        dtype=float,
    )
    VoronoiDelaunay(xy_of_node)
    # VoronoiDelaunayToGraph(xy_of_node)


def test_of_array_is_float(hex_graph, of_property):
    """Test that _of_ properties are arrays of float."""
    xy_of_node = np.array(
        [[0, 0], [2, 0], [4, 0], [1, 1], [3, 1], [5, 1], [0, 2], [2, 2], [4, 2]],
        dtype=float,
    )
    # hex_graph = VoronoiDelaunayToGraph(xy_of_node)
    hex_graph = VoronoiDelaunay(xy_of_node)
    assert getattr(hex_graph, of_property).dtype == float


@pytest.mark.parametrize(
    "element,expected",
    [
        ("nodes", 9),
        ("links", 17),
        ("patches", 9),
        ("corners", 9),
        ("faces", 10),
        ("cells", 2),
    ],
)
def test_element_count_without_perimeter_nodes(hex_graph, element, expected):
    assert getattr(hex_graph, f"number_of_{element}") == expected


@pytest.mark.parametrize(
    "element,expected",
    [
        ("nodes", 9),
        ("links", 16),
        ("patches", 8),
        ("corners", 8),
        ("faces", 8),
        ("cells", 1),
    ],
)
def test_element_count_with_perimeter_nodes(xy_of_hex, element, expected):
    perimeter_links = [[0, 1], [1, 2], [2, 5], [5, 8], [8, 7], [7, 6], [6, 3], [3, 0]]
    graph = VoronoiDelaunayToGraph(xy_of_hex, perimeter_links=perimeter_links)
    assert getattr(graph, f"number_of_{element}") == expected


@pytest.mark.parametrize("at", ("node", "link", "cell", "corner", "face", "cell"))
def test_compact_ids_without_perimeter_nodes(hex_graph, at):
    ids = []
    for name in hex_graph.ids_with_prefix(at):
        ids.append(np.sort(np.unique(getattr(hex_graph, name).reshape((-1,)))))
    ids = np.sort(np.unique(np.concatenate(ids)))
    ids = ids[ids >= 0]

    assert ids[0] >= 0
    assert ids[-1] <= hex_graph._mesh.sizes[at]


@pytest.mark.parametrize("at", ("node", "link", "cell", "corner", "face", "cell"))
def test_compact_ids_with_perimeter_nodes(xy_of_hex, at):
    perimeter_links = [[0, 1], [1, 2], [2, 5], [5, 8], [8, 7], [7, 6], [6, 3], [3, 0]]
    graph = VoronoiDelaunayToGraph(xy_of_hex, perimeter_links=perimeter_links)

    ids = []
    for name in graph.ids_with_prefix(at):
        ids.append(np.sort(np.unique(getattr(graph, name).reshape((-1,)))))
    ids = np.sort(np.unique(np.concatenate(ids)))
    ids = ids[ids >= 0]

    assert ids[0] >= 0
    assert ids[-1] <= graph._mesh.sizes[at]


@pytest.mark.parametrize("at", ["node", "link", "patch", "corner", "face", "cell"])
def test_has_prefix(hex_graph, at):
    expected = {
        "node": ("nodes_at_patch", "nodes_at_face", "node_at_cell", "nodes_at_link"),
        "link": ("links_at_patch",),
        "patch": (),
        "corner": ("corners_at_face", "corners_at_cell"),
        "face": ("faces_at_cell",),
        "cell": ("cell_at_node",),
    }

    assert hex_graph.ids_with_prefix(at) == set(expected[at])


@pytest.mark.parametrize("at", ["node", "link", "patch", "corner", "face", "cell"])
def test_has_suffix(hex_graph, at):
    expected = {
        "node": ("cell_at_node",),
        "link": ("nodes_at_link",),
        "patch": ("nodes_at_patch", "links_at_patch"),
        "corner": (),
        "face": ("corners_at_face", "nodes_at_face"),
        "cell": (
            "n_corners_at_cell",
            "faces_at_cell",
            "node_at_cell",
            "corners_at_cell",
        ),
    }

    assert hex_graph.ids_with_suffix(at) == set(expected[at])


@pytest.mark.parametrize(
    "n_nodes",
    [2**10, 2**11, 2**12, 2**13, 2**14, 2**15],  # , 2 ** 16, 2 ** 20]
)
def test_big_graph(n_nodes):
    xy_of_node = np.random.rand(2 * n_nodes).reshape((-1, 2))
    graph = VoronoiDelaunayToGraph(xy_of_node)
    assert graph.number_of_nodes == n_nodes
