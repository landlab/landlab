"""Test StructuredQuadGraph."""
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal

from landlab.graph import DualRadialGraph


def test_create():
    """Test creating a quad graph."""
    graph = DualRadialGraph((1, 4))

    assert graph.number_of_nodes == 5
    assert graph.number_of_links == 8
    assert graph.number_of_patches == 4

    assert graph.number_of_corners == 4
    assert graph.number_of_faces == 4
    assert graph.number_of_cells == 1


def test_spacing():
    """Test the spacing keyword for raster."""
    graph = DualRadialGraph((1, 4))

    assert_array_almost_equal(
        graph.xy_of_node, [[0.0, -1.0], [-1.0, 0.0], [0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
    )

    assert_array_almost_equal(
        graph.xy_of_corner, [[-0.5, -0.5], [0.5, -0.5], [-0.5, 0.5], [0.5, 0.5]]
    )


def test_spacing_keyword():
    """Test the spacing keyword for raster."""
    graph = DualRadialGraph((1, 4), spacing=2.0)

    assert_array_almost_equal(
        graph.xy_of_node, [[0.0, -2.0], [-2.0, 0.0], [0.0, 0.0], [2.0, 0.0], [0.0, 2.0]]
    )

    assert_array_almost_equal(
        graph.xy_of_corner, [[-1.0, -1.0], [1.0, -1.0], [-1.0, 1.0], [1.0, 1.0]]
    )


def test_origin():
    """Test the origin keyword for raster."""
    graph = DualRadialGraph((1, 4), spacing=2.0, origin=(-1.0, 2))

    assert_array_almost_equal(
        graph.xy_of_node,
        [[2.0, -3.0], [0.0, -1.0], [2.0, -1.0], [4.0, -1.0], [2.0, 1.0]],
    )

    assert_array_almost_equal(
        graph.xy_of_corner, [[1.0, -2.0], [3.0, -2.0], [1.0, 0.0], [3.0, 0.0]]
    )


# def test_perimeter_corners():
#     """Test the perimeter corners."""
#     y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
#     x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
#     graph = DualStructuredQuadGraph((y, x), shape=(3, 3))
#     assert_array_equal(graph.perimeter_corners, [1, 3, 2, 0])


def test_length_of_face_and_link():
    """Test length of faces and links."""
    ROOT_2 = np.sqrt(2.0)

    graph = DualRadialGraph((1, 4))

    assert_array_almost_equal(
        graph.length_of_link, [ROOT_2, 1.0, ROOT_2, 1.0, 1.0, ROOT_2, 1.0, ROOT_2]
    )
    assert_array_almost_equal(graph.length_of_face, [1.0, 1.0, 1.0, 1.0])


def test_area_of_cell_and_patch():
    """Test areas of cells patches."""
    graph = DualRadialGraph((1, 4))

    assert_array_almost_equal(graph.area_of_patch, [0.5, 0.5, 0.5, 0.5])
    assert_array_almost_equal(graph.area_of_cell, [1.0])


def test_corners_at_cell():
    """Test corners of cells."""
    graph = DualRadialGraph((1, 4))

    assert_array_equal(
        graph.nodes_at_patch, [[2, 1, 0], [2, 0, 3], [4, 1, 2], [3, 4, 2]]
    )
    assert_array_equal(graph.corners_at_cell, [[3, 2, 0, 1]])


def test_cells_at_corner():
    """Test areas of patches."""
    graph = DualRadialGraph((1, 4))

    assert_array_equal(
        graph.patches_at_node,
        [[0, 1, -1, -1], [0, 2, -1, -1], [0, 1, 2, 3], [1, 3, -1, -1], [2, 3, -1, -1]],
    )
    assert_array_equal(graph.cells_at_corner, [[0], [0], [0], [0]])


def test_cells_at_face():
    """Test cells on either side of faces."""
    graph = DualRadialGraph((1, 4))

    assert_array_equal(
        graph.patches_at_link,
        [[0, -1], [0, 1], [1, -1], [0, 2], [1, 3], [2, -1], [2, 3], [3, -1]],
    )
    assert_array_equal(graph.cells_at_face, [[0, -1], [0, -1], [0, -1], [0, -1]])


def test_faces_at_cell():
    """Test faces that form cells."""
    graph = DualRadialGraph((1, 4))
    assert_array_equal(
        graph.links_at_patch, [[1, 3, 0], [4, 1, 2], [6, 5, 3], [4, 7, 6]]
    )
    assert_array_equal(graph.faces_at_cell, [[2, 3, 1, 0]])


def test_corners_at_face():
    """Test corners at face tail and head."""
    graph = DualRadialGraph((1, 4))

    assert_array_equal(graph.corners_at_face, [[0, 1], [0, 2], [1, 3], [2, 3]])
    assert_array_equal(
        graph.nodes_at_link,
        [[1, 0], [0, 2], [0, 3], [1, 2], [2, 3], [1, 4], [2, 4], [3, 4]],
    )


def test_faces_at_corner():
    """Test faces around corners."""
    graph = DualRadialGraph((1, 4))

    assert_array_equal(graph.faces_at_corner, [[0, 1], [2, 0], [3, 1], [3, 2]])
    assert_array_equal(
        graph.links_at_node,
        [[2, 1, 0, -1], [3, 5, 0, -1], [4, 6, 3, 1], [7, 4, 2, -1], [5, 6, 7, -1]],
    )


def test_face_dirs_at_corner():
    """Test face directions at corners."""
    graph = DualRadialGraph((1, 4))

    assert_array_equal(graph.face_dirs_at_corner, [[-1, -1], [-1, 1], [-1, 1], [1, 1]])
    assert_array_equal(
        graph.link_dirs_at_node,
        [[-1, -1, 1, 0], [-1, -1, -1, 0], [-1, -1, 1, 1], [-1, 1, 1, 0], [1, 1, 1, 0]],
    )
