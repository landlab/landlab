"""Test StructuredQuadGraph."""

from numpy.testing import assert_array_equal
from pytest import approx

from landlab.graph import DualStructuredQuadGraph
from landlab.graph import DualUniformRectilinearGraph


def test_create():
    """Test creating a quad graph."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = DualStructuredQuadGraph((y, x), shape=(3, 3))

    assert graph.number_of_nodes == 9
    assert graph.number_of_links == 12
    assert graph.number_of_patches == 4

    assert graph.number_of_corners == 4
    assert graph.number_of_faces == 4
    assert graph.number_of_cells == 1


def test_create_raster():
    """Test creating a quad graph."""
    graph = DualUniformRectilinearGraph((3, 4), spacing=(2.0, 3.0))

    assert graph.number_of_nodes == 12
    assert graph.number_of_links == 17
    assert graph.number_of_patches == 6

    assert graph.number_of_corners == 6
    assert graph.number_of_faces == 7
    assert graph.number_of_cells == 2


def test_raster_spacing():
    """Test the spacing keyword for raster."""
    graph = DualUniformRectilinearGraph((3, 4), spacing=(2.0, 3.0))

    assert_array_equal(
        graph.length_of_link,
        [
            3.0,
            3.0,
            3.0,
            2.0,
            2.0,
            2.0,
            2.0,
            3.0,
            3.0,
            3.0,
            2.0,
            2.0,
            2.0,
            2.0,
            3.0,
            3.0,
            3.0,
        ],
    )
    assert_array_equal(graph.length_of_face, [3.0, 3.0, 2.0, 2.0, 2.0, 3.0, 3.0])


def test_raster_spacing_as_scalar():
    """Test the spacing keyword as a scalar for raster."""
    graph = DualUniformRectilinearGraph((3, 4), spacing=2.0)

    assert_array_equal(
        graph.length_of_link,
        [
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
        ],
    )
    assert_array_equal(graph.length_of_face, [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0])


def test_raster_origin():
    """Test the origin keyword for raster."""
    graph = DualUniformRectilinearGraph((3, 4), origin=(-1.0, 10.0))

    assert_array_equal(
        graph.xy_of_node[:, 0],
        [10.0, 11.0, 12.0, 13.0, 10.0, 11.0, 12.0, 13.0, 10.0, 11.0, 12.0, 13.0],
    )
    assert_array_equal(
        graph.xy_of_node[:, 1],
        [-1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
    )

    assert_array_equal(graph.xy_of_corner[:, 0], [10.5, 11.5, 12.5, 10.5, 11.5, 12.5])
    assert_array_equal(graph.xy_of_corner[:, 1], [-0.5, -0.5, -0.5, 0.5, 0.5, 0.5])


def test_raster_origin_as_scalar():
    """Test the origin keyword as a scalar for raster."""
    graph = DualUniformRectilinearGraph((3, 4), origin=-1.0)

    assert_array_equal(
        graph.xy_of_node[:, 0],
        [-1.0, 0.0, 1.0, 2.0, -1.0, 0.0, 1.0, 2.0, -1.0, 0.0, 1.0, 2.0],
    )
    assert_array_equal(
        graph.xy_of_node[:, 1],
        [-1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
    )

    assert_array_equal(graph.xy_of_corner[:, 0], [-0.5, 0.5, 1.5, -0.5, 0.5, 1.5])
    assert_array_equal(graph.xy_of_corner[:, 1], [-0.5, -0.5, -0.5, 0.5, 0.5, 0.5])


def test_perimeter_corners():
    """Test the perimeter corners."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = DualStructuredQuadGraph((y, x), shape=(3, 3))
    assert_array_equal(graph.perimeter_corners, [1, 3, 2, 0])


def test_length_of_face():
    """Test length of faces."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = DualStructuredQuadGraph((y, x), shape=(3, 3))

    assert graph.length_of_face == approx(1.5)


def test_area_of_cell():
    """Test areas of patches."""
    y = [0, 0, 0, 1, 1, 1, 3, 3, 3]
    x = [3, 4, 6, 3, 4, 6, 3, 4, 6]
    graph = DualStructuredQuadGraph((y, x), shape=(3, 3))
    assert graph.area_of_cell == approx([2.25])


def test_corners_at_cell():
    """Test corners of cells."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = DualStructuredQuadGraph((y, x), shape=(3, 3), sort=True)
    assert_array_equal(graph.corners_at_cell, [[3, 2, 0, 1]])


def test_cells_at_corner():
    """Test areas of patches."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = DualStructuredQuadGraph((y, x), shape=(3, 3))
    assert_array_equal(
        graph.cells_at_corner,
        [[0, -1, -1, -1], [-1, 0, -1, -1], [-1, -1, -1, 0], [-1, -1, 0, -1]],
    )


def test_cells_at_face():
    """Test cells on either side of faces."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = DualStructuredQuadGraph((y, x), shape=(3, 3))
    assert_array_equal(graph.cells_at_face, [[-1, 0], [0, -1], [-1, 0], [0, -1]])


def test_faces_at_cell():
    """Test faces that form cells."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = DualStructuredQuadGraph((y, x), shape=(3, 3))
    assert_array_equal(graph.faces_at_cell, [[2, 3, 1, 0]])


def test_corners_at_face():
    """Test corners at face tail and head."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = DualStructuredQuadGraph((y, x), shape=(3, 3))
    assert_array_equal(graph.corners_at_face, [[0, 1], [0, 2], [1, 3], [2, 3]])


def test_faces_at_corner():
    """Test faces around corners."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = DualStructuredQuadGraph((y, x), shape=(3, 3))
    assert_array_equal(
        graph.faces_at_corner,
        [[0, 1, -1, -1], [-1, 2, 0, -1], [3, -1, -1, 1], [-1, -1, 3, 2]],
    )


def test_face_dirs_at_corner():
    """Test face directions at corners."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = DualStructuredQuadGraph((y, x), shape=(3, 3))
    assert_array_equal(
        graph.face_dirs_at_corner,
        [[-1, -1, 0, 0], [0, -1, 1, 0], [-1, 0, 0, 1], [0, 0, 1, 1]],
    )


def test_cell_at_node():
    """Test cell-node connectivity."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = DualStructuredQuadGraph((y, x), shape=(3, 3))
    assert_array_equal(graph.cell_at_node, [-1, -1, -1, -1, 0, -1, -1, -1, -1])

    graph = DualUniformRectilinearGraph((3, 4))
    assert_array_equal(
        graph.cell_at_node, [-1, -1, -1, -1, -1, 0, 1, -1, -1, -1, -1, -1]
    )


def test_link_at_face():
    """Test link-face connectivity."""
    graph = DualUniformRectilinearGraph((3, 4))
    assert_array_equal(graph.link_at_face, [4, 5, 7, 8, 9, 11, 12])
    assert_array_equal(
        graph.face_at_link,
        [-1, -1, -1, -1, 0, 1, -1, 2, 3, 4, -1, 5, 6, -1, -1, -1, -1],
    )


def test_corner_at_face():
    """Test corner-face connectivity."""
    graph = DualUniformRectilinearGraph((3, 4))
    assert_array_equal(
        graph.corners_at_face, [[0, 1], [1, 2], [0, 3], [1, 4], [2, 5], [3, 4], [4, 5]]
    )
    assert_array_equal(graph.corner_at_face_tail, [0, 1, 0, 1, 2, 3, 4])
    assert_array_equal(graph.corner_at_face_head, [1, 2, 3, 4, 5, 4, 5])
