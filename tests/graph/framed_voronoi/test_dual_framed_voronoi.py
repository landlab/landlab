from numpy.testing import assert_array_equal

from landlab.graph import DualFramedVoronoiGraph


def test_rect_create():
    """Test creating a dual hex graph with rectangular layout."""
    graph = DualFramedVoronoiGraph((4, 3))

    assert graph.number_of_nodes == 12
    assert graph.number_of_links == 23
    assert graph.number_of_patches == 12

    assert graph.number_of_corners == 12
    assert graph.number_of_faces == 13
    assert graph.number_of_cells == 2


def test_adjacent_corners_at_corner():
    graph = DualFramedVoronoiGraph(
        (3, 3),
        sort=True,
    )
    assert_array_equal(
        graph.adjacent_corners_at_corner,
        [
            [2, -1, -1],
            [3, 2, -1],
            [4, 0, 1],
            [6, 1, -1],
            [5, 2, -1],
            [6, 4, -1],
            [7, 5, 3],
            [6, -1, -1],
        ],
    )
