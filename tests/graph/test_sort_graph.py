import numpy as np
import pytest
from numpy.testing import assert_allclose
from numpy.testing import assert_array_equal

from landlab.graph._reindexer import Reindexer
from landlab.graph._sort_graph import SortConnections
from landlab.graph._sort_graph import SortGraph
from landlab.graph._sort_graph import SortReindexers
from landlab.graph._sort_graph import calc_xy_of_edge
from landlab.graph._sort_graph import calc_xy_of_polygon
from landlab.graph._sort_graph import redirect_edges
from landlab.graph._sort_graph import remap_graph_connections
from landlab.graph._sort_graph import sort_connected_graphs
from landlab.graph._sort_graph import sort_graph_by_ids
from landlab.graph._sort_graph import sort_graph_by_xy


@pytest.fixture
def triangle():
    """A triangle (3 points, 3 edges, 1 polygon) in scrambled order.

    Physical layout::

        2 (0,1)
        |\\
        | \\
        |  \\
        0---1
      (0,0) (1,0)

    Points are stored as (0,1), (0,0), (1,0) — not in sorted order.
    """
    return SortGraph(
        xy_of_point=np.array([[0.0, 1.0], [0.0, 0.0], [1.0, 0.0]]),
        points_at_edge=np.array([[0, 1], [1, 2], [0, 2]]),
        edges_at_polygon=np.array([[0, 1, 2]]),
    )


@pytest.fixture
def triangle_sorted():
    """The same triangle with points in sorted (y-then-x) order."""
    return SortGraph(
        xy_of_point=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
        points_at_edge=np.array([[0, 1], [2, 0], [2, 1]]),
        edges_at_polygon=np.array([[1, 0, 2]]),
    )


def test_calc_xy_of_edge_single():
    xy = np.array([[0.0, 0.0], [2.0, 4.0]])
    result = calc_xy_of_edge(xy, np.array([[0, 1]]))
    assert_allclose(result, [[1.0, 2.0]])


def test_calc_xy_of_edge_multiple():
    xy = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]])
    edges = np.array([[0, 1], [0, 2], [1, 2]])
    result = calc_xy_of_edge(xy, edges)
    assert_allclose(result, [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])


def test_calc_xy_of_polygon_triangle():
    # Midpoints of a right triangle's edges
    xy_of_edge = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    edges_at_polygon = np.array([[0, 1, 2]])
    result = calc_xy_of_polygon(xy_of_edge, edges_at_polygon)
    assert_allclose(result, [[2 / 3, 2 / 3]])


def test_calc_xy_of_polygon_padded():
    # Two polygons: a square (4 edges) and a triangle (3 edges) padded to width 4
    xy_of_edge = np.array(
        [
            [0.5, 0.0],  # 0: square bottom
            [1.0, 0.5],  # 1: square right
            [0.5, 1.0],  # 2: square top
            [0.0, 0.5],  # 3: square left
            [1.5, 0.0],  # 4: triangle bottom
            [1.75, 0.5],  # 5: triangle right
            [1.25, 0.5],  # 6: triangle left
        ]
    )
    edges_at_polygon = np.array(
        [
            [0, 1, 2, 3],
            [4, 5, 6, -1],
        ]
    )
    result = calc_xy_of_polygon(xy_of_edge, edges_at_polygon)
    assert_allclose(result[0], [0.5, 0.5])
    assert_allclose(result[1], [1.5, 1 / 3])


@pytest.mark.parametrize(
    "tail, head, expected_tail, expected_head",
    [
        # Right (+x): keep
        ([0.0, 0.0], [1.0, 0.0], [0.0, 0.0], [1.0, 0.0]),
        # Left (-x): flip
        ([1.0, 0.0], [0.0, 0.0], [0.0, 0.0], [1.0, 0.0]),
        # Up (+y): keep
        ([0.0, 0.0], [0.0, 1.0], [0.0, 0.0], [0.0, 1.0]),
        # Down (-y): flip
        ([0.0, 1.0], [0.0, 0.0], [0.0, 0.0], [0.0, 1.0]),
        # NE diagonal (45°): keep
        ([0.0, 0.0], [1.0, 1.0], [0.0, 0.0], [1.0, 1.0]),
        # SW diagonal (225°): flip
        ([1.0, 1.0], [0.0, 0.0], [0.0, 0.0], [1.0, 1.0]),
        # SE diagonal (-45°, boundary, dx>0): keep
        ([0.0, 0.0], [1.0, -1.0], [0.0, 0.0], [1.0, -1.0]),
        # NW diagonal (135°, boundary, dx<0): flip
        ([0.0, 0.0], [-1.0, 1.0], [-1.0, 1.0], [0.0, 0.0]),
    ],
)
def test_redirect_edges_directions(tail, head, expected_tail, expected_head):
    xy_of_point = np.array([tail, head])
    points_at_edge = np.array([[0, 1]])
    result = redirect_edges(xy_of_point, points_at_edge)
    et = np.array(expected_tail)
    eh = np.array(expected_head)
    assert_allclose(xy_of_point[result[0, 0]], et)
    assert_allclose(xy_of_point[result[0, 1]], eh)


def test_redirect_edges_docstring_example():
    xy_of_point = np.array(
        [
            [0, 0],
            [1, 0],
            [1, 1],
            [0, 1],
            [-1, 1],
            [-1, 0],
            [-1, 1],
            [0, -1],
            [1, -1],
        ],
        dtype=float,
    )
    points_at_edge = np.array(
        [
            [0, 1],
            [0, 2],
            [0, 3],
            [0, 4],
            [0, 5],
            [0, 6],
            [0, 7],
            [0, 8],
        ]
    )
    result = redirect_edges(xy_of_point, points_at_edge)
    assert_array_equal(
        result,
        [
            [0, 1],
            [0, 2],
            [0, 3],
            [4, 0],
            [5, 0],
            [6, 0],
            [7, 0],
            [0, 8],
        ],
    )


def test_redirect_edges_into_out():
    xy = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    pts = np.array([[1, 0], [0, 2]])  # left-pointing, up-pointing
    out = np.empty_like(pts)
    result = redirect_edges(xy, pts, out=out)
    assert result is out


def test_redirect_edges_wrong_shape_raises():
    xy = np.array([[0.0, 0.0], [1.0, 0.0]])
    with pytest.raises(ValueError, match="shape"):
        redirect_edges(xy, np.array([[0, 1, 0]]))


@pytest.mark.parametrize("shape", ((-1,), (2, -1)))
def test_redirect_edges_wrong_xy_shape_raises(shape):
    xy = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]).reshape(shape)
    with pytest.raises(ValueError, match="^xy_of_point must have shape"):
        redirect_edges(xy, [[1, 0], [0, 2]])


def test_redirect_edges_out_wrong_shape_raises():
    xy = np.array([[0.0, 0.0], [1.0, 0.0]])
    pts = np.array([[0, 1]])
    with pytest.raises(ValueError, match="shape"):
        redirect_edges(xy, pts, out=np.empty((2, 2), dtype=int))


def test_sort_graph_by_xy_identity(triangle_sorted):
    sorted_graph, _ = sort_graph_by_xy(triangle_sorted)
    assert_allclose(sorted_graph.xy_of_point, triangle_sorted.xy_of_point)


def test_sort_graph_by_xy_points_are_sorted(triangle):
    sorted_graph, _ = sort_graph_by_xy(triangle)
    xy = sorted_graph.xy_of_point
    # Primary key: y (col 1), secondary: x (col 0)
    assert np.all(xy[1:, 1] >= xy[:-1, 1])


def test_sort_graph_by_xy_geometric_consistency(triangle):
    sorted_graph, _ = sort_graph_by_xy(triangle)
    xy = sorted_graph.xy_of_point
    pts = sorted_graph.points_at_edge
    # Midpoints should be sorted (y then x) after the edge reindexing
    midpoints = xy[pts].mean(axis=1)
    assert np.all(midpoints[1:, 1] >= midpoints[:-1, 1])


def test_sort_graph_by_xy_physical_coords_preserved(triangle):
    # Each edge should connect the same pair of physical coordinates before/after sort
    orig = triangle
    sorted_graph, _ = sort_graph_by_xy(orig)

    def edge_coord_sets(g):
        return {frozenset(map(tuple, g.xy_of_point[edge])) for edge in g.points_at_edge}

    assert edge_coord_sets(sorted_graph) == edge_coord_sets(orig)


def test_sort_graph_by_xy_returns_reindexers(triangle):
    _, reindex = sort_graph_by_xy(triangle)
    assert len(reindex.point) == 3
    assert len(reindex.edge) == 3
    assert len(reindex.polygon) == 1


def test_sort_graph_by_ids_orders_polygons_by_point_at_polygon():
    # 2 polygons: point_at_polygon=[2,0] → after sort → [0,2] (polygon at pt 0 first)
    secondary = SortGraph(
        xy_of_point=np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]]),
        points_at_edge=np.array([[0, 1], [1, 2], [0, 2]]),
        edges_at_polygon=np.array([[0, 1, 2], [0, 1, 2]]),
    )
    connected = SortConnections(
        point_at_polygon=np.array([2, 0]),
        edge_at_edge=np.array([0, 1, 2]),
        polygon_at_point=np.array([0, 0, 0]),
    )
    _, sorted_conn = sort_graph_by_ids(secondary, connected)
    assert_array_equal(sorted_conn.point_at_polygon, [0, 2])


def test_sort_graph_by_ids_orders_edges_by_edge_at_edge():
    secondary = SortGraph(
        xy_of_point=np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]]),
        points_at_edge=np.array([[0, 1], [1, 2], [0, 2]]),
        edges_at_polygon=np.array([[0, 1, 2]]),
    )
    connected = SortConnections(
        point_at_polygon=np.array([0]),
        edge_at_edge=np.array([2, 0, 1]),
        polygon_at_point=np.array([0, 0, 0]),
    )
    _, sorted_conn = sort_graph_by_ids(secondary, connected)
    assert_array_equal(sorted_conn.edge_at_edge, [0, 1, 2])


def test_sort_graph_by_ids_orders_points_by_polygon_at_point():
    secondary = SortGraph(
        xy_of_point=np.array([[0.0, 2.0], [0.0, 0.0], [0.0, 1.0]]),
        points_at_edge=np.array([[0, 1], [1, 2]]),
        edges_at_polygon=np.array([[0, 1, -1]]),
    )
    connected = SortConnections(
        point_at_polygon=np.array([0]),
        edge_at_edge=np.array([0, 1]),
        polygon_at_point=np.array([2, 0, 1]),
    )
    sorted_secondary, sorted_conn = sort_graph_by_ids(secondary, connected)
    assert_array_equal(sorted_conn.polygon_at_point, [0, 1, 2])
    # Point that was at polygon 0 should now be first
    assert sorted_secondary.xy_of_point[0, 1] == 0.0


@pytest.mark.parametrize(
    "point_order, edge_order, polygon_order, conn, expected",
    [
        (
            [0, 1, 2],
            [0, 1, 2],
            [0, 1, 2],
            ([0, 1, 2], [0, 1, 2], [0, 1, 2]),
            ([0, 1, 2], [0, 1, 2], [0, 1, 2]),
        ),
        (
            [1, 0, 2],
            [2, 1, 0],
            [0, 2, 1],
            ([0, 1, 2], [0, 1, 2], [0, 1, 2]),
            ([1, 0, 2], [2, 1, 0], [0, 2, 1]),
        ),
    ],
)
def test_remap_graph_connections(
    point_order, edge_order, polygon_order, conn, expected
):
    reindex = SortReindexers(
        point=Reindexer(point_order),
        edge=Reindexer(edge_order),
        polygon=Reindexer(polygon_order),
    )
    connected = SortConnections(
        point_at_polygon=np.array(conn[0]),
        edge_at_edge=np.array(conn[1]),
        polygon_at_point=np.array(conn[2]),
    )
    result = remap_graph_connections(connected, reindex)
    assert_array_equal(result.point_at_polygon, expected[0])
    assert_array_equal(result.edge_at_edge, expected[1])
    assert_array_equal(result.polygon_at_point, expected[2])


def test_remap_graph_connections_fill_value_preserved():
    reindex = SortReindexers(
        point=Reindexer([1, 0]),
        edge=Reindexer([1, 0]),
        polygon=Reindexer([1, 0]),
    )
    connected = SortConnections(
        point_at_polygon=np.array([-1, 0]),
        edge_at_edge=np.array([1, -1]),
        polygon_at_point=np.array([-1, -1]),
    )
    result = remap_graph_connections(connected, reindex)
    assert result.point_at_polygon[0] == -1
    assert result.edge_at_edge[1] == -1
    assert_array_equal(result.polygon_at_point, [-1, -1])


def test_sort_connected_graphs_primary_is_sorted():
    # Primary: triangle with points in scrambled order.
    # Secondary has no polygons or edges, just one corner ↔ the one primary polygon.
    primary = SortGraph(
        xy_of_point=np.array([[0.0, 1.0], [0.0, 0.0], [1.0, 0.0]]),
        points_at_edge=np.array([[0, 1], [1, 2], [0, 2]]),
        edges_at_polygon=np.array([[0, 1, 2]]),
    )
    secondary = SortGraph(
        xy_of_point=np.array([[0.33, 0.33]]),
        points_at_edge=np.zeros((0, 2), dtype=int),
        edges_at_polygon=np.zeros((0, 0), dtype=int),
    )
    connected = SortConnections(
        point_at_polygon=np.zeros(0, dtype=int),
        edge_at_edge=np.zeros(0, dtype=int),
        polygon_at_point=np.array([0]),
    )
    sorted_primary, _, _ = sort_connected_graphs(primary, connected, secondary)
    xy = sorted_primary.xy_of_point
    assert np.all(xy[1:, 1] >= xy[:-1, 1])


def test_sort_connected_graphs_connection_consistency():
    # Primary: triangle with scrambled points.
    # Secondary: one corner ↔ the one primary polygon (patch_at_corner).
    # After sorting the primary, polygon_at_point must refer to the same
    # physical polygon (still polygon 0, since there is only one).
    #
    # Primary points (scrambled): (1,0)=0, (0,0)=1, (0.5,1)=2
    primary = SortGraph(
        xy_of_point=np.array([[1.0, 0.0], [0.0, 0.0], [0.5, 1.0]]),
        points_at_edge=np.array([[0, 1], [1, 2], [0, 2]]),
        edges_at_polygon=np.array([[0, 1, 2]]),
    )
    secondary = SortGraph(
        xy_of_point=np.array([[0.5, 0.33]]),
        points_at_edge=np.zeros((0, 2), dtype=int),
        edges_at_polygon=np.zeros((0, 0), dtype=int),
    )
    connected = SortConnections(
        point_at_polygon=np.zeros(0, dtype=int),
        edge_at_edge=np.zeros(0, dtype=int),
        polygon_at_point=np.array([0]),
    )
    sorted_primary, sorted_conn, sorted_secondary = sort_connected_graphs(
        primary, connected, secondary
    )
    # Primary points should be sorted (y then x): (0,0), (1,0), (0.5,1)
    assert_allclose(sorted_primary.xy_of_point[0], [0.0, 0.0])
    assert_allclose(sorted_primary.xy_of_point[1], [1.0, 0.0])
    # Only 1 primary polygon → polygon_at_point[0] is still 0
    assert sorted_conn.polygon_at_point[0] == 0
    # Secondary point (the corner) should remain at the same physical location
    assert_allclose(sorted_secondary.xy_of_point[0], [0.5, 0.33])
