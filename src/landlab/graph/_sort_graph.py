"""Sorting machinery for connected graph pairs.

This module provides utilities for sorting a *primary* graph geometrically
and propagating that ordering to a *secondary* (dual) graph that is connected
to it through one-to-one correspondence arrays.

Typical usage
-------------
Build a :class:`GraphTopology` for each graph and a :class:`GraphConnections`
describing the correspondences between them, then call
:func:`sort_connected_graphs`::

    primary = GraphTopology(
        xy_of_point=...,
        points_at_edge=...,
        edges_at_polygon=...,
    )
    secondary = GraphTopology(
        xy_of_point=...,
        points_at_edge=...,
        edges_at_polygon=...,
    )
    connections = GraphConnections(
        point_at_polygon=...,  # primary point id for each secondary polygon
        edge_at_edge=...,  # primary edge id for each secondary edge
        polygon_at_point=...,  # primary polygon id for each secondary point
    )
    sorted_primary, sorted_connections, sorted_secondary = sort_connected_graphs(
        primary, connections, secondary
    )

The primary graph is sorted by point coordinates (y then x by default).
Edges are then sorted by their midpoint coordinates and polygons by their
centroid coordinates.  The connection arrays are remapped into the new
primary indexing, and the secondary graph is sorted so that its elements
correspond to the already-sorted primary elements.

Relationship to Landlab graph elements
---------------------------------------
This module uses generic terminology (points, edges, polygons) that maps
onto Landlab's dual-graph vocabulary as follows:

+--------------------+-----------------------+-----------------------+
| Generic            | Primary graph         | Dual graph            |
+====================+=======================+=======================+
| point              | node                  | corner                |
+--------------------+-----------------------+-----------------------+
| edge               | link                  | face                  |
+--------------------+-----------------------+-----------------------+
| polygon            | patch                 | cell                  |
+--------------------+-----------------------+-----------------------+

The :class:`GraphConnections` arrays bridge the two graphs:

* ``point_at_polygon`` — node id for each cell (``node_at_cell``)
* ``edge_at_edge``     — link id for each face (``link_at_face``)
* ``polygon_at_point`` — patch id for each corner (``patch_at_corner``)

Examples
--------
Sort a triangle (primary) whose points are stored out of order. The
secondary graph has one corner at the triangle's centroid; no edge or
polygon connections are present in this minimal example.

>>> import numpy as np
>>> primary = GraphTopology(
...     xy_of_point=np.array([[0.0, 1.0], [0.0, 0.0], [1.0, 0.0]]),
...     points_at_edge=np.array([[0, 1], [1, 2], [0, 2]]),
...     edges_at_polygon=np.array([[0, 1, 2]]),
... )
>>> secondary = GraphTopology(
...     xy_of_point=np.array([[1 / 3, 1 / 3]]),
...     points_at_edge=np.zeros((0, 2), dtype=int),
...     edges_at_polygon=np.zeros((0, 0), dtype=int),
... )
>>> connections = GraphConnections(
...     point_at_polygon=np.zeros(0, dtype=int),
...     edge_at_edge=np.zeros(0, dtype=int),
...     polygon_at_point=np.array([0]),
... )
>>> sorted_primary, _, _ = sort_connected_graphs(primary, connections, secondary)
>>> sorted_primary.xy_of_point
array([[0., 0.],
       [1., 0.],
       [0., 1.]])
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from requireit import require_array

from landlab.core.utils import require_id_array
from landlab.graph._reindexer import Reindexer


@dataclass
class GraphTopology:
    """A graph made of points, lines, and polygons.

    Parameters
    ----------
    xy_of_point : ndarray of float, shape (n_points, 2)
        Coordinates of graph points.
    points_at_edge : ndarray of int, shape (n_edges, 2)
        Points for each edge.
    edges_at_polygon : ndarray of int, shape (n_polygons, max_polygon_edges)
        Edges for each polygon.
    """

    xy_of_point: NDArray[np.floating]
    points_at_edge: NDArray[np.integer]
    edges_at_polygon: NDArray[np.integer]


@dataclass(frozen=True, slots=True)
class GraphConnections:
    """One-to-one correspondences between two connected graphs.

    Parameters
    ----------
    point_at_polygon : ndarray of int, shape (n_polygons,)
        Primary-point id for each secondary polygon.
    edge_at_edge : ndarray of int, shape (n_edges,)
        Primary-edge id for each secondary edge.
    polygon_at_point : ndarray of int, shape (n_points,)
        Primary-polygon id for each secondary point.
    """

    point_at_polygon: NDArray[np.integer]
    edge_at_edge: NDArray[np.integer]
    polygon_at_point: NDArray[np.integer]


@dataclass(frozen=True, slots=True)
class GraphReindexers:
    """Reindexers for graph points, edges, and polygons."""

    point: Reindexer
    edge: Reindexer
    polygon: Reindexer


def sort_connected_graphs(
    primary: GraphTopology,
    connections: GraphConnections,
    secondary: GraphTopology,
    *,
    axes: tuple[int, ...] | None = None,
) -> tuple[GraphTopology, GraphConnections, GraphTopology]:
    """Sort two connected graphs.

    The primary graph is sorted geometrically from point coordinates. The
    connecting arrays are then remapped into the new primary indexing.
    The secondary graph is sorted to preserve its correspondence with the
    already-sorted primary graph.
    """
    sorted_primary, reindex = sort_graph_by_xy(primary, axes=axes)
    sorted_connections = remap_graph_connections(connections, reindex)

    sorted_secondary, sorted_connections = sort_graph_by_ids(
        secondary, sorted_connections
    )

    return sorted_primary, sorted_connections, sorted_secondary


def sort_graph_by_xy(
    graph: GraphTopology,
    axes: tuple[int, ...] | None = None,
) -> tuple[GraphTopology, GraphReindexers]:
    """Sort a graph from point coordinates.

    Points are sorted by ``graph.xy_of_point``. Edges are then sorted by their
    derived coordinates, and polygons are sorted by their derived coordinates.
    """
    _validate_graph_topology(graph, name="graph")

    reindex_points = Reindexer.from_xy(graph.xy_of_point, axes=axes)
    xy_of_point = reindex_points.reorder(graph.xy_of_point)
    points_at_edge = reindex_points.remap(graph.points_at_edge, fill_value=-1)

    xy_of_edge = calc_xy_of_edge(xy_of_point, points_at_edge)
    reindex_edges = Reindexer.from_xy(xy_of_edge, axes=axes)
    xy_of_edge = reindex_edges.reorder(xy_of_edge)
    points_at_edge = reindex_edges.reorder(points_at_edge)
    edges_at_polygon = reindex_edges.remap(graph.edges_at_polygon, fill_value=-1)

    xy_of_polygon = calc_xy_of_polygon(xy_of_edge, edges_at_polygon)
    reindex_polygons = Reindexer.from_xy(xy_of_polygon, axes=axes)
    edges_at_polygon = reindex_polygons.reorder(edges_at_polygon)

    sorted_graph = GraphTopology(
        xy_of_point=xy_of_point,
        points_at_edge=points_at_edge,
        edges_at_polygon=edges_at_polygon,
    )
    reindex = GraphReindexers(
        point=reindex_points,
        edge=reindex_edges,
        polygon=reindex_polygons,
    )

    return sorted_graph, reindex


def sort_graph_by_ids(
    secondary: GraphTopology, connections: GraphConnections
) -> tuple[GraphTopology, GraphConnections]:
    """Sort a graph to preserve correspondence with another graph.

    The graph is sorted from ids supplied by ``connections``:

    * polygons from ``point_at_polygon``
    * edges from ``edge_at_edge``
    * points from ``polygon_at_point``

    This is intended for the secondary graph in a connected graph pair, after
    the connecting arrays have already been remapped into the sorted primary
    indexing.
    """
    _validate_graph_topology(secondary, name="secondary")
    _validate_graph_connections(connections, name="connections")

    reindex_polygons = Reindexer.from_ids(connections.point_at_polygon)
    point_at_polygon = reindex_polygons.reorder(connections.point_at_polygon)
    edges_at_polygon = reindex_polygons.reorder(secondary.edges_at_polygon)

    reindex_edges = Reindexer.from_ids(connections.edge_at_edge)
    edge_at_edge = reindex_edges.reorder(connections.edge_at_edge)
    points_at_edge = reindex_edges.reorder(secondary.points_at_edge)
    edges_at_polygon = reindex_edges.remap(edges_at_polygon, fill_value=-1)

    reindex_points = Reindexer.from_ids(connections.polygon_at_point)
    polygon_at_point = reindex_points.reorder(connections.polygon_at_point)
    xy_of_point = reindex_points.reorder(secondary.xy_of_point)
    points_at_edge = reindex_points.remap(points_at_edge, fill_value=-1)

    sorted_connections = GraphConnections(
        point_at_polygon=point_at_polygon,
        edge_at_edge=edge_at_edge,
        polygon_at_point=polygon_at_point,
    )
    sorted_secondary = GraphTopology(
        xy_of_point=xy_of_point,
        points_at_edge=points_at_edge,
        edges_at_polygon=edges_at_polygon,
    )

    return sorted_secondary, sorted_connections


def remap_graph_connections(
    connections: GraphConnections,
    reindex: GraphReindexers,
) -> GraphConnections:
    """Remap graph-connection arrays into a sorted primary indexing."""
    _validate_graph_connections(connections, name="connections")

    point_at_polygon = reindex.point.remap(connections.point_at_polygon, fill_value=-1)
    edge_at_edge = reindex.edge.remap(connections.edge_at_edge, fill_value=-1)
    polygon_at_point = reindex.polygon.remap(
        connections.polygon_at_point, fill_value=-1
    )

    return GraphConnections(
        point_at_polygon=point_at_polygon,
        edge_at_edge=edge_at_edge,
        polygon_at_point=polygon_at_point,
    )


def calc_xy_of_edge(
    xy_of_point: NDArray[np.floating],
    points_at_edge: NDArray[np.integer],
) -> NDArray[np.floating]:
    """Compute edge midpoints."""
    return xy_of_point[points_at_edge].mean(axis=1)


def calc_xy_of_polygon(
    xy_of_edge: NDArray[np.floating],
    edges_at_polygon: NDArray[np.integer],
) -> NDArray[np.floating]:
    """Compute polygon 'centers'."""
    is_edge = edges_at_polygon != -1

    edges = np.where(is_edge, edges_at_polygon, 0)
    coords = xy_of_edge[edges]

    coords[~is_edge, :] = 0

    edges_per_polygon = is_edge.sum(axis=1)
    return coords.sum(axis=1) / edges_per_polygon[:, None]


def redirect_edges(
    xy_of_point: NDArray[np.floating],
    points_at_edge: NDArray[np.integer],
    *,
    out: NDArray[np.integer] | None = None,
) -> NDArray[np.integer]:
    """Direct edges toward the upper-right half-space.

    Parameters
    ----------
    xy_of_point : ndarray of float, shape (n_points, 2)
        Coordinates of points.
    points_at_edge : ndarray of int, shape (n_edges, 2)
        Tail and head point for each edge.
    out : ndarray of int, optional
        Output buffer. Must have the same shape as `points_at_edge`.

    Returns
    -------
    ndarray of int
        Edge endpoints with each edge directed so that its vector angle,
        measured from the positive x-axis, lies in ``[-45°, 135°)``.

    Examples
    --------
    >>> xy_of_point = [
    ...     [0, 0],
    ...     [1, 0],
    ...     [1, 1],
    ...     [0, 1],
    ...     [-1, 1],
    ...     [-1, 0],
    ...     [-1, 1],
    ...     [0, -1],
    ...     [1, -1],
    ... ]
    >>> points_at_edge = [
    ...     [0, 1],
    ...     [0, 2],
    ...     [0, 3],
    ...     [0, 4],
    ...     [0, 5],
    ...     [0, 6],
    ...     [0, 7],
    ...     [0, 8],
    ... ]
    >>> redirect_edges(xy_of_point, points_at_edge)
    array([[0, 1],
           [0, 2],
           [0, 3],
           [4, 0],
           [5, 0],
           [6, 0],
           [7, 0],
           [0, 8]])
    """
    points_at_edge = np.asarray(points_at_edge)
    xy_of_point = np.asarray(xy_of_point)

    if points_at_edge.ndim != 2 or points_at_edge.shape[1] != 2:
        raise ValueError("points_at_edge must have shape (n_edges, 2)")
    if xy_of_point.ndim != 2 or xy_of_point.shape[1] != 2:
        raise ValueError("xy_of_point must have shape (n_points, 2)")

    if out is None:
        out = np.empty_like(points_at_edge)
    elif out.shape != points_at_edge.shape:
        raise ValueError("out must have the same shape as points_at_edge")

    xy = xy_of_point[points_at_edge]
    dxy = xy[:, 1] - xy[:, 0]

    dx, dy = dxy[:, 0], dxy[:, 1]

    eps = 1e-12
    dot_prod = dx + dy
    is_positive = dot_prod > eps
    is_about_zero = ~(is_positive | (dot_prod < -eps))

    keep = is_positive | (is_about_zero & (dx > 0))
    flip = ~keep

    out[keep] = points_at_edge[keep]
    out[flip, 0] = points_at_edge[flip, 1]
    out[flip, 1] = points_at_edge[flip, 0]

    return out


def _validate_graph_topology(
    graph: GraphTopology,
    name: str | None = None,
) -> GraphTopology:
    name = name or "graph"

    require_array(
        graph.xy_of_point,
        shape=("n_points", 2),
        dtype=np.floating,
        name=f"{name}.xy_of_point",
    )
    n_points = len(graph.xy_of_point)
    if not np.all(np.isfinite(graph.xy_of_point)):
        raise ValueError(f"{name}.xy_of_point must be finite")

    require_id_array(
        graph.points_at_edge,
        shape=("n_edges", 2),
        max_id=n_points - 1,
        bad_value=None,
        name=f"{name}.points_at_edge",
    )
    n_edges = len(graph.points_at_edge)

    require_id_array(
        graph.edges_at_polygon,
        shape=("n_polygons", None),
        max_id=n_edges - 1,
        bad_value=-1,
        name=f"{name}.edges_at_polygon",
    )
    edges_per_polygon = np.sum(graph.edges_at_polygon >= 0, axis=1)
    if np.any(edges_per_polygon < 3):
        raise ValueError(
            f"{name}.edges_at_polygon must have at least three edges for every polygon"
        )

    return graph


def _validate_graph_connections(
    connections: GraphConnections,
    name: str | None = None,
) -> GraphConnections:
    name = name or "connections"

    for attr in ("point_at_polygon", "edge_at_edge", "polygon_at_point"):
        require_id_array(
            getattr(connections, attr),
            shape=("n",),
            bad_value=None,
            name=f"{name}.{attr}",
        )

    return connections
