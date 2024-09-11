import numpy as np

from ..dual import DualGraph
from ..graph import Graph
from .voronoi import DelaunayGraph
from .voronoi_to_graph import VoronoiDelaunayToGraph


class DualVoronoiGraph(DualGraph, DelaunayGraph):
    def __init__(
        self, node_y_and_x, max_node_spacing=None, sort=False, perimeter_links=None
    ):
        """Create a voronoi grid.

        Parameters
        ----------
        nodes : tuple of array_like
            Coordinates of every node. First *y*, then *x*.

        Examples
        --------
        >>> from landlab.graph import DualVoronoiGraph
        >>> node_x = [0, 1, 2, 3, 0.2, 1.2, 2.2, 3.2, 0.4, 1.4, 2.4, 3.4]
        >>> node_y = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2]
        >>> graph = DualVoronoiGraph((node_y, node_x), sort=True)
        >>> graph.x_of_corner
        array([0.5,  1.5,  2.5,  0.7,  1.7,  2.7,  0.7,  1.7,  2.7,  0.9,  1.9,
               2.9])
        >>> graph.y_of_corner
        array([0.42,  0.42,  0.42,  0.58,  0.58,  0.58,  1.42,  1.42,  1.42,
               1.58,  1.58,  1.58])
        >>> graph.corners_at_face
        array([[ 0,  3], [ 3,  1], [ 1,  4], [ 4,  2], [ 2,  5],
               [ 3,  6], [ 4,  7], [ 5,  8],
               [ 6,  9], [ 9,  7], [ 7, 10], [10,  8], [ 8, 11]])
        >>> graph.faces_at_corner
        array([[ 0, -1, -1], [ 2,  1, -1], [ 4,  3, -1],
               [ 5,  0,  1], [ 6,  2,  3], [ 7,  4, -1],
               [ 8,  5, -1], [10,  9,  6], [12, 11,  7],
               [ 8,  9, -1], [10, 11, -1], [12, -1, -1]])
        >>> graph.node_at_cell
        array([5, 6])
        """
        mesh = VoronoiDelaunayToGraph(
            np.vstack((node_y_and_x[1], node_y_and_x[0])).T,
            perimeter_links=perimeter_links,
        )

        Graph.__init__(
            self,
            node_y_and_x,
            links=mesh.nodes_at_link,
            patches=mesh.links_at_patch,
            sort=False,
        )
        dual_graph = Graph(
            (mesh.y_of_corner, mesh.x_of_corner),
            links=mesh.corners_at_face,
            patches=mesh.faces_at_cell,
            sort=False,
        )

        self.merge(
            dual_graph, node_at_cell=mesh.node_at_cell, nodes_at_face=mesh.nodes_at_face
        )

        if sort:
            self.sort()
