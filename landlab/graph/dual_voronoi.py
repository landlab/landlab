import numpy as np
from scipy.spatial import Voronoi

from .graph import Graph
from .dual import DualGraphMixIn
from .voronoi import VoronoiGraph
from .voronoi_helpers import setup_voronoi_connectivity


class DualVoronoiGraph(VoronoiGraph, DualGraphMixIn):

    """Dual graph of a voronoi grid."""

    def __init__(self, nodes):
        """Create a voronoi grid.

        Parameters
        ----------
        nodes : tuple of array_like
            Coordinates of every node. First *y*, then *x*.

        Examples
        --------
        >>> from landlab.graph import DualVoronoiGraph
        >>> node_x = [0, 1, 2, 3,
        ...           0.2, 1.2, 2.2, 3.2,
        ...           0.4, 1.4, 2.4, 3.4]
        >>> node_y = [0, 0, 0, 0,
        ...           1, 1, 1, 1,
        ...           2, 2, 2, 2]
        >>> graph = DualVoronoiGraph((node_y, node_x))
        >>> graph.x_of_corner
        array([ 0.9,  1.7,  1.9,  0.7,  2.7,  2.7,  1.7,  2.5,  1.5,  0.7])
        >>> graph.y_of_corner # doctest: +NORMALIZE_WHITESPACE
        array([ 1.58,  1.42,  1.58,  1.42,  0.58,  1.42,  0.58,  0.42,  0.42,
                0.58])
        >>> graph.corners_at_face # doctest: +NORMALIZE_WHITESPACE
        array([[0, 1], [0, 3], [1, 6], [8, 9], [6, 8], [3, 9], [1, 2], [4, 5],
               [2, 5], [4, 7], [6, 7]])
        >>> graph.faces_at_corner # doctest: +NORMALIZE_WHITESPACE
        array([[ 0,  1, -1], [ 0,  2,  6], [ 6,  8, -1], [ 1,  5, -1],
               [ 7,  9, -1], [ 7,  8, -1], [ 2,  4, 10], [ 9, 10, -1],
               [ 3,  4, -1], [ 3,  5, -1]])
        >>> graph.node_at_cell
        array([5, 6])
        """
        super(DualVoronoiGraph, self).__init__(nodes)

        voronoi = Voronoi(list(zip(self.x_of_node, self.y_of_node)))

        (faces_at_cell,
         corners_at_face,
         xy_at_corner,
         node_at_cell) = setup_voronoi_connectivity(voronoi)

        node_y = xy_at_corner[:, 1]
        node_x = xy_at_corner[:, 0]

        self._dual = Graph((node_y, node_x), links=corners_at_face,
                           patches=faces_at_cell)
        self._node_at_cell = node_at_cell
