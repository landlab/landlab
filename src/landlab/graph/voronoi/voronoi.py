import numpy as np

from ..graph import Graph
from .voronoi_to_graph import VoronoiDelaunayToGraph


class DelaunayGraph(Graph):
    """Graph of a voronoi grid.

    Examples
    --------
    >>> from landlab.graph import DelaunayGraph
    """

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
        >>> from landlab.graph import DelaunayGraph
        >>> node_x = [0.0, 1.0, 2.0, 0.9, 1.9, 2.9]
        >>> node_y = [0, 0, 0, 2, 2, 2]
        >>> graph = DelaunayGraph((node_y, node_x), sort=True)
        >>> graph.x_of_node
        array([0. ,  1. ,  2. ,  0.9,  1.9,  2.9])
        >>> graph.y_of_node
        array([0.,  0.,  0.,  2.,  2.,  2.])
        >>> graph.nodes_at_link
        array([[0, 1], [1, 2],
               [0, 3], [1, 3], [1, 4], [2, 4], [2, 5],
               [3, 4], [4, 5]])
        >>> graph.links_at_node
        array([[ 0,  2, -1, -1], [ 1,  4,  3,  0], [ 6,  5,  1, -1],
               [ 7,  2,  3, -1], [ 8,  7,  4,  5], [ 8,  6, -1, -1]])
        >>> graph.links_at_patch
        array([[3, 2, 0], [5, 4, 1], [7, 3, 4], [8, 5, 6]])

        >>> graph.nodes_at_patch
        array([[3, 0, 1], [4, 1, 2], [4, 3, 1], [5, 4, 2]])
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
            sort=sort,
        )
