""" Implement the DualFramedVoronoiGraph

@author sebastien lenard
@date 2022, Aug
"""

from ..dual import DualGraph
from ..voronoi.dual_voronoi import DualVoronoiGraph
from .framed_voronoi import FramedVoronoiGraph


class DualFramedVoronoiGraph(DualGraph, FramedVoronoiGraph):
    """Graph of a unstructured grid of Voronoi Delaunay cells and
    irregular patches. It is a special type of VoronoiDelaunay graph in which
    the initial set of points is arranged in a fixed lattice (e.g. like a rectangular
    raster grid) named here "layout" and the core points are then moved aroung their
    initial position by a random distance, lower than a certain threshold.

    Examples
    --------
    >>> from landlab.graph import DualFramedVoronoiGraph

    >>> graph = DualFramedVoronoiGraph((3, 3), seed=200)
    >>> graph.number_of_nodes
    9

    >>> graph.x_of_node[2:4]
    array([2., 0.])
    >>> graph.y_of_node[2:4]
    array([0.   , 0.749])
    >>> graph.y_of_node[5]
    1.251
    """

    def __init__(
        self,
        shape,
        xy_spacing=(1.0, 1.0),
        xy_of_lower_left=(0.0, 0.0),
        sort=False,
        xy_min_spacing=(0.5, 0.5),
        seed=200,
    ):
        """Create the graph.

        Parameters
        ----------
        shape : tuple of int
            Number of rows and columns of nodes.
        xy_spacing : float or tuple of float, optional
            Node spacing along x and y coordinates. If float, same spacing at x and y.
        xy_of_lower_left : tuple, optional
            Minimum *x*-of-node and *y*-of-node values. Depending on the grid,
            there may not be a node present at this location.
        sort: bool
            If ``True``, nodes, links and patches are re-numbered according to
            their position.
        xy_min_spacing: float or tuple of float, optional
            Final minimal spacing between nodes. Random moves of the core nodes
            around their position cannot be above this threshold:
            ``(xy_spacing - xy_min_spacing) / 2``
            If ``float``, same minimal spacing for *x* and *y*.
        seed: int, optional
            Seed used to generate the random *x* and *y* moves. When set,
            controls a pseudo-randomness of moves to ensure reproducibility.
            When ``None``, seed is random and the moves of coordinates are
            completely random.

        Returns
        -------
        DualFramedVoronoiGraph
            A newly-created graph.

        Examples
        --------
        Create a grid with 3 rows and 2 columns of nodes.

        >>> from landlab.graph import DualFramedVoronoiGraph
        >>> graph = DualFramedVoronoiGraph((3, 2), xy_spacing=1.0)
        >>> graph.number_of_nodes
        6
        """
        FramedVoronoiGraph.__init__(
            self,
            shape,
            xy_spacing=xy_spacing,
            xy_of_lower_left=xy_of_lower_left,
            sort=True,
            xy_min_spacing=xy_min_spacing,
            seed=seed,
        )

        DualVoronoiGraph.__init__(
            self,
            (self._y_of_node, self._x_of_node),
            perimeter_links=self._perimeter_links,
            sort=False,
        )

        if sort:
            self.sort()
