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

    >>> graph = DualFramedVoronoiGraph((3, 3), node_layout='rect', random_seed=False)
    >>> graph.number_of_nodes
    9

    >>> graph.x_of_node[2:5]    # doctest: +NORMALIZE_WHITESPACE
    array([ 2.        ,  0.        ,  1.07341707])
    >>> graph.y_of_node[2:6]    # doctest: +NORMALIZE_WHITESPACE
    array([ 0.        ,  0.749     ,  1.03337157,  1.251     ])
    """

    def __init__(
        self,
        shape,
        xy_spacing=(1.0, 1.0),
        xy_of_lower_left=(0.0, 0.0),
        orientation="horizontal",
        node_layout="rect",
        sort=False,
        xy_min_spacing=(0.5, 0.5),
        random_seed=True,
        seed=(200, 500),
    ):
        """Create the graph.

        Parameters
        ----------
        shape : int or tuple of int
            For a rectangular layout, number of rows and columns of nodes.
            If int, rows number = columns number = value
        xy_spacing : float or tuple of float, optional
            Node spacing along x and y coordinates. If float, same spacing at x and y.
        xy_of_lower_left : tuple, optional
            Minimum x-of-node and y-of-node values. Depending on the grid
            no node may be present at this coordinate. Default is (0., 0.).
        orientation : string, optional
            'horizontal' only
        node_layout : string, optional
            'rect' only. The grid layout of nodes.
        sort: bool
            If true, nodes, links and patches are re-numbered according certain criterias of position
        xy_min_spacing: float or tuple of float, optional
            Final minimal spacing between nodes. Random moves of the core nodes
            around their position cannot be above this threshold:
            (xy_spacing - xy_min_spacing) /2
            If float, same minimal spacing for x and y.
        random_seed: bool, optional
            If True, the moves of coordinates are completely random. False is used
            when reproducibility of moves is needed. Move pseudo-randomness is
            then controlled by the parameter seed.
        seed: int or tuple of int, optional
            Seeds used to generate the random x and y moves. This parameter is unused
            when random_seed = True.

        Returns
        -------
        DualFramedVoronoiGraph
            A newly-created graph.

        Examples
        --------
        Create a grid with 2 rows and 3 columns of nodes.

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
            orientation=orientation,
            node_layout=node_layout,
            sort=True,
            xy_min_spacing=xy_min_spacing,
            random_seed=random_seed,
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