#! /usr/env/python
"""Python implementation of FramedVoronoiGrid, a grid class used to create and
manage unstructured Voronoi-Delaunay grids for 2D numerical models, with a structured
perimeter layout

Do NOT add new documentation here. Grid documentation is now built in a
semi- automated fashion. To modify the text seen on the web, edit the
files `docs/text_for_[gridfile].py.txt`.

@author sebastien lenard
@date 2022, Aug
"""

import numpy
from ..graph import DualFramedVoronoiGraph
from .base import ModelGrid


class FramedVoronoiGrid(DualFramedVoronoiGraph, ModelGrid):
    """A grid of Voronoi Delaunay cells with a structured perimeter layout.

    This inherited class implements a irregular 2D grid with Voronoi Delaunay cells and
    irregular patches. It is a special type of VoronoiDelaunay grid in which
    the initial set of points is arranged in a fixed lattice (e.g. like a rectangular
    raster grid) named here "layout" and the core points are then moved aroung their
    initial position by a random distance, lower than a certain threshold.

    Inheritance diagram
    *******************
                       ModelGrid                                                                       DelaunayGraph
    FramedVoronoiGrid /                                                                               /
                      |                        FramedVoronoiGraph (Layout: HorizontalRectVoronoiGraph)
                       DualFramedVoronoiGraph /
                                              \
                                               DualGraph
                                                        ~ use of static Graph.sort()
    Examples
    --------
    Create a grid with 2 rows and 3 columns of nodes.

    >>> from landlab import FramedVoronoiGrid
    >>> grid = FramedVoronoiGrid((3, 2), xy_spacing=1.0)
    >>> grid.number_of_nodes
    6

    >>> grid = FramedVoronoiGrid((4, 3), orientation="horizontal", node_layout="rect", xy_spacing=(10., 10.), xy_min_spacing=(5., 5.), random_seed=False, seed=(200, 500))
    >>> grid.status_at_node
    array([1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1], dtype=uint8)
    >>> grid.x_of_node[3:6]        # doctest: +NORMALIZE_WHITESPACE
    array([  0.        ,  10.73417072,  20.        ])
    >>> grid.y_of_node[0::3]       # doctest: +NORMALIZE_WHITESPACE
    array([  0.   ,   7.499,  17.499,  30.   ])

    >>> grid = FramedVoronoiGrid((3, 5), orientation="horizontal", node_layout="rect", xy_spacing=(10., 10.), xy_min_spacing=5., random_seed=True)
    >>> grid.boundary_nodes
    array([ 0,  1,  2,  3,  4,  5,  9, 10, 11, 12, 13, 14])
    """

    def __init__(
        self,
        shape,
        xy_spacing=(1.0, 1.0),
        xy_of_lower_left=(0.0, 0.0),
        orientation="horizontal",
        node_layout="rect",
        xy_min_spacing=(0.5, 0.5),
        random_seed=True,
        seed=(200, 500),
        xy_of_reference=(0.0, 0.0),
        xy_axis_name=("x", "y"),
        xy_axis_units="-",
    ):
        """Create a grid of voronoi cells with a structured perimeter.

        Create a irregular 2D grid with voronoi cells and triangular patches.
        It is a special type of VoronoiDelaunay grid in which the initial set
        of points is arranged in a regular lattice determined by the parameters:
        shape, xy_spacing, orientation, node_layout. The coordinates of
        the core points are then randomly moved while the perimeter points
        remaining fixed, in a way determined by the parameters: xy_min_spacing,
        random_seed and seed.

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
        xy_of_reference : tuple, optional
            Coordinate value in projected space of the reference point,
            `xy_of_lower_left`. Default is (0., 0.)
        xy_axis_name: tuple of str, optional
            x y axis names.
        xy_axis_units: str, optional
            x y axis units.

        Returns
        -------
        FramedVoronoiGrid
            A newly-created grid.

        Examples
        --------
        Create a grid with 2 rows and 3 columns of nodes.

        >>> from landlab import FramedVoronoiGrid
        >>> grid = FramedVoronoiGrid((3, 2), xy_spacing=1.0)
        >>> grid.number_of_nodes
        6
        """
        DualFramedVoronoiGraph.__init__(
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
        ModelGrid.__init__(
            self,
            xy_axis_name=xy_axis_name,
            xy_axis_units=xy_axis_units,
            xy_of_reference=xy_of_reference,
        )

        self._node_status = numpy.full(
            self.number_of_nodes, self.BC_NODE_IS_CORE, dtype=numpy.uint8
        )
        self._node_status[self.perimeter_nodes] = self.BC_NODE_IS_FIXED_VALUE

    @classmethod
    def from_dict(cls, kwds):
        args = ()
        return cls(*args, **kwds)