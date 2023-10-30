#! /usr/env/python
"""Python implementation of :class:`~.FramedVoronoiGrid`, a grid class used to create and
manage unstructured Voronoi-Delaunay grids for 2D numerical models, with a structured
perimeter layout

Do NOT add new documentation here. Grid documentation is now built in a
semi- automated fashion. To modify the text seen on the web, edit the
files `docs/text_for_[gridfile].py.txt`.

.. codeauthor:: sebastien lenard
"""

import numpy

from ..graph import DualFramedVoronoiGraph
from .base import ModelGrid


class FramedVoronoiGrid(DualFramedVoronoiGraph, ModelGrid):
    """A grid of Voronoi Delaunay cells with a structured perimeter layout.

    This inherited class implements a irregular 2D grid with Voronoi Delaunay cells and
    irregular patches. It is a special type of :class:`~.VoronoiDelaunayGrid` grid in which
    the initial set of points is arranged in a fixed lattice (e.g. like a
    :class:`~.RasterModelGrid`), named here "layout", and the core points are
    then moved a random distance from their initial positions, bounded by a user-supplied
    threshold.

    Examples
    --------
    Create a grid with 3 rows and 2 columns of nodes.

    >>> from landlab import FramedVoronoiGrid
    >>> grid = FramedVoronoiGrid((3, 2), xy_spacing=1.0)
    >>> grid.number_of_nodes
    6

    >>> grid = FramedVoronoiGrid(
    ...     (4, 3), xy_spacing=(10.0, 10.0), xy_min_spacing=(5.0, 5.0), seed=200
    ... )
    >>> grid.status_at_node.reshape(grid.shape)
    array([[1, 1, 1],
           [1, 0, 1],
           [1, 0, 1],
           [1, 1, 1]], dtype=uint8)
    >>> grid.x_of_node[3]
    0.0
    >>> grid.x_of_node[5]
    20.0
    >>> grid.y_of_node[0::3]
    array([  0.   ,   7.499,  17.499,  30.   ])

    >>> grid = FramedVoronoiGrid(
    ...     (3, 5), xy_spacing=(10.0, 10.0), xy_min_spacing=5.0, seed=None
    ... )
    >>> grid.boundary_nodes
    array([ 0,  1,  2,  3,  4,  5,  9, 10, 11, 12, 13, 14])
    """

    # Inheritance diagram:
    #
    # FramedVoronoiGrid
    # |-- ModelGrid
    # `-- DualFramedVoronoiGraph
    #     |-- FramedVoronoiGraph (Layout: HorizontalRectVoronoiGraph)
    #     |   `-- DelaunayGraph
    #     `-- DualGraph (use of static Graph.sort())

    def __init__(
        self,
        shape,
        xy_spacing=(1.0, 1.0),
        xy_of_lower_left=(0.0, 0.0),
        xy_min_spacing=(0.5, 0.5),
        seed=200,
        xy_of_reference=(0.0, 0.0),
        xy_axis_name=("x", "y"),
        xy_axis_units="-",
    ):
        """Create a grid of voronoi cells with a structured perimeter.

        Create an irregular 2D grid with voronoi cells and triangular patches.
        It is a special type of :class:`~.VoronoiDelaunayGrid` in which the initial set
        of points is arranged in a regular lattice determined by the parameters
        *shape*, and *xy_spacing*. The coordinates of
        the core points are then randomly moved while the perimeter points
        remaining fixed, in a way determined by the parameters *xy_min_spacing*, and
        *seed*.

        Parameters
        ----------
        shape : tuple of int
            Number of rows and columns of nodes.
        xy_spacing : float or tuple of float, optional
            Node spacing along x and y coordinates. If ``float``, same spacing at *x* and *y*.
        xy_of_lower_left : tuple, optional
            Minimum *x*-of-node and *y*-of-node values. Depending on the grid,
            there may not be a node at this coordinate.
        xy_min_spacing: float or tuple of float, optional
            Final minimal spacing between nodes. Random moves of the core nodes
            from their initial positions cannot be above this threshold:
            ``(xy_spacing - xy_min_spacing) / 2``
            If ``float``, same minimal spacing for *x* and *y*.
        seed: int, optional
            Seed used to generate the random *x* and *y* moves.
            When set, controls the pseudo-randomness of moves to ensure
            reproducibility.
            When ``None``, the seed is random and the moves of coordinates are
            completely random.
        xy_of_reference : tuple, optional
            Coordinate value in projected space of the reference point,
            *xy_of_lower_left*.
        xy_axis_name: tuple of str, optional
            *x* and *y* axis names.
        xy_axis_units: str, optional
            *x* and *y* axis units.

        Returns
        -------
        FramedVoronoiGrid
            A newly-created grid.

        Examples
        --------
        Create a grid with 3 rows and 2 columns of nodes.

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
            sort=True,
            xy_min_spacing=xy_min_spacing,
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
