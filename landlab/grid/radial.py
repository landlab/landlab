#! /usr/env/python
"""
Python implementation of RadialModelGrid, a grid class used to create and
manage structured Voronoi-Delaunay grids for 2D numerical models.

Do NOT add new documentation here. Grid documentation is now built in a semi-
automated fashion. To modify the text seen on the web, edit the files
`docs/text_for_[gridfile].py.txt`.
"""
from __future__ import absolute_import

from warnings import warn

import numpy as np

from landlab.utils.decorators import deprecated

from ..graph import DualRadialGraph
from .base import ModelGrid


class RadialModelGrid(DualRadialGraph, ModelGrid):

    """Grid of concentric circles.

    This inherited class implements a circular grid in which grid nodes are
    placed at regular radial and semi-regular arc-wise intervals. That is,
    if the radial spacing between *shells* is *dr*, the nodes are placed around
    the circular shell at regular intervals that get as close as possible to
    *dr*. The points are then arranged in a Delaunay triangulation with Voronoi
    cells. Within each ring, nodes are numbered according to Landlab
    convention, from the first node counterclockwise of east. Numbering
    begins at the centermost node and works outwards through the rings.

    Parameters
    ----------
    num_shells : int
        Number of rings in the grid.
    dr : float, optional
        Radial interval for rings.
    xy_of_center : tuple, optional
        (x, y) coordinate of center point. Default
        is (0., 0.)
    xy_of_reference : tuple, optional
        Coordinate value in projected space of the reference point,
        `xy_of_lower_left`. Default is (0., 0.)

    Returns
    -------
    RadialModelGrid
        A newly-created grid.

    Examples
    --------
    A grid with just one ring will have a node at the origin surrounded
    by six other nodes.

    >>> from landlab import RadialModelGrid
    >>> omg = RadialModelGrid(num_shells=1, dr=1., xy_of_center=(0., 0.))
    >>> omg.number_of_nodes
    7
    >>> omg.number_of_cells
    1

    A second rings will have 13 nodes.

    >>> omg = RadialModelGrid(2)
    >>> omg.number_of_nodes
    20

    >>> omg.radius_at_node
    array([ 2.,  2.,  2.,  2.,  1.,  2.,  1.,  1.,  2.,  2.,  0.,  2.,  1.,
            1.,  2.,  1.,  2.,  2.,  2.,  2.])
    """

    def __init__(self, num_shells=0, dr=1.0, xy_of_center=(0.0, 0.0), **kwds):
        """Create a circular grid.

        Create a circular grid in which grid nodes are placed at regular
        radial and semi-regular arc-wise intervals. That is, if the radial
        spacing between *shells* is *dr*, the nodes are placed around the
        circular shell at regular intervals that get as close as possible to
        *dr*.  The points are then arranged in a Delaunay triangulation with
        Voronoi cells.

        Parameters
        ----------
        num_shells : int
            Number of rings in the grid.
        dr : float, optional
            Radial interval for rings.
        xy_of_center : tuple, optional
            (x, y) coordinate of center point. Default
            is (0., 0.)
        xy_of_reference : tuple, optional
            Coordinate value in projected space of the reference point,
            `xy_of_lower_left`. Default is (0., 0.)

        Returns
        -------
        RadialModelGrid
            A newly-created grid.

        Examples
        --------
        A grid with just one ring will have a node at the origin surrounded
        by six other nodes.

        >>> from landlab import RadialModelGrid
        >>> omg = RadialModelGrid(num_shells=1, dr=1., xy_of_center=(0., 0.))
        >>> omg.number_of_nodes
        7
        >>> omg.number_of_cells
        1

        A grid with two rings will have 13 nodes.

        >>> omg = RadialModelGrid(2)
        >>> omg.number_of_nodes
        19
        """
        # shape = (num_shells + 1, 6)
        # spacing = dr
        # origin = origin_y, origin_x

        if "xy_of_center" not in kwds:
            DualRadialGraph.__init__(
                self,
                (num_shells, int(2. * np.pi / dr)),
                **kwds,
                xy_of_center=(origin_x, origin_y)
            )
        else:
            DualRadialGraph.__init__(self, (num_shells, int(2. * np.pi / dr)), **kwds)
        ModelGrid.__init__(self, **kwds)

        # DualRadialGraph.__init__(self, shape, spacing=spacing, origin=origin)
        # ModelGrid.__init__(self, **kwds)

        self._node_status = np.full(self.number_of_nodes,
                                    self.BC_NODE_IS_CORE, dtype=np.uint8)
        self._node_status[self.perimeter_nodes] = self.BC_NODE_IS_FIXED_VALUE
        # Set number of nodes, and initialize if caller has given dimensions
        if "origin_x" in kwds:
            msg = "The origin_x keyword is deprecated. Use xy_of_center."
            warn(msg, DeprecationWarning)

        if "origin_y" in kwds:
            msg = "The origin_y keyword is deprecated. Use xy_of_center."
            warn(msg, DeprecationWarning)

        # return cls(num_shells=num_shells, dr=dr, origin_x=origin[0], origin_y=origin[1])
        xy_of_center = (
            kwds.get("origin_x", xy_of_center[0]),
            kwds.get("origin_y", xy_of_center[1]),
        )
        xy_of_center = tuple(xy_of_center)

        if num_shells > 0:
            self._initialize(num_shells, dr, xy_of_center)
        super(RadialModelGrid, self).__init__(**kwds)
        self._xy_of_center = xy_of_center
