#! /usr/env/python
"""
Python implementation of RadialModelGrid, a grid class used to create and
manage structured Voronoi-Delaunay grids for 2D numerical models.

Do NOT add new documentation here. Grid documentation is now built in a semi-
automated fashion. To modify the text seen on the web, edit the files
`docs/text_for_[gridfile].py.txt`.
"""

import numpy
from six.moves import range

from .base import ModelGrid
from ..graph import DualRadialGraph
from landlab.utils.decorators import deprecated


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
    origin_x : float, optional
        x-coordinate of origin node.
    origin_y : float, optional
        y-coordinate of origin node.

    Returns
    -------
    RadialModelGrid
        A newly-created grid.

    Examples
    --------
    A grid with just one ring will have a node at the origin surrounded
    by six other nodes.

    >>> from landlab import RadialModelGrid
    >>> omg = RadialModelGrid(num_shells=1, dr=1., origin_x=0., origin_y=0.)
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

    def __init__(self, num_shells=0, dr=1.0, origin_x=0.0, origin_y=0.0,
                 **kwds):
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
        origin_x : float, optional
            x-coordinate of origin node.
        origin_y : float, optional
            y-coordinate of origin node.

        Returns
        -------
        RadialModelGrid
            A newly-created grid.

        Examples
        --------
        A grid with just one ring will have a node at the origin surrounded
        by six other nodes.

        >>> from landlab import RadialModelGrid
        >>> omg = RadialModelGrid(num_shells=1, dr=1., origin_x=0.,
        ...                       origin_y=0.)
        >>> omg.number_of_nodes
        7
        >>> omg.number_of_cells
        1

        A second rings will have 13 nodes.

        >>> omg = RadialModelGrid(2)
        >>> omg.number_of_nodes
        20
        """
        # Set number of nodes, and initialize if caller has given dimensions
        # self._origin_x = origin_x
        # self._origin_y = origin_y
        # if num_shells > 0:
        #     self._initialize(num_shells, dr, origin_x, origin_y)
        # super(RadialModelGrid, self).__init__(**kwds)
        shape = (num_shells + 1, 6)
        spacing = dr
        origin = origin_y, origin_x

        DualRadialGraph.__init__(self, shape, spacing=spacing, origin=origin)
        ModelGrid.__init__(self, **kwds)

        self._node_status = numpy.full(self.number_of_nodes,
                                       self.BC_NODE_IS_CORE, dtype=numpy.uint8)
        self._node_status[self.perimeter_nodes] = self.BC_NODE_IS_FIXED_VALUE

    @classmethod
    def from_dict(cls, params):
        """
        LLCATS: GINF
        """
        num_shells = params['num_shells']
        dr = params.get('dr', 1.)
        origin = params.get('origin', (0., 0.))

        return cls(num_shells=num_shells, dr=dr, origin_x=origin[0],
                   origin_y=origin[1])

    def _initialize(self, num_shells, dr, origin_x=0.0, origin_y=0.0):
        [pts, npts] = self._create_radial_points(num_shells, dr)
        self._n_shells = int(num_shells)
        self._dr = dr
        super(RadialModelGrid, self)._initialize(pts[:, 0], pts[:, 1])

    def _create_radial_points(self, num_shells, dr, origin_x=0.0,
                              origin_y=0.0):
        """Create a set of points on concentric circles.

        Creates and returns a set of (x,y) points placed in a series of
        concentric circles around the origin.
        """
        shells = numpy.arange(0, num_shells) + 1
        twopi = 2 * numpy.pi
        # number of points in each shell
        n_pts_in_shell = numpy.round(twopi * shells)
        dtheta = twopi / n_pts_in_shell
        npts = int(sum(n_pts_in_shell) + 1)
        pts = numpy.zeros((npts, 2))
        r = shells * dr
        startpt = 1
        for i in numpy.arange(0, num_shells):
            theta = (dtheta[i] * numpy.arange(0, n_pts_in_shell[i]) +
                     dtheta[i] / (i + 1))
            ycoord = r[i] * numpy.sin(theta)
            if numpy.isclose(ycoord[-1], 0.):
                # this modification necessary to force the first ring to
                # follow our new CCW from E numbering convention (DEJH, Nov15)
                ycoord[-1] = 0.
                pts[startpt:(startpt + int(n_pts_in_shell[i])),
                    0] = numpy.roll(r[i] * numpy.cos(theta), 1)
                pts[startpt:(startpt + int(n_pts_in_shell[i])),
                    1] = numpy.roll(ycoord, 1)
            else:
                pts[startpt:(startpt + int(n_pts_in_shell[i])),
                    0] = r[i] * numpy.cos(theta)
                pts[startpt:(startpt + int(n_pts_in_shell[i])),
                    1] = ycoord
            etartpt += int(n_pts_in_shell[i])
        pts[:, 0] += origin_x
        pts[:, 1] += origin_y

        return pts, npts
