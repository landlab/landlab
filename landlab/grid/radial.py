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

import numpy

from landlab.utils.decorators import deprecated

from .voronoi import VoronoiDelaunayGrid


class RadialModelGrid(VoronoiDelaunayGrid):

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
        20
        """
        # Set number of nodes, and initialize if caller has given dimensions
        if "origin_x" in kwds:
            msg = "The origin_x keyword is deprecated. Use xy_of_center."
            warn(msg, DeprecationWarning)

        if "origin_y" in kwds:
            msg = "The origin_y keyword is deprecated. Use xy_of_center."
            warn(msg, DeprecationWarning)

        xy_of_center = (
            kwds.get("origin_x", xy_of_center[0]),
            kwds.get("origin_y", xy_of_center[1]),
        )
        xy_of_center = tuple(xy_of_center)

        if num_shells > 0:
            self._initialize(num_shells, dr, xy_of_center)
        super(RadialModelGrid, self).__init__(**kwds)
        self._xy_of_center = xy_of_center

    @property
    def xy_of_center(self):
        """Return (x, y) of the reference point."""
        return self._xy_of_center

    @xy_of_center.setter
    def xy_of_center(self, xy_of_center):
        """Set a new value for the xy_of_lower_left."""
        dx = self._xy_of_center[0] - xy_of_center[0]
        dy = self._xy_of_center[1] - xy_of_center[1]
        self._xy_of_node -= (dx, dy)
        self._xy_of_center = xy_of_center

    def _initialize(self, num_shells, dr, xy_of_center):
        [pts, npts] = self._create_radial_points(num_shells, dr, xy_of_center)
        self._n_shells = int(num_shells)
        self._dr = dr
        super(RadialModelGrid, self)._initialize(pts[:, 0], pts[:, 1])

    def _create_radial_points(self, num_shells, dr, xy_of_center=(0.0, 0.0)):
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
            theta = (dtheta[i] * numpy.arange(0, n_pts_in_shell[i])) + dtheta[i] / (
                i + 1
            )
            ycoord = r[i] * numpy.sin(theta)
            if numpy.isclose(ycoord[-1], 0.0):
                # this modification necessary to force the first ring to
                # follow our new CCW from E numbering convention (DEJH, Nov15)
                ycoord[-1] = 0.0
                pts[startpt : (startpt + int(n_pts_in_shell[i])), 0] = numpy.roll(
                    r[i] * numpy.cos(theta), 1
                )
                pts[startpt : (startpt + int(n_pts_in_shell[i])), 1] = numpy.roll(
                    ycoord, 1
                )
            else:
                pts[startpt : (startpt + int(n_pts_in_shell[i])), 0] = r[i] * numpy.cos(
                    theta
                )
                pts[startpt : (startpt + int(n_pts_in_shell[i])), 1] = ycoord
            startpt += int(n_pts_in_shell[i])

        pts[:, 0] += xy_of_center[0]
        pts[:, 1] += xy_of_center[1]

        return pts, npts

    @property
    def number_of_shells(self):
        """Number of node shells in grid.

        Returns
        -------
        int
            The number of node shells in the radial grid (not counting the
            center node).

        LLCATS: GINF
        """
        return self._n_shells

    @property
    @deprecated(use="spacing_of_shells", version=1.0)
    def shell_spacing(self):
        """Fixed distance between shells.

        LLCATS: DEPR GINF MEAS
        """
        return self._dr

    @property
    def spacing_of_shells(self):
        """Fixed distance between shells.

        LLCATS: GINF MEAS
        """
        return self._dr

    @property
    def number_of_nodes_in_shell(self):
        """Number of nodes in each shell.

        Returns
        -------
        int
            Number of nodes in each shell, excluding the center node.

        LLCATS: GINF NINF
        """
        try:
            return self._nnodes_inshell
        except AttributeError:
            n_pts_in_shell = numpy.round(
                2.0
                * numpy.pi
                * (numpy.arange(self.number_of_shells, dtype=float) + 1.0)
            )
            self._nnodes_inshell = n_pts_in_shell.astype(int)
            return self._nnodes_inshell

    @property
    def radius_at_node(self):
        """Distance for center node to each node.

        Returns
        -------
        ndarray of float
            The distance from the center node of each node.

        >>> mg = RadialModelGrid(num_shells=2)
        >>> mg.radius_at_node
        array([ 2.,  2.,  2.,  2.,  2.,  1.,  1.,  2.,  0.,  1.,  1.,  2.,  2.,
                1.,  1.,  2.,  2.,  2.,  2.,  2.])

        LLCATS: NINF MEAS
        """
        try:
            return self._node_radii
        except AttributeError:
            self._node_radii = numpy.sqrt(
                numpy.square(self.node_x - self._xy_of_center[0])
                + numpy.square(self.node_y - self._xy_of_center[1])
            )
            return self._node_radii
