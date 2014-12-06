#! /usr/bin/env python

import numpy

from .voronoi import VoronoiDelaunayGrid


class RadialModelGrid(VoronoiDelaunayGrid):
    """
    This inherited class implements a circular grid in which grid nodes are
    placed at regular radial and semi-regular arc-wise intervals. That is, if
    the radial spacing between "shells" is dr, the nodes are placed around the
    circular shell at regular intervals that get as close as possible to dr.
    The points are then arranged in a Delaunay triangulation with Voronoi
    cells.
    
    Examples
    --------
    >>> from landlab import RadialModelGrid
    >>> omg = RadialModelGrid(num_shells=1, dr=1., origin_x=0., origin_y=0.)
    >>> omg.number_of_nodes
    7
    """
    
    def __init__(self, num_shells=0, dr=1.0, origin_x=0.0, origin_y=0.0, **kwds):
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
        >>> omg = RadialModelGrid(num_shells=1, dr=1., origin_x=0., origin_y=0.)
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
        #self._num_nodes = num_rows * num_cols
        if num_shells > 0:
            self._initialize(num_shells, dr, origin_x, origin_y)
        super(RadialModelGrid, self).__init__(**kwds)

    def _initialize( self, num_shells, dr, origin_x=0.0, origin_y=0.0):
        if self._DEBUG_TRACK_METHODS:
            print 'RadialModelGrid._initialize('+str(num_shells)+', '+str(dr)+')'
        
        [pts, npts] = self.make_radial_points(num_shells, dr)
        super(RadialModelGrid, self)._initialize(pts[:,0], pts[:,1])
        
        
    def make_radial_points(self, num_shells, dr, origin_x=0.0, origin_y=0.0):
        """
        Creates and returns a set of (x,y) points placed in a series of
        concentric circles around the origin.
        """
        shells = numpy.arange(0, num_shells) + 1
        twopi = 2 * numpy.pi
        n_pts_in_shell = numpy.round(twopi * shells) # number of points in each shell
        dtheta = twopi / n_pts_in_shell
        npts = int(sum(n_pts_in_shell) + 1)
        pts = numpy.zeros((npts, 2))
        r = shells * dr
        startpt = 1
        for i in numpy.arange(0, num_shells):
            theta = (dtheta[i] * numpy.arange(0, n_pts_in_shell[i]) +
                     dtheta[i] / (i + 1))
            pts[startpt:(startpt + int(n_pts_in_shell[i])), 0] = r[i] * numpy.cos(theta)
            pts[startpt:(startpt + int(n_pts_in_shell[i])), 1] = r[i] * numpy.sin(theta)
            startpt += int(n_pts_in_shell[i])
        pts[:,0] += origin_x
        pts[:,1] += origin_y
        
        return pts, npts
