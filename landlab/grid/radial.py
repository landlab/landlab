#! /usr/bin/env python

import numpy

from landlab.grid.voronoi import VoronoiDelaunayGrid


class RadialModelGrid(VoronoiDelaunayGrid):
    """
    This inherited class implements a circular grid in which grid nodes are
    placed at regular radial and semi-regular arc-wise intervals. That is, if
    the radial spacing between "shells" is dr, the nodes are placed around the
    circular shell at regular intervals that get as close as possible to dr.
    The points are then arranged in a Delaunay triangulation with Voronoi cells.
    
   # Examples:
   #     
   #     >>> hmg = HexModelGrid(3, 2, 1.0)
   #     >>> hmg.num_nodes
   #     7
    """
    
    def __init__(self, num_shells=0, dr=1.0, origin_x=0.0, origin_y=0.0):
        """
        Optionally takes number of shells and radial distance between shells. 
        If this are given, calls initialize() to set up the grid.
        
        """
        #print 'RadialModelGrid.init'
        
        # Set number of nodes, and initialize if caller has given dimensions
        #self.num_nodes = num_rows * num_cols
        if num_shells > 0:
            self.initialize(num_shells, dr, origin_x, origin_y)


    def initialize( self, num_shells, dr, origin_x=0.0, origin_y=0.0):
        """
        """
        #if self.DEBUG_TRACK_METHODS:
        print 'RadialModelGrid.initialize('+str(num_shells)+', '+str(dr)+')'
        
        [pts, npts] = self.make_radial_points(num_shells, dr)
        super(RadialModelGrid, self).initialize(pts[:,0], pts[:,1])
        
        
    def make_radial_points(self, num_shells, dr, origin_x=0.0, origin_y=0.0):
        """
        Creates and returns a set of (x,y) points placed in a series of
        concentric circles around the origin.
        """
        shells = numpy.arange(0, num_shells)+1
        twopi = 2*numpy.pi
        n_pts_in_shell = numpy.round( twopi*shells ) # number of points in each shell
        dtheta = twopi / n_pts_in_shell
        npts = int(sum(n_pts_in_shell)+1)
        pts = numpy.zeros((npts,2))
        r =  shells*dr
        startpt = 1
        for i in numpy.arange(0, num_shells):
            theta = dtheta[i]*numpy.arange(0, n_pts_in_shell[i])+dtheta[i]/(i+1)
            pts[startpt:(startpt+n_pts_in_shell[i]),0] = r[i]*numpy.cos(theta)
            pts[startpt:(startpt+n_pts_in_shell[i]),1] = r[i]*numpy.sin(theta)
            startpt += n_pts_in_shell[i]
        pts[:,0] += origin_x
        pts[:,1] += origin_y
        
        return pts, npts
