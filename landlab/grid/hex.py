#! /usr/bin/env python

import numpy

from landlab.grid.voronoi import VoronoiDelaunayGrid


class HexModelGrid(VoronoiDelaunayGrid):
    """
    This inherited class implements a regular 2D grid with hexagonal cells and
    triangular patches. It is a special type of VoronoiDelaunay grid in which
    the initial set of points is arranged in a triangular/hexagonal lattice.
    
    Examples:
        
        >>> hmg = HexModelGrid(3, 2, 1.0)
        >>> hmg.num_nodes
        7
    """
    
    def __init__(self, num_rows=0, base_num_cols=0, dx=1.0):
        """
        Optionally takes numbers of rows and columns and cell size as
        inputs. If this are given, calls initialize() to set up the grid.
        
        """
        #print 'HexModelGrid.init'
        
        # Set number of nodes, and initialize if caller has given dimensions
        #self.num_nodes = num_rows * num_cols
        if num_rows*base_num_cols > 0:
            self.initialize( num_rows, base_num_cols, dx )


    def initialize( self, num_rows, base_num_cols, dx ):
        """
        Sets up a num_rows by num_cols grid with cell spacing dx and
        (by default) regular boundaries (that is, all perimeter cells are
        boundaries and all interior cells are active).

        To be consistent with unstructured grids, the hex grid is
        managed not as a 2D array but rather as a set of arrays that
        describe connectivity information between nodes, links, cells, faces,
        patches, corners, and junctions.
        """
        if self.DEBUG_TRACK_METHODS:
            print 'HexModelGrid.initialize('+str(num_rows)+', ' \
                   +str(base_num_cols)+', '+str(dx)+')'
        
        [pts, self.num_nodes] = self.make_hex_points(num_rows, base_num_cols, dx)
        super(HexModelGrid, self).initialize(pts[:,0], pts[:,1])
        

    def make_hex_points(self, num_rows, base_num_cols, dxh):
        """
        Creates and returns a set of (x,y) points in a staggered grid in which the 
        points represent the centers of regular hexagonal cells, and the points
        could be connected to form equilateral triangles. The overall shape of the
        lattice is hexagonal.
        
        Inputs: num_rows = number of rows in lattice
                base_num_cols = number of columns in the bottom and top rows
                                (middle rows have more)
                dxh = horizontal and diagonal spacing between points
                
        Return: 2D numpy array containing point (x,y) coordinates, and total number
                of points.
                
        Example:
            
            >>> hmg = HexModelGrid()
            >>> [p, npt] = hmg.make_hex_points(3, 2, 1.0)
            >>> npt
            7
            >>> p[1,:]
            array([ 1.,  0.])
            >>> p[:3,0]
            array([ 0. ,  1. , -0.5])
        """

        dxv = dxh*numpy.sqrt(3.)/2.
        half_dxh = dxh/2.

        if numpy.mod(num_rows, 2)==0:  # even number of rows
            npts = num_rows*base_num_cols+(num_rows*num_rows)/4
            #print 'even # rows, npts=', npts
        else:  # odd number of rows
            npts = num_rows*base_num_cols + ((num_rows-1)/2)*((num_rows-1)/2)
            #print 'odd # rows, npts=', npts
        pts = numpy.zeros((npts,2))
        middle_row = num_rows/2
        extra_cols = 0
        xshift = 0.
        i=0
        for r in range(num_rows):
            for c in range(base_num_cols+extra_cols):
                pts[i,0] = c*dxh+xshift
                pts[i,1] = r*dxv
                i += 1
            if r<middle_row:
                extra_cols += 1
            else:
                extra_cols -= 1
            xshift = -half_dxh*extra_cols
        
        return pts, npts


def from_dict(param_dict):
    """
    Create a HexModelGrid from the dictionary-like object, *param_dict*.
    Required keys of the dictionary are NUM_ROWS, NUM_COLS. Raises a KeyError
    if either of these are missing.  If GRID_SPACING is given, use it as the
    HexModelGrid *dx* parameter, otherwise default to unit spacing.
    """
    # Read and create a basic HexModelGrid
    try:
        n_rows = int(param_dict['NUM_ROWS'])
        n_cols = int(param_dict['NUM_COLS'])
        dx = float(param_dict.get('GRID_SPACING', 1.))
    except KeyError as e:
        raise
    except ValueError as e:
        raise
    else:
        hg = HexModelGrid(n_rows, n_cols, dx)
        
    return hg
