#! /usr/bin/env python

import numpy

from landlab.grid.voronoi import VoronoiDelaunayGrid


class HexModelGrid(VoronoiDelaunayGrid):
    """A grid of hexagonal cells.

    This inherited class implements a regular 2D grid with hexagonal cells and
    triangular patches. It is a special type of VoronoiDelaunay grid in which
    the initial set of points is arranged in a triangular/hexagonal lattice.
    
    Examples
    --------
    >>> hmg = HexModelGrid(3, 2, 1.0)
    >>> hmg.number_of_nodes
    7
    """
    
    def __init__(self, base_num_rows=0, base_num_cols=0, dx=1.0, 
                 orientation='horizontal', **kwds):
        """Create a grid of hexagonal cells.

        Create a regular 2D grid with hexagonal cells and triangular patches.
        It is a special type of VoronoiDelaunay grid in which the initial set
        of points is arranged in a triangular/hexagonal lattice.

        Parameters
        ----------
        base_num_rows : int
            Number of rows of nodes in the left column.
        base_num_cols : int
            Number of nodes on the first row.
        dx : float, optional
            Node spacing.
        orientation : string, optional
            One of the 3 cardinal directions in the grid, either 'horizontal' 
            (default) or 'vertical'

        Returns
        -------
        HexModelGrid
            A newly-created grid.

        Examples
        --------
        Create a hex grid with 2 rows of nodes. The first and third rows will
        have 2 nodes, and the second nodes.

        >>> hmg = HexModelGrid(3, 2, 1.0)
        >>> hmg.number_of_nodes
        7
        """
        # Set number of nodes, and initialize if caller has given dimensions
        #self._num_nodes = num_rows * num_cols
        if base_num_rows * base_num_cols > 0:
            self._initialize(base_num_rows, base_num_cols, dx, orientation)
        super(HexModelGrid, self).__init__(**kwds)

    def _initialize(self, base_num_rows, base_num_cols, dx, orientation):
        """
        Sets up a hexagonal grid with cell spacing dx and
        (by default) regular boundaries (that is, all perimeter cells are
        boundaries and all interior cells are active).
        
        Parameters
        ----------
        base_num_rows : int
            Number of rows along left side of grid
        base_num_cols : int
            Number of columns along bottom side of grid
        dx : float
            Distance between nodes
        orientation : string
            Either 'horizontal' (default in __init__) or 'vertical'

        Returns
        -------
        (none)
        
        Creates/modifies
        ----------------
        Creates and initializes self._num_nodes and self._dx
        
        Notes
        -----
        To be consistent with unstructured grids, the hex grid is
        managed not as a 2D array but rather as a set of arrays that
        describe connectivity information between nodes, links, cells, faces,
        patches, corners, and junctions.
        
        'Horizontal' orientation means that one of the 3 axes of the grid is
        horizontal, whereas the other two are at 30 degree angles to the 
        horizontal, like:
            
            \/
           ----
            /\
            
        'Vertical' means that one axis is vertical, with the other
        two at 30 degree angles to the vertical, more like:
            
           \   |   /
             \ | /
             / | \
           /   |   \
           
        (of course, these keyboard characters don't represent the angles quite
        right)
        """
        if self._DEBUG_TRACK_METHODS:
            print 'HexModelGrid._initialize('+str(base_num_rows)+', ' \
                   +str(base_num_cols)+', '+str(dx)+')'
                   
        # Make sure the parameter *orientation* is correct
        assert (orientation[0].lower()=='h' or orientation[0].lower()=='v'), \
               'orientation must be either "horizontal" (default) or "vertical"'
        
        # Create a set of hexagonally arranged points. These will be our nodes.
        if orientation=='horizontal':
            [pts, self._num_nodes] = HexModelGrid.make_hex_points_horizontal(base_num_rows, base_num_cols, dx)
        else:
            [pts, self._num_nodes] = HexModelGrid.make_hex_points_vertical(base_num_rows, base_num_cols, dx)
        
        # Call the VoronoiDelaunayGrid constructor to triangulate/Voronoi
        # the nodes into a grid.
        super(HexModelGrid, self)._initialize(pts[:,0], pts[:,1])
        
        # Remember grid spacing
        self._dx = dx

    def _setup_cell_areas_array(self):
        """
        Creates and returns an array containing the surface areas of the 
        hexagonal (Voronoi) cells.
        
        These cells are perfect hexagons in which the apothem is dx/2. The
        formula for area is:
        
        .. math::
        
            A = 3 dx^2 / 2 \sqrt{3} \approx 0.866 dx^2
        """
        self._cell_areas = 0.8660254*self._dx**2 + numpy.zeros(self.number_of_cells)
        return self._cell_areas
                           
    @staticmethod
    def make_hex_points_horizontal(num_rows, base_num_cols, dxh):
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
            
            >>> [p, npt] = HexModelGrid.make_hex_points_horizontal(3, 2, 1.0)
            >>> npt
            7
            >>> p[1,:]
            array([ 1.,  0.])
            >>> p[:3,0]
            array([ 0. ,  1. , -0.5])
        """

        dxv = dxh * numpy.sqrt(3.) / 2.
        half_dxh = dxh / 2.

        if numpy.mod(num_rows, 2) == 0:  # even number of rows
            npts = num_rows * base_num_cols + (num_rows * num_rows) // 4
            #print 'even # rows, npts=', npts
        else:  # odd number of rows
            npts = num_rows * base_num_cols + ((num_rows - 1) // 2) * ((num_rows - 1) // 2)
            #print 'odd # rows, npts=', npts
        pts = numpy.zeros((npts, 2))
        middle_row = num_rows // 2
        extra_cols = 0
        xshift = 0.
        i = 0
        for r in range(num_rows):
            for c in range(base_num_cols + extra_cols):
                pts[i,0] = c * dxh + xshift
                pts[i,1] = r * dxv
                i += 1
            if r < middle_row:
                extra_cols += 1
            else:
                extra_cols -= 1
            xshift = - half_dxh * extra_cols
        
        return pts, npts


    @staticmethod
    def make_hex_points_vertical(base_num_rows, num_cols, dxv):
        """
        Creates and returns a set of (x,y) points in a staggered grid in which the 
        points represent the centers of regular hexagonal cells, and the points
        could be connected to form equilateral triangles. The overall shape of the
        lattice is hexagonal.
        
        Inputs: base_num_rows = number of columns in the left and right columns
                                (middle columns have more)
                num_cols = number of columns in lattice
                dxv = vertical and diagonal spacing between points
                
        Return: 2D numpy array containing point (x,y) coordinates, and total number
                of points.
                
        Example:
            
            >>> [p, npt] = HexModelGrid.make_hex_points_vertical(2, 3, 1.0)
            >>> npt
            7
            >>> p[1,:]
            array([ 0.,  1.])
            >>> p[:3,1]
            array([ 0. ,  1. , -0.5])
        """

        dxh = dxv * numpy.sqrt(3.) / 2.
        half_dxv = dxv / 2.

        if numpy.mod(num_cols, 2) == 0:  # even number of columns
            npts = base_num_rows * num_cols + (num_cols * num_cols) // 4
            #print 'even # cols, npts=', npts
        else:  # odd number of columns
            npts = base_num_rows * num_cols + ((num_cols - 1) // 2) * ((num_cols - 1) // 2)
            #print 'odd # columns, npts=', npts
        pts = numpy.zeros((npts, 2))
        middle_col = num_cols // 2
        extra_rows = 0
        yshift = 0.
        i = 0
        for c in range(num_cols):
            for r in range(base_num_rows + extra_rows):
                pts[i,1] = r * dxv + yshift
                pts[i,0] = c * dxh
                i += 1
            if c < middle_col:
                extra_rows += 1
            else:
                extra_rows -= 1
            yshift = - half_dxv * extra_rows
        
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
