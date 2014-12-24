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
    >>> from landlab import HexModelGrid
    >>> hmg = HexModelGrid(3, 2, 1.0)
    >>> hmg.number_of_nodes
    7
    """
    
    def __init__(self, base_num_rows=0, base_num_cols=0, dx=1.0, 
                 orientation='horizontal', shape='hex', reorient_links=False, 
                 **kwds):
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

        >>> from landlab import HexModelGrid
        >>> hmg = HexModelGrid(3, 2, 1.0)
        >>> hmg.number_of_nodes
        7
        """
        # Set number of nodes, and initialize if caller has given dimensions
        #self._num_nodes = num_rows * num_cols
        if base_num_rows * base_num_cols > 0:
            self._initialize(base_num_rows, base_num_cols, dx, orientation,
                             shape, reorient_links)
        super(HexModelGrid, self).__init__(**kwds)

    def _initialize(self, base_num_rows, base_num_cols, dx, orientation,
                    shape, reorient_links=False):
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
        shape : string
            Either 'hex' (default in __init__) or 'rect'
        reorient_links : bool
            Whether or not to re-orient all links to point between -45 deg
            and +135 deg clockwise from "north" (i.e., along y axis)

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
            
            \ /
           -----
            / \
            
        'Vertical' means that one axis is vertical, with the other
        two at 30 degree angles to the vertical, more like:
            
           \   |   /
             \ | /
             / | \
           /   |   \
           
        (of course, these keyboard characters don't represent the angles quite
        right)
        
        Numbers of rows and columns: a hex grid with a rectangular shape will
        have a fixed number of rows and columns, and so for rectangular shaped
        grids we record this information in self._nrows and self._ncols. With
        a hex-shaped grid, either the number of columns (if 'horizontal') or
        the number of rows (if 'vertical') will vary across the grid. Therefore,
        for hex-shaped grids we record only self._nrows for 'horizontal' grids,
        and only self._ncols for 'vertical' grids.
        """
        if self._DEBUG_TRACK_METHODS:
            print 'HexModelGrid._initialize('+str(base_num_rows)+', ' \
                   +str(base_num_cols)+', '+str(dx)+')'
                   
        # Make sure the parameter *orientation* is correct
        assert (orientation[0].lower()=='h' or orientation[0].lower()=='v'), \
               'orientation must be either "horizontal" (default) or "vertical"'
        
        # Make sure the parameter *shape* is correct
        assert (shape[0].lower()=='h' or shape[0].lower()=='r'), \
               'shape must be either "hex" (default) or "rect"'
               
        # Create a set of hexagonally arranged points. These will be our nodes.
        if orientation=='horizontal' and shape=='hex':
            [pts, self._num_nodes] = HexModelGrid.make_hex_points_horizontal_hex(base_num_rows, base_num_cols, dx)
            self.orientation = 'horizontal'
            self._nrows = base_num_rows
        elif orientation=='horizontal' and shape=='rect':
            [pts, self._num_nodes] = HexModelGrid.make_hex_points_horizontal_rect(base_num_rows, base_num_cols, dx)
            self.orientation = 'horizontal'
            self._nrows = base_num_rows
            self._ncols = base_num_cols
        elif orientation=='vertical' and shape=='hex':
            [pts, self._num_nodes] = HexModelGrid.make_hex_points_vertical_hex(base_num_rows, base_num_cols, dx)
            self.orientation = 'vertical'
            self._ncols = base_num_cols
        else:
            [pts, self._num_nodes] = HexModelGrid.make_hex_points_vertical_rect(base_num_rows, base_num_cols, dx)
            self.orientation = 'vertical'
            self._nrows = base_num_rows
            self._ncols = base_num_cols
        
        # Call the VoronoiDelaunayGrid constructor to triangulate/Voronoi
        # the nodes into a grid.
        super(HexModelGrid, self)._initialize(pts[:,0], pts[:,1], reorient_links)
        
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
    def make_hex_points_horizontal_hex(num_rows, base_num_cols, dxh):
        """
        Creates and returns a set of (x,y) points in a staggered grid in which the 
        points represent the centers of regular hexagonal cells, and the points
        could be connected to form equilateral triangles. The overall shape of the
        lattice is hexagonal, and one of the 3 axes is horizontal.
        
        Inputs: num_rows = number of rows in lattice
                base_num_cols = number of columns in the bottom and top rows
                                (middle rows have more)
                dxh = horizontal and diagonal spacing between points
                
        Return: 2D numpy array containing point (x,y) coordinates, and total number
                of points.
                
        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> [p, npt] = HexModelGrid.make_hex_points_horizontal_hex(3, 2, 1.0)
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
    def make_hex_points_horizontal_rect(num_rows, num_cols, dxh):
        """
        Creates and returns a set of (x,y) points in a staggered grid in which the 
        points represent the centers of regular hexagonal cells, and the points
        could be connected to form equilateral triangles. The overall shape of the
        lattice is rectangular, and one of the 3 axes is horizontal.
        
        Inputs: num_rows = number of rows in lattice
                num_cols = number of columns in lattice
                dxh = horizontal and diagonal spacing between points
                
        Return: 2D numpy array containing point (x,y) coordinates, and total number
                of points.
                
        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> [p, npt] = HexModelGrid.make_hex_points_horizontal_rect(3, 3, 1.0)
        >>> npt
        9
        >>> p[1,:]
        array([ 1.,  0.])
        >>> p[:3,0]
        array([ 0.,  1.,  2.])
        """

        dxv = dxh * numpy.sqrt(3.) / 2.
        half_dxh = dxh / 2.

        npts = num_rows * num_cols
        pts = numpy.zeros((npts, 2))
        xshift = 0.
        i = 0
        for r in range(num_rows):
            for c in range(num_cols):
                xshift = half_dxh * (r%2)
                pts[i,0] = c * dxh + xshift
                pts[i,1] = r * dxv
                i += 1
        
        return pts, npts


    @staticmethod
    def make_hex_points_vertical_hex(base_num_rows, num_cols, dxv):
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
                
        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> [p, npt] = HexModelGrid.make_hex_points_vertical_hex(2, 3, 1.0)
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
        
        
    @staticmethod
    def make_hex_points_vertical_rect(num_rows, num_cols, dxv):
        """
        Creates and returns a set of (x,y) points in a staggered grid in which the 
        points represent the centers of regular hexagonal cells, and the points
        could be connected to form equilateral triangles. The overall shape of the
        lattice is rectangular.
        
        Inputs: num_rows = number of columns in lattice
                num_cols = number of columns in lattice
                dxv = vertical and diagonal spacing between points
                
        Return: 2D numpy array containing point (x,y) coordinates, and total number
                of points.
                
        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> [p, npt] = HexModelGrid.make_hex_points_vertical_rect(3, 3, 1.0)
        >>> npt
        9
        >>> p[1,:]
        array([ 0.,  1.])
        >>> p[:3,1]
        array([ 0.,  1.,  2.])
        """

        dxh = dxv * numpy.sqrt(3.) / 2.
        half_dxv = dxv / 2.

        npts = num_rows * num_cols
        pts = numpy.zeros((npts, 2))
        yshift = 0.
        i = 0
        for c in range(num_cols):
            for r in range(num_rows):
                yshift = half_dxv * (c%2)
                pts[i,1] = r * dxv + yshift
                pts[i,0] = c * dxh
                i += 1
        
        return pts, npts
        
        
    @property
    def number_of_node_columns(self):
        """Number of node columns in a rectangular-shaped and/or
        vertically oriented hex grid.

        Returns the number of columns, including boundaries.
        
        Notes
        -----
        Will generate an error if called with a hex-shaped, horizontally
        aligned grid.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> grid = HexModelGrid(5, 5, shape='rect')
        >>> grid.number_of_node_columns
        5
        """
        return self._ncols


    @property
    def number_of_node_rows(self):
        """Number of node rows in a rectangular-shaped and/or
        horizontally oriented hex grid.

        Returns the number of rows, including boundaries.
        
        Notes
        -----
        Will generate an error if called with a hex-shaped, vertically
        aligned grid.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> grid = HexModelGrid(5, 5, shape='rect')
        >>> grid.number_of_node_rows
        5
        """
        return self._nrows


    def configure_hexplot(self, data, data_label=None, color_map=None):
        """
        Sets up necessary information for making plots of the hexagonal grid
        colored by a given data element.
        
        Parameters
        ----------
        data : str OR node array (1d numpy array with number_of_nodes entries)
            Data field to be colored
        data_label : str, optional
            Label for colorbar
        color_map : matplotlib colormap object, None
            Color map to apply (defaults to "jet")
        
        Returns
        -------
        (none)
        
        Notes
        -----
        Creates and stores a PatchCollection representing the hexagons. Also 
        stores a handle to the current plotting axis. Both of these are then
        used by hexplot().
        """        
        from numpy import array, sqrt, zeros
        import matplotlib
        from matplotlib.patches import Polygon
        from matplotlib.collections import PatchCollection
        #import matplotlib.pyplot as plt
        
        # color
        if color_map is None:
            color_map = matplotlib.cm.jet

        # geometry
        apothem = self._dx/2.0
        radius = 2.0*apothem / sqrt(3.0)  # distance from node to each hexagon cell vertex
        
        # offsets from node x,y position
        offsets = zeros((6,2))
        poly_verts = zeros((6,2))
        
        # Figure out whether the orientation is horizontal or vertical
        if self.node_y[0]==self.node_y[1]:   # horizontal
            offsets[:,0] = array([0., apothem, apothem, 0., -apothem, -apothem])
            offsets[:,1] = array([radius, radius/2.0, -radius/2.0, -radius, -radius/2.0, radius/2.0])
        else:   # vertical
            offsets[:,0] = array([radius/2.0, radius, radius/2.0, -radius/2.0, -radius, -radius/2.0])
            offsets[:,1] = array([apothem, 0., -apothem, -apothem, 0., apothem])
        
        patches = []
        for i in range(self.number_of_nodes):
            poly_verts[:,0] = self.node_x[i]+offsets[:,0]
            poly_verts[:,1] = self.node_y[i]+offsets[:,1]
            p = Polygon(poly_verts, True)
            patches.append(p)
        
        self._hexplot_pc = PatchCollection(patches, cmap=color_map, edgecolor='none', linewidth=0.0)
        
        self._hexplot_configured=True


    def hexplot(self, data, data_label=None, color_map=None):
        """
        Creates a plot of the grid and one node-data field, showing hexagonal
        cells colored by values in the field.
        
        Parameters
        ----------
        data : str OR node array (1d numpy array with number_of_nodes entries)
            Data field to be colored
        data_label : str, optional
            Label for colorbar
        
        Returns
        -------
        (none)
        """
        from numpy import array, amin, amax
        import matplotlib.pyplot as plt
        
        try:
            self._hexplot_configured is True
        except:
            self.configure_hexplot(data, data_label, color_map)

        # Handle *data*: if it's a numpy array, then we consider it the 
        # data to be plotted. If it's a string, we consider it the name of the 
        # node-field to plot, and we fetch it.
        if type(data) is str:
            data_label = data
            data = self.at_node[data]
            
        ax = plt.gca()
        self._hexplot_pc.set_array(array(data))
        ax.add_collection(self._hexplot_pc)
        plt.xlim([amin(self.node_x)-self._dx, amax(self.node_x)+self._dx])
        plt.ylim([amin(self.node_y)-self._dx, amax(self.node_y)+self._dx])
        #cb = plt.colorbar(self._hexplot_pc)
        #if data_label is not None:
        #    cb.set_label(data_label)
        
        #plt.show()


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
