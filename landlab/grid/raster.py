#! /usr/bin/env python

import numpy

from landlab.grid.base import ModelGrid
import landlab.utils.structured_grid as sgrid
from landlab.utils import count_repeated_values
from landlab.grid.base import (INTERIOR_NODE, FIXED_VALUE_BOUNDARY,
                               FIXED_GRADIENT_BOUNDARY, TRACKS_CELL_BOUNDARY,
                               INACTIVE_BOUNDARY, BAD_INDEX_VALUE,
                              )


def node_has_boundary_neighbor(mg, id):
    for neighbor in mg.get_neighbor_list(id):
        if mg.node_status[neighbor] != INTERIOR_NODE:
            return True
    for neighbor in mg.get_diagonal_list(id):
        if mg.node_status[neighbor] != INTERIOR_NODE:
            return True
    return False


has_boundary_neighbor = numpy.vectorize(node_has_boundary_neighbor,
                                        excluded=['mg'])


class RasterModelGrid(ModelGrid):
    """
    This inherited class implements a regular, raster 2D grid with uniform
    cell dimensions.
    
    Examples:
        
        >>> rmg = RasterModelGrid()
        >>> rmg.num_nodes
        0
        >>> rmg = RasterModelGrid(4, 5, 1.0) # rows, columns, spacing
        >>> rmg.num_nodes
        20
    """

    def __init__(self, num_rows=0, num_cols=0, dx=1.0):
        """
        Optionally takes numbers of rows and columns and cell size as
        inputs. If this are given, calls initialize() to set up the grid.
        
        """
        #print 'RasterModelGrid.init'
        
        # Set number of nodes, and initialize if caller has given dimensions
        self.ncells = num_rows * num_cols   #TBX
        self.num_nodes = num_rows * num_cols
        if self.num_nodes > 0:
            self.initialize( num_rows, num_cols, dx )

    def initialize( self, num_rows, num_cols, dx ):
        """
        Sets up a num_rows by num_cols grid with cell spacing dx and
        (by default) regular boundaries (that is, all perimeter cells are
        boundaries and all interior cells are active).

        To be consistent with unstructured grids, the raster grid is
        managed not as a 2D array but rather as a set of vectors that
        describe connectivity information between nodes, links, active links,
        cells, active cells, faces, patches, junctions, and corners.
        
        Examples and doctests:

            >>> import landlab as ll
            >>> rmg = RasterModelGrid()
            >>> numrows = 20          # number of rows in the grid
            >>> numcols = 30          # number of columns in the grid
            >>> dx = 10.0             # grid cell spacing
            >>> rmg.initialize(numrows, numcols, dx)
            >>> rmg.num_nodes,rmg.num_cells,rmg.num_links,rmg.num_active_links
            (600, 504, 1150, 1054)
            >>> rmg = RasterModelGrid(4, 5)
            >>> rmg.num_nodes,rmg.num_cells,rmg.num_links,rmg.num_active_links
            (20, 6, 31, 17)
            >>> rmg.node_status
            array([1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1], dtype=int8)
            >>> rmg.node_activecell[3] == ll.BAD_INDEX_VALUE
            True
            >>> rmg.node_activecell[8]
            2
            >>> rmg.node_numinlink
            array([0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 2, 2])
            >>> rmg.node_inlink_matrix
            array([[-1, 15, 16, 17, 18,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,
                    12, 13, 14],
                   [-1, -1, -1, -1, -1, -1, 19, 20, 21, 22, -1, 23, 24, 25, 26, -1, 27,
                    28, 29, 30]])
            >>> rmg.node_numoutlink
            array([2, 2, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0])
            >>> rmg.node_outlink_matrix[0]
            array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 27, 28,
                   29, 30, -1])
            >>> rmg.node_numactiveinlink
            array([0, 0, 0, 0, 0, 0, 2, 2, 2, 1, 0, 2, 2, 2, 1, 0, 1, 1, 1, 0])
            >>> rmg.node_active_inlink_matrix
            array([[-1, -1, -1, -1, -1, -1,  0,  1,  2, 12, -1,  3,  4,  5, 16, -1,  6,
                     7,  8, -1],
                   [-1, -1, -1, -1, -1, -1,  9, 10, 11, -1, -1, 13, 14, 15, -1, -1, -1,
                    -1, -1, -1]])
            >>> rmg.node_numactiveoutlink
            array([0, 1, 1, 1, 0, 1, 2, 2, 2, 0, 1, 2, 2, 2, 0, 0, 0, 0, 0, 0])
            >>> rmg.node_active_outlink_matrix
            array([[-1,  0,  1,  2, -1,  9,  3,  4,  5, -1, 13,  6,  7,  8, -1, -1, -1,
                    -1, -1, -1],
                   [-1, -1, -1, -1, -1, -1, 10, 11, 12, -1, -1, 14, 15, 16, -1, -1, -1,
                    -1, -1, -1]])
            >>> rmg.cell_node
            array([ 6,  7,  8, 11, 12, 13])
            >>> rmg.link_fromnode
            array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,  0,  1,
                    2,  3,  5,  6,  7,  8, 10, 11, 12, 13, 15, 16, 17, 18])
            >>> rmg.link_tonode
            array([ 5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,  1,  2,
                    3,  4,  6,  7,  8,  9, 11, 12, 13, 14, 16, 17, 18, 19])
            >>> rmg.link_face[20]
            10
            >>> rmg.active_links
            array([ 1,  2,  3,  6,  7,  8, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26])
        """
        
        if self.DEBUG_TRACK_METHODS:
            print 'RasterModelGrid.initialize('+str(num_rows)+', ' \
                   +str(num_cols)+', '+str(dx)+')'
        
        # Basic info about raster size and shape
        self.nrows = num_rows
        self.ncols = num_cols
        self.ncells = num_rows * num_cols # soon to be deprecated, June 2013
        self._dx = dx
        self.cellarea = dx*dx
        self.num_nodes = num_rows * num_cols
        self.num_cells = (num_rows-2) * (num_cols-2)
        self.num_active_cells = self.num_cells
        self.num_links = num_cols*(num_rows-1)+num_rows*(num_cols-1)
        self.num_active_links = self.num_links-(2*(num_cols-1)+2*(num_rows-1))
        
        # We need at least one row or column of boundary cells on each
        # side, so the grid has to be at least 3x3
        assert(numpy.min((num_rows, num_cols)) >= 3)

        # Record number of boundary and interior cells and the number
        # of interior faces. Ultimately, this info could be overridden
        # if using an irregular geometry of "interior" cells within the
        # rectangular domain. Note that we don't include any faces
        # between boundary cells.
        # NG Do we still have boundary cells?  I thought that cells were only
        # defined on interior nodes.
        self.n_boundary_cells = 2 * ( num_rows - 2 ) + 2 * ( num_cols - 2 ) + 4
        self.n_interior_cells = self.ncells - self.n_boundary_cells
        self.num_faces = ( num_rows - 1 ) * ( num_cols - 2 ) + \
                      ( num_rows - 2 ) * ( num_cols - 1 )
        self.nfaces = self.num_faces # TBX; for backward compatibility
        if self.DEBUG_VERBOSE:
            print self.num_faces
        
        # Assign and store node (x,y,z) coordinates.
        #
        # The relation between node (x,y) coordinates and position is
        # illustrated here for a five-column, four-row grid. The numbers show
        # node positions, and the - and | symbols show the links connecting
        # the nodes.
        #
        # 15------16------17------18------19
        #  |       |       |       |       |
        #  |       |       |       |       |
        #  |       |       |       |       |
        # 10------11------12------13------14
        #  |       |       |       |       |
        #  |       |       |       |       |   
        #  |       |       |       |       |
        #  5-------6-------7-------8-------9
        #  |       |       |       |       |
        #  |       |       |       |       |
        #  |       |       |       |       |
        #  0-------1-------2-------3-------4
        #
        (self._node_x, self._node_y) = sgrid.node_coords(
            (num_rows, num_cols), (self._dx, self._dx), (0., 0.))

        # Node boundary/active status:
        # Next, we set up an array of "node status" values, which indicate 
        # whether a given node is an active, non-boundary node, or some type of 
        # boundary. Here we default to having all perimeter nodes be active
        # fixed-value boundaries.
        self.node_status = sgrid.node_status(
            self.shape, boundary_status=FIXED_VALUE_BOUNDARY)
        
        # Cell lists:
        # For all cells, we create a list of the corresponding node ID for 
        # each cell.
        # We also have a list of the cell IDs of all active cells. By default,
        # all cells are active, so for example if there are six cells, the
        # self.active_cells list reads: 0, 1, 2, 3, 4, 5
        # 
        # Cells and faces in a five-column, four-row grid look like this
        # (where the numbers are cell IDs and lines show faces):
        #
        # |-------|-------|-------|
        # |       |       |       |
        # |   3   |   4   |   5   |
        # |       |       |       |
        # |-------|-------|-------|
        # |       |       |       |
        # |   0   |   1   |   2   |
        # |       |       |       |
        # |-------|-------|-------|
        #
        # While we're at it, we will also build the node_activecell list. This
        # list records, for each node, the ID of its associated active cell, 
        # or None if it has no associated active cell (i.e., it is a boundary)
        self.cell_node = sgrid.cell_node_index((num_rows, num_cols))
        self.node_activecell = sgrid.node_active_cell((num_rows, num_cols))
        self.active_cells = sgrid.active_cells((num_rows, num_cols))
        self.activecell_node = self.cell_node.copy()

        # Link lists:
        # For all links, we encode the "from" and "to" nodes, and the face
        # (if any) associated with the link. If the link does not intersect a
        # face, then face is assigned None.
        # For active links, we store the corresponding link ID.
        #
        # The numbering scheme for links in RasterModelGrid is illustrated with
        # the example of a five-column by four-row grid (each * is a node,
        # the lines show links, and the ^ and > symbols indicate the direction
        # of each link: up for vertical links, and right for horizontal ones):
        #
        #  *--27-->*--28-->*--29-->*--30-->*
        #  ^       ^       ^       ^       ^
        # 10      11      12      13      14
        #  |       |       |       |       |
        #  *--23-->*--24-->*--25-->*--26-->*
        #  ^       ^       ^       ^       ^
        #  5       6       7       8       9   
        #  |       |       |       |       |
        #  *--19-->*--20-->*--21-->*--22-->*
        #  ^       ^       ^       ^       ^
        #  0       1       2       3       4
        #  |       |       |       |       |
        #  *--15-->*--16-->*--17-->*--18-->*
        #
        #   create the fromnode and tonode lists
        (self.link_fromnode,
         self.link_tonode) = sgrid.node_link_index((num_rows, num_cols))

        #   set up in-link and out-link matrices and numbers
        self.setup_inlink_and_outlink_matrices()
        
        #   set up the list of active links
        self.reset_list_of_active_links()

        #   set up link faces
        #
        #   Here we assume that we've already created a list of active links
        # in which all 4 boundaries are "open", such that each boundary node
        # (except the 4 corners) is connected to an adjacent interior node. In
        # this case, there will be the same number of faces as active links,
        # and the numbering of faces will be the same as the corresponding
        # active links. We start off creating a list of all None values. Only
        # those links that cross a face will have this None value replaced with
        # a face ID.
        self.link_face = sgrid.link_faces((num_rows, num_cols),
                                          actives=self.active_links)

        # List of neighbors for each cell: we will start off with no
        # list. If a caller requests it via get_neighbor_list or
        # create_neighbor_list, we'll create it if necessary.
        self.neighbor_list_created = False

        # List of diagonal neighbors. As with the neighbor list, we'll only
        # create it if requested.
        self.diagonal_list_created = False

    def _setup_cell_areas_array(self):
        self.active_cell_areas = numpy.empty(self.num_active_cells)
        self.active_cell_areas.fill(self._dx ** 2)
        return self.active_cell_areas

    @property
    def shape(self):
        return (self.nrows, self.ncols)

    @property
    def dx(self):
        return self._dx

    @property
    def shape(self):
        return (self.nrows, self.ncols)

    def setup_inlink_and_outlink_matrices(self):
        """
        Creates data structures to record the numbers of inlinks and outlinks
        for each node. An inlink of a node is simply a link that has the node as
        its "to" node, and an outlink is a link that has the node as its "from".
        
        We store the inlinks in a 2-row by num_nodes-column matrix called
        node_inlink_matrix. It has two rows because we know that the nodes in
        our raster grid will never have more than two inlinks an two outlinks
        each (a given node could also have zero or one of either). The outlinks
        are stored in a similar matrix.
        
        We also keep track of the total number of inlinks and outlinks at each
        node in the num_inlinks and num_outlinks arrays.
        
        The inlink and outlink matrices are useful in numerical calculations.
        Each row of each matrix contains one inlink or outlink per node. So, if
        you have a corresponding "flux" matrix, you can map incoming or
        outgoing fluxes onto the appropriate nodes. More information on this is
        in the various calculate_flux_divergence... functions.
        
        What happens if a given node does not have two inlinks or outlinks? We
        simply put the default value -1 in this case. This allows us to use a 
        cute little trick when computing inflows and outflows. We make our 
        "flux" array one element longer than the number of links, with the last
        element containing the value 0. Thus, any time we add an influx from 
        link number -1, Python takes the value of the last element in the array,
        which is zero. By doing it this way, we maintain the efficiency that 
        comes with the use of numpy. Again, more info can be found in the 
        description of the flux divergence functions.
        
        Example:
            
            >>> rmg = RasterModelGrid(4, 5, 1.0)
        """

        (self.node_inlink_matrix,
         self.node_numinlink) = sgrid.setup_inlink_matrix(
             self.shape, tonodes=self.link_tonode)

        (self.node_outlink_matrix,
         self.node_numoutlink) = sgrid.setup_outlink_matrix(
             self.shape, fromnodes=self.link_fromnode)
        
    def setup_active_inlink_and_outlink_matrices(self):
        """
        Creates data structures to record the numbers of active inlinks and 
        active outlinks for each node. These data structures are equivalent to
        the "regular" inlink and outlink matrices, except that it uses the IDs
        of active links (only).
        """
        # Create active in-link and out-link matrices.
        self.node_active_inlink_matrix = - numpy.ones((2, self.num_nodes),
                                                       dtype=numpy.int)
        self.node_active_outlink_matrix = - numpy.ones((2, self.num_nodes),
                                                        dtype=numpy.int)
        # Set up the inlink arrays
        tonodes = self.activelink_tonode
        self.node_numactiveinlink = numpy.bincount(tonodes,
                                                   minlength=self.num_nodes)

        counts = count_repeated_values(self.activelink_tonode)
        for (count, (tonodes, active_link_ids)) in enumerate(counts):
            self.node_active_inlink_matrix[count][tonodes] = active_link_ids

        # Set up the outlink arrays
        fromnodes = self.activelink_fromnode
        self.node_numactiveoutlink = numpy.bincount(fromnodes,
                                                    minlength=self.num_nodes)
        counts = count_repeated_values(self.activelink_fromnode)
        for (count, (fromnodes, active_link_ids)) in enumerate(counts):
            self.node_active_outlink_matrix[count][fromnodes] = active_link_ids

    def cell_faces(self, cell_id):
        """
        Returns an array of the face IDs for the faces of a cell with ID,
        *cell_id*. The faces are listed clockwise, starting with the bottom
        face. *cell_id* can be either a scalar or an array. If an array,
        return the faces for each cell of the array.
        
        >>> mg = RasterModelGrid(4, 5)
        >>> mg.cell_faces(0)
        array([ 0,  9,  3, 10])

        >>> mg.cell_faces([0, 5])
        array([[ 0,  9,  3, 10],
               [ 5, 15,  8, 16]])
        """
        node_id = self.cell_node[cell_id]
        inlinks = self.node_inlink_matrix[:, node_id].T
        outlinks = self.node_outlink_matrix[:, node_id].T
        return numpy.concatenate(
            (self.link_face[inlinks], self.link_face[outlinks]), axis=1)

    def get_grid_xdimension(self):
        '''
        Returns the x dimension of the grid. Method added 5/1/13 by DEJH.
        '''
        return (self.ncols * self._dx)
    
    def get_grid_ydimension(self):
        '''
        Returns the y dimension of the grid. Method added 5/1/13 by DEJH.
        '''
        return (self.nrows * self._dx)
        
    def get_count_of_interior_nodes(self):
        """
        Returns the number of interior nodes on the grid.
        """
        return sgrid.interior_node_count(self.shape)
        
    def get_count_of_all_nodes(self):
        """
        Returns total number of nodes, including boundaries.
        """
        return(self.num_nodes)
    
    def get_count_of_cols(self):
        """
        Returns the number of columns, including boundaries.  
        NG, June 2013
        """
        return(self.ncols)
        
    def get_count_of_rows(self):
        """
        Returns the number of rows, including boundaries.  
        NG, June 2013
        """
        return(self.nrows)
    
    def get_grid_spacing(self):
        """
        Returns the spacing between grid nodes.
        DEJH July 2013
        """
        return(self._dx)

    def get_nodes_around_point(self, xcoord, ycoord):
        """
        This method takes an x,y coordinate within the grid, then returns the 
        IDs of the four nodes of the area (enclosure?) around that point as a 
        4 item list. Because the geometry of this grid is so simple, it works 
        purely by counting the number of squares left and below the point. 
        Method added 4/29/13 by DEJH.
        """
        ID = int(ycoord//self._dx * self.ncols + xcoord//self._dx)
        return numpy.array([ID, ID+self.ncols, ID+self.ncols+1, ID+1])
    
    def get_minimum_active_link_length(self):
        """
        Returns the horizontal length of the shortest active link in the grid.
        Overrides ModelGrid.get_minimum_active_link_length().
        """
        return self._dx

    def calculate_max_gradient_across_node(self, u, cell_id):
        '''
            This method calculates the gradients in u across all 4 faces of the 
            cell with ID cell_id, and across the four diagonals. It then returns 
            the steepest (most negative) of these values, followed by its dip 
            direction (e.g.: 0.12, 225). i.e., this is a D8 algorithm. Slopes 
            downward from the cell are reported as positive.
            
            This code is actually calculating slopes, not gradients.  
            The max gradient is the most negative, but the max slope is the most
            positive.  So, this was updated to return the max value, not the 
            min.
            
        GT: Might be possible to speed this up using inlink_matrix and 
        outlink_matrix.
        '''
        #We have poor functionality if these are edge cells! Needs an exception
        neighbor_cells = self.get_neighbor_list(cell_id)
        neighbor_cells.sort()
        #print 'Node is internal: ', self.is_interior(cell_id)
        #print 'Neighbor cells: ', neighbor_cells
        diagonal_cells = []
        if neighbor_cells[0]!=-1:
            diagonal_cells.extend([neighbor_cells[0]-1, neighbor_cells[0]+1])
        if neighbor_cells[3]!=-1:
            diagonal_cells.extend([neighbor_cells[3]-1, neighbor_cells[3]+1])
        slopes = []
        diagonal_dx = numpy.sqrt(2.)
        for a in neighbor_cells:
            #ng I think this is actually slope as defined by a geomorphologist,
            #that is -dz/dx and not the gradient (dz/dx)
            single_slope = (u[cell_id] - u[a])/self._dx
            #print 'cell id: ', cell_id
            #print 'neighbor id: ', a
            #print 'cell, neighbor are internal: ', self.is_interior(cell_id), self.is_interior(a)
            #print 'cell elev: ', u[cell_id]
            #print 'neighbor elev: ', u[a]
            #print single_slope
            if not numpy.isnan(single_slope): #This should no longer be necessary, but retained in case
                slopes.append(single_slope)
            else:
                print 'NaNs present in the grid!'
        for a in diagonal_cells:
            single_slope = (u[cell_id] - u[a])/diagonal_dx
            #print single_slope
            if not numpy.isnan(single_slope):
                slopes.append(single_slope)
            else:
                print 'NaNs present in the grid!'
        #print 'Slopes list: ', slopes
        #ng thinks that the maximum slope should be found here, not the 
        #minimum slope, old code commented out.  New code below it.
        #if slopes:
        #    min_slope, index_min = min((min_slope, index_min) for (index_min, min_slope) in enumerate(slopes))
        #else:
        #    print u
        #    print 'Returning NaN angle and direction...'
        #    min_slope = numpy.nan
        #    index_min = 8
        if slopes:
            max_slope, index_max = max((max_slope, index_max) for (index_max, max_slope) in enumerate(slopes))
        else:
            print u
            print 'Returning NaN angle and direction...'
            max_slope = numpy.nan
            index_max = 8
            
        angles = [180., 270., 90., 0., 225., 135., 315., 45., numpy.nan] #This is inefficient
        
        #ng commented out old code
        #return min_slope, angles[index_min]
        return max_slope, angles[index_max]
        
    def find_node_in_direction_of_max_slope(self, u, node_id):
        '''
            This method calculates the slopes (-dz/dx) in u across all 4 faces of 
            the cell with ID node_id, and across the four diagonals. 
            It then returns the node ID in the direction of the steepest 
            (most positive) of these values,  i.e., this is a 
            D8 algorithm. Slopes downward from the cell are reported as positive.
            Based on code from DH, modified by NG, 6/2013
            
            This doesn't deal with the fixed gradient boundary condition.  
            NG is still confused about that one.
            
            NMG Update.  This is super clumsy. 
            
            DEJH update: Gets confused for the lowest node if w/i grid
            (i.e., closed)- will return a higher neighbour, when it should
            return itself. ->  Now returns itself, but a flagging mechanism
            may be preferable. (i.e., unthinking application will result in a
            doubling of discharge at that node.
        '''
        #We have poor functionality if these are closed boundary nodes! 
        neighbor_nodes = self.get_neighbor_list(node_id)
        neighbor_nodes.sort()
        #print 'Node is internal: ', self.is_interior(cell_id)
        #print 'Neighbor cells: ', neighbor_cells
        diagonal_nodes = []
        #NG also think that this won't happen if you are always sending this 
        #function an id of an interior node.  But maybe there is a case where 
        #this would happen?
        if neighbor_nodes[0]!=-1:
            diagonal_nodes.extend([neighbor_nodes[0]-1, neighbor_nodes[0]+1])
        #ng, if neighbor_nodes is sorted, how could [3] be -1?
        #try commenting out.
        #if neighbor_cells[3]!=-1:
        diagonal_nodes.extend([neighbor_nodes[3]-1, neighbor_nodes[3]+1])
        slopes = []
        diagonal_dx = numpy.sqrt(2.)
        for a in neighbor_nodes:
            if self.node_status[a] != INACTIVE_BOUNDARY:
                single_slope = (u[node_id] - u[a])/self._dx
            else:
                single_slope = -9999
            #print 'cell id: ', cell_id
            #print 'neighbor id: ', a
            #print 'status: ', self.node_status[a]
            #print 'cell, neighbor are internal: ', self.is_interior(cell_id), self.is_interior(a)
            #print 'cell elev: ', u[cell_id]
            #print 'neighbor elev: ', u[a]
            #print single_slope
            if not numpy.isnan(single_slope): #This should no longer be necessary, but retained in case
                slopes.append(single_slope)
            else:
                print 'NaNs present in the grid!'
        for a in diagonal_nodes:
            if self.node_status[a] != INACTIVE_BOUNDARY:
                single_slope = (u[node_id] - u[a])/diagonal_dx
            else:
                single_slope = -9999
            #print single_slope
            if not numpy.isnan(single_slope):
                slopes.append(single_slope)
            else:
                print 'NaNs present in the grid!'
        #print 'Slopes list: ', slopes
        if slopes:
            max_slope, index_max = max((max_slope, index_max) for (index_max, max_slope) in enumerate(slopes))
        else:
            print u
            print 'Returning NaN angle and direction...'
            max_slope = numpy.nan
            index_max = 8
        
        all_neighbor_nodes=numpy.concatenate((neighbor_nodes,diagonal_nodes))
        #print 'all_neighbor_cells ', all_neighbor_cells
        
        #Final check to  allow correct handling of internally draining nodes; DEJH Aug 2013.
        #This remains extremely ad-hoc. An internal node points to itself, but this should never
        #be used to actually route flow. In flow_accumulation, there is an explicit check that flow
        #is not routed to yourself.
        steepest_node = all_neighbor_nodes[index_max]
        #...now if a node is the lowest thing, this method returns that node, not a neighbor:
        if u[steepest_node] > u[node_id]:
            steepest_node=node_id
            #When accumulating flow, it is ESSENTIAL to check if a node routes to itself!!
            #An alternative, creating a list of the internally draining cells
            #try:
            #    self.internal_draining_nodes.append(node_id)
            #except:
            #    self.internal_draining_nodes = [node_id]
            #self.node_status[node_id] = INTERNAL_DRAINING_NODE #==5? ...this doesn't exist yet.
        
        return steepest_node
        
    def set_inactive_boundaries(self, bottom_is_inactive, right_is_inactive, 
                                top_is_inactive, left_is_inactive):
        """
        Handles boundary conditions by setting each of the four sides of the 
        rectangular grid to either 'inactive' or 'active (fixed value)' status.
        Arguments are booleans indicating whether the bottom, right, top, and
        left are inactive (True) or not (False).
        
        For an inactive boundary:
            - the nodes are flagged INACTIVE_BOUNDARY
            - the links between them and the adjacent interior nodes are
              inactive (so they appear on link-based lists, but not
              active_link-based lists)
              
        This means that if you call the calculate_gradients_at_active_links
        method, the inactive boundaries will be ignored: there can be no
        gradients or fluxes calculated, because the links that connect to that
        edge of the grid are not included in the calculation. So, setting a
        grid edge to INACTIVE_BOUNDARY is a convenient way to impose a no-flux
        boundary condition. Note, however, that this applies to the grid as a
        whole, rather than a particular variable that you might use in your
        application. In other words, if you want a no-flux boundary in one
        variable but a different boundary condition for another, then use 
        another method.
        
        The following example sets the top and left boundaries as inactive in a
        four-row by five-column grid that initially has all boundaries active
        and all boundary nodes coded as FIXED_VALUE_BOUNDARY (=1):
        
        >>> rmg = RasterModelGrid(4, 5, 1.0) # rows, columns, spacing
        >>> rmg.num_active_links
        17
        >>> rmg.node_status
        array([1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1], dtype=int8)
        >>> rmg.set_inactive_boundaries(False, False, True, True)
        >>> rmg.num_active_links
        12
        >>> rmg.node_status
        array([1, 1, 1, 1, 1, 4, 0, 0, 0, 1, 4, 0, 0, 0, 1, 4, 4, 4, 4, 4], dtype=int8)
        
        Note that the four corners are treated as follows:
            bottom left = BOTTOM
            bottom right = RIGHT
            top right = TOP
            top left = LEFT
        """
        if self.DEBUG_TRACK_METHODS:
            print 'ModelGrid.set_inactive_boundaries'
            
        bottom_edge = range(0,self.ncols-1)
        right_edge = range(self.ncols-1,self.num_nodes-1,self.ncols)
        top_edge = range((self.nrows-1)*self.ncols+1,self.num_nodes)
        left_edge = range(self.ncols,self.num_nodes,self.ncols)
            
        if bottom_is_inactive:
            self.node_status[bottom_edge] = INACTIVE_BOUNDARY
        else:
            self.node_status[bottom_edge] = FIXED_VALUE_BOUNDARY

        if right_is_inactive:
            self.node_status[right_edge] = INACTIVE_BOUNDARY
        else:
            self.node_status[right_edge] = FIXED_VALUE_BOUNDARY
            
        if top_is_inactive:
            self.node_status[top_edge] = INACTIVE_BOUNDARY
        else:
            self.node_status[top_edge] = FIXED_VALUE_BOUNDARY

        if left_is_inactive:
            self.node_status[left_edge] = INACTIVE_BOUNDARY
        else:
            self.node_status[left_edge] = FIXED_VALUE_BOUNDARY
        
        self.reset_list_of_active_links()

                
    def set_noflux_boundaries( self, bottom, right, top, left,
                               bc = None ):
        """
        Assigns "no flux" status to one or more sides of the rectangular
        domain, for BoundaryCondition "bc", which default's to ModelGrid's
        self.default_bc (i.e., the default BoundaryCondition used when
        the user doesn't specify another one).

        Boundary cells are either "fixed value" (Dirichlet), which is
        the default, "fixed gradient" (Neumann with a zero derivative), or
        "tracks cell" (they track a cell in the interior on the opposite
        side of the grid, e.g., for periodic). Here we implement no flux
        by mirroring adjacent cells in the interior.

        Boundary status is recorded in the bc's TRACKS_CELL vector, 
        defined in bc's initialize() and N_BOUNDARY_CELLS long. For
        no-flux cells, this vector contains the CID of the neighboring
        cell whose value it will mirror in order to maintain a zero
        gradient. For periodic boundary cells, this vector contains the
        CID of the cell on the opposite side whose value it will track.
        Fixed value cells are indicated by a -1 for TRACKS_CELL.

        For no-flux boundaries, the corner cells (which don't really
        matter much anyway), the neighbor is arbitrarily chosen as either 
        the righthand or lefthand cell.
        
        GT: I believe this is now obsolete, because we can do no-flux simply by 
        setting closed boundaries (Aug 2013).
        """
        
        if bc==None:
            bc = self.default_bc
        
        # For no-flux boundaries, we need to know which interior
        # cells to mirror.
        #self.boundary_nbrs = zeros( self.n_boundary_cells, dtype=numpy.int )
        lower_left = 0
        lower_right = self.ncols-1
        upper_right = self.ncols+self.nrows-2
        upper_left = 2*self.ncols+self.nrows-3

        if bottom:
            for id in xrange( 1, self.ncols-1 ):   # Bottom
                bc.boundary_code[id] = bc.TRACKS_CELL_BOUNDARY
                bc.tracks_cell[id] = id+self.ncols
            bc.boundary_code[lower_left] = bc.TRACKS_CELL_BOUNDARY
            bc.tracks_cell[lower_left] = 1
            bc.boundary_code[lower_right] = bc.TRACKS_CELL_BOUNDARY
            bc.tracks_cell[lower_right] = lower_right-1
        if right:
            nbr = 2*self.ncols-2
            for id in xrange( lower_right+1, upper_right ):   # Right
                bc.boundary_code[id] = bc.TRACKS_CELL_BOUNDARY
                bc.tracks_cell[id] = nbr
                nbr += self.ncols
        if top:
            nbr = self.ncells - (self.ncols+2)
            for id in xrange( upper_right+1, upper_left ):   # Top
                bc.boundary_code[id] = bc.TRACKS_CELL_BOUNDARY
                bc.tracks_cell[id] = nbr
                nbr = nbr - 1
            bc.boundary_code[upper_right] = bc.TRACKS_CELL_BOUNDARY
            bc.tracks_cell[upper_right] = self.ncells-2
            bc.boundary_code[upper_left] = bc.TRACKS_CELL_BOUNDARY
            bc.tracks_cell[upper_left] = self.ncells+1-self.ncols
        if left:
            nbr = (self.nrows-2)*self.ncols + 1
            for id in xrange( upper_left+1, self.n_boundary_cells ):   # Left
                bc.boundary_code[id] = bc.TRACKS_CELL_BOUNDARY
                bc.tracks_cell[id] = nbr
                nbr = nbr - self.ncols
        
        if self.DEBUG_VERBOSE:
            print 'tracks_cell:',bc.tracks_cell
    
    def calculate_gradients_at_active_links(self, s, gradient=None):
        """
        Calculates the gradient in quantity s at each active link in the grid.
        This is nearly identical to the method of the same name in ModelGrid,
        except that it uses self._dx for link length to improve efficiency.
        
        Example:
        
            >>> rmg = RasterModelGrid(4, 5, 1.0)
            >>> u = [0., 1., 2., 3., 0.,
            ...     1., 2., 3., 2., 3.,
            ...     0., 1., 2., 1., 2.,
            ...     0., 0., 2., 2., 0.]
            >>> u = numpy.array(u)
            >>> u
            array([ 0.,  1.,  2.,  3.,  0.,  1.,  2.,  3.,  2.,  3.,  0.,  1.,  2.,
                    1.,  2.,  0.,  0.,  2.,  2.,  0.])
            >>> grad = rmg.calculate_gradients_at_active_links(u)
            >>> grad
            array([ 1.,  1., -1., -1., -1., -1., -1.,  0.,  1.,  1.,  1., -1.,  1.,
                    1.,  1., -1.,  1.])
            
        For greater speed, sending a pre-created numpy array as an argument
        avoids having to create a new one with each call:
            
            >>> grad = numpy.zeros(rmg.num_active_links)
            >>> u = u*10
            >>> grad = rmg.calculate_gradients_at_active_links(u, grad)
            >>> grad
            array([ 10.,  10., -10., -10., -10., -10., -10.,   0.,  10.,  10.,  10.,
                   -10.,  10.,  10.,  10., -10.,  10.])
        """
        if self.DEBUG_TRACK_METHODS:
            print 'RasterModelGrid.calculate_gradients_at_active_links'

        if gradient==None:
            gradient = numpy.zeros(self.num_active_links)
            
        assert (len(gradient)==self.num_active_links), \
                "len(gradient)!=num_active_links"
     
        gradient = (s[self.activelink_tonode]-s[self.activelink_fromnode])/self._dx
        
        return gradient
        
    def calculate_flux_divergence_at_active_cells(self, active_link_flux, 
                                                  net_unit_flux=False):
        """
        Given an array of fluxes along links, computes the net total flux
        within each cell, divides by cell area, and stores the result in
        net_unit_flux. Overrides method of the same name in ModelGrid (nearly
        identical, but uses scalars dx and cellarea instead of variable
        link length and active cell area, respectively).
        
        The function works by calling calculate_flux_divergence_at_nodes, then
        slicing out only the values at active cells. Therefore, it is slower
        than calculate_flux_divergence_at_nodes, even though it returns a
        shorter list of numbers.
        
        The input active_link_flux should be flux of
        something (e.g., mass, momentum, energy) per unit face width, positive
        if flowing in the same direction as its link, and negative otherwise.
        There should be one value per active link. Returns an array of net
        total flux per unit area, one value per active cell (creates this
        array if it is not given as an argument).
          By convention, divergence is positive for net outflow, and negative 
        for net outflow. That's why we *add* outgoing flux and *subtract* 
        incoming flux. This makes net_unit_flux have the same sign and 
        dimensions as a typical divergence term in a conservation equation.

        In general, for a polygonal cell with $N$ sides of lengths
        Li and with surface area A, the net influx divided by cell
        area would be:
            .. math::
                {Q_{net} \over A} = {1 \over A} \sum{q_i L_i}

        For a square cell, which is what we have in RasterModelGrid,
        the sum is over 4 sides of length dx, and
        :math:`A = dx^2`, so:
            .. math::
                {Q_{net} \over A} = {1 \over dx} \sum{q_i}

        .. note::
            The net flux is defined as positive outward, negative
            inward. In a diffusion problem, for example, one would use:
                .. math::
                    {du \over dt} = \\text{source} - \\text{fd}
            where fd is "flux divergence".
            
        Example:
            
            >>> rmg = RasterModelGrid(4, 5, 1.0)
            >>> u = [0., 1., 2., 3., 0.,
            ...      1., 2., 3., 2., 3.,
            ...      0., 1., 2., 1., 2.,
            ...      0., 0., 2., 2., 0.]
            >>> u = numpy.array(u)
            >>> grad = rmg.calculate_gradients_at_active_links(u)
            >>> grad
            array([ 1.,  1., -1., -1., -1., -1., -1.,  0.,  1.,  1.,  1., -1.,  1.,
                    1.,  1., -1.,  1.])
            >>> flux = -grad    # downhill flux proportional to gradient
            >>> divflux = rmg.calculate_flux_divergence_at_active_cells(flux)
            >>> divflux
            array([ 2.,  4., -2.,  0.,  1., -4.])
            
        If calculate_gradients_at_active_links is called inside a loop, you can
        improve speed slightly by creating an array outside the loop. For 
        example, do this once, before the loop:
            
            >>> divflux = rmg.create_active_cell_dvector() # outside loop
            
        Then do this inside the loop:
            
            >>> divflux = rmg.calculate_flux_divergence_at_active_cells(flux, divflux)
            
        In this case, the function will not have to create the divflux array.
            
        """
        
        if self.DEBUG_TRACK_METHODS:
            print 'RasterModelGrid.calculate_flux_divergence_at_active_cells'
            
        assert (len(active_link_flux)==self.num_active_links), \
               "incorrect length of active_link_flux array"
            
        # If needed, create net_unit_flux array
        if net_unit_flux is False:
            net_unit_flux = numpy.zeros(self.num_active_cells)
        else:
            net_unit_flux[:] = 0.
            
        assert (len(net_unit_flux))==self.num_active_cells
        
        node_net_unit_flux = self.calculate_flux_divergence_at_nodes(active_link_flux)
                
        net_unit_flux = node_net_unit_flux[self.activecell_node]
        
        return net_unit_flux

    def calculate_max_gradients_at_nodes(self, link_gradients, max_gradient=False):
        """
        Created DEJH Sept 2013. Needs proper documentation. Based on approach of calc_flux_divergence..., below.
        NOT YET FUNCTIONAL.
        """
        
        if self.DEBUG_TRACK_METHODS:
            print 'RasterModelGrid.calculate_max_gradients_at_nodes'
            
        assert (len(link_gradients)==self.num_links), \
               "incorrect length of active_link_flux array"

        # If needed, create max_gradient array
        if max_gradient is False:
            max_gradient = numpy.zeros(self.num_nodes)
        else:
            max_gradient[:] = 0.

        assert(len(max_gradient) == self.num_nodes)

        gradients = numpy.zeros(len(link_gradients)+1)
        gradients[:-1] = link_gradients
        
        #Extract the values for each node_inlink/outlink_matrix, remembering to make the inlinks negative.
        #Search for the maxima.
        #Introduce a way of working on the diagonals!!!
        
        gradients[self.node_active_outlink_matrix[0][:]]


    def calculate_flux_divergence_at_nodes(self, active_link_flux, 
                                           net_unit_flux=False):
        """
        Same as calculate_flux_divergence_at_active_cells, but works with and
        returns a list of net unit fluxes that corresponds to all nodes, rather
        than just active cells.
        
        Note that we DO compute net unit fluxes at boundary nodes (even though
        these don't have active cells associated with them, and often don't have 
        cells of any kind, because they are on the perimeter). It's up to the 
        user to decide what to do with these boundary values.
        
        Example:
            
            >>> rmg = RasterModelGrid(4, 5, 1.0)
            >>> u = [0., 1., 2., 3., 0.,
            ...      1., 2., 3., 2., 3.,
            ...      0., 1., 2., 1., 2.,
            ...      0., 0., 2., 2., 0.]
            >>> u = numpy.array(u)
            >>> grad = rmg.calculate_gradients_at_active_links(u)
            >>> grad
            array([ 1.,  1., -1., -1., -1., -1., -1.,  0.,  1.,  1.,  1., -1.,  1.,
                    1.,  1., -1.,  1.])
            >>> flux = -grad    # downhill flux proportional to gradient
            >>> df = rmg.calculate_flux_divergence_at_nodes(flux)
            >>> df
            array([ 0., -1., -1.,  1.,  0., -1.,  2.,  4., -2.,  1., -1.,  0.,  1.,
                   -4.,  1.,  0., -1.,  0.,  1.,  0.])
            
        If calculate_gradients_at_nodes is called inside a loop, you can
        improve speed by creating an array outside the loop. For example, do
        this once, before the loop:
            
            >>> df = rmg.create_node_dvector() # outside loop
            >>> rmg.num_nodes
            20
            
        Then do this inside the loop:
            
            >>> df = rmg.calculate_flux_divergence_at_nodes(flux, df)
            
        In this case, the function will not have to create the df array.
        """
        
        if self.DEBUG_TRACK_METHODS:
            print 'RasterModelGrid.calculate_flux_divergence_at_nodes'
            
        assert (len(active_link_flux)==self.num_active_links), \
               "incorrect length of active_link_flux array"
            
        # If needed, create net_unit_flux array
        if net_unit_flux is False:
            net_unit_flux = numpy.zeros(self.num_nodes)
        else:
            net_unit_flux[:] = 0.
            
        assert(len(net_unit_flux) == self.num_nodes)
        
        flux = numpy.zeros(len(active_link_flux)+1)
        flux[:len(active_link_flux)] = active_link_flux * self._dx
        
        net_unit_flux = ((flux[self.node_active_outlink_matrix[0][:]] + \
                          flux[self.node_active_outlink_matrix[1][:]]) - \
                         (flux[self.node_active_inlink_matrix[0][:]] + \
                          flux[self.node_active_inlink_matrix[1][:]])) / self.cellarea

        return net_unit_flux
        
    def calculate_flux_divergence( self, q, id ):
        """
        ..todo: UPDATE THIS TO USE NEW DATA STRUCTURES!
        
        This is like calculate_flux_divergences (plural!), but only does
        it for cell "id".
        """
    
        if self.DEBUG_TRACK_METHODS:
            print 'RasterModelGrid.calculate_flux_divergence here with cell',id
            print 'q:',q[self.faces[id,0:4]]
        fd = ( -( q[self.faces[id,2]]   # left face (positive=in)
            + q[self.faces[id,3]] )           # bottom face (positive=in)
            + q[self.faces[id,0]]             # right face (positive=out)
            + q[self.faces[id,1]]             # top face (positive=out)
              ) / self._dx
        return fd
        
    def update_noflux_boundaries( self, u, bc = None ):
        """
        Sets the value of u at all noflux boundary cells equal to the
        value of their interior neighbors, as recorded in the
        "boundary_nbrs" array.

        .. deprecated:: 0.1
            Use `update_boundaries` instead
        """
    
        if bc==None:
            bc = self.default_bc

        inds = (bc.boundary_code[id] == bc.TRACKS_CELL_BOUNDARY)
        u[self.boundary_cells[inds]] = u[bc.tracks_cell[inds]]

        return u

    def update_boundaries( self, u, bc = None ):
        """
        Updates FIXED_GRADIENT and TRACKS_CELL boundaries in 
        BoundaryCondition "bc" (by default ModelGrid's default_bc), so
        that TRACKS_CELL values in u are assigned the value at the 
        corresponding cell-to-track, and FIXED_GRADIENT cells get a value
        equal to the cell-to-track's value plus gradient times distance.

        .. note:: Does NOT change fixed-value boundary cells.

        .. todo::
            use the cell-local distance rather than dx, for use in the base
            class!
            probably now obsolete (GT Aug 2013)
            
        NG Wondering if this should be boundary nodes, not cells.    
        """
    
        if bc == None:
            bc = self.default_bc

        inds = (bc.boundary_code == bc.TRACKS_CELL_BOUNDARY)
        u[self.boundary_cells[inds]] = u[bc.tracks_cell[inds]]

        inds = (bc.boundary_code == bc.FIXED_GRADIENT_BOUNDARY)
        u[self.boundary_cells[inds]] = (u[bc.tracks_cell[id]] +
                                        bc.gradient[id] * self._dx)

        return u

    def node_vector_to_raster(self, u, flip_vertically=False):
        """
        Converts node vector *u* to a 2D array and returns it, so that it
        can be plotted, output, etc.
        
        If the *flip_vertically* keyword is True, this function returns an 
        array that has the rows in reverse order. This is useful for use in
        plot commands (such as the image display functions) that puts the 
        first row at the top of the image. In the landlab coordinate system,
        the first row is thought to be at the bottom. Thus, a flipped matrix
        will plot in the landlab style with the first row at the bottom.

        The returned array is a view of *u*, not a copy.
        
        Example:
            
        >>> rmg = RasterModelGrid(4, 5, 1.0)
        >>> u = rmg.create_node_dvector()
        >>> u = u + range(0, len(u))
        >>> u
        array([  0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.,
                11.,  12.,  13.,  14.,  15.,  16.,  17.,  18.,  19.])
        >>> ur = rmg.node_vector_to_raster(u)
        >>> ur
        array([[  0.,   1.,   2.,   3.,   4.],
               [  5.,   6.,   7.,   8.,   9.],
               [ 10.,  11.,  12.,  13.,  14.],
               [ 15.,  16.,  17.,  18.,  19.]])
        >>> ur = rmg.node_vector_to_raster(u, flip_vertically=True)        
        >>> ur
        array([[ 15.,  16.,  17.,  18.,  19.],
               [ 10.,  11.,  12.,  13.,  14.],
               [  5.,   6.,   7.,   8.,   9.],
               [  0.,   1.,   2.,   3.,   4.]])
        """
        return sgrid.reshape_array(self.shape, u,
                                   flip_vertically=flip_vertically)

    def cell_vector_to_raster(self, u, flip_vertically=False):
        """
        Converts cell (i.e., interior node) vector u to a 2D array and returns it, 
        so that it can be plotted, output, etc.
        
        If the optional argument flip_vertically=True, the function returns an 
        array that has the rows in reverse order, for use in plot commands (such
        as the image display functions) that put the (0,0) axis at the top left 
        instead of the bottom left.
        
        Example:
            
            >>> rmg = RasterModelGrid(4, 5, 1.0)
            >>> u = rmg.create_cell_dvector()
            >>> u = u + range(0, len(u))
            >>> u
            array([ 0.,  1.,  2.,  3.,  4.,  5.])
            >>> ur = rmg.cell_vector_to_raster(u)
            >>> ur
            array([[ 0.,  1.,  2.],
                   [ 3.,  4.,  5.]])
            >>> ur = rmg.cell_vector_to_raster(u, flip_vertically=True)        
            >>> ur
            array([[ 3.,  4.,  5.],
                   [ 0.,  1.,  2.]])
        """
        return sgrid.reshape_array((self.shape[0] - 2, self.shape[1] - 2),
                                   u, flip_vertically=flip_vertically)

    def get_neighbor_list(self, *args):
        """get_neighbor_list([ids])

        Return lists of neighbor nodes for nodes with given *ids*. If *ids*
        is not given, return the neighbors for all of the nodes in the grid.
        For each node, the list gives neighbor ids as [right, top, left,
        bottom]. Set all neighbors of boundary nodes to -1.
        
        >>> mg = RasterModelGrid(4, 5)
        >>> mg.get_neighbor_list([-1, 6])
        array([[-1, -1, -1, -1],
               [ 7, 11,  5,  1]])
        >>> mg.get_neighbor_list(7)
        array([ 8, 12,  6,  2])

        ..todo: could use inlink_matrix, outlink_matrix
        """
        if self.neighbor_list_created == False:
            self.create_neighbor_list()

        if len(args) == 0:
            return self.neighbor_nodes
        elif len(args) == 1:
            return self.neighbor_nodes[args[0], :]
        else:
            raise ValueError('only zero or one arguments accepted')

    def create_neighbor_list( self ):
        """
        Creates a list of IDs of neighbor nodes for each node, as a
        2D array. Only interior nodes are assigned neighbors; boundary
        nodes get -1 for each neighbor. The order of the neighbors is [right,
        top, left, bottom].

        .. note:: This is equivalent to the neighbors of all cells,
            and setting the neighbors of boundary-node cells to -1. In such a
            case, each node has one cell and each node-cell pair have the
            same ID. However, this is the old-style grid structure as
            boundary nodes no longer have associated cells.

        .. todo: could use inlink_matrix, outlink_matrix

        .. todo: Change to use BAD_INDEX_VALUE instead of -1.
        """
        # DH created this.  NG only changed labels.
        assert(self.neighbor_list_created == False)

        self.neighbor_list_created = True 
        self.neighbor_nodes = sgrid.neighbor_node_array(
            self.shape, out_of_bounds=-1, boundary_node_mask=-1)
                
    def has_boundary_neighbor(self, ids):
        """
        Checks to see if one of the eight neighbor nodes of node(s) with
        *id* has a boundary node.  Returns True if a node has a boundary node,
        False if all neighbors are interior.

        >>> mg = RasterModelGrid(5, 5)
        >>> mg.has_boundary_neighbor(6)
        True
        >>> mg.has_boundary_neighbor(12)
        False
        >>> mg.has_boundary_neighbor([12, -1])
        array([False,  True], dtype=bool)

        >>> mg.has_boundary_neighbor(25)
        Traceback (most recent call last):
            ...
        IndexError: index 25 is out of bounds for axis 0 with size 25
        """
        ans = has_boundary_neighbor(self, ids)
        if ans.ndim == 0:
            return bool(ans)
        else:
            return ans

    def get_diagonal_list(self, *args):
        """get_diagonal_list([ids])

        Return lists of diagonals nodes for nodes with given *ids*. If *ids*
        is not given, return the diagonals for all of the nodes in the grid.
        For each node, the list gives diagonal ids as [topright, topleft,
        bottomleft, bottomright]. Set all diagonals for boundary nodes to -1.
        
        >>> mg = RasterModelGrid(4, 5)
        >>> mg.get_diagonal_list([-1, 6])
        array([[-1, -1, -1, -1],
               [12, 10,  0,  2]])
        >>> mg.get_diagonal_list(7)
        array([13, 11,  1,  3])

        ..todo: could use inlink_matrix, outlink_matrix
        """
        #Added DEJH 051513
    
        if self.diagonal_list_created==False:
            self.create_diagonal_list()
        
        if len(args) == 0:
            return self.diagonal_cells
        elif len(args) == 1:
            return self.diagonal_cells[args[0], :]
        else:
            raise ValueError('only zero or one arguments accepted')

    def create_diagonal_list(self):
        """
        Creates a list of IDs of the diagonal nodes to each node, as a 2D
        array.  Only interior nodes are assigned diagonal neighbors; boundary
        nodes get -1 for each neighbor. The order of the diagonal nodes is
        [topright, topleft, bottomleft, bottomright].
        
        .. note:: This is equivalent to the diagonals of all cells,
            and setting the neighbors of boundary-node cells to -1. In such a
            case, each node has one cell and each node-cell pair have the
            same ID. However, this is the old-style grid structure as
            boundary nodes no longer have associated cells.

        .. todo: Change to use BAD_INDEX_VALUE instead of -1.
        """
        assert(self.diagonal_list_created == False)

        self.diagonal_list_created = True
        self.diagonal_cells = sgrid.diagonal_node_array(
            self.shape, out_of_bounds=-1, boundary_node_mask=-1)

    def is_interior( self, id ):
        """
        Returns True if the node is an interior node, False otherwise. 
        Interior status is indicated by a value of 0 in node_status.
        """
        # NG changed this.
        #ng thinks there may be a problem here.
        #return self.boundary_ids[id] < 0
        #return self.node_status[id] < 1
        return self.node_status[id] == INTERIOR_NODE
        
    def get_boundary_code( self, id ):
        """
        Returns the boundary status of a node.
        """
        # ng june 2013
        return self.node_status[id] 
        
    def get_face_connecting_cell_pair(self, cell_a, cell_b):
        """
        Returns an array of face indices that *cell_a* and *cell_b* share.
        If the cells do not share any faces, returns an empty array.
        """
        cell_faces = self.cell_faces([cell_a, cell_b])
        return numpy.intersect1d(cell_faces[0], cell_faces[1],
                                 assume_unique=True)

    def top_edge_node_ids(self):
        """
        Returns a 1D numpy integer array containing the node ID numbers of the 
        nodes along the top (y=ymax) grid edge.
        
        Example:
            
            >>> rmg = RasterModelGrid(4, 5, 1.0)
            >>> rmg.top_edge_node_ids()
            array([15, 16, 17, 18, 19])
        """
        return sgrid.top_edge_node_ids(self.shape)
        
    def bottom_edge_node_ids(self):
        """
        Returns a 1D numpy integer array containing the node ID numbers of the 
        nodes along the bottom (y=0) grid edge.
        
        Example:
            
            >>> rmg = RasterModelGrid(4, 5, 1.0)
            >>> rmg.bottom_edge_node_ids()
            array([0, 1, 2, 3, 4])
        """
        return sgrid.bottom_edge_node_ids(self.shape)
        
    def left_edge_node_ids(self):
        """
        Returns a 1D numpy integer array containing the node ID numbers of the 
        nodes along the left (x=0) grid edge.
        
        Example:
            
            >>> rmg = RasterModelGrid(4, 5, 1.0)
            >>> rmg.left_edge_node_ids()
            array([ 0,  5, 10, 15])
        """
        return sgrid.left_edge_node_ids(self.shape)
        
    def right_edge_node_ids(self):
        """
        Returns a 1D numpy integer array containing the node ID numbers of the 
        nodes along the right (x=xmax) grid edge.
        
        Example:
            
            >>> rmg = RasterModelGrid(4, 5, 1.0)
            >>> rmg.right_edge_node_ids()
            array([ 4,  9, 14, 19])
        """
        return sgrid.right_edge_node_ids(self.shape)
        
    def grid_coords_to_node_id(self, row, col, **kwds):
        """
        Returns the ID of the node at the specified *row* and *col* of the
        raster grid. Since this is a wrapper for the numpy ravel_multi_index
        function, the keyword arguments are the same as that function. In
        addition, *row* and *col* can both be either scalars or arrays (of the
        same length) to get multiple ids.

        As with ravel_multi_index use the *mode* keyword to change the
        behavior of the method when passed an out-of-range *row* or *col*.
        The default is to raise ValueError (not IndexError, as you might
        expect).
        
        .. note::
            The syntax assumes that first row and column are 0,
            so max entry for a mg with 4 rows and 5 cols is row=3, col=4
        
        Example:
            
            >>> mg = RasterModelGrid(4, 5)
            >>> mg.grid_coords_to_node_id(2, 3)
            13

            >>> mg.grid_coords_to_node_id([2, 0], [3, 4])
            array([13,  4])
        """
        return numpy.ravel_multi_index((row, col), self.shape, **kwds)
        
    def unit_test( self ):
        """
        This is just scratch space for testing while developing. More proper
        tests are in the doctests for each function, and in  
        test_raster_model_grid.py.
        """
        
        print 'Performing some tests for RasterModelGrid ...'
        print
        
        num_rows_for_unit_test = 4
        num_cols_for_unit_test = 5
        
        print 'Initializing ...'
        self.initialize( num_rows_for_unit_test, 
                         num_cols_for_unit_test,
                         1.0 )
        print 'done.'
        print
        
        #print 'Testing node lists:'
        #print 'ID   X    Y    Z    Status  Active_cell  #in in1 in2 #out out1 out2'
        #for node in range( 0, self.num_nodes ):
        #    print(str(node)+'    '+str(self._node_x[node])+'  '
        #          +str(self._node_y[node])+'  '
        #          +str(self._node_z[node])+'  '
        #          +str(self.node_status[node])+'  '
        #          +str(self.node_activecell[node])+'  '
        #          +str(self.node_numinlink[node])+'  '
        #          +str(self.node_inlink_matrix[0][node])+'  '
        #          +str(self.node_inlink_matrix[1][node])+'  '
        #          +str(self.node_numoutlink[node])+'  '
        #          +str(self.node_outlink_matrix[0][node])+'  '
        #          +str(self.node_outlink_matrix[1][node])+'  '
        #          +str(self.node_numactiveinlink[node])+'  '
        #          +str(self.node_active_inlink_matrix[0][node])+'  '
        #          +str(self.node_active_inlink_matrix[1][node])+'  '
        #          +str(self.node_numactiveoutlink[node])+'  '
        #          +str(self.node_active_outlink_matrix[0][node])+'  '
        #          +str(self.node_active_outlink_matrix[1][node]))
        #print
        #
        
        print 'Testing fluxes and in/out links:'
        flux2 = numpy.zeros(len(flux)+1)
        flux2[:len(flux)] = flux
        print flux2
        for n in range(0, self.num_nodes):
            mysum = -(flux2[self.node_active_inlink_matrix[0][n]] + \
                      flux2[self.node_active_inlink_matrix[1][n]]) + \
                     (flux2[self.node_active_outlink_matrix[0][n]] + \
                      flux2[self.node_active_outlink_matrix[1][n]])
            print(str(n)+' '+str(flux2[self.node_active_inlink_matrix[0][n]])
                        +' '+str(flux2[self.node_active_inlink_matrix[1][n]])
                        +' '+str(flux2[self.node_active_outlink_matrix[0][n]])
                        +' '+str(flux2[self.node_active_outlink_matrix[1][n]])
                        +' '+str(mysum))
        divg2 = self.calculate_flux_divergence_at_nodes(flux)
        print divg2


def _is_closed_boundary(boundary_string):
    
    return boundary_string.lower() == 'closed'


def from_dict(param_dict):
    """
    Create a RasterModelGrid from the dictionary-like object, *param_dict*.
    Required keys of the dictionary are NUM_ROWS, NUM_COLS. Raises a KeyError
    if either of these are missing.  If GRID_SPACING is given, use it as the
    HexModelGrid *dx* parameter, otherwise default to unit spacing.
    """
    # Read and create basic raster grid
    try:
        nrows = int(param_dict['NUM_ROWS'])
        ncols = int(param_dict['NUM_COLS'])
        dx = float(param_dict.get('GRID_SPACING', 1.))
    except KeyError:
        raise
    except ValueError:
        raise
    else:
        mg = RasterModelGrid(nrows, ncols, dx)
        
    # Set boundaries
    left_boundary_type = param_dict.get('LEFT_BOUNDARY', 'open')
    right_boundary_type = param_dict.get('RIGHT_BOUNDARY', 'open')
    top_boundary_type = param_dict.get('TOP_BOUNDARY', 'open')
    bottom_boundary_type = param_dict.get('BOTTOM_BOUNDARY', 'open')
    mg.set_inactive_boundaries(_is_closed_boundary(bottom_boundary_type), 
                               _is_closed_boundary(right_boundary_type),
                               _is_closed_boundary(top_boundary_type),
                               _is_closed_boundary(left_boundary_type))

    # Return the created and initialized grid
    return mg
