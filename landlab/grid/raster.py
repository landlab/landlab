#! /usr/bin/env python

import numpy
import numpy as np

from landlab.testing.decorators import track_this_method
from landlab.utils import structured_grid as sgrid
from landlab.utils import count_repeated_values

from .base import ModelGrid
from . import grid_funcs as gfuncs
from .base import (INTERIOR_NODE, FIXED_VALUE_BOUNDARY,
                   FIXED_GRADIENT_BOUNDARY, TRACKS_CELL_BOUNDARY,
                   INACTIVE_BOUNDARY, BAD_INDEX_VALUE, )
from . import raster_funcs as rfuncs


def node_has_boundary_neighbor(mg, id):
    for neighbor in mg.get_neighbor_list(id):
        if mg.node_status[neighbor] != INTERIOR_NODE:
            return True
    for neighbor in mg.get_diagonal_list(id):
        if mg.node_status[neighbor] != INTERIOR_NODE:
            return True
    return False


def make_arg_into_array(arg):
    ids = arg
    if not isinstance(ids, list) or not isinstance(ids, numpy.ndarray):
        try:
            ids = list(ids)
        except TypeError:
            ids = [ids]
    return ids

has_boundary_neighbor = numpy.vectorize(node_has_boundary_neighbor,
                                        excluded=['mg'])


class RasterModelGridPlotter(object):
    def imshow(self, group, var_name, **kwds):
        from landlab.plot import imshow_field
        kwds['values_at'] = group
        imshow_field(self, var_name, **kwds)


class RasterModelGrid(ModelGrid, RasterModelGridPlotter):
    """
    This inherited class implements a regular, raster 2D grid with uniform
    cell dimensions.
    
    Examples:
        
        >>> rmg = RasterModelGrid(4, 5, 1.0) # rows, columns, spacing
        >>> rmg.number_of_nodes
        20
    """

    def __init__(self, num_rows=0, num_cols=0, dx=1.0, **kwds):
        """
        Optionally takes numbers of rows and columns and cell size as
        inputs. If this are given, calls initialize() to set up the grid.
        
        ..todo: 
            the option for NOT giving rows, cols, and dx no longer works, 
            because the *field* init requires num_active_cells, etc., to be
            defined. Either we force users to give arguments on instantiation,
            or set it up such that one can create a zero-node grid.
        """
        # Set number of nodes, and initialize if caller has given dimensions
        self._num_nodes = num_rows * num_cols
        if self.number_of_nodes > 0:
            self._initialize(num_rows, num_cols, dx)
        super(RasterModelGrid, self).__init__(**kwds)

    def _initialize(self, num_rows, num_cols, dx):
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
            >>> numrows = 20          # number of rows in the grid
            >>> numcols = 30          # number of columns in the grid
            >>> dx = 10.0             # grid cell spacing
            >>> rmg = RasterModelGrid(numrows, numcols, dx)
            >>> rmg.number_of_nodes, rmg.number_of_cells, rmg.number_of_links, rmg.number_of_active_links
            (600, 504, 1150, 1054)
            >>> rmg = RasterModelGrid(4, 5)
            >>> rmg.number_of_nodes,rmg.number_of_cells,rmg.number_of_links,rmg.number_of_active_links
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
            array([[-1, -1, -1, -1, -1,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,
                    12, 13, 14],
                   [-1, 15, 16, 17, 18, -1, 19, 20, 21, 22, -1, 23, 24, 25, 26, -1, 27,
                    28, 29, 30]])
            >>> rmg.node_numoutlink
            array([2, 2, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0])
            >>> rmg.node_outlink_matrix[0]
            array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, -1, -1,
                   -1, -1, -1])
            >>> rmg.node_numactiveinlink
            array([0, 0, 0, 0, 0, 0, 2, 2, 2, 1, 0, 2, 2, 2, 1, 0, 1, 1, 1, 0])
            >>> rmg.node_active_inlink_matrix
            array([[-1, -1, -1, -1, -1, -1,  0,  1,  2, -1, -1,  3,  4,  5, -1, -1,  6,
                     7,  8, -1],
                   [-1, -1, -1, -1, -1, -1,  9, 10, 11, 12, -1, 13, 14, 15, 16, -1, -1,
                    -1, -1, -1]])
            >>> rmg.node_numactiveoutlink
            array([0, 1, 1, 1, 0, 1, 2, 2, 2, 0, 1, 2, 2, 2, 0, 0, 0, 0, 0, 0])
            >>> rmg.node_active_outlink_matrix
            array([[-1,  0,  1,  2, -1, -1,  3,  4,  5, -1, -1,  6,  7,  8, -1, -1, -1,
                    -1, -1, -1],
                   [-1, -1, -1, -1, -1,  9, 10, 11, 12, -1, 13, 14, 15, 16, -1, -1, -1,
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
            print 'RasterModelGrid._initialize('+str(num_rows)+', ' \
                   +str(num_cols)+', '+str(dx)+')'
        
        # Basic info about raster size and shape
        self._nrows = num_rows
        self._ncols = num_cols

        self._dx = dx
        self.cellarea = dx * dx

        self._num_nodes = sgrid.node_count(self.shape)
        self._num_active_nodes = self.number_of_nodes

        self._num_cells = sgrid.cell_count(self.shape)
        self._num_active_cells = self.number_of_cells

        self._num_links = sgrid.link_count(self.shape)
        self._num_active_links = sgrid.active_link_count(self.shape)
        
        self._num_faces = sgrid.face_count(self.shape)
        self._num_active_faces = sgrid.active_face_count(self.shape)

        # We need at least one row or column of boundary cells on each
        # side, so the grid has to be at least 3x3
        assert(numpy.min((num_rows, num_cols)) >= 3)

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
        self.cell_node = sgrid.node_index_at_cells(self.shape)
        self.node_activecell = sgrid.active_cell_index_at_nodes(self.shape)
        self.active_cells = sgrid.active_cell_index(self.shape)
        self.activecell_node = self.cell_node.copy()
        #self.active_faces = sgrid.active_face_index(self.shape)

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
         self.link_tonode) = sgrid.node_index_at_link_ends(self.shape)

        #   set up in-link and out-link matrices and numbers
        self._setup_inlink_and_outlink_matrices()
        
        # Flag indicating whether we have created diagonal links.
        self._diagonal_links_created = False
        
        #   set up the list of active links
        self._reset_list_of_active_links()

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
        self.link_face = sgrid.face_index_at_links(self.shape,
                                                   actives=self.active_links)

        # List of neighbors for each cell: we will start off with no
        # list. If a caller requests it via get_neighbor_list or
        # create_neighbor_list, we'll create it if necessary.
        self.neighbor_list_created = False

        # List of diagonal neighbors. As with the neighbor list, we'll only
        # create it if requested.
        self.diagonal_list_created = False
        

    def _setup_cell_areas_array(self):
        self.active_cell_areas = numpy.empty(self.number_of_active_cells)
        self.active_cell_areas.fill(self._dx ** 2)
        return self.active_cell_areas

    @property
    def shape(self):
        return (self.number_of_node_rows, self.number_of_node_columns)

    @property
    def dx(self):
        """
        Returns the node spacing of a raster grid. Same as node_spacing.
        Example: my_grid_spacing = my_raster_grid.dx 
        (no parentheses, because it is a property rather than a method)
        """
        return self._dx

    @property
    def node_spacing(self):
        """
        Returns the node spacing of a raster grid. Same as dx.
        Example: my_grid_spacing = my_raster_grid.node_spacing 
        (no parentheses, because it is a property rather than a method)
        """
        return self._dx

    def node_links(self, *args):
        """node_links([node_ids])

        Returns the ids of links attached to grid nodes with *node_ids*. If
        *node_ids* is not given, return links for all of the nodes in the
        grid. Link ids are listed in clockwise order starting with the south
        link. (i.e., [S,W,N,E])
        """
        if len(args) > 1:
            raise ValueError('only zero or one arguments accepted')

        try:
            node_ids = args[0]
        except IndexError: #return all nodes
            return numpy.vstack((self.node_inlink_matrix,
                                 self.node_outlink_matrix))
        else:
            try:
                float(node_ids) #single number (need to reshape the matrices)
            except:
                return numpy.vstack((self.node_inlink_matrix[:, node_ids],
                                 self.node_outlink_matrix[:, node_ids]))
            else:
                return numpy.hstack((self.node_inlink_matrix[:, node_ids],
                                 self.node_outlink_matrix[:, node_ids]))

    def active_node_links(self, *args):
        """active_node_links([node_ids])

        Returns the ids of links attached to active grid nodes with
        *node_ids*. If *node_ids* is not given, return links for all of the
        nodes in the grid. Link ids are listed in clockwise order starting
        with the south link.
        """
        if len(args) > 1:
            raise ValueError('only zero or one arguments accepted')

        try:
            node_ids = args[0]
        except IndexError: #return all nodes
            return numpy.vstack((self.node_active_inlink_matrix,
                                 self.node_active_outlink_matrix))
        else:
            try:
                float(node_ids) #single number (need to alter stacking)
            except:
                return numpy.vstack(
                    (self.node_active_inlink_matrix.take(node_ids, axis=1),
                    self.node_active_outlink_matrix.take(node_ids, axis=1)))
            else:
                return numpy.hstack(
                    (self.node_active_inlink_matrix.take(node_ids, axis=1),
                    self.node_active_outlink_matrix.take(node_ids, axis=1)))

    def _setup_inlink_and_outlink_matrices(self):
        """
        Creates data structures to record the numbers of inlinks and outlinks
        for each node. An inlink of a node is simply a link that has the node as
        its "to" node, and an outlink is a link that has the node as its "from".
        
        We store the inlinks in a 2-row by num_nodes-column matrix called
        node_inlink_matrix. It has two rows because we know that the nodes in
        our raster grid will never have more than two inlinks an two outlinks
        each (a given node could also have zero or one of either). The outlinks
        are stored in a similar matrix.
        
        The order of inlinks is [SOUTH, WEST].
        
        The order of outlinks is [NORTH, EAST].
        
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
         self.node_numinlink) = sgrid.setup_inlink_matrix(self.shape)

        (self.node_outlink_matrix,
         self.node_numoutlink) = sgrid.setup_outlink_matrix(self.shape)
        
    def _setup_active_inlink_and_outlink_matrices(self):
        """
        Creates data structures to record the numbers of active inlinks and 
        active outlinks for each node. These data structures are equivalent to
        the "regular" inlink and outlink matrices, except that it uses the IDs
        of active links (only).
        """ 

        node_status = self.node_status != INACTIVE_BOUNDARY

        (self.node_active_inlink_matrix,
         self.node_numactiveinlink) = sgrid.setup_active_inlink_matrix(
             self.shape, node_status=node_status)

        (self.node_active_outlink_matrix,
         self.node_numactiveoutlink) = sgrid.setup_active_outlink_matrix(
             self.shape, node_status=node_status)
             
             
    def _reset_list_of_active_diagonal_links(self):
        
        assert(self._diagonal_links_created), 'Diagonal links not created'
        
        self._diagonal_active_links = []
        self._diag_activelink_fromnode = []
        self._diag_activelink_tonode = []

        diag_fromnode_status = self.node_status[self._diag_link_fromnode]
        diag_tonode_status = self.node_status[self._diag_link_tonode]
        
        diag_active_links = (((diag_fromnode_status == INTERIOR_NODE) & ~
                              (diag_tonode_status == INACTIVE_BOUNDARY)) |
                             ((diag_tonode_status == INTERIOR_NODE) & ~
                              (diag_fromnode_status == INACTIVE_BOUNDARY)))

        (self._diag_active_links, ) = numpy.where(diag_active_links)

        self._num_diag_active_links = len(self._diag_active_links)
        self._diag_activelink_fromnode = self._diag_link_fromnode[self._diag_active_links]
        self._diag_activelink_tonode = self._diag_link_tonode[self._diag_active_links]


    def _reset_list_of_active_links(self):
        
        super(RasterModelGrid, self)._reset_list_of_active_links()
        if self._diagonal_links_created:
            self._reset_list_of_active_diagonal_links()
       

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
        cell_ids = make_arg_into_array(cell_id)
        node_ids = self.cell_node[cell_ids]
        inlinks = self.node_inlink_matrix[:, node_ids].T
        outlinks = self.node_outlink_matrix[:, node_ids].T
        return numpy.squeeze(numpy.concatenate(
            (self.link_face[inlinks], self.link_face[outlinks]), axis=1))

    def get_grid_xdimension(self):
        '''
        Returns the x dimension of the grid. Method added 5/1/13 by DEJH.
        '''
        return (self.number_of_node_columns * self._dx)
    
    def get_grid_ydimension(self):
        '''
        Returns the y dimension of the grid. Method added 5/1/13 by DEJH.
        '''
        return (self.number_of_node_rows * self._dx)

    @property
    def number_of_interior_nodes(self):
        """
        Returns the number of interior nodes on the grid.
        """
        return sgrid.interior_node_count(self.shape)

    @property
    def number_of_nodes(self):
        """
        Returns total number of nodes, including boundaries.
        """
        return self._num_nodes

    @property
    def number_of_node_columns(self):
        """
        Returns the number of columns, including boundaries.  
        """
        return self._ncols
        
    @property
    def number_of_node_rows(self):
        """
        Returns the number of rows, including boundaries.  
        """
        return self._nrows

    @property
    def node_spacing(self):
        """
        Returns the spacing between grid nodes.
        """
        return self._dx
    
    @property
    def corner_nodes(self):
        """
        Returns an array of the grid corner node IDs.
        """
        return sgrid.corners((self._nrows,self._ncols))

    def is_point_on_grid(self, xcoord, ycoord):
        """
        This method takes x,y coordinates and tests whether they lie within the
        grid. The limits of the grid are taken to be links connecting the 
        boundary nodes. We perform a special test to detect looped boundaries.
        
        Coordinates can be ints or arrays of ints. If arrays, will return an
        array of the same length of truth values.
        """
        x_condition = numpy.logical_and(numpy.less(0.,xcoord), numpy.less(xcoord,(self.get_grid_xdimension()-self._dx)))
        y_condition = numpy.logical_and(numpy.less(0.,ycoord), numpy.less(ycoord,(self.get_grid_ydimension()-self._dx)))
        if numpy.all(self.node_status[sgrid.left_edge_node_ids(self.shape)]==3) or numpy.all(self.node_status[sgrid.right_edge_node_ids(self.shape)]==3):
            try:
                x_condition[:] = 1
            except:
                x_condition = 1
        if numpy.all(self.node_status[sgrid.top_edge_node_ids(self.shape)]==3) or numpy.all(self.node_status[sgrid.bottom_edge_node_ids(self.shape)]==3):
            try:
                y_condition[:] = 1
            except:
                y_condition = 1
        return numpy.logical_and(x_condition, y_condition)

    def get_nodes_around_point(self, xcoord, ycoord):
        """
        This method takes an x,y coordinate within the grid, then returns the 
        IDs of the four nodes of the area (enclosure?) around that point as a 
        4 item array, in order [SOUTHWEST,SOUTHEAST,NORTHWEST,NORTHEAST], i.e.,
        ID order. Because the geometry of this grid is so simple, it works 
        purely by counting the number of squares left and below the point.
        If xcoord and ycoord are arrays, returns a (4,len_array) array. 
        Method added 4/29/13 by DEJH, modified 9/24/13.
        """   
        #try:
        #    assert (len(xcoord) == len(ycoord))
        #except:
        #    assert type(xcoord) == float and type(ycoord) == float
        ID = ycoord//self._dx * self.number_of_node_columns + xcoord//self._dx
        try:
            ID = int(ID)
        except:
            ID = ID.astype(int)
        return numpy.array([ID, ID + self.number_of_node_columns,
                            ID + self.number_of_node_columns + 1, ID + 1])
        
    def snap_coords_to_grid(self, xcoord, ycoord):
        '''
        This method takes existing coordinates, inside the grid, and returns
        the ID of the closest grid node. That node can be active or inactive.
        DEJH, 9/24/13.
        '''
        #This testing suppressed for speed. While suppressed, coordinates provided MUST be within the grid or silent instability will occur.
        #if type(xcoord) == int:
        #    if not self.is_point_on_grid(xcoord, ycoord):
        #        raise LookupError('Coordinates specified are outside the grid area')
        #else: #it's an array
        #    if not numpy.all(self.is_point_on_grid(xcoord, ycoord)):
        #        raise LookupError('One or more pairs of coordinates specified are outside the grid area')
        vertices_array = self.get_nodes_around_point(xcoord, ycoord)
        #print vertices_array
        #vertices_array.reshape((4,-1))
        #print self.node_x[vertices_array], self.node_y[vertices_array]
        xdir_displacement = numpy.tile(xcoord,(4,1)) - self.node_x[vertices_array]
        ydir_displacement = numpy.tile(ycoord,(4,1)) - self.node_y[vertices_array]
        distances_to_vertices = numpy.sqrt(xdir_displacement*xdir_displacement + ydir_displacement*ydir_displacement)
        try:
            return vertices_array[(numpy.argmin(distances_to_vertices, axis=0), xrange(distances_to_vertices.shape[1]))]
        except:
            return vertices_array[numpy.argmin(distances_to_vertices)]
        #...per fancy indexing
    
    def min_active_link_length(self):
        """
        Returns the horizontal length of the shortest active link in the grid.
        Overrides ModelGrid.min_active_link_length().
        """
        return self._dx

    def max_active_link_length(self):
        return self._dx

    def calculate_gradient_across_cell_faces(self, *args, **kwds):
        return rfuncs.calculate_gradient_across_cell_faces(
            self, *args, **kwds)

    def calculate_gradient_across_cell_corners(self, *args, **kwds):
        return rfuncs.calculate_gradient_across_cell_corners(
            self, *args, **kwds)
            
            
    def _setup_diagonal_links(self):
        """
        Creates lists of from and to nodes for diagonal links. A diagonal link
        is a special type of link that connects the diagonal of two raster cells.
        One use for diagonal links has to do with raster digital elevation 
        models: diagonal links allow you to implement "D8" drainage-routing
        algorithms, in which each node is considered to have 8 rather than 4
        neighbors, and flow will go toward whichever of these neighbors lies in
        the steepest downslope direction.
        """
        n_diagonal_links = 2*(self._nrows-1)*(self._ncols-1)
        self._diag_link_fromnode = numpy.zeros(n_diagonal_links, dtype=int)
        self._diag_link_tonode = numpy.zeros(n_diagonal_links, dtype=int)
        i = 0
        for r in range(self._nrows-1):
            for c in range(self._ncols-1):
                self._diag_link_fromnode[i] = c+r*self._ncols
                self._diag_link_tonode[i] = (c+1)+(r+1)*self._ncols
                i += 1
                self._diag_link_fromnode[i] = (c+1)+r*self._ncols
                self._diag_link_tonode[i] = c+(r+1)*self._ncols
                i += 1
                
        self._diagonal_links_created = True
        
        self._reset_list_of_active_diagonal_links()
        
        
    def d8_active_links(self):
        """
        Returns a set of active links that include diagonal connections between
        grid cells, for use with link-based water-routing schemes.
        """
        
        if not self._diagonal_links_created:
            self._setup_diagonal_links()
            
        return self._diag_active_links, \
               numpy.concatenate((self.activelink_fromnode, 
                                  self._diag_activelink_fromnode)), \
               numpy.concatenate((self.activelink_tonode,
                                  self._diag_activelink_tonode))
        

    def calculate_max_gradient_across_cell_faces(self, *args, **kwds):
        """rmg.calculate_max_gradient_across_cell_faces(node_values, [cell_ids], return_face=False)

        Return the maximum gradients across cell faces.

        Refer to :func:`~landlab.grid.raster_funcs.calculate_max_gradient_across_cell_faces`
        for full documentation.

        .. seealso::

            :func:`~landlab.grid.raster_funcs.calculate_max_gradient_across_cell_faces` : equivalent function
        """
        return rfuncs.calculate_max_gradient_across_cell_faces(self, *args,
                                                               **kwds)

    def calculate_max_gradient_across_cell_corners(self, *args, **kwds):
        """
        Return the maximum gradients across diagonal cells.

        Refer to :func:`~landlab.grid.raster_funcs.calculate_max_gradient_across_cell_corners`
        for full documentation.

        .. seealso::

            :func:`~landlab.grid.raster_funcs.calculate_max_gradient_across_cell_corners` : equivalent function
        """
        return rfuncs.calculate_max_gradient_across_cell_corners(self, *args,
                                                                 **kwds)

    def calculate_max_gradient_across_adjacent_cells(self, node_values, *args,
                                                     **kwds):
        """
        Return the maximum gradients across adjacent cells.

        Refer to :func:`~landlab.grid.raster_funcs.calculate_max_gradient_across_adjacent_cells`
        for full documentation.

        .. seealso::

            :func:`~landlab.grid.raster_funcs.calculate_max_gradient_across_adjacent_cells` : equivalent function
        """
        return rfuncs.calculate_max_gradient_across_adjacent_cells(
            self, node_values, *args, **kwds)

    def calculate_max_gradient_across_node(self, u, cell_id):
        """
        .. deprecated:: 0.1
            Use :func:`calculate_max_gradient_across_adjacent_cells`

        This method calculates the gradients in u across all 4 faces of the 
        cell with ID cell_id, and across the four diagonals. It then returns 
        the steepest (most negative) of these values, followed by its dip 
        direction (e.g.: 0.12, 225). i.e., this is a D8 algorithm. Slopes 
        downward from the cell are reported as positive.
            
        This code is actually calculating slopes, not gradients.  
        The max gradient is the most negative, but the max slope is the most
        positive.  So, this was updated to return the max value, not the 
        min.

        """
        return rfuncs.calculate_max_gradient_across_node(self, u, cell_id)
        
    def calculate_max_gradient_across_node_d4(self, u, cell_id):
        """
        .. deprecated:: 0.1
            Use :func:`calculate_max_gradient_across_cell_faces` instead

        This method calculates the gradients in u across all 4 faces of the 
        cell with ID cell_id. It then returns 
        the steepest (most negative) of these values, followed by its dip 
        direction (e.g.: 90 180). i.e., this is a D4 algorithm. Slopes 
        downward from the cell are reported as positive.
            
        Note that this is exactly the same as calculate_max_gradient_across_node
        except that this is d4, and the other is d8.
            
        This code is actually calculating slopes, not gradients.  
        The max gradient is the most negative, but the max slope is the most
        positive.  So, this was updated to return the max value, not the 
        min.
        """
        return rfuncs.calculate_max_gradient_across_node_d4(self, u, cell_id)
        
    def find_node_in_direction_of_max_slope(self, u, node_id):
        """
        .. deprecated:: 0.1
            Use :func:`calculate_max_gradient_across_adjacent_cells` instead

        This method calculates the slopes (-dz/dx) in u across all 4 faces of 
        the cell with ID node_id, and across the four diagonals. 
        It then returns the node ID in the direction of the steepest 
        (most positive) of these values,  i.e., this is a 
        D8 algorithm. Slopes downward from the cell are reported as positive.
            
        This doesn't deal with the fixed gradient boundary condition.  
        """

        # NMG Update.  This is super clumsy. 
            
        # DEJH update: Gets confused for the lowest node if w/i grid
        # (i.e., closed)- will return a higher neighbour, when it should
        # return a null index ->  Now returns -1.

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
                single_slope = (u[node_id] - u[a])/self.dx
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
        #...now if a node is the lowest thing, this method returns -1, not a neighbor:
        if u[steepest_node] > u[node_id]:
            steepest_node=-1
        
        return steepest_node
    
    def find_node_in_direction_of_max_slope_d4(self, u, node_id):
        """
        .. deprecated:: 0.1
            Use :func:`calculate_max_gradient_across_adjacent_cells` instead
        
        This method is exactly the same as find_node_in_direction_of_max_slope
        except that this method only considers nodes that are connected by links,
        or in otherwords, in the 0, 90, 180 and 270 directions.
        
        This method calculates the slopes (-dz/dx) in u across all 4 faces of 
        the cell with ID node_id. 
        It then returns the node ID in the direction of the steepest 
        (most positive) of these values,  i.e., this is a 
        D8 algorithm. Slopes downward from the cell are reported as positive.
            
        This doesn't deal with the fixed gradient boundary condition.  
        """

        # NMG Update.  This is super clumsy. 
            
        # DEJH update: Gets confused for the lowest node if w/i grid
        # (i.e., closed)- will return a higher neighbour, when it should
        # return a null index ->  Now returns -1.

        #We have poor functionality if these are closed boundary nodes! 
        neighbor_nodes = self.get_neighbor_list(node_id)
        neighbor_nodes.sort()
        #print 'Node is internal: ', self.is_interior(cell_id)
        #print 'Neighbor cells: ', neighbor_cells
        slopes = []
        for a in neighbor_nodes:
            if self.node_status[a] != INACTIVE_BOUNDARY:
                single_slope = (u[node_id] - u[a])/self.dx
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

        #print 'Slopes list: ', slopes
        if slopes:
            max_slope, index_max = max((max_slope, index_max) for (index_max, max_slope) in enumerate(slopes))
        else:
            print u
            print 'Returning NaN angle and direction...'
            max_slope = numpy.nan
            index_max = 4
        
        #all_neighbor_nodes=numpy.concatenate((neighbor_nodes,diagonal_nodes))
        #print 'all_neighbor_cells ', all_neighbor_cells
        
        #Final check to  allow correct handling of internally draining nodes; DEJH Aug 2013.
        #This remains extremely ad-hoc. An internal node points to itself, but this should never
        #be used to actually route flow. In flow_accumulation, there is an explicit check that flow
        #is not routed to yourself.
        steepest_node = neighbor_nodes[index_max]
        #...now if a node is the lowest thing, this method returns -1, not a neighbor:
        if u[steepest_node] > u[node_id]:
            steepest_node=-1
        
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
        >>> rmg.number_of_active_links
        17
        >>> rmg.node_status
        array([1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1], dtype=int8)
        >>> rmg.set_inactive_boundaries(False, False, True, True)
        >>> rmg.number_of_active_links
        12
        >>> rmg.node_status
        array([1, 1, 1, 1, 1, 4, 0, 0, 0, 1, 4, 0, 0, 0, 1, 4, 4, 4, 4, 4], dtype=int8)
        
        Note that the four corners are treated as follows:
            bottom left = BOTTOM
            bottom right = BOTTOM
            top right = TOP
            top left = TOP
        This scheme is necessary for internal consistency with looped boundaries.
        """
        if self.DEBUG_TRACK_METHODS:
            print 'ModelGrid.set_inactive_boundaries'
            
        bottom_edge = range(0, self.number_of_node_columns)
        right_edge = range(2*self.number_of_node_columns - 1,
                           self.number_of_nodes - 1,
                           self.number_of_node_columns)
        top_edge = range((self.number_of_node_rows - 1) *
                         self.number_of_node_columns, self.number_of_nodes)
        left_edge = range(self.number_of_node_columns, 
                         self.number_of_nodes-self.number_of_node_columns,
                          self.number_of_node_columns)
            
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
        
        self._reset_list_of_active_links()
        
        
    def set_looped_boundaries(self, top_bottom_are_looped,sides_are_looped):
        """
        Added DEJH Feb 2014
        Handles boundary conditions by setting corresponding parallel grid edges
        as looped "tracks_cell" (==3) status, linked to each other. If top_bottom_are_looped 
        is True, the top and bottom edges will link to each other. If sides_are_
        looped is True, the left and right edges will link to each other.
        
        Note that because of the symmetries this BC implies, the corner nodes
        are all paired with the bottom/top edges, not the sides.
        
        >>> rmg = RasterModelGrid(4, 5, 1.0) # rows, columns, spacing
        >>> rmg.number_of_active_links
        17
        >>> rmg.node_status
        array([1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1], dtype=int8)
        >>> rmg.create_node_array_zeros('planet_surface__elevation')
        array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
                0.,  0.,  0.,  0.,  0.,  0.,  0.])
        >>> rmg.set_looped_boundaries(True, True)
        >>> rmg.looped_node_properties['boundary_node_IDs']
        array([ 0,  1,  2,  3,  4,  5,  9, 10, 14, 15, 16, 17, 18, 19])
        >>> rmg.looped_node_properties['linked_node_IDs']
        array([10, 11, 12, 13, 14,  8,  6, 13, 11,  5,  6,  7,  8,  9])
        
            15  16  17  18  19
            10  11  12  13  14
             5   6   7   8   9
             0   1   2   3   4

        TODO: Assign BC_statuses also to *links*
        """
        
        bottom_edge = numpy.array(range(0, self.number_of_node_columns))
        right_edge = numpy.array(range(2 * self.number_of_node_columns  - 1,
                           self.number_of_nodes - 1,
                           self.number_of_node_columns))
        top_edge = numpy.array(range((self.number_of_node_rows - 1) *
                          self.number_of_node_columns, self.number_of_nodes))
        left_edge = numpy.array(range(self.number_of_node_columns,
                          (self.number_of_nodes - self.number_of_node_columns),
                          self.number_of_node_columns))
        these_boundary_IDs = numpy.array([])
        these_linked_nodes = numpy.array([])
        
        if top_bottom_are_looped:
            self.node_status[bottom_edge] = TRACKS_CELL_BOUNDARY
            self.node_status[top_edge] = TRACKS_CELL_BOUNDARY
            these_boundary_IDs = numpy.concatenate((these_boundary_IDs,
                                   bottom_edge, top_edge))
            these_linked_nodes = numpy.concatenate((
                                   these_linked_nodes,
                                   top_edge-self.number_of_node_columns,
                                   bottom_edge+self.number_of_node_columns))

        if sides_are_looped:
            self.node_status[right_edge] = TRACKS_CELL_BOUNDARY
            self.node_status[left_edge] = TRACKS_CELL_BOUNDARY
            these_boundary_IDs = numpy.concatenate((these_boundary_IDs,
                                   left_edge, right_edge))
            these_linked_nodes = numpy.concatenate((
                                   these_linked_nodes,
                                   right_edge-1, left_edge+1))

        self._reset_list_of_active_links()
        
        try:
            type(self.looped_node_properties)
        except:
            existing_IDs = numpy.array([])
            existing_links = numpy.array([])
        else:
            unrepeated_node_entries = numpy.logical_not(numpy.in1d(self.looped_node_properties['boundary_node_IDs'], these_linked_nodes))
            existing_IDs = self.looped_node_properties['boundary_node_IDs'][unrepeated_node_entries]
            existing_links = self.looped_node_properties['linked_node_IDs'][unrepeated_node_entries]
                        
        self.looped_node_properties = {}
        all_the_IDs = numpy.concatenate((these_boundary_IDs, existing_IDs))
        ID_ordering = numpy.argsort(all_the_IDs)
        self.looped_node_properties['boundary_node_IDs'] = all_the_IDs[ID_ordering].astype(int)
        self.looped_node_properties['linked_node_IDs'] = numpy.concatenate((these_linked_nodes,existing_links))[ID_ordering].astype(int)

        if numpy.any(self.node_status[self.looped_node_properties['boundary_node_IDs']] == 2):
            raise AttributeError('Switching a boundary between fixed gradient and looped will result in bad BC handling! Bailing out...')        


    def set_fixed_gradient_boundaries(self, bottom_is_fixed,
            right_is_fixed, top_is_fixed, left_is_fixed, gradient_in=numpy.nan,
            gradient_of='planet_surface__elevation'):
        """
        Added DEJH Jan 2014
        Handles boundary conditions by setting each of the four sides of the 
        rectangular grid to 'active (fixed gradient)' (==2) status.
        Arguments are booleans indicating whether the bottom, right, top, and
        left are fixed gradient (True) or fixed value (False).
        
        For an fixed gradient boundary:
            - the nodes on the specified edges are flagged
              FIXED_GRADIENT_BOUNDARY (== 2). Other edges are ignored, and
              presumed to be set elsewhere.
              
            - the links between them and the adjacent interior nodes are
              active, but the links between each other are not. Corners are an
              awkward special case; they do not have any active links, but when
              boundary conditions are updated, each corner has a "pseudo-active"
              link that connects to one of its edge neighbors that lets it 
              update (see examples below).
              
            - the gradient is assumed by default to be the surface elevation,
              and this is assumed to be named "planet_surface__elevation" in the
              grid. If the gradient is in another surface, or the elevation
              surface is named differently, you need to set 'gradient_of' equal
              to the relevant string. self.fixed_gradient_of stores this string
              for access elsewhere.
              
            - The critical IDs and values relevant to the boundary conditions
              are stored in two special dictionaries,
                grid.fixed_gradient_link_properties, &
                grid.fixed_gradient_node_properties.
              The link dictionary stores the fixed gradients and the IDs of the
              links these are defined on: 'boundary_link_gradients', 
              'boundary_link_IDs'.
              The node dictionary stores the nodes on the boundary which are
              fixed gradient, the nodes at the other end of the links from each
              of these nodes, and the (fixed) value differences between these 
              two pairs of nodes: 'boundary_node_IDs', 'anchor_node_IDs',
              'values_to_add'. These can be used to update the boundaries fast,
              without needing to interrogate the links.
              
            - if *gradient* is provided, either as a float or an as a iterable
              of length number_of_boundary_nodes, then 'boundary_link_gradients'
              is set equal to *gradient*, and all the other properties updated 
              using these values. If it is not, then this method will
              attempt to access the link gradients and/or node elevations which
              were already in the grid when the method was called (i.e., the
              initial conditions), and use these to set the properties.
              Remember, gradient is as a fractional slope, not radians or
              degrees, and downslope gradients are positive!
              If gradient is a float, this method will assume you mean downslope
              flow out of all the edges of the grid. If you want some edges
              pointing in and some out, you'll need to call the function more
              than once, or provide an array of values.
              
            - If initial conditions are present in the grid ::and:: *gradient*
              is set, *gradient* will override the initial conditions provided.
              
            - if *gradient* is not provided (or is the wrong length), and
              initial conditions have not yet been set, the method will raise an
              exception.
            
            - Note that the four corners are treated as follows:
                bottom left = BOTTOM
                bottom right = BOTTOM
                top right = TOP
                top left = TOP,
              ...and the gradient on the link (if supplied) corresponds to the 
              link which points in the same direction as the rest of its edge 
              (i.e., the fixed gradient links of the bottom left and right 
              corners point up). This handling is necessary for internal
              consistency with looped BCs.
                      
        The following example sets all boundaries as fixed gradient in a
        four-row by five-column grid, but does so three times. The first time,
        initial conditions are allowed to set the fixed value. The second time,
        this is overridden by setting *gradient* in the function call as a
        constant. The third time, values are specified in an array:
        
        >>> rmg = RasterModelGrid(4, 5, 1.0) # rows, columns, spacing
        >>> rmg.number_of_active_links
        17
        >>> rmg.node_status
        array([1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1], dtype=int8)
        >>> rmg.create_node_array_zeros('planet_surface__elevation')
        array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
                0.,  0.,  0.,  0.,  0.,  0.,  0.])
        >>> rmg['node']['planet_surface__elevation'] += 1.
        >>> rmg['node']['planet_surface__elevation'][sgrid.boundary_nodes(rmg.shape)] = 0.8
        >>> rmg.set_fixed_gradient_boundaries(True, True, True, True) #first case
        Fixed gradients will be set according to existing data in the grid...
        >>> rmg.node_status
        array([2, 2, 2, 2, 2, 2, 0, 0, 0, 2, 2, 0, 0, 0, 2, 2, 2, 2, 2, 2], dtype=int8)
        >>> rmg.fixed_gradient_of
        'planet_surface__elevation'
        >>> rmg.fixed_gradient_node_properties['boundary_node_IDs']
        array([ 0,  1,  2,  3,  4,  9, 14, 15, 16, 17, 18, 19,  5, 10])
        >>> rmg.fixed_gradient_link_properties['boundary_link_IDs']
        array([ 0,  1,  2,  3,  4, 22, 26, 10, 11, 12, 13, 14, 19, 23])
        >>> rmg.fixed_gradient_link_properties['boundary_link_gradients']
        array([-0. , -0.2, -0.2, -0.2, -0. ,  0.2,  0.2, -0. ,  0.2,  0.2,  0.2,
               -0. , -0.2, -0.2])
        >>> rmg.set_fixed_gradient_boundaries(True, True, True, True, 0.1, gradient_of='planet_surface__elevation') #second case
        >>> rmg.fixed_gradient_link_properties['boundary_link_gradients']
        array([-0. , -0.1, -0.1, -0.1, -0. ,  0.1,  0.1, -0. ,  0.1,  0.1,  0.1,
               -0. , -0.1, -0.1])
        >>> rmg['node']['planet_surface__elevation']
        array([ 0.9,  0.9,  0.9,  0.9,  0.9,  0.9,  1. ,  1. ,  1. ,  0.9,  0.9,
                1. ,  1. ,  1. ,  0.9,  0.9,  0.9,  0.9,  0.9,  0.9])
                
        -> All boundaries end up with the same dip outwards. Note that the
        corners have been automatically set with "true" gradients of 0., so they
        mimic their edge neighbor. This is almost always what you want to
        happen.
        
        >>> my_gradients = numpy.array([0.5,0.5,0.5,0.5,]) #remember these are in edge, then ID order, with the corners attached to the other edges
        >>> rmg.set_fixed_gradient_boundaries(False, True, False, True, my_gradients) #third case
        >>> rmg.fixed_gradient_link_properties['boundary_link_gradients']
        array([ 0.5,  0.5,  0.5,  0.5, -0.6, -0.1, -0.1, -0.1,  0.4,  0.6,  0.1,
                0.1,  0.1, -0.4])
        >>> rmg.fixed_gradient_node_properties['boundary_node_IDs']
        array([ 9, 14,  5, 10,  0,  1,  2,  3,  4, 15, 16, 17, 18, 19])
        >>> rmg.fixed_gradient_node_properties['anchor_node_IDs']
        array([ 8, 13,  6, 11,  6,  6,  7,  8,  8, 11, 11, 12, 13, 13])
        >>> rmg.fixed_gradient_node_properties['values_to_add']
        array([-0.5, -0.5,  0.5,  0.5, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1,
               -0.1, -0.1, -0.1])
        >>> rmg['node']['planet_surface__elevation']
        array([ 0.9,  0.9,  0.9,  0.9,  0.9,  1.5,  1. ,  1. ,  1. ,  0.5,  1.5,
                1. ,  1. ,  1. ,  0.5,  0.9,  0.9,  0.9,  0.9,  0.9])
        
        ...i.e.,    0.9  0.9  0.9  0.9  0.9
                    1.5  1.   1.   1.   0.5
                    1.5  1.   1.   1.   0.5
                    0.9  0.9  0.9  0.9  0.9
        
        Now note we can easily update these boundary conditions much faster:
        
        >>> elevs = rmg['node']['planet_surface__elevation']
        >>> updated_elevs = elevs
        >>> updated_elevs[rmg.fixed_gradient_node_properties['boundary_node_IDs']] = updated_elevs[rmg.fixed_gradient_node_properties['anchor_node_IDs']] + rmg.fixed_gradient_node_properties['values_to_add']
        >>> numpy.all(numpy.equal(elevs, updated_elevs))
        True
        
        """
        
        bottom_edge = range(0, self.number_of_node_columns)
        if type(bottom_edge) != list:
            bottom_edge = numpy.array([bottom_edge])
        else:
            bottom_edge = numpy.array(bottom_edge)
        right_edge = range(2*self.number_of_node_columns - 1,
                           self.number_of_nodes - 1,
                           self.number_of_node_columns)
        if type(right_edge) != list:
            right_edge = numpy.array([right_edge])
        else:
            right_edge = numpy.array(right_edge)
        top_edge = range((self.number_of_node_rows - 1) *
                         self.number_of_node_columns, self.number_of_nodes)
        if type(top_edge) != list:
            top_edge = numpy.array([top_edge])
        else:
            top_edge = numpy.array(top_edge)
        left_edge = range(self.number_of_node_columns, 
                          self.number_of_nodes - self.number_of_node_columns,
                          self.number_of_node_columns)
        if type(left_edge) != list:
            left_edge = numpy.array([left_edge])
        else:
            left_edge = numpy.array(left_edge)

        fixed_gradient_nodes = numpy.array([], dtype=int)
        fixed_gradient_linked_nodes = numpy.array([], dtype=int)
        boundary_links = numpy.array([], dtype=int)
        #fixed_gradient_values_to_add = numpy.array([], dtype=float)
        
        if bottom_is_fixed:
            self.node_status[bottom_edge] = FIXED_GRADIENT_BOUNDARY
            fixed_gradient_nodes = numpy.concatenate((fixed_gradient_nodes, bottom_edge))
            bottom_anchor_nodes = bottom_edge+self.number_of_node_columns
            bottom_anchor_nodes[0] = bottom_edge[1]+self.number_of_node_columns
            bottom_anchor_nodes[-1] = bottom_edge[-2]+self.number_of_node_columns
            fixed_gradient_linked_nodes = numpy.concatenate((fixed_gradient_linked_nodes, bottom_anchor_nodes))
            bottom_links = self.node_links(bottom_edge)[2,:]
            boundary_links = numpy.concatenate((boundary_links, bottom_links))
        if right_is_fixed:
            self.node_status[right_edge] = FIXED_GRADIENT_BOUNDARY
            fixed_gradient_nodes = numpy.concatenate((fixed_gradient_nodes, right_edge))
            right_anchor_nodes = right_edge-1
            fixed_gradient_linked_nodes = numpy.concatenate((fixed_gradient_linked_nodes, right_anchor_nodes))
            right_links = self.node_links(right_edge)[1,:]
            boundary_links = numpy.concatenate((boundary_links, right_links))
        if top_is_fixed:
            self.node_status[top_edge] = FIXED_GRADIENT_BOUNDARY
            fixed_gradient_nodes = numpy.concatenate((fixed_gradient_nodes, top_edge))
            top_anchor_nodes = top_edge-self.number_of_node_columns
            top_anchor_nodes[0] = top_edge[1]-self.number_of_node_columns
            top_anchor_nodes[-1] = top_edge[-2]-self.number_of_node_columns
            fixed_gradient_linked_nodes = numpy.concatenate((fixed_gradient_linked_nodes, top_anchor_nodes))
            top_links = self.node_links(top_edge)[0,:]
            boundary_links = numpy.concatenate((boundary_links, top_links))
        if left_is_fixed:
            self.node_status[left_edge] = FIXED_GRADIENT_BOUNDARY
            fixed_gradient_nodes = numpy.concatenate((fixed_gradient_nodes, left_edge))
            left_anchor_nodes = left_edge+1
            fixed_gradient_linked_nodes = numpy.concatenate((fixed_gradient_linked_nodes, left_anchor_nodes))
            left_links = self.node_links(left_edge)[3,:]
            boundary_links = numpy.concatenate((boundary_links, left_links))
        
        self._reset_list_of_active_links()
        
        try:
            no_val_provided = numpy.all(numpy.isnan(gradient_in))
        except:
            no_val_provided = False
        if no_val_provided:
            #Set the gradients by reference to existing data on grid
            print 'Fixed gradients will be set according to existing data in the grid...'
            fixed_gradient_array = -self.calculate_gradients_at_links(self['node'][gradient_of])[boundary_links] #this grid func gives slopes UP as positive
            fixed_gradient_values_to_add = self['node'][gradient_of][fixed_gradient_nodes] - self['node'][gradient_of][fixed_gradient_linked_nodes]
            
        else:
            try:
                fixed_gradient = float(gradient_in)
            except:
                try:
                    fixed_gradient_array = numpy.array(gradient_in) #an iterable of gradients was supplied
                except TypeError:
                    raise TypeError('The supplied gradient parameter must be a single number or iterable of length number_of_fixed_gradient_nodes')
                self.force_boundaries_from_gradients(boundary_links, fixed_gradient_array, gradient_of)
                fixed_gradient_values_to_add = self['node'][gradient_of][fixed_gradient_nodes] - self['node'][gradient_of][fixed_gradient_linked_nodes]
            else:
                #the supplied gradient was a single number
                fixed_gradient_array = numpy.array([], dtype=float)
                if bottom_is_fixed:
                    bottom_fixed_gradients = numpy.ones(bottom_edge.size, dtype=float)*-fixed_gradient
                    #force the corner gradient:
                    bottom_fixed_gradients[0] = 0.
                    bottom_fixed_gradients[-1] = 0.
                    fixed_gradient_array = numpy.concatenate((fixed_gradient_array, bottom_fixed_gradients))
                if right_is_fixed:
                    right_fixed_gradients = numpy.ones(right_edge.size,dtype=float)*fixed_gradient
                    #right_fixed_gradients[0] = 0.
                    fixed_gradient_array = numpy.concatenate((fixed_gradient_array, right_fixed_gradients))
                if top_is_fixed:
                    top_fixed_gradients = numpy.ones(top_edge.size,dtype=float)*fixed_gradient
                    top_fixed_gradients[0] = 0.
                    top_fixed_gradients[-1] = 0.
                    fixed_gradient_array = numpy.concatenate((fixed_gradient_array, top_fixed_gradients))
                if left_is_fixed:
                    left_fixed_gradients = numpy.ones(left_edge.size,dtype=float)*-fixed_gradient
                    #left_fixed_gradients[-1] = 0.
                    fixed_gradient_array = numpy.concatenate((fixed_gradient_array, left_fixed_gradients))
                self.force_boundaries_from_gradients(boundary_links, fixed_gradient_array, gradient_of)
                fixed_gradient_values_to_add = self['node'][gradient_of][fixed_gradient_nodes] - self['node'][gradient_of][fixed_gradient_linked_nodes]

        #Now we need to save the various fixed_gradient property arrays to the
        #grid, but making sure we don't duplicate any entries that might already
        #be in there from a previous (diff constant gradient?) run of this
        #method...
        try:
            self.fixed_gradient_node_properties['boundary_node_IDs']
        except AttributeError:
            #easy case; there's nothing there already
            self.fixed_gradient_of = gradient_of
            self.fixed_gradient_node_properties = {}
            self.fixed_gradient_link_properties = {}
            self.fixed_gradient_node_properties['boundary_node_IDs'] = fixed_gradient_nodes.astype(int)
            self.fixed_gradient_node_properties['anchor_node_IDs'] = fixed_gradient_linked_nodes.astype(int)
            self.fixed_gradient_node_properties['values_to_add'] = fixed_gradient_values_to_add
            self.fixed_gradient_link_properties['boundary_link_IDs'] = boundary_links.astype(int)
            #Update the link gradients over whole grid, as if there's values in the grid already, there could be compatibility issues...
            self.fixed_gradient_link_properties['boundary_link_gradients'] = -self.calculate_gradients_at_links(self['node'][gradient_of])[boundary_links]

        else:
            #there's something in there, which we need to merge with, and
            #overwrite some entries of if appropriate.
            if self.fixed_gradient_of != gradient_of:
                raise ValueError('At the moment, you have to define all your boundaries on the same set of values!') #need to sort this ASAP...
                #...We probably want the syntax to be rmg.BCs['process_module']['node'][gradient_of] as AN OBJECT, to which we can pin these properties
            #The fixed_gradient_nodes should be uniquely defined...
            unrepeated_node_entries = numpy.logical_not(numpy.in1d(self.fixed_gradient_node_properties['boundary_node_IDs'], fixed_gradient_nodes))
            unrepeated_link_entries = numpy.logical_not(numpy.in1d(self.fixed_gradient_link_properties['boundary_link_IDs'], boundary_links))
            fixed_gradient_array = numpy.concatenate((fixed_gradient_array, self.fixed_gradient_link_properties['boundary_link_gradients'][unrepeated_link_entries]))
            boundary_links = numpy.concatenate((boundary_links, self.fixed_gradient_link_properties['boundary_link_IDs'][unrepeated_link_entries]))
            fixed_gradient_nodes = numpy.concatenate((fixed_gradient_nodes, self.fixed_gradient_node_properties['boundary_node_IDs'][unrepeated_node_entries]))
            fixed_gradient_linked_nodes = numpy.concatenate((fixed_gradient_linked_nodes, self.fixed_gradient_node_properties['anchor_node_IDs'][unrepeated_node_entries]))
            fixed_gradient_values_to_add = numpy.concatenate((fixed_gradient_values_to_add, self.fixed_gradient_node_properties['values_to_add'][unrepeated_node_entries]))

            if numpy.any(self.node_status[fixed_gradient_nodes] == 3):
                raise AttributeError('Switching a boundary between fixed gradient and looped will result in bad BC handling! Bailing out...')

            self.fixed_gradient_node_properties = {}
            self.fixed_gradient_link_properties = {}
            self.fixed_gradient_node_properties['boundary_node_IDs'] = fixed_gradient_nodes.astype(int)
            self.fixed_gradient_node_properties['anchor_node_IDs'] = fixed_gradient_linked_nodes.astype(int)
            self.fixed_gradient_node_properties['values_to_add'] = fixed_gradient_values_to_add
            self.fixed_gradient_link_properties['boundary_link_IDs'] = boundary_links.astype(int)
            #Update the link gradients over whole grid, as if there's values in the grid already, there could be compatibility issues...
            self.fixed_gradient_link_properties['boundary_link_gradients'] = -self.calculate_gradients_at_links(self['node'][gradient_of])[boundary_links]


###This method depreciated in favor of force_boundaries_from_gradients(), below.
###It works fine, but its handling of corners is very problematic and counterintuitive.
###DEJH, Feb 14.
    #def forced_values_from_gradients(self, floating_nodes, anchored_nodes,
    #        link_gradients, value='planet_surface__elevation',
    #        run_stability_test=True):
    #    """
    #    Returns the new value at the "floating" end of a link,
    #    where the anchored nodes at the other end of the link and the gradient
    #    on the link are specified. The value is set by the "value" flag, and
    #    defaults to the surface elevation, 'planet_surface__elevation'. If you 
    #    want to work with some other, or differently named, set of values, 
    #    change the flag.
    #    Gradients are, as always, positive down.
    #    Note that the true "sense" of the link between the nodes is preserved,
    #    and is not set by the order of anchor and floating nodes.
    #    
    #    If any of the node pairs cannot be connected by a single link, an
    #    exception is raised.
    #    
    #    If some nodes appear in both the from and to lists, instability may
    #    result! A test is applied to raise an error if this is the case, unless 
    #    suppressed by the input flag.
    #    
    #    >>> rmg = RasterModelGrid(3, 4, 1.0) # rows, columns, spacing
    #    >>> rmg.create_node_array_zeros('planet_surface__elevation')
    #    array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])
    #    >>> rmg['node']['planet_surface__elevation'] += 1.
    #    >>> rmg['node']['planet_surface__elevation'][[3,4,6]] = rmg.forced_values_from_gradients(numpy.array([3,4,6]),numpy.array([7,5,2]),numpy.array([0.5,-0.5,1.0]))
    #    >>> rmg['node']['planet_surface__elevation']
    #    array([ 1. ,  1. ,  1. ,  1.5,  0.5,  1. ,  0. ,  1. ,  1. ,  1. ,  1. ,
    #            1. ])
    #    
    #    ...Remember this printed array is inverted!, i.e., we recover the
    #    raster:
    #            1.0  1.0  1.0  1.0
    #            0.5  1.0  0.0  1.0
    #            1.0  1.0  1.0  1.5
    #    
    #    FOR THE FUTURE: This should transfer cleanly across to base.py, once 
    #    equivalent function to node_links is available there.
    #    """
    #    #determine the grid corners. This method has to be clever enough to realize when it's been given them!
    #    corner_nodes = self.corner_nodes
    #    if run_stability_test:
    #        if numpy.intersect1d(floating_nodes,anchored_nodes).size != 0:
    #            raise ValueError('Some nodes appear in both node lists! Aborting to prevent silent instability...')
    #    floating_links = self.node_links(floating_nodes)
    #    anchored_links = self.node_links(anchored_nodes)
    #    try:
    #        all_links = numpy.sort(numpy.concatenate((floating_links, anchored_links), axis=0), axis=0)
    #    except ValueError:
    #        raise ValueError('The two lists of nodes were not the same length!')
    #    shared_links_incl_null = numpy.where(numpy.diff(all_links, axis=0)==0,all_links[:-1,:],-1).flatten('F')
    #    shared_links = shared_links_incl_null[shared_links_incl_null>=0]
    #    try:
    #        anchored_node_is_tonode = numpy.equal(anchored_nodes, self.link_tonode[shared_links])
    #    except ValueError:
    #        raise ValueError('One or more of the node pairs were not connected by a link!')
    #    floating_node_is_tonode = numpy.logical_not(anchored_node_is_tonode)
    #    final_elevs = numpy.empty(floating_nodes.size, dtype='float')
    #    final_elevs[anchored_node_is_tonode] = self['node'][value][anchored_nodes[anchored_node_is_tonode]] + link_gradients[anchored_node_is_tonode]*self.dx
    #    final_elevs[floating_node_is_tonode] = self['node'][value][anchored_nodes[floating_node_is_tonode]] - link_gradients[floating_node_is_tonode]*self.dx
    #    
    #    return final_elevs

    def force_boundaries_from_gradients(self, link_IDs, link_gradients,
                value='planet_surface__elevation'):
        """
        Calculates and updates new values at the boundary nodes of a grid, when 
        provided with a list of fixed gradient link IDs, and the fixed values of 
        the gradients on these links.
        
        The "value" flag specifies which kind of data (e.g., elevation) the
        gradients refer to.
        
        This method follows the convention POSITIVE GRADIENT IS DOWN.
        
        The routine will automatically test to ensure the provided links are
        boundary links, and will raise an exception if they aren't. It is clever
        enough to distinguish for itself if any of the links provided are corner
        links (i.e., joining an edge node to a corner node). In such cases, the
        values of the corner nodes are updated *last*, such that the edge nodes
        they refer to have already been updated.
        
        Some examples:
                
        >>> rmg = RasterModelGrid(3, 4, 1.0) # rows, columns, spacing
        >>> rmg.create_node_array_zeros('planet_surface__elevation')
        array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])
        >>> rmg['node']['planet_surface__elevation'] += 2.
        >>> rmg.force_boundaries_from_gradients(numpy.array([  0,  1,  2,  5,  6,  7, 10, 11, 13, 14]),numpy.array([ 0.,-1.,-1., 1., 1., 0., 0.,-1., 1., 0.]))
        >>> rmg['node']['planet_surface__elevation']
        array([ 1.,  1.,  1.,  1.,  1.,  2.,  2.,  1.,  1.,  1.,  1.,  1.])
        >>> rmg.force_boundaries_from_gradients(numpy.array([ 11, 13]),numpy.array([ 2., 2.]))
        >>> rmg['node']['planet_surface__elevation']
        array([ 1.,  1.,  1.,  1.,  4.,  2.,  2.,  0.,  1.,  1.,  1.,  1.])
        
        ...and now we demonstrate an exception if an interior link is included:
        >>> rmg.force_boundaries_from_gradients(numpy.array([ 12, 13]),numpy.array([ 2., 2.]))
        Traceback (most recent call last):
            File "/Applications/Canopy.app/appdata/canopy-1.3.0.1715.macosx-x86_64/Canopy.app/Contents/lib/python2.7/doctest.py", line 1289, in __run
                compileflags, 1) in test.globs
            File "<doctest landlab.grid.raster.RasterModelGrid.force_boundaries_from_gradients[7]>", line 1, in <module>
                rmg.force_boundaries_from_gradients(numpy.array([ 12, 13]),numpy.array([ 2., 2.]))
            File "/Users/danhobley/Xcodesvn/PyLL/trunk/landlab/grid/raster.py", line 1381, in force_boundaries_from_gradients
                raise ValueError('One or more of the supplied links was neither an edge link, nor a link to a corner!')
        ValueError: One or more of the supplied links was neither an edge link, nor a link to a corner!
        
        """
        #determine the grid corners. This method has to be clever enough to realize when it's been given them!
        corner_nodes = self.corner_nodes
        tonodes = self.link_tonode[link_IDs]
        fromnodes = self.link_fromnode[link_IDs]
        tonode_boundaries = self.node_status[tonodes] != 0
        fromnode_boundaries = self.node_status[fromnodes] != 0
        edge_links = numpy.logical_xor(tonode_boundaries, fromnode_boundaries)
        edge_tonode_boundaries = numpy.logical_and(tonode_boundaries,edge_links)
        edge_fromnode_boundaries = numpy.logical_and(fromnode_boundaries,edge_links)
        self['node'][value][tonodes[edge_tonode_boundaries]] = self['node'][value][fromnodes[edge_tonode_boundaries]] - link_gradients[edge_tonode_boundaries]*self.dx
        self['node'][value][fromnodes[edge_fromnode_boundaries]] = self['node'][value][tonodes[edge_fromnode_boundaries]] + link_gradients[edge_fromnode_boundaries]*self.dx

        if not numpy.all(edge_links):
            corner_links = numpy.logical_not(edge_links)
            tonode_is_corner = numpy.in1d(tonodes[corner_links], corner_nodes)
            fromnode_is_corner = numpy.in1d(fromnodes[numpy.logical_not(edge_links)], corner_nodes)
            if not numpy.all(numpy.logical_xor(tonode_is_corner, fromnode_is_corner)):
                raise ValueError('One or more of the supplied links was neither an edge link, nor a link to a corner!')
            self['node'][value][(tonodes[corner_links])[tonode_is_corner]] = self['node'][value][(fromnodes[corner_links])[tonode_is_corner]] - (link_gradients[corner_links])[tonode_is_corner]*self.dx
            self['node'][value][(fromnodes[corner_links])[fromnode_is_corner]] = self['node'][value][(tonodes[corner_links])[fromnode_is_corner]] + (link_gradients[corner_links])[fromnode_is_corner]*self.dx

                
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
        lower_right = self.number_of_node_columns - 1
        upper_right = (self.number_of_node_columns +
                       self.number_of_node_rows - 2)
        upper_left = (2 * self.number_of_node_columns +
                      self.number_of_node_rows - 3)

        if bottom:
            for id in xrange(1, self.number_of_node_columns - 1):   # Bottom
                bc.boundary_code[id] = bc.TRACKS_CELL_BOUNDARY
                bc.tracks_cell[id] = id + self.number_of_node_columns
            bc.boundary_code[lower_left] = bc.TRACKS_CELL_BOUNDARY
            bc.tracks_cell[lower_left] = 1
            bc.boundary_code[lower_right] = bc.TRACKS_CELL_BOUNDARY
            bc.tracks_cell[lower_right] = lower_right-1
        if right:
            nbr = 2 * self.number_of_node_columns - 2
            for id in xrange( lower_right+1, upper_right ):   # Right
                bc.boundary_code[id] = bc.TRACKS_CELL_BOUNDARY
                bc.tracks_cell[id] = nbr
                nbr += self.number_of_node_columns
        if top:
            ncells = self.number_of_node_rows * self.number_of_node_columns
            nbr = ncells - (self.number_of_node_columns + 2)
            for id in xrange( upper_right+1, upper_left ):   # Top
                bc.boundary_code[id] = bc.TRACKS_CELL_BOUNDARY
                bc.tracks_cell[id] = nbr
                nbr = nbr - 1
            bc.boundary_code[upper_right] = bc.TRACKS_CELL_BOUNDARY
            bc.tracks_cell[upper_right] = ncells - 2
            bc.boundary_code[upper_left] = bc.TRACKS_CELL_BOUNDARY
            bc.tracks_cell[upper_left] = ncells + 1 - self.number_of_node_columns
        if left:
            n_boundary_cells = (2 * (self.num_rows - 2) +
                                2 * (self.num_cols - 2) + 4)
            nbr = (self.number_of_node_rows - 2) * self.number_of_node_columns + 1
            for id in xrange( upper_left+1, n_boundary_cells ):   # Left
                bc.boundary_code[id] = bc.TRACKS_CELL_BOUNDARY
                bc.tracks_cell[id] = nbr
                nbr = nbr - self.number_of_node_columns
        
        if self.DEBUG_VERBOSE:
            print 'tracks_cell:',bc.tracks_cell
    
    def calculate_gradients_at_links(self, node_values, out=None):
        diffs = gfuncs.calculate_diff_at_links(self, node_values, out=out)
        return numpy.divide(diffs, self._dx, out=diffs)

    @track_this_method
    def calculate_gradients_at_active_links(self, node_values, out=None):
        """
        Calculates the gradient in quantity s at each active link in the grid.
        This is nearly identical to the method of the same name in ModelGrid,
        except that it uses self._dx for link length to improve efficiency.
        
        Note that a negative gradient corresponds to a lower node in the direction
        of the link.
        
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
            
            >>> grad = numpy.zeros(rmg.number_of_active_links)
            >>> u = u*10
            >>> grad = rmg.calculate_gradients_at_active_links(u, grad)
            >>> grad
            array([ 10.,  10., -10., -10., -10., -10., -10.,   0.,  10.,  10.,  10.,
                   -10.,  10.,  10.,  10., -10.,  10.])
        """
        diffs = gfuncs.calculate_diff_at_active_links(self, node_values,
                                                      out=out)
        return numpy.divide(diffs, self._dx, out=diffs)


    def calculate_gradients_at_d8_active_links(self, node_values, out=None):
        """
        
        """
        
        diag_dist = 1.4142*self._dx
        straight_link_slopes = (node_values[self.activelink_tonode] - 
                                node_values[self.activelink_fromnode]) / \
                                self._dx
        diagonal_link_slopes = (node_values[self._diag_activelink_tonode] - 
                                node_values[self._diag_activelink_fromnode]) / \
                                diag_dist
        return numpy.concatenate((straight_link_slopes, diagonal_link_slopes))
        
        
    def calculate_steepest_descent_on_nodes(self, elevs_in, link_gradients, max_slope=False, dstr_node_ids=False):
        """
        Created DEJH Sept 2013. Based on approach of calc_flux_divergence..., below.
        Takes the elevations across a raster and the active link_gradients between those nodes, and returns the
        magnitude of the most downward slope from each node, and the direction of that cell. In the case where a node
        is a local minimum, method returns the lowest upward slope *as a negative slope*, and returns the downslope 
        node id as -1. i.e., downhill slopes are returned as positive values.
        
        At the moment, handling when equal gradients are present is poor. Method will currently preferentially route 
        flow according the priority scheme [N, E, S, W, NE, NW, SW, SE] for the equal height nodes in these cases.
        """
        
        if self.DEBUG_TRACK_METHODS:
            print 'RasterModelGrid.calculate_steepest_descent_on_nodes'
            
        assert (len(link_gradients)==self.number_of_active_links), \
               "incorrect length of active_link_gradients array"

        # If needed, create max_gradient array and the ID array
        if max_slope is False:
            max_slope = numpy.zeros(self.number_of_nodes)
        else:
            max_slope[:] = 0.
        
        if dstr_node_ids is False:
            dstr_node_ids = - numpy.ones(self.number_of_nodes, dtype=int)
        else:
            dstr_node_ids[:] = -1

        assert(len(max_slope) == self.number_of_nodes)
        assert(len(dstr_node_ids) == self.number_of_nodes)

        gradients = numpy.zeros(len(link_gradients)+1)
        gradients[:-1] = link_gradients
        
        #Make a matrix of the links. Need to append to this the gradients *on the diagonals*.
        node_links = numpy.vstack((gradients[self.node_active_outlink_matrix[0][:]], gradients[self.node_active_outlink_matrix[1][:]], -gradients[self.node_active_inlink_matrix[0][:]], -gradients[self.node_active_inlink_matrix[1][:]]))

        #calc the gradients on the diagonals:
        diagonal_nodes = (sgrid.diagonal_node_array(self.shape, out_of_bounds=-1)).T
        #Set the diagonals pointing to inactive nodes as inactive
        diagonal_nodes[numpy.where(self.node_status[diagonal_nodes] == 4)] = -1
        #Repeat the -1 indexing trick from above:
        elevs = numpy.zeros(len(elevs_in)+1)
        elevs[-1] = 9999999999999. #...as we want the gradients to inhibit flow
        elevs[:-1] = elevs_in
        
        slopes_diagonal_nodes = ((elevs[diagonal_nodes])-numpy.tile(elevs_in, (4,1)))/(1.41421356*self.dx)
        #Debug:
        #print 'Shape of node_links: ', node_links.shape
        #print 'Shape of diagonal array: ', slopes_diagonal_nodes.shape
        gradients_all_nodes = numpy.vstack((node_links, slopes_diagonal_nodes))
        #The ordering of this array is now [N, E, S, W, NE, NW, SW, SE][:].
        max_slope_indices = numpy.argmin(gradients_all_nodes, axis=0)
        max_slope = -gradients_all_nodes[max_slope_indices, xrange(gradients_all_nodes.shape[1])]
        #...per fancy indexing
        #print gradients_all_nodes
        #print max_slope_indices
        #print max_slope.reshape((5,5))
        #Assemble a node index array which corresponds to this gradients array from which to draw the dstr IDs:
        neighbors_ENWS = (self.get_neighbor_list()).T
        #print neighbors_ENWS.shape
        #print diagonal_nodes.shape
        dstr_id_source_array = numpy.vstack((neighbors_ENWS[1][:], neighbors_ENWS[0][:], neighbors_ENWS[3][:], neighbors_ENWS[2][:], diagonal_nodes))
        #print dstr_id_source_array.shape
        most_negative_gradient_node_ids = dstr_id_source_array[max_slope_indices, xrange(dstr_id_source_array.shape[1])]
        #print numpy.array([range(24),])
        #print dstr_id_source_array[:,:-1]
        #print most_negative_gradient_node_ids.reshape((5,5))
        #But we only want to return an id if the "dstr" node is actually downstream! So ->
        downslope_nodes = numpy.where(max_slope > 0)
        #print downslope_nodes
        dstr_node_ids[downslope_nodes] = most_negative_gradient_node_ids[downslope_nodes]
        #Local topo lows will retain the -1 index in dstr_node_ids
        #print 'dstr node ids: '
        #print dstr_node_ids.reshape((5,5))
        
        return max_slope, dstr_node_ids 
        

    @track_this_method
    def calculate_flux_divergence_at_nodes(self, active_link_flux, out=None):
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
            
            >>> df = rmg.zeros(centering='node') # outside loop
            >>> rmg.number_of_nodes
            20
            
        Then do this inside the loop:
            
            >>> df = rmg.calculate_flux_divergence_at_nodes(flux, df)
            
        In this case, the function will not have to create the df array.
        """
        return rfuncs.calculate_flux_divergence_at_nodes(
            self, active_link_flux, out=out)
        
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

#    def update_boundaries( self, u, bc = None ):
#        """
#        Updates FIXED_GRADIENT and TRACKS_CELL boundaries in 
#        BoundaryCondition "bc" (by default ModelGrid's default_bc), so
#        that TRACKS_CELL values in u are assigned the value at the 
#        corresponding cell-to-track, and FIXED_GRADIENT cells get a value
#        equal to the cell-to-track's value plus gradient times distance.
#
#        .. note:: Does NOT change fixed-value boundary cells.
#
#        .. todo::
#            use the cell-local distance rather than dx, for use in the base
#            class!
#            probably now obsolete (GT Aug 2013)
#            -->Agree DEJH, Jan 2014, and commented out
#            
#        NG Wondering if this should be boundary nodes, not cells.    
#        """
#    
#        if bc == None:
#            bc = self.default_bc
#
#        inds = (bc.boundary_code == bc.TRACKS_CELL_BOUNDARY)
#        u[self.boundary_cells[inds]] = u[bc.tracks_cell[inds]]
#
#        inds = (bc.boundary_code == bc.FIXED_GRADIENT_BOUNDARY)
#        u[self.boundary_cells[inds]] = (u[bc.tracks_cell[id]] +
#                                        bc.gradient[id] * self._dx)
#
#        return u

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
        >>> u = rmg.zeros(centering='node')
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
            >>> u = rmg.zeros(centering='cell')
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
        bottom]. Boundary nodes receive their actual neighbors (see example
        below); references to positions which are off the grid from boundary
        nodes receive -1.
        
        >>> mg = RasterModelGrid(4, 5)
        >>> mg.get_neighbor_list([-1, 6, 2])
        array([[-1, -1, 18, 14],
               [ 7, 11,  5,  1],
               [ 3,  7,  1, -1]])
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
        2D array. Formerly, only interior nodes were assigned neighbors;
        boundary nodes got -1 for each neighbor. As of 01/30/14, DEJH has
        modified the method so that neighboring nodes at the other end of active
        links are returned, but otherwise -1. The order of the neighbors is
        [right, top, left, bottom].

        .. note:: This is equivalent to the neighbors of all cells, getting the
            links between them, and setting neighbors at the other end of
            inactive links to -1. 

        .. todo: remove the loop by creating functions to return inactive links
            in sgrid.

        .. todo: Change to use BAD_INDEX_VALUE instead of -1.
        """
        # DH created this.  NG only changed labels.
        assert(self.neighbor_list_created == False)

        self.neighbor_list_created = True 
        #Former method which assigned -1 to ALL neighbors of boundary nodes:
        #self.neighbor_nodes = sgrid.neighbor_node_array(
        #    self.shape, out_of_bounds=-1, boundary_node_mask=-1)
        
        #DEJH new method 01/30/14. Returns neighbors of boundary nodes if there
        #are active links between them:
        self.neighbor_nodes = sgrid.neighbor_node_array(
            self.shape, out_of_bounds=-1)
        #get the boundary nodes
        boundary_nodes = sgrid.boundary_nodes(self.shape)
        #get the link and active link arrays for these nodes, in order SWNE
        links_from_each_BN = self.node_links(boundary_nodes)
        active_links_from_each_BN = self.active_node_links(boundary_nodes)
        #make a mask of links which AREN'T active (not v efficient)
        #...but OK as this method is seldom called
        mask = numpy.empty_like(links_from_each_BN)
        for i in xrange(links_from_each_BN.shape[0]):
            mask[i,:] = numpy.in1d(links_from_each_BN[i,:],
                active_links_from_each_BN[i,:], invert=True)
        #make sure to align the element orders, and blank the nodes at the end
        #of inactive links
        self.neighbor_nodes[numpy.fliplr(mask)] = -1
                
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
        self._initialize(num_rows_for_unit_test, num_cols_for_unit_test, 1.0)
        print 'done.'
        print
        
        print 'Testing fluxes and in/out links:'
        flux2 = numpy.zeros(len(flux)+1)
        flux2[:len(flux)] = flux
        print flux2
        for n in range(0, self.number_of_nodes):
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
        
    def calculate_aspect_at_nodes_bestFitPlane(self, id, val):
        """
        .. codeauthor:: Katy Barnhart <katherine.barnhart@colorado.edu>

        Calculates the aspect at each node based on the elevation of 
        the node and its neighbors using a best fit plane calculated
        using single value decomposition. 

        requires:
        id: id of nodes at which to calculate the aspect
        val: elevation at all nodes
        
        returns:
        a: the aspect at the nodes given by id
        """
        # additional note, KRB has written three codes in raster.py
        # one to calculate slope, one to calculate aspect, and one
        # to calculate both
        
        # get the list of neighboring nodes for the nodes given by id
        n=self.get_neighbor_list(id)
        a=[]
        
        # for each node in id make a list with the node id and the ids of
        # its neighbors. 
        
        # determine the values for the x, y, and z coordinates of each node, 
        # pass these to rfuncs.calculate_slope_aspect_BFP to calculate the
        # slope and aspect. 
        
        for i in range(len(id)):
            ns=list(n[i])
            ns.append(id[i])
            x=self.node_x[ns]
            y=self.node_y[ns]
            z=val[ns]
            slope, aspect = rfuncs.calculate_slope_aspect_BFP(x, y, z)
            a.append(aspect)
            del ns        
        # return aspect alone    
        return a
        
    def calculate_slope_at_nodes_bestFitPlane(self, id, val):
        """
        .. codeauthor:: Katy Barnhart <katherine.barnhart@colorado.edu>

        Calculates the aspect at each node based on the elevation of 
        the node and its neighbors using a best fit plane calculated
        using single value decomposition. 

        requires:
        id: id of nodes at which to calculate the aspect
        val: elevation at all nodes

        returns:
        a: the aspect at the nodes given by id
        """
        #
        # additional note, KRB has written three codes in raster.py
        # one to calculate slope, one to calculate aspect, and one
        # to calculate both
        
        # get the list of neighboring nodes for the nodes given by id
        n=self.get_neighbor_list(id)
        s=[]
        
        # for each node in id make a list with the node id and the ids of
        # its neighbors. 
        
        # determine the values for the x, y, and z coordinates of each node, 
        # pass these to rfuncs.calculate_slope_aspect_BFP to calculate the
        # slope and aspect. 
        
        for i in range(len(id)):
            ns=list(n[i])
            ns.append(id[i])
            x=self.node_x[ns]
            y=self.node_y[ns]
            z=val[ns]
            slope, aspect = rfuncs.calculate_slope_aspect_BFP(x, y, z)
            s.append(slope)
            del ns    
        # return slope alone    
        return s
        
    def calculate_slope_aspect_at_nodes_bestFitPlane(self, id, val):
        """
        .. codeauthor:: Katy Barnhart <katherine.barnhart@colorado.edu>

        Calculates the aspect at each node based on the elevation of 
        the node and its neighbors using a best fit plane calculated
        using single value decomposition. 

        requires:
        id: id of nodes at which to calculate the aspect
        val: elevation at all nodes

        returns:
        a: the aspect at the nodes given by id
        """

        # additional note, KRB has written three codes in raster.py
        # one to calculate slope, one to calculate aspect, and one
        # to calculate both
        
        # get the list of neighboring nodes for the nodes given by id
        n=self.get_neighbor_list(id)
        a=[]
        s=[]
        
        # for each node in id make a list with the node id and the ids of
        # its neighbors. 
        
        # determine the values for the x, y, and z coordinates of each node, 
        # pass these to rfuncs.calculate_slope_aspect_BFP to calculate the
        # slope and aspect. 
    
        for i in range(len(id)):
            ns=list(n[i])
            ns.append(id[i])
            x=self.node_x[ns]
            y=self.node_y[ns]
            z=val[ns]
            slope, aspect = rfuncs.calculate_slope_aspect_BFP(x, y, z)
            a.append(aspect)
            s.append(slope)
            
            del ns
        # return slope and aspect            
        return s, a
        
    def hillshade(self, alt, az, slp, asp):
        """
        .. codeauthor:: Katy Barnhart <katherine.barnhart@colorado.edu>
        
        code taken from GeospatialPython.com example from December 14th, 2014
        
        az = sun azimuth (degrees from north)
        alt = sun altitude (degrees up from horizon)
        slp = slope of cells at surface (degrees)
        asp = aspect of cells at surface (degrees from north)
        """
        
        # krb note: I don't know where this code is best put, probably not in
        # raster.py, but I'll set someone else (Dan, Greg, Eric) who knows the
        # whole landlab structure better put it where you think it should go.
        # it doesn't need self as an input, though if at some point, slope and
        # aspect were properties of the nodes, then they wouldn't need to be
        # passed around

        (alt, az, slp, asp) = (np.radians(alt), np.radians(az),
                               np.radians(slp), np.radians(asp))
        
        shaded = (
            np.sin(alt) * np.sin(slp) +
            np.cos(alt) * np.cos(slp) * np.cos(az - asp)
        )
    
        return shaded
        
        
def _is_closed_boundary(boundary_string):
    
    return boundary_string.lower() == 'closed'


def from_dict(param_dict):
    """
    Create a RasterModelGrid from the dictionary-like object, *param_dict*.
    Required keys of the dictionary are NUM_ROWS, NUM_COLS. Raises a KeyError
    if either of these are missing is given, use it as the
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


