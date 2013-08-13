#! /usr/env/python
"""
Python implementation of ModelGrid, a class used to
create and manage grids for 2D numerical models.

First version GT, July 2010
Last modified July 2013
"""

import numpy
import warnings

from landlab import model_parameter_dictionary as mpd
import landlab.utils.structured_grid as sgrid
from landlab.utils import count_repeated_values

BAD_INDEX_VALUE = numpy.iinfo(numpy.int).max

# Define the boundary-type codes
INTERIOR_NODE = 0
FIXED_VALUE_BOUNDARY = 1
FIXED_GRADIENT_BOUNDARY = 2
TRACKS_CELL_BOUNDARY = 3
INACTIVE_BOUNDARY = 4

_SLOW = False

if _SLOW:
    BAD_INDEX_VALUE = None


def _is_closed_boundary(boundary_string):
    
    return boundary_string.lower()=='closed'

    

def create_and_initialize_grid(input_source):
    """
    Creates, initializes, and returns a new grid object using parametes 
    specified either in a ModelParameterDictionary (param_dict) or a named input
    file (input_file_name).
    
    Example:
        
    >>> from StringIO import StringIO
    >>> test_file = StringIO(\"\"\"
    ... GRID_TYPE:
    ... raster
    ... NUM_ROWS:
    ... 4
    ... NUM_COLS:
    ... 5
    ... GRID_SPACING: 
    ... 2.5
    ... \"\"\")
    >>> from landlab import model_grid
    >>> mg = model_grid.create_and_initialize_grid(test_file)
    >>> mg.num_nodes
    20
        
    """    
    # Handle input source. 
    #In this code block, we do the following:
    #   - handle the case in which caller provides neither a parameter
    #     dictionary nor an input file name.
    #   - if we're given an input file name, create a parameter dictionary
    #     object that reads the specified file name
    if type(input_source) is mpd.ModelParameterDictionary:
        param_dict = input_source
    else:
        param_dict = mpd.ModelParameterDictionary(from_file=input_source)
        
    # Find out what type of grid the user wants
    #
    # Dev note: could handle defaults like: param_dict.get('GRID_TYPE','raster')
    # so if no GRID_TYPE is specified you get the second arg as default. If no
    # second arg, then exception
    grid_type = param_dict.read_string('GRID_TYPE')
    grid_type.strip().lower()   # make LC w/o leading/trailing spaces
    
    # Read parameters appropriate to that type, create it, and initialize it
    if grid_type=='raster':
        
        # Read and create basic raster grid
        nrows = param_dict.read_int('NUM_ROWS')
        ncols = param_dict.read_int('NUM_COLS')
        dx = param_dict.read_float('GRID_SPACING')
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

    else:
        print 'Non-raster grids not yet available in ModelGrid'
        mg = None
    
    # Return the created and initialized grid
    return mg


class BoundaryCondition(object):
    """
    The BoundaryCondition class stores boundary condition information for
    a particular variable. It includes:

        * A list of boundary codes

          - 1 = fixed value
          - 2 = fixed gradient
          - 3 = tracks cell
          - 4 = no-flux / inactive
          
        * A list of boundary gradients
        * A list of cell IDs to track
        
        NOTE: we may not need this anymore ... GT 6/13
    """
    def __init__( self, n_boundary_cells = 0 ):
    
        # Create the 3 vectors
        self.boundary_code = numpy.zeros(n_boundary_cells, dtype=numpy.short)
        self.boundary_gradient = numpy.zeros(n_boundary_cells, dtype=numpy.double)
        self.tracks_cell = numpy.zeros(n_boundary_cells, dtype=numpy.long)
        
        # Define the boundary-type codes
        self.INTERIOR_NODE = 0
        self.FIXED_VALUE_BOUNDARY = 1
        self.FIXED_GRADIENT_BOUNDARY = 2
        self.TRACKS_CELL_BOUNDARY = 3
        self.INACTIVE_BOUNDARY = 4

        # Set defaults
        self.boundary_code[:] = self.FIXED_VALUE_BOUNDARY
        self.boundary_gradient[:] = 0.0
        self.tracks_cell[:] = -1
        #DEJH, 7/17/13 - I don't think we should be using -1. This way, accidentally using one of these cells will silently just reference the final cell in the dvector, rather than throwing an exception. " = numpy.nan" would probably work better.


class ModelGrid(object):
    """
    Base class for creating and manipulating 2D structured or
    unstructured grids for numerical models.
    
    The idea is to have at least two inherited
    classes, RasterModelGrid and DelaunayModelGrid, that can create and
    manage grids. To this might be added a GenericModelGrid, which would
    be an unstructured polygonal grid that doesn't necessarily obey or
    understand the Delaunay triangulation, but rather simply accepts
    an input grid from the user. Also a HexModelGrid for hexagonal.
    """

    # Define the boundary-type codes
    INTERIOR_NODE = 0
    FIXED_VALUE_BOUNDARY = 1
    FIXED_GRADIENT_BOUNDARY = 2
    NODE_TRACKING_BOUNDARY = 3
    INACTIVE_BOUNDARY = 4
    
    # Debugging flags (if True, activates some output statements)
    DEBUG_VERBOSE = False
    DEBUG_TRACK_METHODS = False

    #-------------------------------------------------------------------
    def __init__( self ):
        
        pass
    
    #-------------------------------------------------------------------
    def initialize( self ):
    
        pass
        
    def get_node_status(self):
        """
        Returns an array of node boundary-status codes.
        """
        return self.node_status

    def create_cell_dvector( self ):
        """
        Returns a vector of floating point numbers the same length as 
        the number of cells.
        """
    
        return numpy.zeros( self.num_cells )

    def create_active_cell_dvector( self ):
        """
        Returns a vector of floating point numbers the same length as 
        the number of active cells.
        """
    
        return numpy.zeros( self.num_active_cells )

    def create_node_dvector( self ):
        """
        Returns a vector of floating point numbers the same length as 
        the number of nodes.
        """
    
        return numpy.zeros( self.num_nodes )
        
    def create_link_dvector(self):
        """
        Returns a 1D numpy array of floating point numbers the same length as 
        the number of links.
        """
        
        return numpy.zeros(self.num_links)

    def create_active_link_dvector(self):
        """
        Returns a 1D numpy array of floating point numbers the same length as 
        the number of active links.
        """
        
        return numpy.zeros(self.num_active_links)

    def create_face_dvector( self ):
        """
        Returns a vector of floating point numbers the same length as 
        the number of interior faces.
        """
    
        return numpy.zeros( self.num_faces )
        
    def set_fixed_value_boundaries(self, node_ids):
        """
        Assignes FIXED_VALUE_BOUNDARY status to specified nodes.
        """
        self.node_status[node_ids] = self.FIXED_VALUE_BOUNDARY
        if _SLOW:
            self.reset_list_of_active_links_slow()
        else:
            self.reset_list_of_active_links()

    def calculate_gradients_at_active_links(self, s, gradient=None):
        """
        Calculates the gradient in quantity s at each active link in the grid.
        """
        if self.DEBUG_TRACK_METHODS:
            print 'ModelGrid.calculate_gradients_at_active_links'
        
        if gradient==None:
            gradient = numpy.zeros(self.num_active_links)
            
        assert (len(gradient)==self.num_active_links), \
                "len(gradient)!=num_active_links"
                
        active_link_id = 0
        for link_id in self.active_links:
            gradient[active_link_id] = (s[self.link_tonode[link_id]]
                                        -s[self.link_fromnode[link_id]]) / \
                                        self.link_length[link_id]
            active_link_id += 1
        
        return gradient
        
    def calculate_flux_divergence_at_active_cells(self, active_link_flux, 
                                                  net_unit_flux=False):
        """
        Given an array of fluxes along links, computes the net total flux
        within each cell, divides by cell area, and stores the result in
        net_unit_flux. 
        
        The input active_link_flux should be flux of
        something (e.g., mass, momentum, energy) per unit face width, positive
        if flowing in the same direction as its link, and negative otherwise.
        There should be one value per active link. Returns an array of net
        total flux per unit area, one value per active cell (creates this
        array if it is not given as an argument).
          By convention, divergence is positive for net outflow, and 
        negative for net outflow. That's why we *add* outgoing flux and
        *subtract* incoming flux. This makes net_unit_flux have the same sign
        and dimensions as a typical divergence term in a conservation equation.

        In general, for a polygonal cell with N sides of lengths
        Li and with surface area A, the net influx divided by cell
        area would be:
            .. math::
                {Q_{net} \over A} = {1 \over A} \sum{q_i L_i}

        .. note::
            The net flux is defined as positive outward, negative
            inward. In a diffusion problem, for example, one would use:
                .. math::
                    {du \over dt} = \\text{source} - \\text{fd}
            where fd is "flux divergence".
            
        .. todo::
            This needs to be re-implemented with the faster algorithm used in
            the RasterModelGrid version.
        """
        
        if self.DEBUG_TRACK_METHODS:
            print 'ModelGrid.calculate_flux_divergence_at_active_cells'
            
        assert (len(active_link_flux)==self.num_active_links), \
               "incorrect length of active_link_flux array"
            
        # If needed, create net_unit_flux array
        if net_unit_flux==False:
            net_unit_flux = numpy.zeros(self.num_active_cells)
        else:
            net_unit_flux[:] = 0.
            
        assert (len(net_unit_flux))==self.num_active_cells
        
        # For each active link, add up the flux out of the "from" cell and 
        # into the "to" cell.
        active_link_id = 0
        for link_id in self.active_links:
            from_cell = self.node_activecell[self.link_fromnode[link_id]]
            to_cell = self.node_activecell[self.link_tonode[link_id]]
            total_flux = active_link_flux[active_link_id] * \
                         self.face_width[self.link_face[link_id]]
            #print('Flux '+str(total_flux)+' from '+str(from_cell) \
            #      +' to '+str(to_cell)+' along link '+str(link_id))
            if from_cell != BAD_INDEX_VALUE:
                net_unit_flux[from_cell] += total_flux
                #print('cell '+str(from_cell)+' net='+str(net_unit_flux[from_cell]))
            if to_cell != BAD_INDEX_VALUE:
                net_unit_flux[to_cell] -= total_flux
                #print('cell '+str(to_cell)+' net='+str(net_unit_flux[to_cell]))
            active_link_id += 1
        
        # Divide by cell area
        net_unit_flux = net_unit_flux / self.active_cell_areas
        
        return net_unit_flux
        
    def calculate_flux_divergence_at_nodes(self, active_link_flux, 
                                           net_unit_flux=False):
        """
        Same as calculate_flux_divergence_at_active_cells, but works with and
        returns a list of net unit fluxes that corresponds to all nodes, rather
        than just active cells. 
        
        Note that we don't compute net unit fluxes at
        boundary nodes (which don't have active cells associated with them, and 
        often don't have cells of any kind, because they are on the perimeter), 
        but simply return zeros for these entries. The advantage is that the 
        caller can work with node-based arrays instead of active-cell-based 
        arrays.

        .. todo::
            This needs to be re-implemented with the faster algorithm used in
            the RasterModelGrid version.
        """
        
        if self.DEBUG_TRACK_METHODS:
            print 'ModelGrid.calculate_flux_divergence_at_nodes'
            
        assert (len(active_link_flux)==self.num_active_links), \
               "incorrect length of active_link_flux array"
            
        # If needed, create net_unit_flux array
        if net_unit_flux==False:
            net_unit_flux = numpy.zeros(self.num_nodes)
        else:
            net_unit_flux[:] = 0.
            
        assert (len(net_unit_flux))==self.num_nodes
        
        # For each active link, add up the flux out of the "from" cell and 
        # into the "to" cell.
        active_link_id = 0
        for link_id in self.active_links:
            from_node = self.link_fromnode[link_id]
            from_cell = self.node_activecell[from_node]
            to_node = self.link_tonode[link_id]
            to_cell = self.node_activecell[to_node]
            total_flux = active_link_flux[active_link_id] * \
                         self.face_width[self.link_face[link_id]]
            #print('Flux '+str(total_flux)+' from '+str(from_node) \
            #      +' to '+str(to_node)+' along link '+str(link_id))
            if from_cell != BAD_INDEX_VALUE:
                net_unit_flux[from_node] += total_flux / \
                                            self.active_cell_areas[from_cell]
                #print('node '+str(from_node)+' net='+str(net_unit_flux[from_node]))
            if to_cell != BAD_INDEX_VALUE:
                net_unit_flux[to_node] -= total_flux / \
                                          self.active_cell_areas[to_cell]
                #print('node '+str(to_node)+' net='+str(net_unit_flux[to_node]))
            active_link_id += 1
        
        return net_unit_flux
        
    def get_node_x( self, id ):
        """
        Returns the x coordinate of node "id".
        """
        return self._node_x[id]
        
    def get_node_y( self, id ):
        """
        Returns the y coordinate of node "id".
        """
        return self._node_y[id]

    #Decorator alternatives to above node getters added DEJH late Jul '13.
    #Should allow either .get_property(ID) or .property[ID]
    @property
    def node_x(self):
        return self._node_x
    
    @property
    def node_y(self):
        return self._node_y
        
    def get_cell_x_coords( self ):
        """
        Returns vector of node x coordinates (same as get_node_x_coords).
        """
        warnings.warn('Use get_node_x_coords instead', DeprecationWarning)
        return self._node_x

    def get_cell_y_coords( self ):
        """
        Returns vector of node y coordinates (same as get_node_y_coords).
        """
        warnings.warn('Use get_node_y_coords instead', DeprecationWarning)
        return self._node_y

    def get_node_x_coords( self ):
        """
        Returns vector of node x coordinates.
        """
        return self._node_x           

    def get_node_y_coords( self ):
        """
        Returns vector of node y coordinates.
        """
        return self._node_y           

    def get_node_coords(self, axis=0):
        assert(axis in (0, 1))

        if axis == 0:
            return self.get_node_y_coords()
        else:
            return self.get_node_x_coords()

    def get_coordinate_units(self, axis=0):
        assert(axis in (0, 1))

        if axis == 0:
            return 'degrees_north'
        else:
            return 'degrees_east'

    def get_coordinate_name(self, axis=0):
        assert(axis in (0, 1))

        if axis == 0:
            return 'latitude'
        else:
            return 'longitude'

    def get_active_cell_node_ids( self ):
        """
        Returns an integer vector of the node IDs of all active cells.
        """
        return self.activecell_node
        
    def get_active_link_connecting_node_pair(self, node1, node2):
        """
        Returns the ID number of the active link that connects the given pair of
        nodes, or None if not found.
        
        Example:
            
            >>> rmg=RasterModelGrid(4, 5)
            >>> rmg.get_active_link_connecting_node_pair(8, 3)
            2
        """
        active_link = None
        for alink in range(0, self.num_active_links):
            if (self.activelink_fromnode[alink]==node1 \
                and self.activelink_tonode[alink]==node2) \
               or (self.activelink_tonode[alink]==node1
                and self.activelink_fromnode[alink]==node2):
                active_link = alink
                break
        return active_link
        
    def get_minimum_active_link_length(self):
        """
        Returns the horizontal length of the shortest active link in the grid.
        """
        return amin(self.link_length[self.active_links])

    def assign_upslope_vals_to_faces( self, u, v=0 ):
        """
        Assigns to each face the values of u at whichever of its
        neighbors has a higher value of v. If v is omitted, uses u for
        both.
        """
        
        fv = numpy.zeros( self.nfaces )
        if len(v) < len(u):
            for i in xrange( 0, self.nfaces ):
                fv[i] = max( u[self.fromcell[i]], u[self.tocell[i]] )
        else:
            for i in xrange( 0, self.nfaces ):
                if v[self.fromcell[i]] > v[self.tocell[i]]:
                    fv[i] = u[self.fromcell[i]]
                else:
                    fv[i] = u[self.tocell[i]]
        return fv
        
    def reset_list_of_active_links(self):
        """
        Creates or resets a list of active links. We do this by sweeping
        through the given lists of from and to nodes, and checking the status
        of these as given in the node_status list. A link is active if both its
        nodes are active interior points, or if one is an active interior and
        the other is an active boundary.
        """
        if self.DEBUG_TRACK_METHODS:
            print 'ModelGrid.reset_list_of_active_links'
            
        fromnode_status = self.node_status[self.link_fromnode]
        tonode_status = self.node_status[self.link_tonode]

        active_links = (((fromnode_status == self.INTERIOR_NODE) & ~
                         (tonode_status == self.INACTIVE_BOUNDARY)) |
                        ((tonode_status == self.INTERIOR_NODE) & ~
                         (fromnode_status == self.INACTIVE_BOUNDARY)))

        (self.active_links, ) = numpy.where(active_links)

        self.num_active_links = len(self.active_links)
        self.activelink_fromnode = self.link_fromnode[self.active_links]
        self.activelink_tonode = self.link_tonode[self.active_links]
        
        # Set up active inlink and outlink matrices
        self.setup_active_inlink_and_outlink_matrices()

    def reset_list_of_active_links_slow(self):
        """
        Creates or resets a list of active links. We do this by sweeping
        through the given lists of from and to nodes, and checking the status
        of these as given in the node_status list. A link is active if both its
        nodes are active interior points, or if one is an active interior and
        the other is an active boundary.
        """
        if self.DEBUG_TRACK_METHODS:
            print 'ModelGrid.reset_list_of_active_links'
            
        # Create or reset empy list of active links (we'll convert this to
        # a numpy array below)
        self.active_links = []
        
        # Create or reset lists of active outlinks and active inlinks for each
        # node. These are lists of lists.
        #self.node_activeoutlinks = [[] for x in xrange(self.num_nodes)]
        #self.node_activeinlinks = [[] for x in xrange(self.num_nodes)]
        
        for link in range(0, len(self.link_fromnode)):
            fromnode_status = self.node_status[self.link_fromnode[link]]
            tonode_status = self.node_status[self.link_tonode[link]]
            if ((fromnode_status==self.INTERIOR_NODE and
                 not tonode_status==self.INACTIVE_BOUNDARY ) or
                (tonode_status==self.INTERIOR_NODE and
                 not fromnode_status==self.INACTIVE_BOUNDARY)):
                self.active_links.append(link)
        
        self.num_active_links = len(self.active_links)
        self.active_links = numpy.array(self.active_links)
        self.activelink_fromnode = numpy.array(self.link_fromnode[self.active_links])
        self.activelink_tonode = numpy.array(self.link_tonode[self.active_links])
        
        # Set up active inlink and outlink matrices
        self.setup_active_inlink_and_outlink_matrices_slow()
        
    def deactivate_nodata_nodes(self, node_data, nodata_value):
        """
        Sets self.node_status to INACTIVE_BOUNDARY for all nodes whose value of
        node_data is equal to the nodata_value.
        
        Example:
            
            >>> mg = RasterModelGrid(3, 4, 1.0)
            >>> mg.node_status
            array([1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1], dtype=int8)
            >>> h = numpy.array([-9999,-9999,-9999,-9999,-9999,-9999,12345.,0.,-9999,0.,0.,0.])
            >>> mg.deactivate_nodata_nodes(h, -9999)
            >>> mg.node_status
            array([4, 4, 4, 4, 4, 4, 0, 1, 4, 1, 1, 1], dtype=int8)
        """
        
        # Find locations where value equals the NODATA code and set these nodes
        # as inactive boundaries.
        nodata_locations = numpy.nonzero(node_data==nodata_value)
        self.node_status[nodata_locations] = self.INACTIVE_BOUNDARY
        
        # Recreate the list of active cell IDs
        node_ids = numpy.array(range(0,self.num_nodes))
        self.activecell_node = node_ids[numpy.where(self.node_status == self.INTERIOR_NODE)]
        
        # Recreate the list of active links
        self.reset_list_of_active_links()
        
    def active_link_max(self, node_data):
        
        """
        For each active link, finds and returns the maximum value of node_data
        at either of the two ends. Use this, for example, if you want to find
        the maximum value of water depth at linked pairs of nodes (by passing
        in an array of water depth values at nodes).
        
        node_data: a 1D numpy array with length = number of nodes
        returns: a 1D numpy array of maximum values, with length = number of
            active links.
        
        Example:
            
            >>> mg = RasterModelGrid(3, 4, 1.0)
            >>> h = numpy.array([2.,2.,8.,0.,8.,0.,3.,0.,5.,6.,8.,3.])
            >>> mg.active_link_max(h)
            array([ 2.,  8.,  6.,  8.,  8.,  3.,  3.])
        """
        return numpy.maximum(node_data[self.activelink_fromnode],
                             node_data[self.activelink_tonode])
        

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
        describe connectivity information between cells and faces. Each
        cell in the grid has four faces. Each face has a "fromcell" and
        a "tocell"; the convention is that these always "point" up or
        right (so a negative flux across a face is either going left or
        down).
        
        Examples and doctests:
        
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
            >>> rmg.node_activecell[3] == BAD_INDEX_VALUE
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
        if 0:
            self.node_status = numpy.zeros( self.num_nodes, numpy.int8 )
            self.node_status[:] = self.INTERIOR_NODE
            bottom = range(0, num_cols)
            top = range(num_cols*(num_rows-1), self.num_nodes) 
            left = range(0, self.num_nodes, num_cols)
            right = range(num_cols-1, self.num_nodes, num_cols)
            self.node_status[bottom] = self.FIXED_VALUE_BOUNDARY
            self.node_status[top] = self.FIXED_VALUE_BOUNDARY
            self.node_status[left] = self.FIXED_VALUE_BOUNDARY
            self.node_status[right] = self.FIXED_VALUE_BOUNDARY
        else:
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
        if _SLOW:
            self.cell_node = []
            self.node_activecell = []
            node_id = 0
            cell_id = 0
            for r in range(0, num_rows):
                for c in range(0, num_cols):
                    if r!=0 and r!=(num_rows-1) and c!=0 and c!=(num_cols-1):
                        self.cell_node.append(node_id)
                        self.node_activecell.append(cell_id)
                        cell_id += 1
                    else:
                        self.node_activecell.append(None)
                    node_id += 1
            self.active_cells = list(range(0, len(self.cell_node)))
            self.activecell_node = self.cell_node   # default to all cells active
        else:
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
        if _SLOW:
            self.link_fromnode = []
            self.link_tonode = []
            
            #   vertical links
            for r in range(0, num_rows-1):
                for c in range(0, num_cols):
                    self.link_fromnode.append(c+r*num_cols)
                    self.link_tonode.append(c+(r+1)*num_cols)
            
            #   horizontal links
            for r in range(0, num_rows):
                for c in range(0, num_cols-1):
                    self.link_fromnode.append(c+r*num_cols)
                    self.link_tonode.append(c+r*num_cols+1)
            
            #   convert to numpy arrays
            self.link_fromnode = numpy.array(self.link_fromnode)
            self.link_tonode = numpy.array(self.link_tonode)
        else:
            (self.link_fromnode,
             self.link_tonode) = sgrid.node_link_index((num_rows, num_cols))

        #   set up in-link and out-link matrices and numbers
        if _SLOW:
            self.setup_inlink_and_outlink_matrices_slow()
        else:
            self.setup_inlink_and_outlink_matrices()
        
        #   set up the list of active links
        if _SLOW:
            self.reset_list_of_active_links_slow()
        else:
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
        if _SLOW:
            self.link_face = [None]*self.num_links  # make the list
            face_id = 0
            for link in self.active_links:
                self.link_face[link] = face_id
                face_id += 1
        else:
            self.link_face = sgrid.link_faces((num_rows, num_cols),
                                              actives=self.active_links)

        # List of neighbors for each cell: we will start off with no
        # list. If a caller requests it via get_neighbor_list or
        # create_neighbor_list, we'll create it if necessary.
        self.neighbor_list_created = False
        if self.DEBUG_VERBOSE:
          print 'Setting nlc flag'

        # List of diagonal neighbors. As with the neighbor list, we'll only
        # create it if requested.
        self.diagonal_list_created = False

        #--------OLDER STUFF BELOW----------
        if 0:
            # Keep track of pairs of cells that lie on either side of
            # each face. Cells are numbered 0=(0,0), 1=(0,1), 2=(0,2), etc.
            # Faces are numbered as follows: first vertical faces, going
            # left to right then bottom to top, then horizontal faces, again
            # left to right then bottom to top.
            #
            # Example, 3-row by 4-column grid, with 12 cells, 2 interior
            # cells, 10 boundary cells, and 7 faces (3 vertical and
            # 4 horizontal):
            #
            # |-------|-------|-------|-------|
            # |       |       |       |       |
            # |   8   |   9   |  10   |  11   |
            # |       |       |       |       |
            # |-------|---5---|---6---|-------|
            # |       |       |       |       |
            # |   4   0   5   1   6   2   7   |
            # |       |       |       |       |
            # |-------|---3---|---4---|-------|
            # |       |       |       |       |
            # |   0   |   1   |   2   |   3   |
            # |       |       |       |       |
            # |-------|-------|-------|-------|
            #
            # The from and to cells of face 0 are 4 and 5, respectively.
            # The faces of cell 5 are 0, 1, 3 and 5.
            #
            # Along the way, we store the x and y coordinates of the center
            # of each face.
            #
            self.fromcell = numpy.zeros(self.nfaces, dtype=numpy.int) #TBX
            self.tocell = numpy.zeros(self.nfaces, dtype=numpy.int) #TBX
            self.facex = numpy.zeros( self.nfaces ) #TBX
            self.facey = numpy.zeros( self.nfaces ) #TBX
            halfdx = self._dx / 2.0
            face_id = 0
            for r in xrange( 1, num_rows-1 ):
                for c in range( 1, num_cols ):
                    self.fromcell[face_id] = r * num_cols + ( c - 1 )
                    self.tocell[face_id] = self.fromcell[face_id] + 1
                    self.facex[face_id] = c*self._dx - halfdx
                    self.facey[face_id] = r*self._dx
                    face_id += 1
            for r in xrange( 1, num_rows ):
                for c in range( 1, num_cols-1 ):
                    self.fromcell[face_id] = ( r - 1 ) * num_cols + c
                    self.tocell[face_id] = self.fromcell[face_id] + num_cols
                    self.facex[face_id] = c*self._dx
                    self.facey[face_id] = r*self._dx - halfdx
                    face_id += 1
            if self.DEBUG_VERBOSE:
                print 'fromcell:',self.fromcell
                print 'tocell:',self.tocell
                print 'facex:',self.facex
                print 'facey:',self.facey

            # Now we find the face IDs connected to each cell.
            # The boundary cells don't have faces between them, so they will be flagged
            # with a -1. However, the bottom row of boundary cells have top faces,
            # the top row has bottom faces, etc.
            # Note that the four faces are numbered counter-clockwise
            # starting from the right face, that is, 0 is right, 1 is top,
            # 2 is left, and 3 is bottom.
            self.faces = - numpy.ones([self.ncells, 4], dtype=numpy.int)
            n_vert_faces = ( num_rows - 2 ) * ( num_cols - 1 )
            for r in xrange( 1, num_rows-1 ):   # Faces for interior cells
                for c in xrange( 1, num_cols-1 ):
                    cell_id = r * num_cols + c
                    self.faces[cell_id,2] = (r-1)*(num_cols-1)+(c-1)   # left
                    self.faces[cell_id,0] = self.faces[cell_id,2] + 1  # right
                    self.faces[cell_id,3] = n_vert_faces+(r-1)*(num_cols-2)+(c-1) # bottom
                    self.faces[cell_id,1] = self.faces[cell_id,3] + (num_cols-2)  # top
            for cell_id in xrange( 1, num_cols-1 ):  # Top faces for bottom row
                self.faces[cell_id,1] = n_vert_faces+(cell_id-1)
            r = num_rows - 1
            for c in xrange( 1, num_cols-1 ): # Bottom faces for top row
                cell_id = r * num_cols + c
                self.faces[cell_id,3] = n_vert_faces+(r-1)*(num_cols-2)+(c-1)
            c = 0
            for r in xrange( 1, num_rows-1 ): # Right faces for left column
                cell_id = r * num_cols + c
                self.faces[cell_id,0] = (r-1)*(num_cols-1)+c
            c = num_cols-1
            for r in xrange( 1, num_rows-1 ): # Left faces for right column
                cell_id = r * num_cols + c
                self.faces[cell_id,2] = (r-1)*(num_cols-1)+(c-1)
            

            if self.DEBUG_VERBOSE:
                for i in xrange( 1, self.ncells ):
                    print i,self.faces[i,:]
                    
            # List of neighbors for each cell: we will start off with no
            # list. If a caller requests it via get_neighbor_list or
            # create_neighbor_list, we'll create it if necessary.
            self.neighbor_list_created = False
            if self.DEBUG_VERBOSE:
              print 'Setting nlc flag'

            # List of diagonal neighbors. As with the neighbor list, we'll only
            # create it if requested.
            self.diagonal_list_created = False

            # For each node, we also need to know which way the face
            # points. A flux of mass across a face is considered positive
            # when it flows toward the "to" cell, and negative otherwise.
            # The "face_sign" matrix records whether the face points
            # toward the cell (+1) or away from it (-1). Note that in this
            # raster grid, the left and bottom faces (2 and 3) point into the cell,
            # while the right and top faces (0 and 1) point away from the cell.
            #    This means that you can calculate the "cell gradient", the
            # gradient of a field from the perspective of a cell, using:
            #        self.face_sign[id,:] * g[self.faces[id,:]]
            # where g is a gradient of a field defined at faces. The
            # "cell gradient" at a face is positive when you have to "walk uphill"
            # into the cell (the cell is higher than its neighbor on the other
            # side of the face), and negative when you "walk downhill" to the
            # cell (the cell is lower than its neighbor).
            self.face_sign = numpy.zeros( [self.ncells, 4], dtype=numpy.short )
            self.face_sign[:,0:2] = -1
            self.face_sign[:,2:] = 1
            if self.DEBUG_VERBOSE:
                print 'face sign:',self.face_sign
            
            # Set up list of interior cells
            # (Note that this could be superceded if you wanted an irregular
            # boundary inside the rectangular grid)
            self.interior_cells = numpy.zeros(self.n_interior_cells, dtype=numpy.int)
            id = 0
            for r in xrange( 1, num_rows-1 ):
                for c in range( 1, num_cols-1 ):
                    self.interior_cells[id] = r * num_cols + c
                    id += 1
            
            if self.DEBUG_VERBOSE:
                print self.interior_cells
            
            #NG is pretty sure that the code below these comments are obsolete.
            #For now it is just commented out.  Maybe we can delete it?
            #
            # Boundary condition handling: 
            # 
            # To handle the boundaries properly, we need to do the following:
            #  1. Find out whether cell J is a boundary, and what type
            #  2. Identify which cells are fixed-value boundaries, so they
            #     can be updated as needed.
            #  3. For boundary cells that track the value of another cell,
            #     update all their values.
            #  4. For these cells, remember the cell IDs of the cells they
            #     are to track.
            #  5. For any single boundary cell J, update its value according
            #     to its boundary type.
            #
            # To do accomplish these, we use three data structures. The first
            # is a vector of the cell ID ("CID") of all
            # boundary cells (self.boundary_cells). The second is
            # a list of the boundary id ("BID") for all cells. The boundary
            # id is the index number (0 to self.n_boundary_cells-1) of the
            # corresponding cell in the self.boundary_cells vector. In this
            # respect, self.boundary_cells and self.boundary_ids point to each
            # other. Obviously, not all cells are boundaries, and so those 
            # cells that are interior cells are simply flagged with a -1.
            # This way, the self.boundary_ids list encodes two pieces of
            # information: (1) whether the cell is a boundary or not, and 
            # (2) where to look for more info if it is a boundary (i.e., 
            # what is its BID).
            #
            # Finally, for those cells (if any) that are no-flux or periodic,
            # we need to know the CID of the cell whose value they mirror.
            # Now, the user may want to have different kinds of boundary 
            # conditions apply to different variables. For example, for
            # modeling floods across an active fault, you might want the
            # land surface elevation boundary condition to be fixed value
            # at some places, while the water depth is fixed gradient.
            # For this reason, we might need multiple version of boundary
            # status and related information, depending on the user's needs.
            # To manage this, we use the BoundaryCondition class. We have
            # BoundaryCondition called self.default_bc that represents the
            # default handling; if the user wants, they can make more 
            # BoundaryCondition objects and treat them differently.
            #
            # Here we work counter-clockwise around the perimeter, starting
            # with CID/BID 0 at the lower left. BID's increase in counter-
            # clockwise order around the edge.
            #
            # (Note that this could be superceded if you wanted an irregular
            # boundary inside the rectangular grid)
            #
            #self.boundary_ids = - numpy.ones( self.ncells, dtype=numpy.int )
            #self.boundary_cells = numpy.zeros( self.n_boundary_cells, dtype=numpy.int )
            #id = 0
            #for r in xrange( 0, num_cols-1 ):       # Bottom
            #    self.boundary_cells[id] = r
            #    self.boundary_ids[r] = id
            #    id += 1
            #for c in xrange( num_cols-1, self.ncells-num_cols, num_cols ):  # Right
            #    self.boundary_cells[id] = c
            #    self.boundary_ids[c] = id
            #    id += 1
            #for r in xrange( self.ncells-1, num_cols*(num_rows-1), -1 ):       # Top
            #    self.boundary_cells[id] = r
            #    self.boundary_ids[r] = id
            #    id += 1
            #for c in xrange( num_cols*(num_rows-1), 0, -num_cols ):  # Left
            #    self.boundary_cells[id] = c
            #    self.boundary_ids[c] = id
            #    id += 1
            #
            #if self.DEBUG_VERBOSE:
            #    print 'Boundary CIDs:',self.boundary_cells
            #    print 'Cell BIDs:',self.boundary_ids
                
            self.default_bc = BoundaryCondition( self.n_boundary_cells )

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

        if 0:
            # Create active in-link and out-link matrices.
            self.node_inlink_matrix = - numpy.ones((2, self.num_nodes), dtype=numpy.int)
            self.node_outlink_matrix = - numpy.ones((2, self.num_nodes), dtype=numpy.int)

            # Set up the inlink arrays
            tonodes = self.link_tonode
            self.node_numinlink = numpy.bincount(tonodes,
                                                 minlength=self.num_nodes)

            counts = count_repeated_values(self.link_tonode)
            for (count, (tonodes, link_ids)) in enumerate(counts):
                self.node_inlink_matrix[count][tonodes] = link_ids

            # Set up the outlink arrays
            fromnodes = self.link_fromnode
            self.node_numoutlink = numpy.bincount(fromnodes,
                                                  minlength=self.num_nodes)
            counts = count_repeated_values(self.link_fromnode)
            for (count, (fromnodes, link_ids)) in enumerate(counts):
                self.node_outlink_matrix[count][fromnodes] = link_ids
        
    def setup_inlink_and_outlink_matrices_slow(self):
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
        # Create in-link and out-link matrices.
        self.node_numinlink = numpy.zeros(self.num_nodes,dtype=numpy.int)
        self.node_numoutlink = numpy.zeros(self.num_nodes,dtype=numpy.int)
        self.node_inlink_matrix = -numpy.ones((2,self.num_nodes),dtype=numpy.int)
        self.node_outlink_matrix = -numpy.ones((2,self.num_nodes),dtype=numpy.int)
        
        # For each link, assign it as an "inlink" of its "to" node and an 
        # "outlink" of its from node. Keep track of the total number of inlinks
        # and outlinks of each node.
        for link_id in range(0, self.num_links):
            
            # enter this link as an inlink of its "to" node, and increment the
            # total number of inlinks for this node
            tonode = self.link_tonode[link_id]
            inlinknum = self.node_numinlink[tonode]
            self.node_inlink_matrix[inlinknum][tonode] = link_id
            self.node_numinlink[tonode] += 1

            # enter this link as an outlink of its "from" node, and increment the
            # total number of outlinks for this node
            fromnode = self.link_fromnode[link_id]
            outlinknum = self.node_numoutlink[fromnode]
            self.node_outlink_matrix[outlinknum][fromnode] = link_id
            self.node_numoutlink[fromnode] += 1
        
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

    def setup_active_inlink_and_outlink_matrices_slow(self):
        """
        Creates data structures to record the numbers of active inlinks and 
        active outlinks for each node. These data structures are equivalent to
        the "regular" inlink and outlink matrices, except that it uses the IDs
        of active links (only).
        """
        
        # Create active in-link and out-link matrices.
        self.node_numactiveinlink = numpy.zeros(self.num_nodes,dtype=numpy.int)
        self.node_numactiveoutlink = numpy.zeros(self.num_nodes,dtype=numpy.int)
        self.node_active_inlink_matrix = -numpy.ones((2,self.num_nodes),dtype=numpy.int)
        self.node_active_outlink_matrix = -numpy.ones((2,self.num_nodes),dtype=numpy.int)
        
        # For each active link, assign it as an "inlink" of its "to" node and an 
        # "outlink" of its from node. Keep track of the total number of active
        # inlinks and outlinks of each node.
        for active_link_id in range(0, self.num_active_links):
            
            # enter this link as an inlink of its "to" node, and increment the
            # total number of inlinks for this node
            tonode = self.activelink_tonode[active_link_id]
            inlinknum = self.node_numactiveinlink[tonode]
            self.node_active_inlink_matrix[inlinknum][tonode] = active_link_id
            self.node_numactiveinlink[tonode] += 1

            # enter this link as an outlink of its "from" node, and increment the
            # total number of outlinks for this node
            fromnode = self.activelink_fromnode[active_link_id]
            outlinknum = self.node_numactiveoutlink[fromnode]
            self.node_active_outlink_matrix[outlinknum][fromnode] = active_link_id
            self.node_numactiveoutlink[fromnode] += 1
        
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
        
    # no such thing as interior cells anymore.  NG got rid of this.    
    #def get_count_of_interior_cells(self):
    #    """
    #    Returns the number of interior cells on the grid.  
    #    NG, June 2013
    #    """
    #    return(self.num_active_cells)
    
    def get_count_of_interior_nodes(self):
        """
        Returns the number of interior nodes on the grid. Functionally identical to get_count_of_interior_cells()
        DEJH, July 2013
        """
        return(self.num_active_cells)
        
    #this was wrong, delete if it doesn't blow-up code.
    #def get_count_of_all_cells(self):
    #    """
    #    Returns total number of nodes, including boundaries. Note this call is misleadingly named; it returns the number of nodes, not cells (i.e., it includes all boundary nodes).
    #    NG, June 2013
    #    """
    #    return(self.num_nodes)
    
    def get_count_of_all_nodes(self):
        """
        Returns total number of nodes, including boundaries. Functionally identical to get_count_of_all_cells().
        DEJH, July 2013
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
        return numpy.array([ID, ID+1, ID+self.ncols, ID+self.ncols+1])
    
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
            if self.node_status[a] != self.INACTIVE_BOUNDARY:
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
            if self.node_status[a] != self.INACTIVE_BOUNDARY:
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
        
        return all_neighbor_nodes[index_max]
        
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
            self.node_status[bottom_edge] = self.INACTIVE_BOUNDARY
        else:
            self.node_status[bottom_edge] = self.FIXED_VALUE_BOUNDARY

        if right_is_inactive:
            self.node_status[right_edge] = self.INACTIVE_BOUNDARY
        else:
            self.node_status[right_edge] = self.FIXED_VALUE_BOUNDARY
            
        if top_is_inactive:
            self.node_status[top_edge] = self.INACTIVE_BOUNDARY
        else:
            self.node_status[top_edge] = self.FIXED_VALUE_BOUNDARY

        if left_is_inactive:
            self.node_status[left_edge] = self.INACTIVE_BOUNDARY
        else:
            self.node_status[left_edge] = self.FIXED_VALUE_BOUNDARY
        
        if _SLOW:
            self.reset_list_of_active_links_slow()
        else:
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
        
    def calculate_gradients_at_active_links_slow(self, s, gradient=None):
        """
        Slow implementation of calculate_gradients_at_active_links, kept here
        as an example of slow vs fast code using numpy arrays.
        """
        if self.DEBUG_TRACK_METHODS:
            print 'RasterModelGrid.calculate_gradients_at_active_links_slow'
        
        if gradient==None:
            gradient = numpy.zeros(self.num_active_links)
            
        assert (len(gradient)==self.num_active_links), \
                "len(gradient)!=num_active_links"
                
        active_link_id = 0
        for link_id in self.active_links:
            gradient[active_link_id] = (s[self.link_tonode[link_id]]
                                        -s[self.link_fromnode[link_id]]) / \
                                        self._dx
            active_link_id += 1
        
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

    def calculate_flux_divergence_at_active_cells_slow(self, active_link_flux, 
                                                  net_unit_flux=False):
        """
        Slower version of calculate_flux_divergence_at_active_cells. Retained
        here to illustrate value of using numpy array operations instead of
        FOR loops.
        """
        
        if self.DEBUG_TRACK_METHODS:
            print 'RasterModelGrid.calculate_flux_divergence_at_active_cells_slow'
            
        assert (len(active_link_flux)==self.num_active_links), \
               "incorrect length of active_link_flux array"
            
        # If needed, create net_unit_flux array
        if net_unit_flux is False:
            net_unit_flux = numpy.zeros(self.num_active_cells)
        else:
            net_unit_flux[:] = 0.
            
        assert (len(net_unit_flux))==self.num_active_cells
        
        # For each active link, add up the flux out of the "from" cell and 
        # into the "to" cell.
        active_link_id = 0
        for link_id in self.active_links:
            from_cell = self.node_activecell[self.link_fromnode[link_id]]
            to_cell = self.node_activecell[self.link_tonode[link_id]]
            total_flux = active_link_flux[active_link_id] * self._dx
            #print('Flux '+str(total_flux)+' from '+str(from_cell) \
            #      +' to '+str(to_cell)+' along link '+str(link_id))
            if from_cell != BAD_INDEX_VALUE:
                net_unit_flux[from_cell] += total_flux
                #print('cell '+str(from_cell)+' net='+str(net_unit_flux[from_cell]))
            if to_cell != BAD_INDEX_VALUE:
                net_unit_flux[to_cell] -= total_flux
                #print('cell '+str(to_cell)+' net='+str(net_unit_flux[to_cell]))
            active_link_id += 1
        
        # Divide by cell area
        net_unit_flux = net_unit_flux / self.cellarea
        
        return net_unit_flux  
        
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
        
    def calculate_flux_divergence_at_nodes_slow(self, active_link_flux, 
                                           net_unit_flux=False):
        """
        Same as calculate_flux_divergence_at_active_cells, but works with and
        returns a list of net unit fluxes that corresponds to all nodes, rather
        than just active cells. 
        
        Note that we don't compute net unit fluxes at
        boundary nodes (which don't have active cells associated with them, and 
        often don't have cells of any kind, because they are on the perimeter), 
        but simply return zeros for these entries. The advantage is that the 
        caller can work with node-based arrays instead of active-cell-based 
        arrays.
        
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
            >>> df = rmg.calculate_flux_divergence_at_nodes_slow(flux)
            >>> df # doctest: +NORMALIZE_WHITESPACE
            array([ 0.,  0.,  0.,  0.,  0.,  0.,  2.,  4., -2.,  0.,  0.,  0.,  1.,
                   -4.,  0.,  0.,  0.,  0.,  0.,  0.])
        """
        
        if self.DEBUG_TRACK_METHODS:
            print 'RasterModelGrid.calculate_flux_divergence_at_nodes_slow'
            
        assert (len(active_link_flux)==self.num_active_links), \
               "incorrect length of active_link_flux array"
            
        # If needed, create net_unit_flux array
        if net_unit_flux is False:
            net_unit_flux = numpy.zeros(self.num_nodes)
        else:
            net_unit_flux[:] = 0.
            
        assert(len(net_unit_flux) == self.num_nodes)
        
        # For each active link, add up the flux out of the "from" cell and 
        # into the "to" cell.
        active_link_id = 0
        for link_id in self.active_links:
            from_node = self.link_fromnode[link_id]
            from_cell = self.node_activecell[from_node]
            to_node = self.link_tonode[link_id]
            to_cell = self.node_activecell[to_node]
            total_flux = active_link_flux[active_link_id] * self._dx
            #print('Flux '+str(total_flux)+' from '+str(from_node) \
            #      +' to '+str(to_node)+' along link '+str(link_id))
            if from_cell != BAD_INDEX_VALUE:
                net_unit_flux[from_node] += total_flux / self.cellarea
                #print('node '+str(from_node)+' net='+str(net_unit_flux[from_node]))
            if to_cell != BAD_INDEX_VALUE:
                net_unit_flux[to_node] -= total_flux / self.cellarea
                #print('node '+str(to_node)+' net='+str(net_unit_flux[to_node]))
            active_link_id += 1

        return net_unit_flux

    def calculate_flux_divergences( self, q ):
        """
        TBX: TO BE REPLACED WITH METHODS ABOVE
        
        Calculates the net flux at each cell by adding up the fluxes
        through all four faces, and divides the total by cell area.

        In general, for a polygonal cell with N sides of lengths
        Li and with surface area A, the net influx divided by cell
        area would be:
            .. math::
                {Q_{net} \over A} = {1 \over A} \sum{q_i L_i}

        For a square cell, the sum is over 4 sides of length dx, and
        :math:`A = dx^2`, so:
            .. math::
                {Q_{net} \over A} = {1 \over dx} \sum{q_i}

        .. note::
            The net flux is defined as positive outward, negative
            inward. In a diffusion problem, for example, one would use:
                .. math::
                    {du \over dt} = \\text{source} - \\text{fd}
            where fd is "flux divergence".
        """
    
        if self.DEBUG_TRACK_METHODS:
            print 'RasterModelGrid.calculate_flux_divergences here'
        
        fd = numpy.zeros( self.ncells )
        for cell in self.interior_cells:
            if self.DEBUG_VERBOSE:
                print 'Cell',cell
                print 'q:',q[self.faces[cell,0:4]]
            fd[cell] = ( -( q[self.faces[cell,2]]   # left face (positive=in)
                + q[self.faces[cell,3]] )           # bottom face (positive=in)
                + q[self.faces[cell,0]]             # right face (positive=out)
                + q[self.faces[cell,1]]             # top face (positive=out)
                  ) / self._dx
        return fd
        
    def calculate_flux_divergence( self, q, id ):
        """
        TODO: UPDATE THIS TO USE NEW DATA STRUCTURES!
        
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

        if _SLOW:
            for id in xrange( 0, self.n_boundary_cells ):
                if bc.boundary_code[id] == bc.TRACKS_CELL_BOUNDARY:
                    u[self.boundary_cells[id]] = u[bc.tracks_cell[id]]
        else:
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
            
        NG Wondering if this should be boundary nodes, not cells.    
        """
    
        if bc == None:
            bc = self.default_bc

        if _SLOW:
            for id in xrange( 0, self.n_boundary_cells ):
                if bc.boundary_code[id] == bc.TRACKS_CELL_BOUNDARY:
                    u[self.boundary_cells[id]] = u[bc.tracks_cell[id]]
                elif bc.boundary_code[id] == bc.FIXED_GRADIENT_BOUNDARY:
                    u[self.boundary_cells[id]] = u[bc.tracks_cell[id]] \
                                                 + bc.gradient[id]*self._dx
        else:
            inds = (bc.boundary_code == bc.TRACKS_CELL_BOUNDARY)
            u[self.boundary_cells[inds]] = u[bc.tracks_cell[inds]]

            inds = (bc.boundary_code == bc.FIXED_GRADIENT_BOUNDARY)
            u[self.boundary_cells[inds]] = (u[bc.tracks_cell[id]] +
                                            bc.gradient[id] * self._dx)

        return u

    def node_vector_to_raster(self, u, flip_vertically=False):
        """
        Converts node vector u to a 2D array and returns it, so that it
        can be plotted, output, etc.
        
        If the optional argument flip_vertically=True, the function returns an 
        array that has the rows in reverse order, for use in plot commands (such
        as the image display functions) that put the (0,0) axis at the top left 
        instead of the bottom left.
        
        NG : The comment above confuses me.  Based on the example below, looks
        (0,0) goes on the top left if it's False.
        
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
        
        assert(len(u)==self.num_nodes), ('u should have '+str(self.num_nodes) \
                                         +' elements')
    
        rast = numpy.zeros( [self.nrows, self.ncols] )
        if flip_vertically==False:
            rows = range(0, self.nrows)
        else:
            rows = range(self.nrows-1, -1, -1)
        id = 0
        for r in rows:
            rast[r,:] = u[id:(id+self.ncols)]
            id += self.ncols
        return rast

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
    
        assert(len(u)==self.num_cells), ('u should have '+str(self.num_cells) \
                                         +' elements')
    
        rast = numpy.zeros( [self.nrows-2, self.ncols-2] )
        if flip_vertically==False:
            rows = range(0, self.nrows-2)
        else:
            rows = range(self.nrows-3, -1, -1)
        id = 0
        for r in rows:
            rast[r,:] = u[id:(id+(self.ncols-2))]
            id += self.ncols-2
        return rast

    #OLD CODE, same code is below, but with cells changed to node because
    #NG is easily confused.
    #def get_neighbor_list( self, id = -1 ):
    #    """
    #    If id is specified, returns a list of neighboring cell IDs for
    #    the node "id". Otherwise, returns lists for all cells as a 2D
    #    array. The list is in the order [right, top, left, bottom].
    #    """
    #
    #    if self.neighbor_list_created==False:
    #        self.create_neighbor_list()
    #    
    #    if id > -1:
    #        return self.neighbor_cells[id,:]
    #    else:
    #        return self.neighbor_cells
    #        
    #def create_neighbor_list( self ):
    #    """
    #    Creates a list of IDs of neighbor cells for each cell, as a
    #    2D array. Only interior cells are assigned neighbors; boundary
    #    cells get -1 for each neighbor. The order of the neighbors is [right, top, left, bottom].
    #    
    #    """
    #
    #    assert self.neighbor_list_created == False
    #    
    #    self.neighbor_list_created = True       
    #    self.neighbor_cells = -ones( [self.ncells, 4], dtype=numpy.int )
    #    for r in xrange( 1, self.nrows-1 ):
    #        for c in xrange( 1, self.ncols-1 ):
    #            cell_id = r * self.ncols + c
    #            self.neighbor_cells[cell_id,2] = cell_id - 1   # left
    #            self.neighbor_cells[cell_id,0] = cell_id + 1  # right
    #            self.neighbor_cells[cell_id,3] = cell_id - self.ncols # bottom
    #            self.neighbor_cells[cell_id,1] = cell_id + self.ncols # top

    def get_neighbor_list( self, id = -1 ):
        """
        If id is specified, returns a list of neighboring node IDs. 
        Otherwise, returns lists for all nodes as a 2D
        array. The list is in the order [right, top, left, bottom].
        DH created this.  NG only changed labels.
        """
    
        if self.neighbor_list_created==False:
            self.create_neighbor_list()
        
        if id > -1:
            return self.neighbor_nodes[id,:]
        else:
            return self.neighbor_nodes
            
    def create_neighbor_list( self ):
        """
        Creates a list of IDs of neighbor nodes for each node, as a
        2D array. Only interior nodes are assigned neighbors; boundary
        nodes get -1 for each neighbor. 
        The order of the neighbors is [right, top, left, bottom].
        DH created this.  NG only changed labels.
        """
    
        assert self.neighbor_list_created == False
        
        self.neighbor_list_created = True 
        #below had ncells instead of num_nodes.  This could be an issue? 
        #looks like it should still work though.  I think ncells was defined
        #as rows*cols before, so same as num_nodes     
        self.neighbor_nodes = - numpy.ones([self.num_nodes, 4], dtype=numpy.int)
        for r in xrange( 1, self.nrows-1 ):
            for c in xrange( 1, self.ncols-1 ):
                node_id = r * self.ncols + c
                self.neighbor_nodes[node_id,2] = node_id - 1   # left
                self.neighbor_nodes[node_id,0] = node_id + 1  # right
                self.neighbor_nodes[node_id,3] = node_id - self.ncols # bottom
                self.neighbor_nodes[node_id,1] = node_id + self.ncols # top
                
    def has_boundary_neighbor( self, id ):
        """
        Checks to see if one of the eight neighbor nodes of node with passed id
        is a boundary node.  Returns true if it has a boundary node, false if
        all neighbors are interior.
        
        Assumes that a valid id is passed.
        
        ng aug 2013
        """
        nbr_nodes=self.get_neighbor_list(id)
        #print "nbr nodes ", nbr_nodes
        diag_nbrs=self.get_diagonal_list(id)
        #print "diag nbrs ",diag_nbrs
        
        i=0
        #print 'id ',nbr_cells[i],' i ',i,' is interior? ',self.is_interior(nbr_cells[i])
        while (i<4 and self.is_interior(nbr_nodes[i]) ) :
            i += 1
            #print 'id ',nbr_nodes[i],' i ',i
        
        if i<4:
            return True
        else:
            r=0
            while (r<4 and self.is_interior(diag_nbrs[r]) ):
                r += 1
                #print 'id ',self.diag_nbrs[r],' r ',r
        
        if r<4 :
            return True
        else:
            return False
              
        
    
    def get_diagonal_list( self, id = -1 ):
        """
        If id is specified, returns a list of IDs for the diagonal cells of the node "id". 
        Otherwise, returns lists for all cells as a 2D array. 
        The list is in the order [topright, topleft, bottomleft, bottomright].
        
        NG didn't touch this, but she thinks this should be nodes, not cells.
        """
        #Added DEJH 051513
    
        if self.diagonal_list_created==False:
            self.create_diagonal_list()
        
        if id > -1:
            return self.diagonal_cells[id,:]
        else:
            return self.diagonal_cells

    def create_diagonal_list( self ):
        """
        Creates a list of IDs of the diagonal cells to each cell, as a 2D array. 
        Only interior cells are assigned neighbors; boundary cells get -1 for each neighbor. 
        The order of the diagonal cells is [topright, topleft, bottomleft, bottomright].
        
        NG didn't touch this, but she thinks this should be nodes, not cells.
        """
        #Added DEJH 051513
        
        assert self.diagonal_list_created == False
        
        self.diagonal_list_created = True
        self.diagonal_cells = - numpy.ones([self.ncells, 4], dtype=numpy.int)
        for r in xrange( 1, self.nrows-1 ):
            for c in xrange( 1, self.ncols-1 ):
                cell_id = r * self.ncols + c
                self.diagonal_cells[cell_id,2] = cell_id - self.ncols - 1   # bottom left
                self.diagonal_cells[cell_id,0] = cell_id + self.ncols + 1  # top right
                self.diagonal_cells[cell_id,3] = cell_id - self.ncols + 1 # bottom right
                self.diagonal_cells[cell_id,1] = cell_id + self.ncols - 1 # top left

    def is_interior( self, id ):
        """
        Returns True if the cell is an interior cell, False otherwise. 
        Interior status is indicated by a value of 0 in node_status.
    
        NG changed this.
        """
        #ng thinks there may be a problem here.
        #return self.boundary_ids[id] < 0
        #return self.node_status[id] < 1
        return self.node_status[id] == INTERIOR_NODE
        
    def get_boundary_code( self, id ):
        """
        Returns the boundary status of a node.
        ng june 2013
        """
        return self.node_status[id] 
        
#    def update_boundary_cell( self, id, u, bc = None ):
#        """
#        If cell ID tracks the value at another cell, this function sets
#        the value of U at cell ID to the value at its tracking cell in
#        BoundaryCondition "bc", which defaults to the self.default_bc.
#        If it is a fixed-gradient boundary, it updates the gradient.
#
#        .. todo::
#            Generalize, or create base-class version, that uses local
#            link length instead of self._dx.
#            
#        NG thinks this method is obsolete.  Maybe we can delete it?
#        """
#
#        if bc == None:
#            bc = self.default_bc
#        bid = self.boundary_ids[id]
#        if bid > -1:
#            if bc.boundary_code[bid] == bc.TRACKS_CELL_BOUNDARY:
#            	u[id] = u[bc.tracks_cell[bid]]
#            elif bc.boundary_code[bid] == bc.FIXED_GRADIENT_BOUNDARY:
#                u[id] = u[bc.tracks_cell[bid]] + bc.gradient[bid]*self._dx
#
#        if self.DEBUG_TRACK_METHODS:
#            print 'In RasterModelGrid.update_boundary_cell with cell',id
#            print 'Its value is now',u[id]
#            if bid < 0:
#                print 'It is NOT a boundary cell!'
#            if ( bid > -1 ) and ( bc.tracks_cell[bid] > -1 ):
#                print 'It tracks cell',bc.tracks_cell[bid]
#            else:
#                print 'It does not track any cell.'

    #-------------------------------------------------------------------
    # RasterModelGrid.calculate_cell_pair_gradient:
    #
    # 
    #-------------------------------------------------------------------
#   def calculate_cell_pair_gradient( self, cid1, cid2 ):
#   
#       my_faces = self.faces[cid1]
        #if self.

    def get_face_connecting_cell_pair( self, cid1, cid2 ):
        """
        TODO: UPDATE FOR NEW DATA STRUCTURES
        Returns the face that connects cells cid1 and cid2, or -1 if
        no such face is found.
        """
        
        #print 'Cells',cid1,'and',cid2
        for i in xrange( 0, 4 ):
            #print 'face num:',i
            fid = self.faces[cid1,i]
            #print 'face id:',fid,'fc:',self.fromcell[fid],'tc:',self.tocell[fid]
            if ( ( fid > -1 ) and ( ( cid1==self.fromcell[fid] and cid2==self.tocell[fid] )
               or ( cid1==self.tocell[fid] and cid2==self.fromcell[fid] ) ) ):
                return fid
        return -1
        
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
        
    def grid_coords_to_node_id(self, row, col):
        """
        Returns the ID of the node at the specified row and column of the raster
        grid.
        
        ng notes:
        note that syntax assumes that first row and column are 0,
        so max entry for a mg with 4 rows and 5 cols is row =3, col=4
        
        Example:
            
            >>> mg = RasterModelGrid(4, 5)
            >>> mg.grid_coords_to_node_id(2, 3)
            13
        """
        return row*self.ncols+col
        
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
        

if __name__ == '__main__':
    import doctest
    doctest.testmod()
