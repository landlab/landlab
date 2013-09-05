#! /usr/env/python
"""
Python implementation of ModelGrid, a class used to
create and manage grids for 2D numerical models.

First version GT, July 2010
Last modified August 2013
"""

import numpy
import warnings

from landlab.utils import count_repeated_values
from landlab.field import ModelDataFields

BAD_INDEX_VALUE = numpy.iinfo(numpy.int).max

# Map names grid elements to the ModelGrid attribute that contains the count
# of that element in the grid.
_ARRAY_LENGTH_ATTRIBUTES = {
    'node': 'num_nodes',
    'cell': 'num_active_cells',
    'link': 'num_active_links',
    'face': 'num_faces',
}

# Define the boundary-type codes
INTERIOR_NODE = 0
FIXED_VALUE_BOUNDARY = 1
FIXED_GRADIENT_BOUNDARY = 2
TRACKS_CELL_BOUNDARY = 3
INACTIVE_BOUNDARY = 4

BOUNDARY_STATUS_FLAGS_LIST = [
    FIXED_VALUE_BOUNDARY,
    FIXED_GRADIENT_BOUNDARY,
    TRACKS_CELL_BOUNDARY,
    INACTIVE_BOUNDARY,
]
BOUNDARY_STATUS_FLAGS = set(BOUNDARY_STATUS_FLAGS_LIST)


class Error(Exception):
    """
    Base class for exceptions from this module.
    """
    pass


def _sort_points_into_quadrants(x, y, nodes):
    """
    Divide points with locations given in the *x*, and *y* arrays into north,
    south, east, and west quadrants. Returns nodes contained in quadrants
    (west, east, north, south).

    >>> x = numpy.array([0, 1, 0, -1])
    >>> y = numpy.array([1, 0, -1, 0])
    >>> nodes = numpy.array([1, 2, 3, 4])
    >>> _sort_points_into_quadrants(x, y, nodes)
    (array([4]), array([2]), array([1]), array([3]))
    """
    above_x_axis = y > 0
    right_of_y_axis = x > 0
    closer_to_y_axis = numpy.abs(y) >= numpy.abs(x)

    north_nodes = nodes[above_x_axis & closer_to_y_axis]
    south_nodes = nodes[(~ above_x_axis) & closer_to_y_axis]
    east_nodes = nodes[right_of_y_axis & (~ closer_to_y_axis)]
    west_nodes = nodes[(~ right_of_y_axis) & (~ closer_to_y_axis)]

    return (west_nodes, east_nodes, north_nodes, south_nodes)


class ModelGrid(ModelDataFields):
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

    # Debugging flags (if True, activates some output statements)
    DEBUG_VERBOSE = False
    DEBUG_TRACK_METHODS = False

    #-------------------------------------------------------------------
    def __init__(self):
        for centering in ['node', 'cell', 'link', 'face', ]:
            try:
                self.add_group(centering, self._array_length(centering))
            except AttributeError:
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

    def _array_length(self, centering):
        try:
            return getattr(self, _ARRAY_LENGTH_ATTRIBUTES[centering])
        except KeyError:
            raise TypeError('centering value not understood')

    def zeros(self, **kwds):
        """
        Returns a numpy array of zeros that is the same length as the number
        of nodes in the grid. Use the *centering* keyword to return an
        array for other elements of the grid. *centering* is a string that is
        one of *node*, *cell*, *link*, or *face*.

        All other keywords are the same as for the numpy zeros function.
        """
        centering = kwds.pop('centering', 'node')
        return numpy.zeros(self._array_length(centering), **kwds)

    def empty(self, **kwds):
        """
        Returns a numpy array of uninitialized values that is the same length
        as the number of nodes in the grid. Use the *centering* keyword to
        return an array for other elements of the grid. *centering* is a
        string that is one of *node*, *cell*, *link*, or *face*.

        All other keywords are the same as for the numpy zeros function.
        """
        centering = kwds.pop('centering', 'node')
        return numpy.empty(self._array_length(centering), **kwds)

    def ones(self, **kwds):
        """
        Returns a numpy array of ones that is the same length as the number
        of nodes in the grid. Use the *centering* keyword to return an
        array for other elements of the grid. *centering* is a string that is
        one of *node*, *cell*, *link*, or *face*.

        All other keywords are the same as for the numpy zeros function.
        """
        centering = kwds.pop('centering', 'node')
        return numpy.ones(self._array_length(centering), **kwds)

    def set_fixed_value_boundaries(self, node_ids):
        """
        Assignes FIXED_VALUE_BOUNDARY status to specified nodes.
        """
        self.node_status[node_ids] = FIXED_VALUE_BOUNDARY
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
        
        #print 'cfdn here'
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
                #if from_node==46:
                #    print('node '+str(from_node)+' net='+str(net_unit_flux[from_node]))
                #    print 'fw=', self.face_width[self.link_face[link_id]]
            if to_cell != BAD_INDEX_VALUE:
                net_unit_flux[to_node] -= total_flux / \
                                          self.active_cell_areas[to_cell]
                #if to_node==46:
                #    print('node '+str(to_node)+' net='+str(net_unit_flux[to_node]))
                #    print 'fw=', self.face_width[self.link_face[link_id]]
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
        
        .. todo: 
            Could be useful to return a Numpy array of x-coords of the cell's
            corners.
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
        """
        Return node coordinates from a given *axis* (defaulting to 0). Axis
        numbering is the same as that for numpy arrays. That is, the zeroth
        axis is along the rows, and the first along the columns.
        """
        assert(axis in (0, 1))

        if axis == 0:
            return self.get_node_y_coords()
        else:
            return self.get_node_x_coords()

    def get_coordinate_units(self, axis=0):
        """
        .. todo:
            GT: coordinate units should be model/component dependent.
        """
        assert(axis in (0, 1))

        if axis == 0:
            return 'degrees_north'
        else:
            return 'degrees_east'

    def get_coordinate_name(self, axis=0):
        """
        .. todo:
            GT: coordinate units should be model/component dependent.
        """
        assert(axis in (0, 1))

        if axis == 0:
            return 'latitude'
        else:
            return 'longitude'

    def get_cell_areas(self, *args):
        """get_cell_areas([ids])

        This function returns the area of an interior cell, or an array of
        all interior cell areas if *ids* is not given. *ids* can be either
        a scalar index, or an array of cell indices.

        .. note::
            It is up to a specific grid class, which has inherited from
            ModelGrid, to construct its own cell_areas array.
        """
        assert(len(args) <= 1)

        if len(args) == 0:
            return self.cell_areas
        else:
            return self.cell_areas[args[0]]
    
    @property
    def cell_areas(self):
        """
        Returns an array of grid-cell areas.

        .. note::
            Sometimes it may make sense for a grid to not always calculate
            its cell areas but, instead, only calculate them once they are
            required. In such cases, the grid class must implement a
            _setup_cell_areas_array method, which will be called the first
            time cell areas are requested.
        """
        try:
            return self.active_cell_areas
        except AttributeError:
            return self._setup_cell_areas_array()

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
            
            >>> import landlab as ll
            >>> rmg = ll.RasterModelGrid(4, 5)
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
        return numpy.amin(self.link_length[self.active_links])

    def assign_upslope_vals_to_active_links( self, u, v=0 ):
        """
        Assigns to each active link the value of u at whichever of its
        neighbors has a higher value of v. If v is omitted, uses u for
        both.
        """
        fv = numpy.zeros( self.num_active_links )
        if len(v) < len(u):
            for i in xrange( 0, self.num_active_links ):
                fv[i] = max( u[self.activelink_fromnode[i]], 
                             u[self.activelink_tonode[i]] )
        else:
            for i in xrange( 0, self.num_active_links ):
                if v[self.activelink_fromnode[i]] > v[self.activelink_tonode[i]]:
                    fv[i] = u[self.activelink_fromnode[i]]
                else:
                    fv[i] = u[self.activelink_tonode[i]]
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

        active_links = (((fromnode_status == INTERIOR_NODE) & ~
                         (tonode_status == INACTIVE_BOUNDARY)) |
                        ((tonode_status == INTERIOR_NODE) & ~
                         (fromnode_status == INACTIVE_BOUNDARY)))

        (self.active_links, ) = numpy.where(active_links)

        self.num_active_links = len(self.active_links)
        self.activelink_fromnode = self.link_fromnode[self.active_links]
        self.activelink_tonode = self.link_tonode[self.active_links]
        
        # Set up active inlink and outlink matrices
        self.setup_active_inlink_and_outlink_matrices()

    def deactivate_nodata_nodes(self, node_data, nodata_value):
        """
        Sets self.node_status to INACTIVE_BOUNDARY for all nodes whose value of
        node_data is equal to the nodata_value.
        
        Example:
            
            >>> import landlab as ll
            >>> mg = ll.RasterModelGrid(3, 4, 1.0)
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
        self.node_status[nodata_locations] = INACTIVE_BOUNDARY
        
        # Recreate the list of active cell IDs
        node_ids = numpy.array(range(0,self.num_nodes))
        self.activecell_node = node_ids[numpy.where(self.node_status == INTERIOR_NODE)]
        
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
            
            >>> import landlab as ll
            >>> mg = ll.RasterModelGrid(3, 4, 1.0)
            >>> h = numpy.array([2.,2.,8.,0.,8.,0.,3.,0.,5.,6.,8.,3.])
            >>> mg.active_link_max(h)
            array([ 2.,  8.,  6.,  8.,  8.,  3.,  3.])
        """
        return numpy.maximum(node_data[self.activelink_fromnode],
                             node_data[self.activelink_tonode])
        
    def calculate_link_lengths(self, pts, link_from, link_to):
        """
        Calculates and returns length of links between nodes.
        
        Inputs: pts - Nx2 numpy array containing (x,y) values
                link_from - 1D numpy array containing index numbers of nodes at 
                            starting point ("from") of links
                link_to - 1D numpy array containing index numbers of nodes at 
                          ending point ("to") of links
                          
        Returns: 1D numpy array containing horizontal length of each link
        
        Example:
            
            >>> pts = numpy.array([[0.,0.],[3.,0.],[3.,4.]]) # 3:4:5 triangle
            >>> lfrom = numpy.array([0,1,2])
            >>> lto = numpy.array([1,2,0])
            >>> mg = ModelGrid()
            >>> ll = mg.calculate_link_lengths(pts, lfrom, lto)
            >>> ll
            array([ 3.,  4.,  5.])
        """
        dx = pts[link_to,0]-pts[link_from,0]
        dy = pts[link_to,1]-pts[link_from,1]
        link_length = numpy.sqrt( dx*dx + dy*dy )
        return link_length
        
        
    def calculate_numbers_of_node_neighbors(self):
        """
        Calculates the number of neighboring nodes for each node, and returns
        the result as a 1D numpy array. Used to find the maximum number of
        neighbors, so that inlink and outlink matrices can be dimensioned
        accordingly. Assumes that self.num_nodes, self.link_fromnode, and
        self.link_tonode have already been set up.
        
        Algorithm works by simply looping through all links; for each, the 
        endpoints are neighbors of one another, so we increment the number of
        neighbors for both the endpoint nodes.
        """
        num_nbrs = numpy.zeros(self.num_nodes, dtype=int)
        for link in range(self.num_links):
            num_nbrs[self.link_fromnode[link]] += 1
            num_nbrs[self.link_tonode[link]] += 1
        return num_nbrs


    def setup_inlink_and_outlink_matrices(self):
        """
        Creates data structures to record the numbers of inlinks and outlinks
        for each node. An inlink of a node is simply a link that has the node as
        its "to" node, and an outlink is a link that has the node as its "from".
        
        We store the inlinks in an NM-row by num_nodes-column matrix called
        node_inlink_matrix. NM is the maximum number of neighbors for any node.
        
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
            
        """
        
        # Find the maximum number of neighbors for any node
        num_nbrs = self.calculate_numbers_of_node_neighbors()
        self.max_num_nbrs = numpy.amax(num_nbrs)

        # Create active in-link and out-link matrices.
        self.node_inlink_matrix = - numpy.ones((self.max_num_nbrs, self.num_nodes), dtype=numpy.int)
        self.node_outlink_matrix = - numpy.ones((self.max_num_nbrs, self.num_nodes), dtype=numpy.int)

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
                
        
    def setup_active_inlink_and_outlink_matrices(self):
        """
        Creates data structures to record the numbers of active inlinks and 
        active outlinks for each node. These data structures are equivalent to
        the "regular" inlink and outlink matrices, except that it uses the IDs
        of active links (only).
        """
        # Create active in-link and out-link matrices.
        self.node_active_inlink_matrix = - numpy.ones((self.max_num_nbrs, self.num_nodes),
                                                       dtype=numpy.int)
        self.node_active_outlink_matrix = - numpy.ones((self.max_num_nbrs, self.num_nodes),
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

    
    def display_grid(self, draw_voronoi=False):
        """
        Displays the grid (mainly for purposes of debugging/testing and
        visual examples).
        """
        import matplotlib.pyplot as plt
        
        # Plot nodes, colored by boundary vs interior
        plt.plot(self._node_x[self.interior_nodes], 
                 self._node_y[self.interior_nodes], 'go')
        plt.plot(self._node_x[self.boundary_nodes], 
                 self._node_y[self.boundary_nodes], 'ro')
                 
        # Draw links
        for i in range(self.num_links):
            plt.plot([self._node_x[self.link_fromnode[i]],
                     self._node_x[self.link_tonode[i]]],
                     [self._node_y[self.link_fromnode[i]],
                     self._node_y[self.link_tonode[i]]], 'k-')
                     
        # Draw active links
        for link in self.active_links:
            plt.plot([self._node_x[self.link_fromnode[link]],
                     self._node_x[self.link_tonode[link]]],
                     [self._node_y[self.link_fromnode[link]],
                     self._node_y[self.link_tonode[link]]], 'g-')
                     
        # If caller asked for a voronoi diagram, draw that too
        if draw_voronoi!=None:
            from scipy.spatial import Voronoi, voronoi_plot_2d
            pts = numpy.zeros((self.num_nodes, 2))
            pts[:,0] = self._node_x
            pts[:,1] = self._node_y
            vor = Voronoi(pts)
            voronoi_plot_2d(vor)
                     
        plt.show()
        
        
    def is_boundary(self, ids, boundary_flag=None):
        """
        Check if nodes at given *ids* are boundary nodes. Use the
        *boundary_flag* to specify a particular boundary type status flag.
        """
        if boundary_flag is None:
            return ~ (self.node_status[ids] == INTERIOR_NODE)
        else:
            return self.node_status[ids] == boundary_flag
    
    def assign_boundary_nodes_to_grid_sides(self):
        """
        For each boundary node, determines whether it belongs to the left, 
        right, top or bottom of the grid, based on its distance from the grid's
        centerpoint (mean (x,y) position). Returns lists of nodes on each of 
        the four grid sides. Assumes self.node_status, self.num_nodes, 
        self.boundary_nodes, self._node_x, and self._node_y have been initialized.
        
        Example:

            >>> import landlab as ll
            >>> m = ll.HexModelGrid(5, 3, 1.0)
            >>> [l,r,t,b] = m.assign_boundary_nodes_to_grid_sides()
            >>> l
            array([ 7, 12,  3], dtype=int32)
            >>> r
            array([11, 15,  6], dtype=int32)
            >>> t
            array([16, 18, 17], dtype=int32)
            >>> b
            array([0, 2, 1], dtype=int32)
        """
        # Calculate x and y distance from centerpoint
        dx = self._node_x[self.boundary_nodes] - numpy.mean(self._node_x)
        dy = self._node_y[self.boundary_nodes] - numpy.mean(self._node_y)

        return _sort_points_into_quadrants(dx, dy, self.boundary_nodes)
        
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
        
        >>> import landlab as ll
        >>> rmg = ll.HexModelGrid(5, 3, 1.0) # rows, columns, spacing
        >>> rmg.num_active_links
        30
        >>> rmg.node_status
        array([1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1], dtype=int8)
        >>> rmg.set_inactive_boundaries(False, False, True, True)
        >>> rmg.num_active_links
        21
        >>> rmg.node_status
        array([1, 1, 1, 4, 0, 0, 1, 4, 0, 0, 0, 1, 4, 0, 0, 1, 4, 4, 4], dtype=int8)
        """
        if self.DEBUG_TRACK_METHODS:
            print 'ModelGrid.set_inactive_boundaries'
            
        [left_edge, right_edge, top_edge, bottom_edge] = \
                self.assign_boundary_nodes_to_grid_sides()
            
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
                

if __name__ == '__main__':
    import doctest
    doctest.testmod()
