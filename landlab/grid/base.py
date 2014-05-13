#! /usr/env/python
"""
Python implementation of ModelGrid, a class used to
create and manage grids for 2D numerical models.

First version GT, July 2010
Last modified May 2014
"""

import numpy
import warnings

from landlab.testing.decorators import track_this_method
from landlab.utils import count_repeated_values
from landlab.utils.decorators import make_return_array_immutable, deprecated
from landlab.field import ModelDataFields
from . import grid_funcs as gfuncs


#: Indicates that an index is, in some way, *bad*.
BAD_INDEX_VALUE = numpy.iinfo(numpy.int).max


# Map names grid elements to the ModelGrid attribute that contains the count
# of that element in the grid.
_ARRAY_LENGTH_ATTRIBUTES = {
    'node': 'number_of_nodes',
    'cell': 'number_of_cells',
    'link': 'number_of_links',
    'face': 'number_of_faces',
    'core_node': 'number_of_core_nodes',
    'core_cell': 'number_of_core_cells',
    'active_link': 'number_of_active_links',
    'active_face': 'number_of_active_faces',
}

# Define the boundary-type codes
CORE_NODE = 0
FIXED_VALUE_BOUNDARY = 1
FIXED_GRADIENT_BOUNDARY = 2
TRACKS_CELL_BOUNDARY = 3

#: Indicates that a boundary node is *inactive*
<<<<<<< HEAD
=======
INACTIVE_BOUNDARY = 4
>>>>>>> FETCH_HEAD
CLOSED_BOUNDARY = 4

BOUNDARY_STATUS_FLAGS_LIST = [
    FIXED_VALUE_BOUNDARY,
    FIXED_GRADIENT_BOUNDARY,
    TRACKS_CELL_BOUNDARY,
    CLOSED_BOUNDARY,
<<<<<<< HEAD
=======
    #INACTIVE_BOUNDARY,
>>>>>>> FETCH_HEAD
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


def default_axis_names(n_dims):
    '''
    Returns a tuple of the default axis names.
    (Helper function)
    '''
    _DEFAULT_NAMES = ('z', 'y', 'x')
    return _DEFAULT_NAMES[- n_dims:]


def default_axis_units(n_dims):
    '''
    Returns a tuple of the default axis units.
    (Helper function)
    '''
    return ('-', ) * n_dims


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

    def __init__(self, **kwds):
        #print 'ModelGrid.__init__'
        super(ModelGrid, self).__init__()
        for element_name in _ARRAY_LENGTH_ATTRIBUTES:
            array_length = self.number_of_elements(element_name)
            try:
                self.new_field_location(element_name, array_length)
            except AttributeError:
                pass

        self.axis_name = kwds.get('axis_name', default_axis_names(self.ndim))
        self.axis_units = kwds.get('axis_units', default_axis_units(self.ndim))

    def _initialize( self ):
        pass

    @property
    def ndim(self):
        """Number of spatial dimensions of the grid"""
        return 2

    @property
    def node_index_at_cells(self):
        """Node ID associated with grid cells"""
        return self.cell_node

    @property
    def active_nodes(self):
        """
        Node IDs of all active (core & open boundary) nodes.
        core_nodes will return just core nodes.
        """
        (active_node_ids, ) = numpy.where(self.node_status != CLOSED_BOUNDARY)
        return active_node_ids

    @property
    def core_nodes(self):
        """
        Node IDs of all core nodes.
        """
        (core_node_ids, ) = numpy.where(self.node_status == CORE_NODE)
        return core_node_ids

    @property
    def node_boundary_status(self):
        """
        Node BC status codes for all nodes:
            0: core (nonboundary) node
            1: fixed value open boundary
            2: fixed gradient open boundary
            3: looped open boundary
            4: closed boundary
        """
        return self.node_status

    @property
    def open_nodes(self):
        """
        .. deprecated:: 0.6
            This terminology is no longer preferred, "active_nodes" is a synonym.
            
        Node id for all nodes not marked as a closed boundary
        """
        (open_node_ids, ) = numpy.where(self.node_status != CLOSED_BOUNDARY)
        return open_node_ids
    
    @property
    def open_boundary_nodes(self):
        """
        Node id of all open boundary nodes.
        """
        (open_boundary_node_ids, ) = numpy.where(
            (self.node_status != CLOSED_BOUNDARY) &
            (self.node_status != CORE_NODE))
        return open_boundary_node_ids
    
    @property
    def closed_boundary_nodes(self):
        """Node id of all closed boundary nodes.
        """
        (closed_boundary_node_ids, ) = numpy.where(
            self.node_status == CLOSED_BOUNDARY)
        return closed_boundary_node_ids
    
    @property
    def active_links(self):
        """Link IDs of all active links"""
        try:
            return self.active_link_ids
        except AttributeError:
            self._reset_list_of_active_links()
            return self.active_link_ids

    @property
    def node_index_at_active_cells(self):
        """
        .. deprecated:: 0.6
            Deprecated due to out-of-date terminology; 
            use :func:`node_index_at_core_cells` for an exact equivalent.
        Node ID associated with active grid cells
        """
        (active_cell_ids, ) = numpy.where(self.node_status == CORE_NODE)
        return active_cell_ids

    @property
    def node_index_at_core_cells(self):
        """
        Node ID associated with core grid cells
        """
        (core_cell_ids, ) = numpy.where(self.node_status == CORE_NODE)
        return core_cell_ids

    @property
    def active_cell_index_at_nodes(self):
        """
        .. deprecated:: 0.6
            "active" terminology now superceded by "core", unless explicitly
            referring to the open boundaries as well as core cells.
            
        Active cell ID associated with grid nodes.
        """
        return self.node_activecell

    @property
    def active_cell_index(self):
        """
        .. deprecated:: 0.6
            "active" terminology now superceded by "core", unless explicitly
            referring to the open boundaries as well as core cells.
        IDs of active cells
        """
        return self.active_cells
    
    @property
    def core_cell_index_at_nodes(self):
        """
        Core cell ID associated with grid nodes.
        """
        return self.node_corecell
        
    @property
    def core_cell_index(self):
        """
        IDs of core cells
        """
        return self.core_cells

    @property
    def node_index_at_link_head(self):
        """Node ID that defines the start of a link"""
        return self.link_fromnode

    @property
    def node_index_at_link_tail(self):
        """Node ID that defines the end of a link"""
        return self.link_tonode

    @property
    def face_index_at_links(self):
        """ID of the face associated with a link between two grid nodes"""
        return self.link_face
        
    @property
    def number_of_nodes(self):
        """Total number of nodes in the grid"""
        return self._num_nodes
    
    @property
    def number_of_cells(self):
        """Total number of cells in the grid"""
        return self._num_cells
    
    @property
    def number_of_links(self):
        """Total number of links in the grid"""
        return self._num_links
    
    @property
    def number_of_faces(self):
        """Total number of faces in the grid"""
        return self._num_faces
    
    @property
    def number_of_active_nodes(self):
        """Number of active nodes in the grid (i.e., core + open boundary)"""
        return self._num_active_nodes
    
    @property
    def number_of_core_nodes(self):
        """Number of core nodes in the grid (i.e., not boundaries)"""
        return self._num_core_nodes

    @property
    def number_of_active_cells(self):
        """
        Number of active cells in the grid (includes any possible
        boundary cells)
        """
        return self._num_active_cells
    
    @property
    def number_of_core_cells(self):
        """
        Number of core cells in the grid (excludes all boundary cells).
        """
        return self._num_core_cells
        

    @property
    def number_of_active_links(self):
        """Number of active links in the grid"""
        return self._num_active_links

    @property
    def number_of_active_faces(self):
        """Number of active faces in the grid"""
        return self._num_active_faces

    def number_of_elements(self, element_name):
        """Return the number of elements, given by the *element_name* string
        in a grid. *element_name* must be one of:
            * node
            * cell
            * link
            * face
            * core_node
            * core_cell
            * active_link
            * active_face
        """
        try:
            return getattr(self, _ARRAY_LENGTH_ATTRIBUTES[element_name])
        except KeyError:
            raise TypeError('element name not understood')

    def get_interior_nodes(self):
        """
        .. deprecated:: 0.6
            Deprecated due to outdated terminology;
            use :func:`get_core_nodes` instead.
            
        Return node IDs of all of a grid's interior nodes. Interior nodes
        are active nodes that are not on a boundary.
        """
        return numpy.where(self.node_status == CORE_NODE)[0]

    def get_core_nodes(self):
        """
        Return node IDs of all of a grid's core nodes.
        """
        return numpy.where(self.node_status == CORE_NODE)[0]

    @make_return_array_immutable
    def get_node_status(self):
        """
        Returns an array of node boundary-status codes.
        """
        return self.node_status

    @property
    @make_return_array_immutable
    def node_x(self):
        """X-coordinates of all nodes."""
        return self._node_x
    
    @property
    @make_return_array_immutable
    def node_y(self):
        """Y-coordinates of all nodes."""
        return self._node_y

    @make_return_array_immutable
    def node_axis_coordinates(self, axis=0):
        """
        Return node coordinates from a given *axis* (defaulting to 0). Axis
        numbering is the same as that for numpy arrays. That is, the zeroth
        axis is along the rows, and the first along the columns.
        """
        AXES = ('node_y', 'node_x')
        try:
            return getattr(self, AXES[axis])
        except IndexError:
            raise ValueError("'axis' entry is out of bounds")

    @property
    def axis_units(self):
        """A tuple of the units (as a string) for each of a grid's
        coordinates."""
        return self._axis_units

    @axis_units.setter
    def axis_units(self, new_units):
        """Set the units for each a grid's coordinates"""
        if len(new_units) != self.ndim:
            raise ValueError('length of units does not match grid dimension')
        self._axis_units = tuple(new_units)

    @property
    def axis_name(self):
        """A tuple of coordinate names for the grid"""
        return self._axis_name

    @axis_name.setter
    def axis_name(self, new_names):
        """Set the names of a grid's coordinates"""
        if len(new_names) != self.ndim:
            raise ValueError('length of names does not match grid dimension')
        self._axis_name = tuple(new_names)

    def create_node_array_zeros( self, name=None ):
        """
        Returns a 1D numpy array the same length as the number of nodes. If
        user gives optional argument 'name', we add this data to the grid with
        the specified name and return a reference to it; otherwise, we just
        create and return a 1D numpy array.
        """
        if name is None:
            return numpy.zeros(self.number_of_nodes)
        else: 
            self.add_zeros('node', name)
            return self.at_node[name]
        
    def create_active_link_array_zeros( self, name=None ):
        """
        Returns a 1D numpy array the same length as the number of nodes. If
        user gives optional argument 'name', we add this data to the grid with
        the specified name and return a reference to it; otherwise, we just
        create and return a 1D numpy array.
        """
        if name is None:
            return numpy.zeros(self.number_of_active_links)
        else: 
            self.add_zeros('link', name)
            return self.at_link[name]

    def zeros(self, **kwds):
        """
        Returns a numpy array of zeros that is the same length as the number
        of nodes in the grid. Use the *centering* keyword to return an
        array for other elements of the grid. *centering* is a string that is
        one of *node*, *cell*, *link*, or *face*.

        All other keywords are the same as for the numpy zeros function.
        """
        centering = kwds.pop('centering', 'node')
        try:
            return numpy.zeros(self.number_of_elements(centering), **kwds)
        except KeyError:
            raise TypeError(centering)

    def empty(self, **kwds):
        """
        Returns a numpy array of uninitialized values that is the same length
        as the number of nodes in the grid. Use the *centering* keyword to
        return an array for other elements of the grid. *centering* is a
        string that is one of *node*, *cell*, *link*, or *face*.

        All other keywords are the same as for the numpy zeros function.
        """
        centering = kwds.pop('centering', 'node')
        try:
            return numpy.empty(self.number_of_elements(centering), **kwds)
        except KeyError:
            raise TypeError(centering)

    def ones(self, **kwds):
        """
        Returns a numpy array of ones that is the same length as the number
        of nodes in the grid. Use the *centering* keyword to return an
        array for other elements of the grid. *centering* is a string that is
        one of *node*, *cell*, *link*, or *face*.

        All other keywords are the same as for the numpy zeros function.
        """
        centering = kwds.pop('centering', 'node')
        try:
            return numpy.ones(self.number_of_elements(centering), **kwds)
        except KeyError:
            raise TypeError(centering)

    def set_fixed_value_boundaries(self, node_ids):
        """
        Assignes FIXED_VALUE_BOUNDARY status to specified nodes.
        """
        self.node_status[node_ids] = FIXED_VALUE_BOUNDARY
        self._reset_list_of_active_links()

    @track_this_method
    def calculate_diff_at_links(self, node_values, out=None):
        """
        Calculates the difference in quantity *node_values* at every link
        in the grid.
        Note that this is tonode-fromnode along links, and is thus equivalent to
        positive gradient up.
        """
        return gfuncs.calculate_diff_at_links(self, node_values, out=out)
        
    @track_this_method
    def calculate_diff_at_active_links(self, node_values, out=None):
        """
        Calculates the difference in quantity *node_values* at each active link
        in the grid.
        Note that this is tonode-fromnode along links, and is thus equivalent to
        positive gradient up.
        """
        return gfuncs.calculate_diff_at_active_links(self, node_values,
                                                     out=out)
        
    @track_this_method
    def calculate_gradients_at_links(self, node_values, out=None):
        """
        Calculates the gradient in quantity *node_values* at every link
        in the grid.
        This method follows the convention POSITIVE UP.
        """
        return gfuncs.calculate_gradients_at_links(self, node_values, out=out)
        
    @track_this_method
    def calculate_gradients_at_active_links(self, node_values, out=None):
        """
        Calculates the gradient in quantity *node_values* at each active link
        in the grid.
        This method follows the convention POSITIVE UP.
        """
        return gfuncs.calculate_gradients_at_active_links(self, node_values,
                                                          out=out)
        
    @track_this_method
    def calculate_gradients_at_active_links_slow(self, s, gradient=None):
        """
        .. deprecated:: 0.1
            Use :func:`calculate_gradients_at_active_links`
        
        Calculates the gradient in quantity s at each active link in the grid.
        """
        if gradient==None:
            gradient = numpy.zeros(self.number_of_active_links)
            
        assert (len(gradient) == self.number_of_active_links), \
                "len(gradient)!=number_of_active_links"
                
        active_link_id = 0
        for link_id in self.active_link_ids:
            gradient[active_link_id] = (s[self.link_tonode[link_id]]
                                        -s[self.link_fromnode[link_id]]) / \
                                        self.link_length[link_id]
            active_link_id += 1
        
        return gradient
        
    def resolve_values_on_links(self, link_values, out=None):
        """
        Resolves values provided defined on links into the x and y directions.
        Returns values_along_x, values_along_y
        """
        return gfuncs.resolve_values_on_links(self, link_values, out=out)

    def resolve_values_on_active_links(self, link_values, out=None):
        """
        Resolves values provided defined on active links into the x and y 
        directions.
        Returns values_along_x, values_along_y
        """
        return gfuncs.resolve_values_on_active_links(self, link_values, out=out)
        
    def calculate_flux_divergence_at_active_cells(self, active_link_flux, 
                                                  net_unit_flux=None):
        """
        .. deprecated:: 0.6
            Uses outdated terminology; use the exact equivalent
            :func:`calculate_flux_divergence_at_core_nodes` instead.
            
        Given an array of fluxes along links, computes the net total flux
        within each cell, divides by cell area, and stores the result in
        net_unit_flux.
        
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
            
            >>> from landlab import RasterModelGrid
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
            
            >>> divflux = rmg.zeros(centering='active_cell') # outside loop
            
        Then do this inside the loop:
            
            >>> divflux = rmg.calculate_flux_divergence_at_active_cells(flux, divflux)
            
        In this case, the function will not have to create the divflux array.
        
        Note this method is untested with looped boundary conditions.
        """
        
        if self.DEBUG_TRACK_METHODS:
            print 'ModelGrid.calculate_flux_divergence_at_active_cells'
            
        assert (len(active_link_flux) == self.number_of_active_links), \
               "incorrect length of active_link_flux array"
            
        # If needed, create net_unit_flux array
        if net_unit_flux is None:
            net_unit_flux = numpy.zeros(self.number_of_active_cells)
        else:
            net_unit_flux[:] = 0.
            
        assert (len(net_unit_flux)) == self.number_of_active_cells
        
        node_net_unit_flux = self.calculate_flux_divergence_at_nodes(active_link_flux)
                
        net_unit_flux = node_net_unit_flux[self.activecell_node]
                
        return net_unit_flux
        
        
    def calculate_flux_divergence_at_core_nodes(self, active_link_flux, 
                                                  net_unit_flux=None):
        """
        Given an array of fluxes along links, computes the net total flux
        within each cell, divides by cell area, and stores the result in
        net_unit_flux.
        
        The function works by calling calculate_flux_divergence_at_nodes, then
        slicing out only the values at core nodes. Therefore, it is slower
        than calculate_flux_divergence_at_nodes, even though it returns a
        shorter list of numbers.
        
        The input active_link_flux should be flux of
        something (e.g., mass, momentum, energy) per unit face width, positive
        if flowing in the same direction as its link, and negative otherwise.
        There should be one value per active link. Returns an array of net
        total flux per unit area, one value per core node (creates this
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
            
            >>> from landlab import RasterModelGrid
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
            >>> divflux = rmg.calculate_flux_divergence_at_core_nodes(flux)
            >>> divflux
            array([ 2.,  4., -2.,  0.,  1., -4.])
            
        If calculate_gradients_at_core_nodes is called inside a loop, you can
        improve speed slightly by creating an array outside the loop. For 
        example, do this once, before the loop:
            
            >>> divflux = rmg.zeros(centering='active_cell') # outside loop
            
        Then do this inside the loop:
            
            >>> divflux = rmg.calculate_flux_divergence_at_core_nodes(flux, divflux)
            
        In this case, the function will not have to create the divflux array.
        
        Note this method is untested with looped boundary conditions.
        """
        
        if self.DEBUG_TRACK_METHODS:
            print 'ModelGrid.calculate_flux_divergence_at_core_nodes'
            
        assert (len(active_link_flux) == self.number_of_active_links), \
               "incorrect length of active_link_flux array"
            
        # If needed, create net_unit_flux array
        if net_unit_flux is None:
            net_unit_flux = numpy.zeros(self.number_of_core_nodes)
        else:
            net_unit_flux[:] = 0.
            
        assert (len(net_unit_flux)) == self.number_of_core_nodes
        
        node_net_unit_flux = self.calculate_flux_divergence_at_nodes(active_link_flux)
                
        net_unit_flux = node_net_unit_flux[self.corecell_node]
                
        return net_unit_flux


    def calculate_flux_divergence_at_active_cells_slow(self, active_link_flux, 
                                                  net_unit_flux=False):
        """
        .. deprecated:: 0.1
            Use :func:`calculate_flux_divergence_at_active_cells`
            
        Original, slower version of calculate_flux_divergence_at_active_cells, 
        using a for-loop instead of simply calling the node-based version of
        the method. Kept here as illustration of what the method is intended
        to do.
        """
        
        if self.DEBUG_TRACK_METHODS:
            print 'ModelGrid.calculate_flux_divergence_at_active_cells'
            
        assert (len(active_link_flux) == self.number_of_active_links), \
               "incorrect length of active_link_flux array"
            
        # If needed, create net_unit_flux array
        if net_unit_flux==False:
            net_unit_flux = numpy.zeros(self.number_of_active_cells)
        else:
            net_unit_flux[:] = 0.
            
        assert (len(net_unit_flux))==self.number_of_active_cells
        
        # For each active link, add up the flux out of the "from" cell and 
        # into the "to" cell.
        active_link_id = 0
        for link_id in self.active_link_ids:
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
        net_unit_flux = net_unit_flux / self._cell_areas
        
        return net_unit_flux

    @track_this_method
    def calculate_flux_divergence_at_nodes(self, active_link_flux, out=None):
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
        
        This method is untested with looped boundary conditions.
        """
        return gfuncs.calculate_flux_divergence_at_nodes(self, active_link_flux,
                                                        out=out)
        
                        
    @property
    @make_return_array_immutable
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
            return self._cell_areas
        except AttributeError:
            return self._setup_cell_areas_array()

    @property
    @make_return_array_immutable    
    def forced_cell_areas(self):
        """
        Returns an array of grid cell areas. In the cases of inactive nodes,
        this method forces the area of those nodes so it can return an nnodes-
        long array. For a raster, it assumes areas are equal to the normal case.
        For a voronoi...?
        """
        try:
            return self.forced_cell_areas
        except AttributeError:
            return self._setup_cell_areas_array_force_inactive()    
            
    def _setup_cell_areas_array_force_inactive(self):
        '''
        Sets up an array of cell areas which is nnodes long. Nodes which have 
        cells receive the area of that cell. Nodes which do not receive
        numpy.nan entries.
        Note this method is typically only required for some raster purposes, 
        and is overridden in raster.py. It is unlikely this parent method will
        ever need to be called.
        '''
        self.forced_cell_areas = numpy.empty(self.number_of_nodes)
        self.forced_cell_areas.fill(numpy.nan)
        cell_node_ids = self.get_active_cell_node_ids()
        self.forced_cell_areas[cell_node_ids] = self.cell_areas

    def get_active_cell_node_ids( self ):
        """
        Returns an integer vector of the node IDs of all active (i.e., core +
        open boundary) cells.
        get_core_cell_node_ids may be preferable.
        """
        return self.activecell_node
        
        
    def get_core_cell_node_ids( self ):
        """
        Returns an integer vector of the node IDs of all core cells.
        """
        return self.corecell_node

        
    def get_active_link_connecting_node_pair(self, node1, node2):
        """
        Returns the ID number of the active link that connects the given pair of
        nodes, or BAD_INDEX_VALUE if not found.
        This method is slow, and can only take single ints as *node1* and 
        *node2*. It should ideally be overridden for optimal functionality in
        more specialized grid modules (e.g., raster).
        
        Example:
            
            >>> import landlab as ll
            >>> rmg = ll.RasterModelGrid(4, 5)
            >>> rmg.get_active_link_connecting_node_pair(8, 3)
            array([2])
        """
        active_link = BAD_INDEX_VALUE
        for alink in xrange(0, self.number_of_active_links):
            link_connects_nodes = (
                (self.activelink_fromnode[alink] == node1 and
                self.activelink_tonode[alink] == node2) or
                (self.activelink_tonode[alink] == node1 and
                self.activelink_fromnode[alink] == node2))
            if link_connects_nodes:
                active_link = alink
                break
        return numpy.array([active_link])
        

    @property
    def active_link_length(self):
        """Returns the lengths of all active links, in ID order"""
        return self.link_length[self.active_link_ids]

    @property
    def link_length(self):
        """Returns the lengths of all links, in ID order"""
        try:
            return self._link_length
        except AttributeError:
            return self.calculate_link_length()

    def min_active_link_length(self):
        """
        Returns the horizontal length of the shortest active link in the grid.
        """
        return numpy.amin(self.link_length[self.active_link_ids])

    def max_active_link_length(self):
        """
        Returns the horizontal length of the longest active link in the grid.
        """
        return numpy.amax(self.link_length[self.active_link_ids])

    def calculate_link_length(self):
        """
        Calculates, returns, and stores as a property of the grid the lengths
        of all the links in the grid.
        """
        if not hasattr(self, '_link_length'):
            self._link_length = self.empty(centering='link')
        dx = (self.node_x[self.node_index_at_link_head] -
              self.node_x[self.node_index_at_link_tail])
        dy = (self.node_y[self.node_index_at_link_head] -
              self.node_y[self.node_index_at_link_tail])
        numpy.sqrt(dx ** 2 + dy **2, out=self._link_length)
        return self.link_length

    def assign_upslope_vals_to_active_links( self, u, v=[0] ):
        """
        Assigns to each active link the value of u at whichever of its
        neighbors has a higher value of v. If v is omitted, uses u for
        both.
        """
        fv = numpy.zeros(self.number_of_active_links)
        if len(v) < len(u):
            for i in xrange(0, self.number_of_active_links):
                fv[i] = max(u[self.activelink_fromnode[i]], 
                            u[self.activelink_tonode[i]] )
        else:
            for i in xrange(0, self.number_of_active_links):
                if v[self.activelink_fromnode[i]] > v[self.activelink_tonode[i]]:
                    fv[i] = u[self.activelink_fromnode[i]]
                else:
                    fv[i] = u[self.activelink_tonode[i]]
        return fv
        
    def _reset_list_of_active_links(self):
        """
        Creates or resets a list of active links. We do this by sweeping
        through the given lists of from and to nodes, and checking the status
        of these as given in the node_status list. A link is active if both its
        nodes are active interior points, or if one is an active interior and
        the other is an active boundary.
        """
        if self.DEBUG_TRACK_METHODS:
            print 'ModelGrid._reset_list_of_active_links'
            
        fromnode_status = self.node_status[self.link_fromnode]
        tonode_status = self.node_status[self.link_tonode]

        active_links = (((fromnode_status == CORE_NODE) & ~
                         (tonode_status == CLOSED_BOUNDARY)) |
                        ((tonode_status == CORE_NODE) & ~
                         (fromnode_status == CLOSED_BOUNDARY)))

        (self.active_link_ids, ) = numpy.where(active_links)

        self._num_active_links = len(self.active_link_ids)
        self._num_active_faces = self._num_active_links
        self.activelink_fromnode = self.link_fromnode[self.active_link_ids]
        self.activelink_tonode = self.link_tonode[self.active_link_ids]
        
        # Set up active inlink and outlink matrices
        self._setup_active_inlink_and_outlink_matrices()

    def set_nodata_nodes_to_inactive(self, node_data, nodata_value):
        """
        Sets self.node_status to CLOSED_BOUNDARY for all nodes whose value of
        node_data is equal to the nodata_value.
        
        Example:
            
            >>> import landlab as ll
            >>> mg = ll.RasterModelGrid(3, 4, 1.0)
            >>> mg.node_status
            array([1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1], dtype=int8)
            >>> h = numpy.array([-9999,-9999,-9999,-9999,-9999,-9999,12345.,0.,-9999,0.,0.,0.])
            >>> mg.set_nodata_nodes_to_inactive(h, -9999)
            >>> mg.node_status
            array([4, 4, 4, 4, 4, 4, 0, 1, 4, 1, 1, 1], dtype=int8)
        """
        
        # Find locations where value equals the NODATA code and set these nodes
        # as inactive boundaries.
        nodata_locations = numpy.nonzero(node_data==nodata_value)
        self.node_status[nodata_locations] = CLOSED_BOUNDARY
        
        # Recreate the list of active cell IDs
        node_ids = numpy.array(range(0, self.number_of_nodes))
        self.activecell_node = node_ids[numpy.where(self.node_status != CLOSED_BOUNDARY)]
        self.corecell_node = node_ids[numpy.where(self.node_status == CORE_NODE)]
        
        # Recreate the list of active links
        self._reset_list_of_active_links()
        
        
    def max_of_link_end_node_values(self, node_data):
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
            >>> mg.max_of_link_end_node_values(h)
            array([ 2.,  8.,  6.,  8.,  8.,  3.,  3.])
        """
        return numpy.maximum(node_data[self.activelink_fromnode],
                             node_data[self.activelink_tonode])
        
    def calculate_numbers_of_node_neighbors(self):
        """
        Calculates the number of neighboring nodes for each node, and returns
        the result as a 1D numpy array. Used to find the maximum number of
        neighbors, so that inlink and outlink matrices can be dimensioned
        accordingly. Assumes that self.number_of_nodes, self.link_fromnode, and
        self.link_tonode have already been set up.
        
        Algorithm works by simply looping through all links; for each, the 
        endpoints are neighbors of one another, so we increment the number of
        neighbors for both the endpoint nodes.
        """
        num_nbrs = numpy.zeros(self.number_of_nodes, dtype=int)
        for link in range(self.number_of_links):
            num_nbrs[self.link_fromnode[link]] += 1
            num_nbrs[self.link_tonode[link]] += 1
        return num_nbrs


    def _setup_inlink_and_outlink_matrices(self):
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
        self.node_inlink_matrix = - numpy.ones(
            (self.max_num_nbrs, self.number_of_nodes), dtype=numpy.int)
        self.node_outlink_matrix = - numpy.ones(
            (self.max_num_nbrs, self.number_of_nodes), dtype=numpy.int)

        # Set up the inlink arrays
        tonodes = self.link_tonode
        self.node_numinlink = numpy.bincount(tonodes,
                                             minlength=self.number_of_nodes)

        counts = count_repeated_values(self.link_tonode)
        for (count, (tonodes, link_ids)) in enumerate(counts):
            self.node_inlink_matrix[count][tonodes] = link_ids

        # Set up the outlink arrays
        fromnodes = self.link_fromnode
        self.node_numoutlink = numpy.bincount(fromnodes,
                                              minlength=self.number_of_nodes)
        counts = count_repeated_values(self.link_fromnode)
        for (count, (fromnodes, link_ids)) in enumerate(counts):
            self.node_outlink_matrix[count][fromnodes] = link_ids
                
        
    def _setup_active_inlink_and_outlink_matrices(self):
        """
        Creates data structures to record the numbers of active inlinks and 
        active outlinks for each node. These data structures are equivalent to
        the "regular" inlink and outlink matrices, except that it uses the IDs
        of active links (only).
        """
        # Create active in-link and out-link matrices.
        self.node_active_inlink_matrix = - numpy.ones(
            (self.max_num_nbrs, self.number_of_nodes), dtype=numpy.int)
        self.node_active_outlink_matrix = - numpy.ones(
            (self.max_num_nbrs, self.number_of_nodes), dtype=numpy.int)

        # Set up the inlink arrays
        tonodes = self.activelink_tonode
        self.node_numactiveinlink = numpy.bincount(
            tonodes, minlength=self.number_of_nodes)

        counts = count_repeated_values(self.activelink_tonode)
        for (count, (tonodes, active_link_ids)) in enumerate(counts):
            self.node_active_inlink_matrix[count][tonodes] = active_link_ids

        # Set up the outlink arrays
        fromnodes = self.activelink_fromnode
        self.node_numactiveoutlink = numpy.bincount(
            fromnodes, minlength=self.number_of_nodes)
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
        plt.plot(self._node_x[self.core_nodes], 
                 self._node_y[self.core_nodes], 'go')
        plt.plot(self._node_x[self.boundary_nodes], 
                 self._node_y[self.boundary_nodes], 'ro')
                 
        # Draw links
        for i in range(self.number_of_links):
            plt.plot([self._node_x[self.link_fromnode[i]],
                     self._node_x[self.link_tonode[i]]],
                     [self._node_y[self.link_fromnode[i]],
                     self._node_y[self.link_tonode[i]]], 'k-')
                     
        # Draw active links
        for link in self.active_link_ids:
            plt.plot([self._node_x[self.link_fromnode[link]],
                     self._node_x[self.link_tonode[link]]],
                     [self._node_y[self.link_fromnode[link]],
                     self._node_y[self.link_tonode[link]]], 'g-')
                     
        # If caller asked for a voronoi diagram, draw that too
        if draw_voronoi!=None:
            from scipy.spatial import Voronoi, voronoi_plot_2d
            pts = numpy.zeros((self.number_of_nodes, 2))
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
            return ~ (self.node_status[ids] == CORE_NODE)
        else:
            return self.node_status[ids] == boundary_flag
    
    def get_boundary_nodes(self):
        """
        Returns ids of all open and closed boundary nodes in the grid.
        """
        return numpy.where(self.node_status != 0)[0]
    
    def _assign_boundary_nodes_to_grid_sides(self):
        """
        For each boundary node, determines whether it belongs to the left, 
        right, top or bottom of the grid, based on its distance from the grid's
        centerpoint (mean (x,y) position). Returns lists of nodes on each of 
        the four grid sides. Assumes self.node_status, self.number_of_nodes, 
        self.boundary_nodes, self._node_x, and self._node_y have been initialized.
        
        Example:

            >>> import landlab as ll
            >>> m = ll.HexModelGrid(5, 3, 1.0)
            >>> [l,r,t,b] = m._assign_boundary_nodes_to_grid_sides()
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
        .. deprecated:: 0.6
            Due to imprecise terminology. Use :func:`set_closed_boundaries`
            instead.
        Handles boundary conditions by setting each of the four sides of the 
        rectangular grid to either 'inactive' or 'active (fixed value)' status.
        Arguments are booleans indicating whether the bottom, right, top, and
        left are inactive (True) or not (False).
        
        For an inactive boundary:
            - the nodes are flagged CLOSED_BOUNDARY
            - the links between them and the adjacent core nodes are
              inactive (so they appear on link-based lists, but not
              active_link-based lists)
              
        This means that if you call the calculate_gradients_at_active_links
        method, the inactive boundaries will be ignored: there can be no
        gradients or fluxes calculated, because the links that connect to that
        edge of the grid are not included in the calculation. So, setting a
        grid edge to CLOSED_BOUNDARY is a convenient way to impose a no-flux
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
        >>> rmg.number_of_active_links
        30
        >>> rmg.node_status
        array([1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1], dtype=int8)
        >>> rmg.set_inactive_boundaries(False, False, True, True)
        >>> rmg.number_of_active_links
        21
        >>> rmg.node_status
        array([1, 1, 1, 4, 0, 0, 1, 4, 0, 0, 0, 1, 4, 0, 0, 1, 4, 4, 4], dtype=int8)
        """
        if self.DEBUG_TRACK_METHODS:
            print 'ModelGrid.set_inactive_boundaries'
            
        [left_edge, right_edge, top_edge, bottom_edge] = \
                self._assign_boundary_nodes_to_grid_sides()
            
        if bottom_is_inactive:
            self.node_status[bottom_edge] = CLOSED_BOUNDARY
        else:
            self.node_status[bottom_edge] = FIXED_VALUE_BOUNDARY

        if right_is_inactive:
            self.node_status[right_edge] = CLOSED_BOUNDARY
        else:
            self.node_status[right_edge] = FIXED_VALUE_BOUNDARY
            
        if top_is_inactive:
            self.node_status[top_edge] = CLOSED_BOUNDARY
        else:
            self.node_status[top_edge] = FIXED_VALUE_BOUNDARY

        if left_is_inactive:
            self.node_status[left_edge] = CLOSED_BOUNDARY
        else:
            self.node_status[left_edge] = FIXED_VALUE_BOUNDARY
        
        self._reset_list_of_active_links()

    def set_inactive_nodes(self, nodes):
        """
        Sets the given nodes' boundary condition statuses to INACTIVE (==4),
        and resets the list of active links to reflect any changes.
        """
        self.node_status[nodes] = CLOSED_BOUNDARY
        node_ids = numpy.array(range(0, self.number_of_nodes))
        self.activecell_node = node_ids[numpy.where(self.node_status != CLOSED_BOUNDARY)]
        self.corecell_node = node_ids[numpy.where(self.node_status == CORE_NODE)]
        self._reset_list_of_active_links()

    def get_distances_of_nodes_to_point(self, tuple_xy, get_az=None, node_subset=numpy.nan, out_distance=None, out_azimuth=None):
        """
        Returns an array of distances for each node to a provided point.
        If "get_az" is set to 'angles', returns both the distance array and an
        array of azimuths from up/north. If it is set to 'displacements', it
        returns the azimuths as a 2xnnodes array of x and y displacements.
        If it is not set, returns just the distance array.
        If "node_subset" is set as an ID, or list/array/etc of IDs method
        returns just the distance (and optionally azimuth) for that node.
        Point is provided as a tuple (x,y).
        If out_distance (& out_azimuth) are provided, these arrays are used to
        store the outputs. This is recommended for memory management reasons if
        you are working with node subsets.
        
        ***Developer's note***
        Once you start working with node subsets in Landlab, which can change
        size between loops, it's quite possible for Python's internal memory
        management to crap out after large numbers of loops (~>10k). This is
        to do with the way it block allocates memory for arrays of differing
        lengths, then cannot free this memory effectively.
        The solution - as implemented here - is to pre-allocate all arrays as
        nnodes long, then only work with the first [len_subset] entries by
        slicing (in a pseudo-C-style). Care has to be taken not to
        "accidentally" allow Python to allocate a new array you don't have
        control over.
        Then, to maintain efficient memory allocation, we create some "dummy"
        nnode-long arrays to store intermediate parts of the solution in.
        """
        assert isinstance(tuple_xy, tuple)
        assert len(tuple_xy) == 2
            
        if numpy.any(numpy.isnan(node_subset)):
            subset_flag = False
        else:
            subset_flag = True
        
        if subset_flag:
            if type(node_subset) == int:
                node_subset = numpy.array([node_subset])
        
        azimuths_as_displacements = numpy.empty((2, self.number_of_nodes))
        dummy_nodes_1 = numpy.empty(self.number_of_nodes)
        dummy_nodes_2 = numpy.empty(self.number_of_nodes)
        dummy_nodes_3 = numpy.empty(self.number_of_nodes)
        dummy_bool = numpy.empty(self.number_of_nodes, dtype=bool)
        
        if out_distance is None:
            try:
                out_distance = numpy.empty(node_subset.size)
            except:
                out_distance = self.empty(centering='node')
        else:
            if subset_flag:
                assert out_distance.size == node_subset.size
            else:
                assert out_distance.size == self.number_of_nodes
        if out_azimuth is None and get_az:
            try:
                out_azimuth = numpy.empty((2,node_subset.size))
            except:
                out_azimuth = numpy.empty((2, self.number_of_nodes))
            #only one of these colums will get used if get_az == 'angles'
        elif out_azimuth is not None:
            if subset_flag:
                if get_az == 'displacements':
                    assert out_azimuth.shape == (2,node_subset.shape[1])
                elif get_az == 'angles':
                    assert out_azimuth.size == node_subset.size
            else:
                assert out_azimuth.size == self.number_of_nodes
            

        try:
            len_subset = node_subset.size
        except:
            len_subset = self.number_of_nodes
        
        try:
            azimuths_as_displacements[0,:len_subset] = self.node_x[node_subset]-tuple_xy[0]
            azimuths_as_displacements[1,:len_subset] = self.node_y[node_subset]-tuple_xy[1]
        except:
            azimuths_as_displacements[0] = (self.node_x-tuple_xy[0])
            azimuths_as_displacements[1] = (self.node_y-tuple_xy[1])

        numpy.square(azimuths_as_displacements[0,:len_subset], out=dummy_nodes_1[:len_subset])
        numpy.square(azimuths_as_displacements[1,:len_subset], out=dummy_nodes_2[:len_subset])
        numpy.add(dummy_nodes_1[:len_subset], dummy_nodes_2[:len_subset], out=dummy_nodes_3[:len_subset])
        numpy.sqrt(dummy_nodes_3[:len_subset], out=out_distance)
        
        if get_az:
            if get_az == 'displacements':
                out_azimuth[:len_subset] = self.azimuths_as_displacements[:len_subset]
                return out_distance, out_azimuth
            elif get_az == 'angles':
                try:
                    numpy.divide(azimuths_as_displacements[1,:len_subset],
                                 azimuths_as_displacements[0,:len_subset],
                                 out=dummy_nodes_1[:len_subset])
                    numpy.arctan(dummy_nodes_1[:len_subset],
                                 out=dummy_nodes_2[:len_subset]) #"angle_to_xaxis"
                except: #These cases have the impact right on a gridline.
                    if len_subset == 1: #this is the single node case, point directly N or S of the node of interest
                        if azimuths_as_displacements[1]<0:
                            out_azimuth[0] = numpy.pi
                        else:
                            out_azimuth[0] = 0.
                    else: #general case with whole array, with the impact right one one of the gridlines
                        num_nonzero_nodes = numpy.count_nonzero(azimuths_as_displacements[0,:len_subset])
                        dummy_nodes_3[:num_nonzero_nodes] = azimuths_as_displacements[0,:len_subset].nonzero() #"nonzero_nodes"
                        nonzero_nodes = dummy_nodes_3[:num_nonzero_nodes]
                        numpy.divide(azimuths_as_displacements[1,:len_subset][nonzero_nodes],
                                     azimuths_as_displacements[0,:len_subset][nonzero_nodes],
                                     out=dummy_nodes_1[:len_subset][nonzero_nodes])
                        numpy.arctan(dummy_nodes_1[:len_subset][nonzero_nodes],
                                     out=dummy_nodes_2[:len_subset][nonzero_nodes]) #"angle_to_xaxis"
                        ##angle_to_xaxis = numpy.arctan(y_displacement[nonzero_nodes]/x_displacement[nonzero_nodes])
                        numpy.less(azimuths_as_displacements[0,:len_subset][nonzero_nodes], 0., out=dummy_bool[:len_subset][nonzero_nodes])
                        out_azimuth[nonzero_nodes][dummy_bool[:len_subset][nonzero_nodes]] = 1.5*numpy.pi-dummy_nodes_2[:len_subset][nonzero_nodes][dummy_bool[:len_subset][nonzero_nodes]]
                        numpy.logical_not(dummy_bool[:len_subset][nonzero_nodes], out=dummy_bool[:len_subset][nonzero_nodes])
                        out_azimuth[nonzero_nodes][dummy_bool[:len_subset][nonzero_nodes]] = 0.5*numpy.pi-dummy_nodes_2[:len_subset][nonzero_nodes][dummy_bool[:len_subset][nonzero_nodes]]
                        #out_azimuth[nonzero_nodes] = numpy.where(azimuths_as_displacements[0,:len_subset][nonzero_nodes]<0,
                        #                                         1.5*numpy.pi-dummy_nodes_2[:len_subset][nonzero_nodes],
                        #                                         0.5*numpy.pi-dummy_nodes_2[:len_subset][nonzero_nodes])
                    num_zero_nodes = len_subset - num_nonzero_nodes
                    #numpy.equal(azimuths_as_displacements[0,:len_subset], 0., out=dummy_bool[:num_zero_nodes]) #not clear if this will work, as output might be 2D
                    dummy_nodes_3[:num_zero_nodes] = numpy.where(azimuths_as_displacements[0,:len_subset]==0.)[0] #"zero_nodes" ##POTENTIAL MEMORY LEAK REMAINS
                    zero_nodes = dummy_nodes_3[:num_zero_nodes]
                    numpy.less(azimuths_as_displacements[1,:len_subset][zero_nodes], 0., out=dummy_bool[:num_zero_nodes][zero_nodes])
                    out_azimuth[zero_nodes][dummy_bool[:num_zero_nodes][zero_nodes]] = numpy.pi
                    numpy.logical_not(dummy_bool[:num_zero_nodes][zero_nodes], out=dummy_bool[:num_zero_nodes][zero_nodes])
                    out_azimuth[zero_nodes][dummy_bool[:num_zero_nodes][zero_nodes]] = 0.        
                    #out_azimuth[zero_nodes] = numpy.where(azimuths_as_displacements[1,:len_subset][zero_nodes]<0.,numpy.pi,0.)
                else: #the normal case
                    numpy.sign(azimuths_as_displacements[0,:len_subset],
                               out=dummy_nodes_1[:len_subset])
                    numpy.subtract(1., dummy_nodes_1[:len_subset],
                                   out=dummy_nodes_3[:len_subset])
                    numpy.multiply(dummy_nodes_3[:len_subset], 0.5*numpy.pi,
                                   out=dummy_nodes_1[:len_subset])
                    numpy.subtract(0.5*numpy.pi, dummy_nodes_2[:len_subset],
                                   out=dummy_nodes_3[:len_subset])
                    if out_azimuth is not None:
                        numpy.add(dummy_nodes_1[:len_subset],
                                  dummy_nodes_3[:len_subset],
                                  out=out_azimuth)
                    else:
                        numpy.add(dummy_nodes_1[:len_subset],
                                  dummy_nodes_3[:len_subset],
                                  out=out_azimuth[0,:])
                    ##azimuth_array = ((1.-numpy.sign(x_displacement))*0.5)*numpy.pi + (0.5*numpy.pi-angle_to_xaxis) #duplicated by the above
                if out_azimuth.shape[0] == 2 and len(out_azimuth.shape) == 2:
                    return out_distance, out_azimuth[0,:]
                else:
                    return out_distance, out_azimuth
            else:
                print "Option set for get_az not recognised. Should be 'displacements' or 'angles'."
        else:
            return out_distance
            
    def build_all_node_distances_azimuths_maps(self):
        """
        This function creates and stores in the grid field two nnodes*nnodes 
        arrays that map the distances and azimuths of all nodes in the grid to 
        all nodes in the grid.
        This is useful if your module needs to make repeated lookups of distances
        between the same nodes, but does potentially use up a lot of memory so
        should be used with caution.
        The map is symmetrical, so it does not matter whether rows are "from" or
        "to".
        The arrays are called:
            self.all_node_distances_map
            self.all_node_azimuths_map
        
        The method returns these two arrays as output.
        """
        
        self.all_node_distances_map = numpy.empty((self.number_of_nodes,
                                                  self.number_of_nodes))
        self.all_node_azimuths_map = numpy.empty((self.number_of_nodes,
                                                 self.number_of_nodes))
        
        node_coords = numpy.empty((self.number_of_nodes, 2))
        node_coords[:,0] = self.node_x
        node_coords[:,1] = self.node_y
        
        for i in xrange(self.number_of_nodes):
            self.all_node_distances_map[i,:], self.all_node_azimuths_map[i,:] = self.get_distances_of_nodes_to_point((node_coords[i,0],node_coords[i,1]), get_az='angles')
        
        assert numpy.all(self.all_node_distances_map >= 0.)
        
        return self.all_node_distances_map, self.all_node_azimuths_map
        

if __name__ == '__main__':
    import doctest
    doctest.testmod()
