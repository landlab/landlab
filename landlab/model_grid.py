#! /usr/env/python
"""
Python implementation of ModelGrid, a class used to
create and manage grids for 2D numerical, flux-conservative models.

GT, July 2010
"""

import numpy
from numpy import *


class BoundaryCondition:
    """
    The BoundaryCondition class stores boundary condition information for
    a particular variable. It includes:
        * A list of boundary codes
          - 1 = fixed value
          - 2 = fixed gradient
          - 3 = tracks cell
        * A list of boundary gradients
        * A list of cell IDs to track
    """
    def __init__( self, n_boundary_cells = 0 ):
    
        # Create the 3 vectors
        self.boundary_code = zeros( n_boundary_cells, dtype = short )
        self.boundary_gradient = zeros( n_boundary_cells, dtype = double )
        self.tracks_cell = zeros( n_boundary_cells, dtype = long )
        
        # Define the boundary-type codes
        self.FIXED_VALUE_BOUNDARY = 1
        self.FIXED_GRADIENT_BOUNDARY = 2
        self.TRACKS_CELL_BOUNDARY = 3

        # Set defaults
        self.boundary_code[:] = self.FIXED_VALUE_BOUNDARY
        self.boundary_gradient[:] = 0.0
        self.tracks_cell[:] = -1


class ModelGrid:
    """
    This is the base class. The idea is to have at least two inherited
    classes, RasterModelGrid and DelaunayModelGrid, that can create and
    manage grids. To this might be added a GenericModelGrid, which would
    be an unstructured polygonal grid that doesn't necessarily obey or
    understand the Delaunay triangulation, but rather simply accepts
    an input grid from the user. Also a HexModelGrid for hexagonal.
    """

    """Base class for creating and manipulating 2D structured or
       unstructured grids for flux-conservative numerical models."""
    
    #-------------------------------------------------------------------
    def __init__( self ):
    
        pass

    #-------------------------------------------------------------------
    def initialize( self ):
    
        pass

    def create_cell_dvector( self ):
        """
        Returns a vector of floating point numbers the same length as 
        the number of cells.
        """
    
        return zeros( self.ncells )

    def create_boundary_condition( self ):
        """
        Creates, initializes, and returns a BoundaryCondition of the 
        proper size.
        """
    
        return BoundaryCondition( self.n_boundary_cells )
        
    def create_face_dvector( self ):
        """
        Returns a vector of floating point numbers the same length as 
        the number of interior faces.
        """
    
        return zeros( self.nfaces )

    def calculate_face_gradients( self, u ):
        """
        Calculates and returns gradients in u across all interior faces.

        .. todo::
            At the moment, we just use self.dx, but for unstructured
            grids this needs to use link length associated with each face!
        """
    
        g = zeros( self.nfaces )
        g = ( u[self.tocell[:]] - u[self.fromcell[:]] ) / self.dx
        #for i in arange( 0, self.nfaces ):
            #g[i] = ( u[self.tocell[i]] - u[self.fromcell[i]] ) / self.dx
            #print 'face',i,'from',self.fromcell[i],'to',self.tocell[i]
        return g
        
    def calculate_flux_divergences( self, q ):
        """
        At the moment, this is just a virtual function that does nothing,
        and in fact I don't understand why it doesn't print its "hello"
        message when called from an inherited class ...
        """
    
        print 'ModelGrid.calculate_flux_divergences here'
        #pass    # this is a virtual function
            
    def x( self, id ):
        """
        Returns the x coordinate of cell "id".
        """
        return self.cellx[id]
        
    def y( self, id ):
        """
        Returns the y coordinate of cell "id".
        """
        return self.celly[id]
        
    def get_cell_x_coords( self ):
        """
        Returns vector of cell x coordinates.
        """
        return self.cellx           

    def get_cell_y_coords( self ):
        """
        Returns vector of cell y coordinates.
        """
        return self.celly           

    def get_face_x_coords( self ):
        """
        Returns vector of cell x coordinates.
        """
        return self.facex           

    def get_face_y_coords( self ):
        """
        Returns vector of face y coordinates.
        """
        return self.facey           

    def get_interior_cells( self ):
        """
        Returns an integer vector of the IDs of all interior cells.
        """
        return self.interior_cells
                
    def assign_upslope_vals_to_faces( self, u, v=0 ):
        """
        Assigns to each face the values of u at whichever of its
        neighbors has a higher value of v. If v is omitted, uses u for
        both.
        """
        
        fv = zeros( self.nfaces )
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
        
        
class RasterModelGrid ( ModelGrid ):
    """
    This inherited class implements a regular, raster 2D grid with uniform
    cell dimensions.
    """

    def __init__( self, num_rows=0, num_cols=0, dx=1.0 ):
        """
        Optionally takes numbers of rows and columns and cell size as
        inputs. If this are given, calls initialize() to set up the grid.
        """
    
        self.ncells = num_rows * num_cols
        if self.ncells > 0:
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
        """
    
        # Debugging output flag
        self.debug = False
        
        # Basic info about raster size and shape
        self.nrows = num_rows
        self.ncols = num_cols
        self.ncells = num_rows * num_cols
        self.dx = dx
        self.cellarea = dx*dx
        
        # We need at least one row or column of boundary cells on each
        # side, so the grid has to be at least 3x3
        assert self.ncells >= 9

        # Record number of boundary and interior cells and the number
        # of interior faces. Ultimately, this info could be overridden
        # if using an irregular geometry of "interior" cells within the
        # rectangular domain. Note that we don't include any faces
        # between boundary cells.
        self.n_boundary_cells = 2 * ( num_rows - 2 ) + 2 * ( num_cols - 2 ) + 4
        self.n_interior_cells = self.ncells - self.n_boundary_cells
        self.nfaces = ( num_rows - 1 ) * ( num_cols - 2 ) + \
                      ( num_rows - 2 ) * ( num_cols - 1 )
        if self.debug:
            print self.nfaces
        
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
        self.fromcell = zeros( self.nfaces, dtype = int )
        self.tocell = zeros( self.nfaces, dtype = int )
        self.facex = zeros( self.nfaces )
        self.facey = zeros( self.nfaces )
        halfdx = self.dx / 2.0
        face_id = 0
        for r in xrange( 1, num_rows-1 ):
            for c in range( 1, num_cols ):
                self.fromcell[face_id] = r * num_cols + ( c - 1 )
                self.tocell[face_id] = self.fromcell[face_id] + 1
                self.facex[face_id] = c*self.dx - halfdx
                self.facey[face_id] = r*self.dx
                face_id += 1
        for r in xrange( 1, num_rows ):
            for c in range( 1, num_cols-1 ):
                self.fromcell[face_id] = ( r - 1 ) * num_cols + c
                self.tocell[face_id] = self.fromcell[face_id] + num_cols
                self.facex[face_id] = c*self.dx
                self.facey[face_id] = r*self.dx - halfdx
                face_id += 1
        if self.debug:
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
        self.faces = -ones( [self.ncells, 4], dtype=int )
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
        

        if self.debug:
            for i in xrange( 1, self.ncells ):
                print i,self.faces[i,:]
                
        # List of neighbors for each cell: we will start off with no
        # list. If a caller requests it via get_neighbor_list or
        # create_neighbor_list, we'll create it if necessary.
        self.neighbor_list_created = False
        if self.debug:
        	print 'Setting nlc flag'

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
        self.face_sign = zeros( [self.ncells, 4], dtype=short )
        self.face_sign[:,0:2] = -1
        self.face_sign[:,2:] = 1
        if self.debug:
            print 'face sign:',self.face_sign
        
        # Set up list of interior cells
        # (Note that this could be superceded if you wanted an irregular
        # boundary inside the rectangular grid)
        self.interior_cells = zeros( self.n_interior_cells, dtype=int )
        id = 0
        for r in xrange( 1, num_rows-1 ):
            for c in range( 1, num_cols-1 ):
                self.interior_cells[id] = r * num_cols + c
                id += 1
        
        if self.debug:
            print self.interior_cells
        
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
        self.boundary_ids = -ones( self.ncells, dtype=int )
        self.boundary_cells = zeros( self.n_boundary_cells, dtype=int )
        id = 0
        for r in xrange( 0, num_cols-1 ):       # Bottom
            self.boundary_cells[id] = r
            self.boundary_ids[r] = id
            id += 1
        for c in xrange( num_cols-1, self.ncells-num_cols, num_cols ):  # Right
            self.boundary_cells[id] = c
            self.boundary_ids[c] = id
            id += 1
        for r in xrange( self.ncells-1, num_cols*(num_rows-1), -1 ):       # Top
            self.boundary_cells[id] = r
            self.boundary_ids[r] = id
            id += 1
        for c in xrange( num_cols*(num_rows-1), 0, -num_cols ):  # Left
            self.boundary_cells[id] = c
            self.boundary_ids[c] = id
            id += 1
        
        if self.debug:
            print 'Boundary CIDs:',self.boundary_cells
            print 'Cell BIDs:',self.boundary_ids
            
        self.default_bc = BoundaryCondition( self.n_boundary_cells )
        
        # Store cell x and y coordinates
        self.cellx = zeros( self.ncells )
        self.celly = zeros( self.ncells )
        id = 0
        for r in range( 0, num_rows ):
            for c in xrange( 0, num_cols ):
                self.cellx[id] = c*self.dx
                self.celly[id] = r*self.dx
                id += 1
                
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
        #self.boundary_nbrs = zeros( self.n_boundary_cells, dtype=int )
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
        
        if self.debug:
            print 'tracks_cell:',bc.tracks_cell
    
    def calculate_flux_divergences( self, q ):
        """
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

        Note: the net flux is defined as positive outward, negative
        inward. In a diffusion problem, for example, one would use:
            .. math::
                {du \over dt} = source - fd
        where fd is "flux divergence".
        """
    
        if self.debug:
            print 'RasterModelGrid.calculate_flux_divergences here'
        
        fd = zeros( self.ncells )
        for cell in self.interior_cells:
            if self.debug:
                print 'Cell',cell
                print 'q:',q[self.faces[cell,0:4]]
            fd[cell] = ( -( q[self.faces[cell,2]]   # left face (positive=in)
                + q[self.faces[cell,3]] )           # bottom face (positive=in)
                + q[self.faces[cell,0]]             # right face (positive=out)
                + q[self.faces[cell,1]]             # top face (positive=out)
                  ) / self.dx
        return fd
        
    def calculate_flux_divergence( self, q, id ):
        """
        This is like calculate_flux_divergences (plural!), but only does
        it for cell "id".
        """
    
        if self.debug:
            print 'RasterModelGrid.calculate_flux_divergence here with cell',id
            print 'q:',q[self.faces[id,0:4]]
        fd = ( -( q[self.faces[id,2]]   # left face (positive=in)
            + q[self.faces[id,3]] )           # bottom face (positive=in)
            + q[self.faces[id,0]]             # right face (positive=out)
            + q[self.faces[id,1]]             # top face (positive=out)
              ) / self.dx
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
        for id in xrange( 0, self.n_boundary_cells ):
            if bc.boundary_code[id] == bc.TRACKS_CELL_BOUNDARY:
                u[self.boundary_cells[id]] = u[bc.tracks_cell[id]]
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
        """
    
        if bc==None:
            bc = self.default_bc
        for id in xrange( 0, self.n_boundary_cells ):
            if bc.boundary_code[id] == bc.TRACKS_CELL_BOUNDARY:
                u[self.boundary_cells[id]] = u[bc.tracks_cell[id]]
            elif bc.boundary_code[id] == bc.FIXED_GRADIENT_BOUNDARY:
                u[self.boundary_cells[id]] = u[bc.tracks_cell[id]] \
                                             + bc.gradient[id]*self.dx
        return u

    def cell_vector_to_raster( self, u ):
        """
        Converts cell vector u to a 2D array and returns it, so that it
        can be plotted, output, etc.
        """
    
        rast = zeros( [self.nrows, self.ncols] )
        id = 0
        for r in xrange( 0, self.nrows ):
            rast[r,:] = u[id:(id+self.ncols)]
            id += self.ncols
        return rast

    def get_neighbor_list( self, id = -1 ):
        """
        If id is specified, returns a list of neighboring cell IDs for
        the node "id". Otherwise, returns lists for all cells as a 2D
        array.
        """
    
        if self.neighbor_list_created==False:
            self.create_neighbor_list()
        
        if id > -1:
            return self.neighbor_cells[id,:]
        else:
            return self.neighbor_cells
            
    def create_neighbor_list( self ):
        """
        Creates a list of IDs of neighbor cells for each cell, as a
        2D array. Only interior cells are assigned neighbors; boundary
        cells get -1 for each neighbor.
        """
    
        assert self.neighbor_list_created == False
        
        self.neighbor_list_created = True       
        self.neighbor_cells = -ones( [self.ncells, 4], dtype=int )
        for r in xrange( 1, self.nrows-1 ):
            for c in xrange( 1, self.ncols-1 ):
                cell_id = r * self.ncols + c
                self.neighbor_cells[cell_id,2] = cell_id - 1   # left
                self.neighbor_cells[cell_id,0] = cell_id + 1  # right
                self.neighbor_cells[cell_id,3] = cell_id - self.ncols # bottom
                self.neighbor_cells[cell_id,1] = cell_id + self.ncols # top

    def is_interior( self, id ):
        """
        Returns True if the cell is an interior cell, False otherwise. 
        Interior status is indicated by a value of -1 in boundary_ids.
        """
    
        return self.boundary_ids[id] < 0
        
    def update_boundary_cell( self, id, u, bc = None ):
        """
        If cell ID tracks the value at another cell, this function sets
        the value of U at cell ID to the value at its tracking cell in
        BoundaryCondition "bc", which defaults to the self.default_bc.
        If it is a fixed-gradient boundary, it updates the gradient.

        .. todo::
            Generalize, or create base-class version, that uses local
            link length instead of self.dx.
        """

        if bc == None:
            bc = self.default_bc
        bid = self.boundary_ids[id]
        if bid > -1:
            if bc.boundary_code[bid] == bc.TRACKS_CELL_BOUNDARY:
            	u[id] = u[bc.tracks_cell[bid]]
            elif bc.boundary_code[bid] == bc.FIXED_GRADIENT_BOUNDARY:
                u[id] = u[bc.tracks_cell[bid]] + bc.gradient[bid]*self.dx

        if self.debug:
            print 'In RasterModelGrid.update_boundary_cell with cell',id
            print 'Its value is now',u[id]
            if bid < 0:
                print 'It is NOT a boundary cell!'
            if ( bid > -1 ) and ( bc.tracks_cell[bid] > -1 ):
                print 'It tracks cell',bc.tracks_cell[bid]
            else:
                print 'It does not track any cell.'

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
