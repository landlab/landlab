#! /usr/bin/env python

import numpy

from landlab.grid.base import ModelGrid, INTERIOR_NODE, BAD_INDEX_VALUE


def simple_poly_area(x, y):
    """
    Calculates and returns the area of a 2-D simple polygon.  Input vertices
    must be in sequence (clockwise or counterclockwise). *x* and *y* are
    arrays that give the x- and y-axis coordinates of the polygon's
    vertices.
        
    >>> import numpy as np
    >>> x = np.array([3., 1., 1., 3.])
    >>> y = np.array([1.5, 1.5, 0.5, 0.5])
    >>> simple_poly_area(x, y)
    2.0

    If the input coordinate arrays are 2D, calculate the area of each polygon.
    Note that when used in this mode, all polygons must have the same
    number of vertices, and polygon vertices are listed column-by-column.

    >>> x = np.array([[ 3.,  1.,  1.,  3.],
    ...               [-2., -2., -1., -1.]]).T
    >>> y = np.array([[1.5, 1.5, 0.5, 0.5],
    ...               [ 0.,  1.,  2.,  0.]]).T
    >>> simple_poly_area(x, y)
    array([ 2. ,  1.5])
    """
    # For short arrays (less than about 100 elements) it seems that the
    # Python sum is faster than the numpy sum. Likewise for the Python
    # built-in abs.
    return .5 * abs(sum(x[:-1] * y[1:] - x[1:] * y[:-1]) +
                    x[-1] * y[0] - x[0] * y[-1])


class VoronoiDelaunayGrid(ModelGrid):
    """
    This inherited class implements an unstructured grid in which cells are
    Voronoi polygons and nodes are connected by a Delaunay triangulation. Uses
    scipy.spatial module to build the triangulation.
    
    """
    def __init__(self):
        pass
        
    def initialize(self, x, y):
        """
        Creates an unstructured grid around the given (x,y) points.
        """
        
        assert type(x)==numpy.ndarray, 'x must be a numpy array'
        assert type(y)==numpy.ndarray, 'y must be a numpy array'
        assert len(x)==len(y), 'x and y arrays must have the same size'
        
        # Make a copy of the points in a 2D array (useful for calls to geometry
        # routines, but takes extra memory space).
        pts = numpy.zeros((len(x), 2))
        pts[:,0] = x
        pts[:,1] = y
        
        # NODES AND CELLS: Set up information pertaining to nodes and cells:
        #   - number of nodes
        #   - node x, y coordinates
        #   - default boundary status 
        #   - interior and boundary nodes
        #   - nodes associated with each cell and active cell
        #   - cells and active cells associated with each node 
        #     (or BAD_VALUE_INDEX if none)
        #
        # Assumptions we make here:
        #   - all interior (non-perimeter) nodes have cells (this should be 
        #       guaranteed in a Delaunay triangulation, but there may be 
        #       special cases)
        #   - all cells are active (later we'll build a mechanism for the user
        #       specify a subset of cells as active)
        #
        self.num_nodes = len(x)
        #print x, y
        self._node_x = x
        self._node_y = y
        [self.node_status, self.interior_nodes, self.boundary_nodes] = \
                self.find_perimeter_nodes(pts)
        self.num_cells = len(self.interior_nodes)
        self.num_activecells = self.num_cells
        [self.node_cell, self.cell_node] = self.setup_node_cell_connectivity( \
                                            self.node_status, self.num_cells)
        self.node_activecell = self.node_cell
        self.activecell_node = self.cell_node
        
        # ACTIVE CELLS: Construct Voronoi diagram and calculate surface area of
        # each active cell.
        from scipy.spatial import Voronoi
        vor = Voronoi(pts)
        self.active_cell_areas = numpy.zeros(self.num_activecells)
        for node in self.activecell_node:
            xv = vor.vertices[vor.regions[vor.point_region[node]],0]
            yv = vor.vertices[vor.regions[vor.point_region[node]],1]
            self.active_cell_areas[self.node_activecell[node]] = simple_poly_area(xv, yv)
        
        # LINKS: Construct Delaunay triangulation and construct lists of link
        # "from" and "to" nodes.
        [self.link_fromnode, self.link_tonode, self.active_links, self.face_width] \
                = self.create_links_and_faces_from_voronoi_diagram(vor)
        self.num_links = len(self.link_fromnode)
                    
        # LINKS: Calculate link lengths
        self.link_length = self.calculate_link_lengths(pts, self.link_fromnode, 
                                                       self.link_tonode)
                                                       
        # LINKS: inlink and outlink matrices
        self.setup_inlink_and_outlink_matrices()
        
        # ACTIVE LINKS: Create list of active links, as well as "from" and "to"
        # nodes of active links.
        self.reset_list_of_active_links()

        # LINKS: ID of corresponding face, if any
        self.link_face = numpy.zeros(self.num_links, dtype=int)+BAD_INDEX_VALUE  # make the list
        face_id = 0
        for link in self.active_links:
            self.link_face[link] = face_id
            face_id += 1
            

    def find_perimeter_nodes(self, pts):
    
        # Calculate the convex hull for the set of points
        from scipy.spatial import ConvexHull
        hull = ConvexHull(pts, qhull_options='Qc') # see below why we use 'Qt'
        
        # The ConvexHull object lists the edges that form the hull. We need to
        # get from this list of edges the unique set of nodes. To do this, we
        # first flatten the list of vertices that make up all the hull edges 
        # ("simplices"), so it becomes a 1D array. With that, we can use the set()
        # function to turn the array into a set, which removes duplicate vertices.
        # Then we turn it back into an array, which now contains the set of IDs for
        # the nodes that make up the convex hull.
        #   The next thing to worry about is the fact that the mesh perimeter 
        # might contain nodes that are co-planar (that is, co-linear in our 2D 
        # world). For example, if you make a set of staggered points for a
        # hexagonal lattice using make_hex_points(), there will be some 
        # co-linear points along the perimeter. The ones of these that don't 
        # form convex corners won't be included in convex_hull_nodes, but they
        # are nonetheless part of the perimeter and need to be included in
        # the list of boundary_nodes. To deal with this, we pass the 'Qt'
        # option to ConvexHull, which makes it generate a list of coplanar
        # points. We include these in our set of boundary nodes.
        convex_hull_nodes = numpy.array(list(set(hull.simplices.flatten())))
        coplanar_nodes = hull.coplanar[:,0]
        boundary_nodes = numpy.concatenate((convex_hull_nodes, coplanar_nodes))
    
        # Now we'll create the "node_status" array, which contains the code
        # indicating whether the node is interior and active (=0) or a
        # boundary (=1). This means that all perimeter (convex hull) nodes are
        # initially flagged as boundary code 1. An application might wish to change
        # this so that, for example, some boundaries are inactive.
        node_status = numpy.zeros(len(pts[:,0]), dtype=numpy.int8)
        node_status[boundary_nodes] = 1
        
        # It's also useful to have a list of interior nodes
        interior_nodes = numpy.where(node_status==0)[0]
    
        # Return the results
        return node_status, interior_nodes, boundary_nodes
        
        
    def setup_node_cell_connectivity(self, node_status, ncells):
        """
        Creates and returns the following arrays:
            1) for each node, the ID of the corresponding cell, or
                BAD_INDEX_VALUE if the node has no cell.
            2) for each cell, the ID of the corresponding node.
            
        Inputs:
            node_status: 1D numpy array containing the boundary status code
                         for each node
            ncells: the number of cells (must equal the number of occurrences of
                    INTERIOR_NODE in node_status)
                    
        Example:
            
            >>> import landlab as ll
            >>> vdmg = VoronoiDelaunayGrid()
            >>> ns = numpy.array([1,0,0,1,0])  # 3 interior, 2 boundary nodes
            >>> [node_cell,cell_node] = vdmg.setup_node_cell_connectivity(ns, 3)
            >>> node_cell[1:3]
            array([0, 1])
            >>> node_cell[0]==BAD_INDEX_VALUE
            True
            >>> cell_node
            array([1, 2, 4])
        """
        assert ncells==numpy.count_nonzero(node_status==INTERIOR_NODE), \
               'ncells must equal number of INTERIOR_NODE values in node_status'

        cell = 0
        node_cell = numpy.ones(len(node_status), dtype=int)*BAD_INDEX_VALUE
        cell_node = numpy.zeros(ncells, dtype=int)
        for node in range(len(node_cell)):
            if node_status[node] == INTERIOR_NODE:
                node_cell[node] = cell
                cell_node[cell] = node
                cell += 1
                
        return node_cell, cell_node
        

    def create_links_from_triangulation(self, tri):
        """
        From a Delaunay Triangulation of a set of points, contained in a
        scipy.spatial.Delaunay object "tri", creates and returns:
            1) a numpy array containing the ID of the "from" node for each link
            2) a numpy array containing the ID of the "to" node for each link
            3) the number of links in the triangulation
        
        Example:
            
            >>> vdmg = VoronoiDelaunayGrid()
            >>> pts = numpy.array([[ 0., 0.],[  1., 0.],[  1., 0.87],[-0.5, 0.87],[ 0.5, 0.87],[  0., 1.73],[  1., 1.73]])
            >>> from scipy.spatial import Delaunay
            >>> dt = Delaunay(pts)
            >>> [myfrom,myto,nl] = vdmg.create_links_from_triangulation(dt)
            >>> print myfrom, myto, nl
            [5 3 4 6 4 3 0 4 1 1 2 6] [3 4 5 5 6 0 4 1 0 2 4 2] 12
        
        """
    
        # Calculate how many links there will be and create the arrays.
        #
        # The number of links equals 3 times the number of triangles minus
        # half the number of shared links. Finding out the number of shared links
        # is easy: for every shared link, there is an entry in the tri.neighbors
        # array that is > -1 (indicating that the triangle has a neighbor opposite
        # a given vertex; in other words, two triangles are sharing an edge).
        #
        num_shared_links = numpy.count_nonzero(tri.neighbors>-1)
        num_links = 3*tri.nsimplex - num_shared_links/2
        link_fromnode = numpy.zeros(num_links, dtype=int)
        link_tonode = numpy.zeros(num_links, dtype=int)
        
        # Sweep through the list of triangles, assigning "from" and "to" nodes to
        # the list of links.
        #
        # The basic algorithm works as follows. For each triangle, we will add its
        # 3 edges as links. However, we have to make sure that each shared edge
        # is added only once. To do this, we keep track of whether or not each
        # triangle has been processed yet using a boolean array called "tridone".
        # When we look at a given triangle, we check each vertex in turn. If there
        # is no neighboring triangle opposite that vertex, then we need to add the
        # corresponding edge. If there is a neighboring triangle but we haven't
        # processed it yet, we also need to add the edge. If neither condition is
        # true, then this edge has already been added, so we skip it.
        link_id = 0
        tridone = numpy.zeros(tri.nsimplex, dtype=bool)    
        for t in range(tri.nsimplex):  # loop over triangles
            for i in range(0, 3):       # loop over vertices & neighbors
                if tri.neighbors[t,i] == -1 or not tridone[tri.neighbors[t,i]]:
                    link_fromnode[link_id] = tri.simplices[t,numpy.mod(i+1,3)]
                    link_tonode[link_id] = tri.simplices[t,numpy.mod(i+2,3)]
                    link_id += 1
            tridone[t] = True
    
        # Return the results
        return link_fromnode, link_tonode, num_links
    

    def is_valid_voronoi_ridge(self, vor, n):
        
        SUSPICIOUSLY_BIG = 40000000.0
        return vor.ridge_vertices[n][0]!=-1 and vor.ridge_vertices[n][1]!=-1 \
                and numpy.amax(numpy.abs(vor.vertices[vor.ridge_vertices[n]]))<SUSPICIOUSLY_BIG

        
        
    def create_links_and_faces_from_voronoi_diagram(self, vor):
        """
        From a Voronoi diagram object created by scipy.spatial.Voronoi(),
        builds and returns:
            1) Arrays of link "from" and "to" nodes
            2) Array of link IDs for each active link
            3) Array containing with of each face
        
        Inputs: vor = a scipy.spatial.Voronoi() object that was initialized
                      with the grid nodes.
                      
        Returns four 1D numpy arrays:
            
            link_fromnode = "from" node for each link (len=num_links)
            link_tonode   = "to" node for each link (len=num_links)
            active_links  = link ID for each active link (len=num_active_links)
            face_width    = width of each face (len=num_active_links
        
        Example:
            
            >>> vdmg = VoronoiDelaunayGrid()
            >>> pts = numpy.array([[ 0., 0.],[  1., 0.],[  1.5, 0.87],[-0.5, 0.87],[ 0.5, 0.87],[  0., 1.73],[  1., 1.73]])
            >>> from scipy.spatial import Voronoi
            >>> vor = Voronoi(pts)
            >>> [fr,to,al,fw] = vdmg.create_links_and_faces_from_voronoi_diagram(vor)
            >>> fr
            array([0, 0, 0, 1, 1, 3, 3, 6, 6, 6, 4, 4])
            >>> to
            array([3, 1, 4, 2, 4, 4, 5, 4, 2, 5, 2, 5])
            >>> al
            array([ 2,  4,  5,  7, 10, 11])
            >>> fw
            array([ 0.57669199,  0.57669199,  0.575973  ,  0.57836419,  0.575973  ,
                    0.57836419])
        """
        # Each Voronoi "ridge" corresponds to a link. The Voronoi object has an
        # attribute ridge_points that contains the IDs of the nodes on either
        # side (including ridges that have one of their endpoints undefined).
        # So, we set the number of links equal to the number of ridges.
        num_links = len(vor.ridge_points)
        
        # Create the arrays for link from and to nodes
        link_fromnode = -numpy.ones(num_links, dtype=int)
        link_tonode = -numpy.ones(num_links, dtype=int)
        
        # Ridges along the perimeter of the grid will have one of their 
        # endpoints undefined. The endpoints of each ridge are contained in
        # vor.ridge_vertices, and an undefined vertex is flagged with -1.
        # Ridges with both vertices defined correspond to faces and active 
        # links, while ridges with an undefined vertex correspond to inactive
        # links. So, to find the number of active links, we subtract from the
        # total number of links the number of occurrences of an undefined
        # vertex.
        num_active_links = num_links \
                    - numpy.count_nonzero(numpy.array(vor.ridge_vertices)==-1)
        #print 'num_links=', num_links,'num_active_links=',num_active_links
        
        # Create arrays for active links and width of faces (which are Voronoi
        # ridges).
        active_links = -numpy.ones(num_active_links, dtype=int)
        face_width = -numpy.ones(num_active_links)
    
        # Loop through the list of ridges. For each ridge, there is a link, and
        # its "from" and "to" nodes are the associated "points". In addition, if
        # the ridge endpoints are defined, we have a face and an active link,
        # so we add them to our arrays as well.
        j = 0
        for i in range(num_links):
            link_fromnode[i] = vor.ridge_points[i,0]
            link_tonode[i] = vor.ridge_points[i,1]
            face_corner1 = vor.ridge_vertices[i][0]
            face_corner2 = vor.ridge_vertices[i][1]
            if self.is_valid_voronoi_ridge(vor, i):  # means it's a valid face
                dx = vor.vertices[face_corner2,0]-vor.vertices[face_corner1,0]
                dy = vor.vertices[face_corner2,1]-vor.vertices[face_corner1,1]
                face_width[j] = numpy.sqrt(dx*dx+dy*dy)
                if abs(face_width[j])>=40000.0:
                    print 'link',i,'from',link_fromnode[i],'to',link_tonode[i],'has face width',face_width[j]
                    print vor.ridge_vertices[i]
                    print vor.vertices[vor.ridge_vertices[i]]
                    from scipy.spatial import voronoi_plot_2d
                    voronoi_plot_2d(vor)
                assert face_width[j] < 40000., 'face width must be less than earth circumference!'
                active_links[j] = i
                j += 1
        #print 'active links:',active_links
        #print 'vor ridge points:',vor.ridge_points
        return link_fromnode, link_tonode, active_links, face_width


