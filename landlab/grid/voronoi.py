#! /usr/env/python
"""
Python implementation of VoronoiDelaunayGrid, a class used to create and manage
unstructured, irregular grids for 2D numerical models.

Do NOT add new documentation here. Grid documentation is now built in a semi-
automated fashion. To modify the text seen on the web, edit the files
`docs/text_for_[gridfile].py.txt`.
"""
import numpy as np
from scipy.spatial import Voronoi
from six.moves import range

from landlab.core.utils import (
    argsort_points_by_x_then_y,
    as_id_array,
    sort_points_by_x_then_y,
)
from landlab.grid.base import BAD_INDEX_VALUE, CORE_NODE, ModelGrid

from .decorators import return_readonly_id_array


def simple_poly_area(x, y):
    """Calculates and returns the area of a 2-D simple polygon.

    Input vertices must be in sequence (clockwise or counterclockwise). *x*
    and *y* are arrays that give the x- and y-axis coordinates of the
    polygon's vertices.

    Parameters
    ----------
    x : ndarray
        x-coordinates of of polygon vertices.
    y : ndarray
        y-coordinates of of polygon vertices.

    Returns
    -------
    out : float
        Area of the polygon

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.voronoi import simple_poly_area
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
    return 0.5 * abs(sum(x[:-1] * y[1:] - x[1:] * y[:-1]) + x[-1] * y[0] - x[0] * y[-1])


class VoronoiDelaunayGrid(ModelGrid):
    """
    This inherited class implements an unstructured grid in which cells are
    Voronoi polygons and nodes are connected by a Delaunay triangulation. Uses
    scipy.spatial module to build the triangulation.

    Create an unstructured grid from points whose coordinates are given
    by the arrays *x*, *y*.

    Parameters
    ----------
    x : array_like
        x-coordinate of points
    y : array_like
        y-coordinate of points
    reorient_links (optional) : bool
        whether to point all links to the upper-right quadrant

    Returns
    -------
    VoronoiDelaunayGrid
        A newly-created grid.

    Examples
    --------
    >>> from numpy.random import rand
    >>> from landlab.grid import VoronoiDelaunayGrid
    >>> x, y = rand(25), rand(25)
    >>> vmg = VoronoiDelaunayGrid(x, y)  # node_x_coords, node_y_coords
    >>> vmg.number_of_nodes
    25

    >>> import numpy as np
    >>> x = [0, 0.1, 0.2, 0.3,
    ...      1, 1.1, 1.2, 1.3,
    ...      2, 2.1, 2.2, 2.3,]
    >>> y = [0, 1, 2, 3,
    ...      0, 1, 2, 3,
    ...      0, 1, 2, 3]
    >>> vmg = VoronoiDelaunayGrid(x, y)
    >>> vmg.node_x # doctest: +NORMALIZE_WHITESPACE
    array([ 0. ,  1. ,  2. ,
            0.1,  1.1,  2.1,
            0.2,  1.2,  2.2,
            0.3,  1.3,  2.3])
    >>> vmg.node_y # doctest: +NORMALIZE_WHITESPACE
    array([ 0.,  0.,  0.,
            1.,  1.,  1.,
            2.,  2.,  2.,
            3.,  3.,  3.])
    >>> vmg.adjacent_nodes_at_node
    array([[ 1,  3, -1, -1, -1, -1],
           [ 2,  4,  3,  0, -1, -1],
           [ 5,  4,  1, -1, -1, -1],
           [ 4,  6,  0,  1, -1, -1],
           [ 5,  7,  6,  3,  1,  2],
           [ 8,  7,  4,  2, -1, -1],
           [ 7,  9,  3,  4, -1, -1],
           [ 8, 10,  9,  6,  4,  5],
           [11, 10,  7,  5, -1, -1],
           [10,  6,  7, -1, -1, -1],
           [11,  9,  7,  8, -1, -1],
           [10,  8, -1, -1, -1, -1]])
    """

    def __init__(self, x=None, y=None, reorient_links=True, **kwds):
        """
        Create a Voronoi Delaunay grid from a set of points.

        Create an unstructured grid from points whose coordinates are given
        by the arrays *x*, *y*.

        Parameters
        ----------
        x : array_like
            x-coordinate of points
        y : array_like
            y-coordinate of points
        reorient_links (optional) : bool
            whether to point all links to the upper-right quadrant
        xy_of_reference : tuple, optional
            Coordinate value in projected space of (0., 0.)
            Default is (0., 0.)

        Returns
        -------
        VoronoiDelaunayGrid
            A newly-created grid.

        Examples
        --------
        >>> from numpy.random import rand
        >>> from landlab.grid import VoronoiDelaunayGrid
        >>> x, y = rand(25), rand(25)
        >>> vmg = VoronoiDelaunayGrid(x, y)  # node_x_coords, node_y_coords
        >>> vmg.number_of_nodes
        25
        """
        if (x is not None) and (y is not None):
            self._initialize(x, y, reorient_links)
        super(VoronoiDelaunayGrid, self).__init__(**kwds)

    def _initialize(self, x, y, reorient_links=True):
        """
        Creates an unstructured grid around the given (x,y) points.
        """
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)

        if x.size != y.size:
            raise ValueError("x and y arrays must have the same size")

        # Make a copy of the points in a 2D array (useful for calls to geometry
        # routines, but takes extra memory space).
        xy_of_node = np.hstack((x.reshape((-1, 1)), y.reshape((-1, 1))))
        self._xy_of_node = sort_points_by_x_then_y(xy_of_node)

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
        self._find_perimeter_nodes_and_BC_set(self._xy_of_node)
        [self._cell_at_node, self._node_at_cell] = self._node_to_cell_connectivity(
            self.status_at_node, self.number_of_cells
        )

        # ACTIVE CELLS: Construct Voronoi diagram and calculate surface area of
        # each active cell.
        vor = Voronoi(self._xy_of_node)
        self.vor = vor
        self._area_of_cell = np.zeros(self.number_of_cells)
        for node in self._node_at_cell:
            xv = vor.vertices[vor.regions[vor.point_region[node]], 0]
            yv = vor.vertices[vor.regions[vor.point_region[node]], 1]
            self._area_of_cell[self.cell_at_node[node]] = simple_poly_area(xv, yv)

        # LINKS: Construct Delaunay triangulation and construct lists of link
        # "from" and "to" nodes.
        (
            node_at_link_tail,
            node_at_link_head,
            _,
            self._face_width,
        ) = self._create_links_and_faces_from_voronoi_diagram(vor)

        self._nodes_at_link = np.hstack(
            (node_at_link_tail.reshape((-1, 1)), node_at_link_head.reshape((-1, 1)))
        )

        # Sort them by midpoint coordinates
        self._sort_links_by_midpoint()

        # Optionally re-orient links so that they all point within upper-right
        # semicircle
        if reorient_links:
            self._reorient_links_upper_right()

        # NODES & LINKS: IDs and directions of links at each node
        self._create_links_and_link_dirs_at_node()

        # LINKS: set up link unit vectors and node unit-vector sums
        self._create_link_unit_vectors()

    @property
    def number_of_patches(self):
        """Number of patches.

        Returns the number of patches over the grid.

        LLCATS: PINF
        """
        try:
            return self._number_of_patches
        except AttributeError:
            self._create_patches_from_delaunay_diagram(self._xy_of_node, self.vor)
            return self._number_of_patches

    @property
    def nodes_at_patch(self):
        """Get the four nodes at the corners of each patch in a regular grid.

        LLCATS: PINF NINF CONN
        """
        try:
            return self._nodes_at_patch
        except AttributeError:
            self._create_patches_from_delaunay_diagram(self._xy_of_node, self.vor)
            return self._nodes_at_patch

    @property
    @return_readonly_id_array
    def patches_at_node(self):
        """
        Return a (nnodes, max_voronoi_polygon_sides) array of patches at nodes.

        The patches are returned in LL standard order (ccw from E), with any
        nonexistent patches recorded after the ids of existing faces.
        Nonexistent patches are ID'ed as -1.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> mg = HexModelGrid(3, 3)
        >>> mg.patches_at_node # doctest: +SKIP
        array([[ 0,  2, -1, -1, -1, -1],
               [ 1,  3,  0, -1, -1, -1],
               [ 4,  1, -1, -1, -1, -1],
               [ 5,  2, -1, -1, -1, -1],
               [ 6,  8,  5,  2,  0,  3],
               [ 7,  9,  6,  3,  1,  4],
               [ 7,  4, -1, -1, -1, -1],
               [ 5,  8, -1, -1, -1, -1],
               [ 8,  6,  9, -1, -1, -1],
               [ 9,  7, -1, -1, -1, -1]])

        LLCATS: NINF PINF CONN
        """
        try:
            return self._patches_at_node
        except AttributeError:
            self._create_patches_from_delaunay_diagram(self._xy_of_node, self.vor)
            return self._patches_at_node

    @property
    @return_readonly_id_array
    def links_at_patch(self):
        """Returns the links forming each patch.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> mg = HexModelGrid(3, 2)
        >>> mg.links_at_patch
        array([[ 3,  2,  0],
               [ 5,  1,  2],
               [ 6,  3,  4],
               [ 8,  7,  5],
               [10,  9,  6],
               [11,  8,  9]])

        LLCATS: LINF PINF CONN
        """
        try:
            return self._links_at_patch
        except AttributeError:
            self._create_patches_from_delaunay_diagram(self._xy_of_node, self.vor)
            return self._links_at_patch

    @property
    @return_readonly_id_array
    def patches_at_link(self):
        """Returns the patches adjoined to each link.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> mg = HexModelGrid(3, 2)
        >>> mg.patches_at_link
        array([[ 0, -1],
               [ 1, -1],
               [ 0,  1],
               [ 0,  2],
               [ 2, -1],
               [ 1,  3],
               [ 2,  4],
               [ 3, -1],
               [ 3,  5],
               [ 4,  5],
               [ 4, -1],
               [ 5, -1]])

        LLCATS: PINF LINF CONN
        """
        try:
            return self._patches_at_link
        except AttributeError:
            self._create_patches_from_delaunay_diagram(self._xy_of_node, self.vor)
            return self._patches_at_link

    def _find_perimeter_nodes_and_BC_set(self, pts):
        """
        Uses a convex hull to locate the perimeter nodes of the Voronoi grid,
        then sets them as fixed value boundary nodes.
        It then sets/updates the various relevant node lists held by the grid,
        and returns *node_status*, *core_nodes*, *boundary_nodes*.
        """

        # Calculate the convex hull for the set of points
        from scipy.spatial import ConvexHull

        hull = ConvexHull(pts, qhull_options="Qc")  # see below why we use 'Qt'

        # The ConvexHull object lists the edges that form the hull. We need to
        # get from this list of edges the unique set of nodes. To do this, we
        # first flatten the list of vertices that make up all the hull edges
        # ("simplices"), so it becomes a 1D array. With that, we can use the
        # set() function to turn the array into a set, which removes duplicate
        # vertices. Then we turn it back into an array, which now contains the
        # set of IDs for the nodes that make up the convex hull.
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
        convex_hull_nodes = np.array(list(set(hull.simplices.flatten())))
        coplanar_nodes = hull.coplanar[:, 0]
        boundary_nodes = as_id_array(
            np.concatenate((convex_hull_nodes, coplanar_nodes))
        )

        # Now we'll create the "node_status" array, which contains the code
        # indicating whether the node is interior and active (=0) or a
        # boundary (=1). This means that all perimeter (convex hull) nodes are
        # initially flagged as boundary code 1. An application might wish to
        # change this so that, for example, some boundaries are inactive.
        node_status = np.zeros(len(pts[:, 0]), dtype=np.uint8)
        node_status[boundary_nodes] = 1

        # It's also useful to have a list of interior nodes
        core_nodes = as_id_array(np.where(node_status == 0)[0])

        # save the arrays and update the properties
        self._node_status = node_status
        self._node_at_cell = core_nodes
        self._boundary_nodes = boundary_nodes

        self.status_at_node = node_status

        # Return the results
        return node_status, core_nodes, boundary_nodes

    def _create_cell_areas_array(self):
        """Set up an array of cell areas."""
        self._cell_areas = self.active_cell_areas
        return self._cell_areas

    @staticmethod
    def _node_to_cell_connectivity(node_status, ncells):
        """Set up node connectivity.

        Creates and returns the following arrays:

        *  For each node, the ID of the corresponding cell, or
           BAD_INDEX_VALUE if the node has no cell.
        *  For each cell, the ID of the corresponding node.

        Parameters
        ----------
        node_status : ndarray of ints
            1D array containing the boundary status code for each node.
        ncells : ndarray of ints
            Number of cells (must equal the number of occurrences of CORE_NODE
            in node_status).

        Examples
        --------
        >>> from landlab import VoronoiDelaunayGrid as vdg
        >>> import numpy as np
        >>> from landlab.grid import BAD_INDEX_VALUE
        >>> ns = np.array([1, 0, 0, 1, 0])  # 3 interior, 2 boundary nodes
        >>> [node_cell, cell_node] = vdg._node_to_cell_connectivity(ns, 3)
        >>> node_cell[1:3]
        array([0, 1])
        >>> node_cell[0] == BAD_INDEX_VALUE
        True
        >>> cell_node
        array([1, 2, 4])
        """
        assert ncells == np.count_nonzero(
            node_status == CORE_NODE
        ), "ncells must equal number of CORE_NODE values in node_status"

        cell = 0
        node_cell = np.ones(len(node_status), dtype=int) * BAD_INDEX_VALUE
        cell_node = np.zeros(ncells, dtype=int)
        for node in range(len(node_cell)):
            if node_status[node] == CORE_NODE:
                node_cell[node] = cell
                cell_node[cell] = node
                cell += 1

        return node_cell, cell_node

    @staticmethod
    def _create_links_from_triangulation(tri):
        """Create links from a Delaunay triangulation.

        From a Delaunay Triangulation of a set of points, contained in a
        scipy.spatial.Delaunay object "tri", creates and returns:

        *  a numpy array containing the ID of the "from" node for each link
        *  a numpy array containing the ID of the "to" node for each link
        *  the number of links in the triangulation

        Examples
        --------
        >>> from scipy.spatial import Delaunay
        >>> import numpy as np
        >>> from landlab.grid import VoronoiDelaunayGrid as vdg
        >>> pts = np.array([[ 0., 0.], [ 1., 0.], [ 1., 0.87],
        ...                 [-0.5, 0.87], [ 0.5, 0.87], [ 0., 1.73],
        ...                 [ 1., 1.73]])
        >>> dt = Delaunay(pts)
        >>> [myfrom,myto,nl] = vdg._create_links_from_triangulation(dt)
        >>> print myfrom, myto, nl # doctest: +SKIP
        [5 3 4 6 4 3 0 4 1 1 2 6] [3 4 5 5 6 0 4 1 0 2 4 2] 12
        """

        # Calculate how many links there will be and create the arrays.
        #
        # The number of links equals 3 times the number of triangles minus
        # half the number of shared links. Finding out the number of shared
        # links is easy: for every shared link, there is an entry in the
        # tri.neighbors array that is > -1 (indicating that the triangle has a
        # neighbor opposite a given vertex; in other words, two triangles are
        # sharing an edge).
        num_shared_links = np.count_nonzero(tri.neighbors > -1)
        num_links = 3 * tri.nsimplex - num_shared_links // 2
        link_fromnode = np.zeros(num_links, dtype=int)
        link_tonode = np.zeros(num_links, dtype=int)

        # Sweep through the list of triangles, assigning "from" and "to" nodes
        # to the list of links.
        #
        # The basic algorithm works as follows. For each triangle, we will add
        # its 3 edges as links. However, we have to make sure that each shared
        # edge is added only once. To do this, we keep track of whether or not
        # each triangle has been processed yet using a boolean array called
        # "tridone". When we look at a given triangle, we check each vertex in
        # turn. If there is no neighboring triangle opposite that vertex, then
        # we need to add the corresponding edge. If there is a neighboring
        # triangle but we haven't processed it yet, we also need to add the
        # edge. If neither condition is true, then this edge has already been
        # added, so we skip it.
        link_id = 0
        tridone = np.zeros(tri.nsimplex, dtype=bool)
        for t in range(tri.nsimplex):  # loop over triangles
            for i in range(0, 3):  # loop over vertices & neighbors
                if tri.neighbors[t, i] == -1 or not tridone[tri.neighbors[t, i]]:
                    link_fromnode[link_id] = tri.simplices[t, np.mod(i + 1, 3)]
                    link_tonode[link_id] = tri.simplices[t, np.mod(i + 2, 3)]
                    link_id += 1
            tridone[t] = True

        # save the results
        # self.node_at_link_tail = link_fromnode
        # self.node_at_link_head = link_tonode

        # Return the results
        return link_fromnode, link_tonode, num_links

    @staticmethod
    def _is_valid_voronoi_ridge(vor, n):

        SUSPICIOUSLY_BIG = 40000000.0
        return (
            vor.ridge_vertices[n][0] != -1
            and vor.ridge_vertices[n][1] != -1
            and np.amax(np.abs(vor.vertices[vor.ridge_vertices[n]])) < SUSPICIOUSLY_BIG
        )

    @staticmethod
    def _create_links_and_faces_from_voronoi_diagram(vor):
        """
        From a Voronoi diagram object created by scipy.spatial.Voronoi(),
        builds and returns:
        1. Arrays of link tail and head nodes
        2. Array of link IDs for each active link
        3. Array containing with of each face

        Parameters
        ----------
        vor : scipy.spatial.Voronoi
            Voronoi object initialized with the grid nodes.

        Returns
        -------
        out : tuple of ndarrays
            - link_fromnode = "from" node for each link (len=num_links)
            - link_tonode   = "to" node for each link (len=num_links)
            - active_links  = link ID for each active link
                                                    (len=num_active_links)
            - face_width    = width of each face (len=num_active_links

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.grid import VoronoiDelaunayGrid as vdg
        >>> pts = np.array([[0., 0.], [1., 0.], [-0.5, 0.87], [0.5, 0.87],
        ...                 [1.5, 0.87], [0., 1.73], [1., 1.73]])
        >>> from scipy.spatial import Voronoi
        >>> vor = Voronoi(pts)
        >>> [tn,hn,al,fw] = vdg._create_links_and_faces_from_voronoi_diagram(
        ...     vor)
        >>> tn
        array([0, 0, 0, 1, 1, 2, 3, 2, 3, 6, 6, 6])
        >>> hn
        array([1, 2, 3, 3, 4, 3, 4, 5, 5, 3, 4, 5])
        >>> al
        array([2, 3, 5, 6, 8, 9])
        >>> fw
        array([ 0.57669199,  0.57669199,  0.575973  ,  0.575973  ,  0.57836419,
                0.57836419])
        """
        # Each Voronoi "ridge" corresponds to a link. The Voronoi object has an
        # attribute ridge_points that contains the IDs of the nodes on either
        # side (including ridges that have one of their endpoints undefined).
        # So, we set the number of links equal to the number of ridges.
        num_links = len(vor.ridge_points)

        # Create the arrays for link from and to nodes
        link_fromnode = -np.ones(num_links, dtype=int)
        link_tonode = -np.ones(num_links, dtype=int)

        # Ridges along the perimeter of the grid will have one of their
        # endpoints undefined. The endpoints of each ridge are contained in
        # vor.ridge_vertices, and an undefined vertex is flagged with -1.
        # Ridges with both vertices defined correspond to faces and active
        # links, while ridges with an undefined vertex correspond to inactive
        # links. So, to find the number of active links, we subtract from the
        # total number of links the number of occurrences of an undefined
        # vertex.
        num_active_links = num_links - np.count_nonzero(
            np.array(vor.ridge_vertices) == -1
        )

        # Create arrays for active links and width of faces (which are Voronoi
        # ridges).
        active_links = -np.ones(num_active_links, dtype=int)
        face_width = -np.ones(num_active_links)

        # Find the order to sort by link midpoints
        link_midpoints = np.zeros((num_links, 2))
        for i in range(num_links):
            link_midpoints[i][:] = (
                vor.points[vor.ridge_points[i, 0]] + vor.points[vor.ridge_points[i, 1]]
            ) / 2.0
        ind = argsort_points_by_x_then_y(link_midpoints)

        # Loop through the list of ridges. For each ridge, there is a link, and
        # its "from" and "to" nodes are the associated "points". In addition,
        # if the ridge endpoints are defined, we have a face and an active
        # link, so we add them to our arrays as well.
        j = 0
        for i in range(num_links):
            link_fromnode[i] = vor.ridge_points[ind[i], 0]
            link_tonode[i] = vor.ridge_points[ind[i], 1]
            face_corner1 = vor.ridge_vertices[ind[i]][0]
            face_corner2 = vor.ridge_vertices[ind[i]][1]
            # means it's a valid face
            if VoronoiDelaunayGrid._is_valid_voronoi_ridge(vor, ind[i]):
                dx = vor.vertices[face_corner2, 0] - vor.vertices[face_corner1, 0]
                dy = vor.vertices[face_corner2, 1] - vor.vertices[face_corner1, 1]
                face_width[j] = np.sqrt(dx * dx + dy * dy)
                active_links[j] = i
                j += 1

        return link_fromnode, link_tonode, active_links, face_width

    def _reorient_links_upper_right(self):
        r"""Reorient links to all point within the upper-right semi-circle.

        Notes
        -----
        "Upper right semi-circle" means that the angle of the link with respect
        to the vertical (measured clockwise) falls between -45 and +135. More
        precisely, if :math:`\theta' is the angle,
        :math:`-45 \ge \theta < 135`.
        For example, the link could point up and left as much as -45, but not
        -46. It could point down and right as much as 134.9999, but not 135. It
        will never point down and left, or up-but-mostly-left, or
        right-but-mostly-down.

        Examples
        --------
        >>> from landlab.grid import HexModelGrid
        >>> hg = HexModelGrid(3, 2, 1., reorient_links=True)
        >>> hg.node_at_link_tail
        array([0, 0, 0, 1, 1, 2, 3, 2, 3, 3, 4, 5])
        >>> hg.node_at_link_head
        array([1, 2, 3, 3, 4, 3, 4, 5, 5, 6, 6, 6])
        """

        # Calculate the horizontal (dx) and vertical (dy) link offsets
        link_dx = (
            self.node_x[self.node_at_link_head] - self.node_x[self.node_at_link_tail]
        )
        link_dy = (
            self.node_y[self.node_at_link_head] - self.node_y[self.node_at_link_tail]
        )

        # Calculate the angle, clockwise, with respect to vertical, then rotate
        # by 45 degrees counter-clockwise (by adding pi/4)
        link_angle = np.arctan2(link_dx, link_dy) + np.pi / 4

        # The range of values should be -180 to +180 degrees (but in radians).
        # It won't be after the above operation, because angles that were
        # > 135 degrees will now have values > 180. To correct this, we
        # subtract 360 (i.e., 2 pi radians) from those that are > 180 (i.e.,
        # > pi radians).
        link_angle -= 2 * np.pi * (link_angle >= np.pi)

        # Find locations where the angle is negative; these are the ones we
        # want to flip
        (flip_locs,) = np.where(link_angle < 0.0)

        # If there are any flip locations, proceed to switch their fromnodes
        # and tonodes; otherwise, we're done
        self._nodes_at_link[flip_locs, :] = self._nodes_at_link[flip_locs, ::-1]

    def _create_patches_from_delaunay_diagram(self, pts, vor):
        """
        Uses a delaunay diagram drawn from the provided points to
        generate an array of patches and patch-node-link connectivity.
        Returns ...
        DEJH, 10/3/14, modified May 16.
        """
        from scipy.spatial import Delaunay
        from landlab.core.utils import anticlockwise_argsort_points_multiline
        from .cfuncs import create_patches_at_element, create_links_at_patch

        tri = Delaunay(pts)
        assert np.array_equal(tri.points, vor.points)
        nodata = -1
        self._nodes_at_patch = as_id_array(tri.simplices)
        # self._nodes_at_patch = np.empty_like(_nodes_at_patch)
        self._number_of_patches = tri.simplices.shape[0]
        # get the patches in order:
        patches_xy = np.empty((self._number_of_patches, 2), dtype=float)
        patches_xy[:, 0] = np.mean(self.node_x[self._nodes_at_patch], axis=1)
        patches_xy[:, 1] = np.mean(self.node_y[self._nodes_at_patch], axis=1)
        orderforsort = argsort_points_by_x_then_y(patches_xy)
        self._nodes_at_patch = self._nodes_at_patch[orderforsort, :]
        patches_xy = patches_xy[orderforsort, :]

        # perform a CCW sort without a line-by-line loop:
        patch_nodes_x = self.node_x[self._nodes_at_patch]
        patch_nodes_y = self.node_y[self._nodes_at_patch]
        anticlockwise_argsort_points_multiline(
            patch_nodes_x, patch_nodes_y, out=self._nodes_at_patch
        )

        # need to build a squared off, masked array of the patches_at_node
        # the max number of patches for a node in the grid is the max sides of
        # the side-iest voronoi region.
        max_dimension = len(max(vor.regions, key=len))

        self._patches_at_node = np.full(
            (self.number_of_nodes, max_dimension), nodata, dtype=int
        )

        self._nodes_at_patch = as_id_array(self._nodes_at_patch)
        self._patches_at_node = as_id_array(self._patches_at_node)

        create_patches_at_element(
            self._nodes_at_patch, self.number_of_nodes, self._patches_at_node
        )

        # build the patch-link connectivity:
        self._links_at_patch = np.empty((self._number_of_patches, 3), dtype=int)
        create_links_at_patch(
            self._nodes_at_patch,
            self._links_at_node,
            self._number_of_patches,
            self._links_at_patch,
        )
        patch_links_x = self.x_of_link[self._links_at_patch]
        patch_links_y = self.y_of_link[self._links_at_patch]
        anticlockwise_argsort_points_multiline(
            patch_links_x, patch_links_y, out=self._links_at_patch
        )

        self._patches_at_link = np.empty((self.number_of_links, 2), dtype=int)
        self._patches_at_link.fill(-1)
        create_patches_at_element(
            self._links_at_patch, self.number_of_links, self._patches_at_link
        )
        # a sort of the links will be performed here once we have corners

        self._patches_created = True

    def save(self, path, clobber=False):
        """Save a grid and fields.

        This method uses cPickle to save a Voronoi grid as a cPickle file.
        At the time of coding, this is the only convenient output format
        for Voronoi grids, but support for netCDF is likely coming.

        All fields will be saved, along with the grid.

        The recommended suffix for the save file is '.grid'. This will
        be added to your save if you don't include it.

        This method is equivalent to
        :py:func:`~landlab.io.native_landlab.save_grid`, and
        :py:func:`~landlab.io.native_landlab.load_grid` can be used to
        load these files.

        Caution: Pickling can be slow, and can produce very large files.
        Caution 2: Future updates to Landlab could potentially render old
        saves unloadable.

        Parameters
        ----------
        path : str
            Path to output file.
        clobber : bool (defaults to false)
            Set to true to allow overwriting

        Examples
        --------
        >>> from landlab import VoronoiDelaunayGrid
        >>> import numpy as np
        >>> import os
        >>> x = np.random.rand(20)
        >>> y = np.random.rand(20)
        >>> vmg = VoronoiDelaunayGrid(x,y)
        >>> vmg.save('./mytestsave.grid')
        >>> os.remove('mytestsave.grid') #to remove traces of this test

        LLCATS: GINF
        """
        import os
        from six.moves import cPickle

        if os.path.exists(path) and not clobber:
            raise ValueError("file exists")

        (base, ext) = os.path.splitext(path)
        if ext != ".grid":
            ext = ext + ".grid"
        path = base + ext

        with open(path, "wb") as fp:
            cPickle.dump(self, fp)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
