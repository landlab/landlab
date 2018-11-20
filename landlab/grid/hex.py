#! /usr/env/python
"""
Python implementation of HexModelGrid, a grid class used to create and manage
structured Voronoi-Delaunay grids for 2D numerical models.

Do NOT add new documentation here. Grid documentation is now built in a semi-
automated fashion. To modify the text seen on the web, edit the files
`docs/text_for_[gridfile].py.txt`.
"""
import warnings

import numpy
import six

from ..core.utils import as_id_array
from ..graph import DualHexGraph
from .base import ModelGrid, BAD_INDEX_VALUE


class HexModelGrid(DualHexGraph, ModelGrid):
    """A grid of hexagonal cells.

    This inherited class implements a regular 2D grid with hexagonal cells and
    triangular patches. It is a special type of VoronoiDelaunay grid in which
    the initial set of points is arranged in a triangular/hexagonal lattice.

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
    shape : string, optional
        Controls the shape of the bounding hull, i.e., are the nodes arranged
        in a hexagon, or a rectangle? Either 'hex' (default) or 'rect'.

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

    def __init__(
        self,
        base_num_rows=0,
        base_num_cols=0,
        dx=1.0,
        origin=(0., 0.),
        orientation="horizontal",
        shape="hex",
        reorient_links=True,
        **kwds
    ):
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
        node_layout = shape
        shape = (base_num_rows, base_num_cols)
        spacing = dx

        # if orientation.startswith('vert'):
        #     warnings.warn(
        #         "orientation vert is deprecated. Use vertical",
        #         category=DeprecationWarning,
        #     )
        #     orientation = 'vertical'
        # elif orientation.startswith('horiz'):
        #     warnings.warn(
        #         "orientation horiz is deprecated. Use horizontal",
        #         category=DeprecationWarning,
        #     )
        #     orientation = 'horizontal'

        DualHexGraph.__init__(
            self,
            shape,
            spacing=spacing,
            origin=origin,
            orientation=orientation,
            node_layout=node_layout,
        )
        ModelGrid.__init__(self, **kwds)

        self._node_status = numpy.full(self.number_of_nodes,
                                       self.BC_NODE_IS_CORE, dtype=numpy.uint8)
        self._node_status[self.perimeter_nodes] = self.BC_NODE_IS_FIXED_VALUE

    @classmethod
    def from_dict(cls, params):
        """
        LLCATS: GINF
        """
        shape = params["shape"]
        spacing = params.get("spacing", 1.)

        return cls(shape[0], shape[1], spacing)

    def _set_boundary_stat_at_rect_grid_ragged_edges(self, orientation, dx):
        """Assign boundary status to all edge nodes along the 'ragged' edges.

        Handle special case of boundary nodes in rectangular grid shape.
        One pair of edges will have the nodes staggered. By default, only the
        outer nodes will be assigned boundary status; we need the inner edge
        nodes on these "ragged" edges also to be flagged as boundary nodes.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> hg = HexModelGrid(3, 3, shape='rect', dx=2.0)
        >>> hg.status_at_node
        array([1, 1, 1, 1, 0, 1, 1, 1, 1], dtype=uint8)
        >>> hg = HexModelGrid(3, 3, shape='rect', orientation="vertical")
        >>> hg.status_at_node
        array([1, 1, 1, 1, 1, 0, 1, 1, 1], dtype=uint8)
        >>> hg = HexModelGrid(4, 4, shape='rect', orientation="vertcal")
        >>> hg.status_at_node
        array([1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1], dtype=uint8)
        >>> hg = HexModelGrid(3, 4, shape='rect')
        >>> hg.status_at_node
        array([1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1], dtype=uint8)
        """
        if orientation[0].lower() == "v":  # vert; top & bottom staggered
            bot_row = numpy.where(self.y_of_node <= 0.5 * dx)[0]
            self.status_at_node[bot_row] = self.status_at_node[0]
            top_row = numpy.where(self.y_of_node >= (self._nrows - 1) * dx)[0]
            self.status_at_node[top_row] = self.status_at_node[0]
        else:  # horizontal orientation; left & right staggered
            left_row = numpy.where(self.x_of_node <= 0.5 * dx)[0]
            self.status_at_node[left_row] = self.status_at_node[0]
            right_row = numpy.where(self.x_of_node >= (self._ncols - 1) * dx)[0]
            self.status_at_node[right_row] = self.status_at_node[0]

    def _create_cell_areas_array(self):
        r"""Create an array of surface areas of hexagonal cells.

        Creates and returns an array containing the surface areas of the
        hexagonal (Voronoi) cells.

        These cells are perfect hexagons in which the apothem is dx/2. The
        formula for area is:

        .. math::
            A = 3 dx^2 / 2 \sqrt{3} \approx 0.866 dx^2
        """
        self._area_of_cell = 0.8660254 * self._dx ** 2 + numpy.zeros(
            self.number_of_cells
        )
        return self._area_of_cell

    @staticmethod
    def _hex_points_with_horizontal_hex(num_rows, base_num_cols, dxh):
        """Create a set of points on a staggered grid.

        Creates and returns a set of (x,y) points in a staggered grid in which
        the points represent the centers of regular hexagonal cells, and the
        points could be connected to form equilateral triangles. The overall
        shape of the lattice is hexagonal, and one of the 3 axes is horizontal.

        Parameters
        ----------
        num_rows : int
            Number of rows in lattice
        base_num_cols : int
            Number of columns in the bottom and top rows (middle rows have
            more)
        dxh : float
            Horizontal and diagonal spacing between points

        Returns
        -------
        poinst : ndarray
            A 2D numpy array containing point (x,y) coordinates, and total
            number of points.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> points = HexModelGrid._hex_points_with_horizontal_hex(3, 2, 1.0)
        >>> len(points)
        7
        >>> points[1, :]
        array([ 1.,  0.])
        >>> points[:3, 0]
        array([ 0. ,  1. , -0.5])
        """
        dxv = dxh * numpy.sqrt(3.) / 2.
        half_dxh = dxh / 2.

        if numpy.mod(num_rows, 2) == 0:  # even number of rows
            npts = num_rows * base_num_cols + (num_rows * num_rows) // 4
        else:  # odd number of rows
            npts = num_rows * base_num_cols + ((num_rows - 1) // 2) * (
                (num_rows - 1) // 2
            )
        pts = numpy.zeros((npts, 2))
        middle_row = num_rows // 2
        extra_cols = 0
        xshift = 0.
        i = 0
        for r in range(num_rows):
            for c in range(base_num_cols + extra_cols):
                pts[i, 0] = c * dxh + xshift
                pts[i, 1] = r * dxv
                i += 1
            if r < middle_row:
                extra_cols += 1
            else:
                extra_cols -= 1
            xshift = -half_dxh * extra_cols

        return pts

    @staticmethod
    def _hex_points_with_horizontal_rect(num_rows, num_cols, dxh):
        """Create a set of points in a taggered grid.
        Creates and returns a set of (x,y) points in a staggered grid in which
        the points represent the centers of regular hexagonal cells, and the
        points could be connected to form equilateral triangles. The overall
        shape of the lattice is rectangular, and one of the 3 axes is
        horizontal.

        Parameters
        ----------
        num_rows : int
            Number of rows in lattice
        num_cols : int
            Number of columns in lattice
        dxh : float
            Horizontal and diagonal spacing between points

        Returns
        -------
        points : ndarray of shape `(n_points, 2)`
            A 2D numpy array containing point (x, y) coordinates, and total
            number of points.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> points = HexModelGrid._hex_points_with_horizontal_rect(3, 3, 1.0)
        >>> len(points)
        9
        >>> points[1, :]
        array([ 1.,  0.])
        >>> points[:3, 0]
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
                xshift = half_dxh * (r % 2)
                pts[i, 0] = c * dxh + xshift
                pts[i, 1] = r * dxv
                i += 1

        return pts

    @staticmethod
    def _hex_points_with_vertical_hex(base_num_rows, num_cols, dxv):
        """
        Creates and returns a set of (x,y) points in a staggered grid in which
        the points represent the centers of regular hexagonal cells, and the
        points could be connected to form equilateral triangles. The overall
        shape of the lattice is hexagonal.

        Parameters
        ----------
        base_num_rows : int
            Number of columns in the left and right columns (middle columns
            have more)
        num_cols : int
            Number of columns in lattice
        dxv : float
            Vertical and diagonal spacing between points

        Returns
        -------
        points : ndarray of shape `(n_points, 2)`
            2D numpy array containing point (x,y) coordinates, and total
            number of points.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> points = HexModelGrid._hex_points_with_vertical_hex(2, 3, 1.0)
        >>> len(points)
        7
        >>> points[1, :]
        array([ 0.,  1.])
        >>> points[:3, 1]
        array([ 0. ,  1. , -0.5])
        """
        dxh = dxv * numpy.sqrt(3.) / 2.
        half_dxv = dxv / 2.

        if numpy.mod(num_cols, 2) == 0:  # even number of columns
            npts = base_num_rows * num_cols + (num_cols * num_cols) // 4
        else:  # odd number of columns
            npts = base_num_rows * num_cols + ((num_cols - 1) // 2) * (
                (num_cols - 1) // 2
            )
        pts = numpy.zeros((npts, 2))
        middle_col = num_cols // 2
        extra_rows = 0
        yshift = 0.
        i = 0
        for c in range(num_cols):
            for r in range(base_num_rows + extra_rows):
                pts[i, 1] = r * dxv + yshift
                pts[i, 0] = c * dxh
                i += 1
            if c < middle_col:
                extra_rows += 1
            else:
                extra_rows -= 1
            yshift = -half_dxv * extra_rows

        return pts

    @staticmethod
    def _hex_points_with_vertical_rect(num_rows, num_cols, dxv):
        """
        Creates and returns a set of (x,y) points in a staggered grid in which
        the points represent the centers of regular hexagonal cells, and the
        points could be connected to form equilateral triangles. The overall
        shape of the lattice is rectangular.

        Parameters
        ----------
        num_rows : int
            Number of columns in lattice
        num_cols : int
            Number of columns in lattice
        dxv : float
            Vertical and diagonal spacing between points

        Returns
        -------
        points : ndarray of shape `(n_points, 2)`
            2D numpy array containing point (x,y) coordinates, and total
            number of points.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> points = HexModelGrid._hex_points_with_vertical_rect(3, 3, 1.0)
        >>> len(points)
        9
        >>> points[1, :]
        array([ 0.,  1.])
        >>> points[:3, 1]
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
                yshift = half_dxv * (c % 2)
                pts[i, 1] = r * dxv + yshift
                pts[i, 0] = c * dxh
                i += 1

        return pts

    @property
    def number_of_node_columns(self):
        """Number of node columns hex grid.

        Number of node columns in a rectangular-shaped and/or
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

        LLCATS: GINF NINF
        """
        return self.shape[1]

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

        LLCATS: GINF NINF
        """
        return self._shape[0]

    def node_row_and_column(self, node_id):
        """Row and column from node ID, FOR VERT RECT CONFIGURATION ONLY.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> grid = HexModelGrid(3, 4, shape='rect', orientation='vertical')
        >>> grid.node_row_and_column(5)
        (1, 2)
        >>> grid = HexModelGrid(3, 5, shape='rect', orientation='vertical')
        >>> grid.node_row_and_column(13)
        (2, 1)
        """
        assert self.orientation[0] == "v", "grid orientation must be vertical"
        try:
            (nr, nc) = self._shape
        except AttributeError:
            raise AttributeError(
                "Only rectangular Hex grids have defined rows and columns."
            )

        row = node_id // nc
        n_mod_nc = node_id % nc
        half_nc = (nc + 1) // 2
        col = 2 * (n_mod_nc % half_nc) + n_mod_nc // half_nc
        return (row, col)

    def _configure_hexplot(self, data, data_label=None, color_map=None):
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

        # color
        if color_map is None:
            color_map = matplotlib.cm.jet

        # geometry
        apothem = self.spacing / 2.0
        # distance from node to each hexagon cell vertex
        radius = 2.0 * apothem / sqrt(3.0)

        # offsets from node x,y position
        offsets = zeros((6, 2))
        poly_verts = zeros((6, 2))

        # Figure out whether the orientation is horizontal or vertical
        if self.orientation[0] == "h":  # horizontal
            offsets[:, 0] = array([0., apothem, apothem, 0., -apothem, -apothem])
            offsets[:, 1] = array(
                [
                    radius,
                    radius / 2.0,
                    -radius / 2.0,
                    -radius,
                    -radius / 2.0,
                    radius / 2.0,
                ]
            )
        else:  # vertical
            offsets[:, 0] = array(
                [
                    radius / 2.0,
                    radius,
                    radius / 2.0,
                    -radius / 2.0,
                    -radius,
                    -radius / 2.0,
                ]
            )
            offsets[:, 1] = array([apothem, 0., -apothem, -apothem, 0., apothem])

        patches = []
        for i in range(self.number_of_nodes):
            poly_verts[:, 0] = self.node_x[i] + offsets[:, 0]
            poly_verts[:, 1] = self.node_y[i] + offsets[:, 1]
            p = Polygon(poly_verts, True)
            patches.append(p)

        self._hexplot_pc = PatchCollection(
            patches, cmap=color_map, edgecolor="none", linewidth=0.0
        )

        self._hexplot_configured = True

    def hexplot(self, data, data_label=None, color_map=None):
        """Create a plot of the grid elements.

        Creates a plot of the grid and one node-data field, showing hexagonal
        cells colored by values in the field.

        Parameters
        ----------
        data : str or node array (1d numpy array with number_of_nodes entries)
            Data field to be colored.
        data_label : str, optional
            Label for colorbar.
        color_map : matplotlib colormap object, None
            Color map to apply (defaults to "jet")

        See also
        --------
        plot.imshow_grid
            Another Landlab function capable of producing hexplots, with a
            fuller-featured set of options.

        LLCATS: GINF
        """
        from numpy import array, amin, amax
        import matplotlib.pyplot as plt
        import copy

        try:
            self._hexplot_configured
        except AttributeError:
            self._configure_hexplot(data, data_label, color_map)
        else:
            if self._hexplot_pc.cmap != color_map:
                self._configure_hexplot(data, data_label, color_map)

        # Handle *data*: if it's a numpy array, then we consider it the
        # data to be plotted. If it's a string, we consider it the name of the
        # node-field to plot, and we fetch it.
        if type(data) is str:
            data_label = data
            data = self.at_node[data]

        ax = plt.gca()
        self._hexplot_pc.set_array(array(data))
        copy_of_pc = copy.copy(self._hexplot_pc)
        ax.add_collection(copy_of_pc)
        plt.xlim([amin(self.node_x) - self.spacing, amax(self.node_x) + self.spacing])
        plt.ylim([amin(self.node_y) - self.spacing, amax(self.node_y) + self.spacing])

        return ax

    def node_has_boundary_neighbor(self, ids):
        """Check if HexModelGrid nodes have neighbors that are boundary nodes.

        Parameters
        ----------
        mg : HexModelGrid
            Source grid
        node_id : int
            ID of node to test.

        Returns
        -------
        boolean
            ``True`` if node has a neighbor with a boundary ID,
            ``False`` otherwise.


        Checks to see if one of the eight neighbor nodes of node(s) with
        *id* has a boundary node.  Returns True if a node has a boundary node,
        False if all neighbors are interior.

                0,  1,  2,  3,
              4,  5,  6,  7,  8,
            9, 10,  11, 12, 13, 14,
              15, 16, 17, 18, 19,
                20, 21, 22, 23

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> hmg = HexModelGrid(5, 4)
        >>> hmg.node_has_boundary_neighbor(6)
        True
        >>> hmg.node_has_boundary_neighbor(12)
        False
        >>> hmg.node_has_boundary_neighbor([12, 0])
        [False, True]

        LLCATS: NINF CONN BC
        """
        ans = []
        for i in numpy.atleast_1d(numpy.asarray(ids)):
            neighbors = self.adjacent_nodes_at_node[i]
            real_neighbors = neighbors[neighbors != BAD_INDEX_VALUE]
            if real_neighbors.size == 0:
                ans.append(True)
            else:
                neighbor_status = self.status_at_node[real_neighbors].astype(bool)
                if numpy.any(neighbor_status != self.BC_NODE_IS_CORE):
                    ans.append(True)
                else:
                    ans.append(False)

        if len(ans) == 1:
            return ans[0]
        else:
            return ans

    def set_watershed_boundary_condition_outlet_id(
        self, outlet_id, node_data, nodata_value=-9999.
    ):
        """Set the boundary conditions for a watershed on a HexModelGrid.

        All nodes with nodata_value are set to CLOSED_BOUNDARY (4).
        All nodes with data values are set to CORE_NODES (0), with the
        exception that the outlet node is set to a FIXED_VALUE_BOUNDARY (1).

        Note that the outer ring of the HexModelGrid is set to CLOSED_BOUNDARY, even
        if there are nodes that have values.  The only exception to this would
        be if the outlet node is on the boundary, which is acceptable.

        Assumes that the id of the outlet is already known.

        This assumes that the grid has a single watershed.  If this is not
        the case this will not work.

        Parameters
        ----------
        outlet_id : integer
            id of the outlet node
        node_data : ndarray
            Data values.
        nodata_value : float, optional
            Value that indicates an invalid value.

        Examples
        --------
        The example will use a HexModelGrid with node data values
        as illustrated:

                1. ,  2. ,  3. ,  4. ,
            0.5,  1.5,  2.5,  3.5,  4.5,
          0. ,  1. ,  2. ,  3. ,  4. ,  5.,
            0.5,  1.5,  2.5,  3.5,  4.5,
                1. ,  2. ,  3. ,  4.

        >>> from landlab import HexModelGrid
        >>> hmg = HexModelGrid(5, 4)
        >>> z = hmg.add_zeros('node', 'topographic__elevation')
        >>> z += hmg.x_of_node + 1.0

        >>> hmg.status_at_node
        array([1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1,
           1], dtype=uint8)

        >>> outlet = hmg.set_watershed_boundary_condition_outlet_id(
        ...          9, z, -9999.)
        >>> hmg.status_at_node
        array([4, 4, 4, 4, 4, 0, 0, 0, 4, 1, 0, 0, 0, 0, 4, 4, 0, 0, 0, 4, 4, 4, 4,
           4], dtype=uint8)

        LLCATS: BC
        """
        # make ring of no data nodes
        self.status_at_node[self.boundary_nodes] = self.BC_NODE_IS_CLOSED
        # set no data nodes to inactive boundaries
        self.set_nodata_nodes_to_closed(node_data, nodata_value)
        # set the boundary condition (fixed value) at the outlet_node
        self.status_at_node[outlet_id] = self.BC_NODE_IS_FIXED_VALUE

    def set_watershed_boundary_condition(
        self, node_data, nodata_value=-9999., return_outlet_id=False
    ):
        """
        Finds the node adjacent to a boundary node with the smallest value.
        This node is set as the outlet.  The outlet node must have a data
        value.  Can return the outlet id as a one element numpy array if
        return_outlet_id is set to True.

        All nodes with nodata_value are set to CLOSED_BOUNDARY
        (grid.status_at_node == 4). All nodes with data values are set to
        CORE_NODES (grid.status_at_node == 0), with the exception that the
        outlet node is set to a FIXED_VALUE_BOUNDARY (grid.status_at_node == 1).

        Note that the outer ring (perimeter) of the grid is set to
        CLOSED_BOUNDARY, even if there are nodes that have values. The only
        exception to this would be if the outlet node is on the perimeter, which
        is acceptable.

        This routine assumes that all of the nodata_values are on the outside of
        the data values. In other words, there are no islands of nodata_values
        surrounded by nodes with data.

        This also assumes that the grid has a single watershed (that is a single
        outlet node).

        Parameters
        ----------
        node_data : ndarray
            Data values.
        nodata_value : float, optional
            Value that indicates an invalid value.
        return_outlet_id : boolean, optional
            Indicates whether or not to return the id of the found outlet

        Examples
        --------
        The example will use a HexModelGrid with node data values
        as illustrated:

                1. ,  2. ,  3. ,  4. ,
            0.5,  1.5,  2.5,  3.5,  4.5,
          0. ,  1. ,  2. ,  3. ,  4. ,  5.,
            0.5,  1.5,  2.5,  3.5,  4.5,
                1. ,  2. ,  3. ,  4.

        >>> from landlab import HexModelGrid
        >>> hmg = HexModelGrid(5, 4)
        >>> z = hmg.add_zeros('node', 'topographic__elevation')
        >>> z += hmg.x_of_node + 1.0
        >>> out_id = hmg.set_watershed_boundary_condition(z, -9999.,
        ...                                               True)
        >>> out_id
        array([9])
        >>> hmg.status_at_node
        array([4, 4, 4, 4, 4, 0, 0, 0, 4, 1, 0, 0, 0, 0, 4, 4, 0, 0, 0, 4, 4, 4, 4,
           4], dtype=uint8)

        LLCATS: BC
        """
        # make ring of no data nodes
        self.status_at_node[self.boundary_nodes] = self.BC_NODE_IS_CLOSED

        # set no data nodes to inactive boundaries
        self.set_nodata_nodes_to_closed(node_data, nodata_value)

        # locs is a list that contains locations where
        # node data is not equal to the nodata value
        locs = numpy.where(node_data != nodata_value)
        if len(locs) < 1:
            raise ValueError("All data values are no_data values")

        # now find minimum of the data values
        min_val = numpy.min(node_data[locs])

        # now find where minimum values are
        min_locs = numpy.where(node_data == min_val)[0]

        # check all the locations with the minimum value to see if one
        not_found = True
        while not_found:
            # now check the min locations to see if any are next to
            # a boundary node
            local_not_found = True
            next_to_boundary = []

            # check all nodes rather than selecting the first node that meets
            # the criteria
            for i in range(len(min_locs)):
                next_to_boundary.append(self.node_has_boundary_neighbor(min_locs[i]))

            # if any of those nodes were adjacent to the boundary, check
            # that  there is only one. If only one, set as outlet loc, else,
            # raise a value error
            if any(next_to_boundary):
                local_not_found = False
                if sum(next_to_boundary) > 1:
                    potential_locs = min_locs[
                        numpy.where(numpy.asarray(next_to_boundary))[0]
                    ]
                    raise ValueError(
                        (
                            "Grid has two potential outlet nodes."
                            "They have the following node IDs: \n"
                            + str(potential_locs)
                            + "\nUse the method set_watershed_boundary_condition_outlet_id "
                            "to explicitly select one of these "
                            "IDs as the outlet node."
                        )
                    )
                else:
                    outlet_loc = min_locs[numpy.where(next_to_boundary)[0][0]]

            # checked all of the min vals, (so done with inner while)
            # and none of the min values were outlet candidates
            if local_not_found:
                # need to find the next largest minimum value
                # first find the locations of all values greater
                # than the old minimum
                # not done with outer while
                locs = numpy.where((node_data > min_val) & (node_data != nodata_value))
                # now find new minimum of these values
                min_val = numpy.min(node_data[locs])
                min_locs = numpy.where(node_data == min_val)[0]
            else:
                # if locally found, it is also globally found
                # so done with outer while
                not_found = False

        # set outlet boundary condition
        self.status_at_node[outlet_loc] = self.BC_NODE_IS_FIXED_VALUE

        if return_outlet_id:
            return as_id_array(numpy.array([outlet_loc]))


def from_dict(param_dict):
    """
    Create a HexModelGrid from the dictionary-like object, *param_dict*.
    Required keys of the dictionary are NUM_ROWS, NUM_COLS. Raises a KeyError
    if either of these are missing.  If GRID_SPACING is given, use it as the
    HexModelGrid *dx* parameter, otherwise default to unit spacing.
    """
    # Read and create a basic HexModelGrid
    try:
        n_rows = int(param_dict["NUM_ROWS"])
        n_cols = int(param_dict["NUM_COLS"])
        dx = float(param_dict.get("GRID_SPACING", 1.))
    except KeyError:
        raise
    except ValueError:
        raise
    else:
        hg = HexModelGrid(n_rows, n_cols, dx)

    return hg
