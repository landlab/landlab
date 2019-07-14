#! /usr/env/python
"""
Python implementation of HexModelGrid, a grid class used to create and manage
structured Voronoi-Delaunay grids for 2D numerical models.

Do NOT add new documentation here. Grid documentation is now built in a semi-
automated fashion. To modify the text seen on the web, edit the files
`docs/text_for_[gridfile].py.txt`.
"""
from __future__ import absolute_import

from warnings import warn

import numpy
import six

from landlab.grid.voronoi import VoronoiDelaunayGrid

from ..core.utils import as_id_array
from .base import CLOSED_BOUNDARY, CORE_NODE, FIXED_VALUE_BOUNDARY


class HexModelGrid(VoronoiDelaunayGrid):
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
    xy_of_lower_left : tuple, optional
        Minimum x-of-node and y-of-node values. Depending on the grid
        no node may be present at this coordinate. Default is (0., 0.).
    xy_of_reference : tuple, optional
        Coordinate value in projected space of the reference point,
        `xy_of_lower_left`. Default is (0., 0.)
    orientation : string, optional
        One of the 3 cardinal directions in the grid, either 'horizontal'
        (default) or 'vertical'
    node_layout : string, optional
        Controls the shape of the bounding hull, i.e., are the nodes arranged
        in a hexagon, or a rectangle? Either 'hex' (default) or 'rect'.
    shape : tuple of 2 int
        Alternative way to specify (base_num_rows, base_num_cols)

    Returns
    -------
    HexModelGrid
        A newly-created grid.

    Examples
    --------
    Create a hex grid with 2 rows of nodes. The first and third rows will
    have 2 nodes, and the second nodes.

    >>> from landlab import HexModelGrid
    >>> hmg = HexModelGrid(shape=(3, 2), dx=1.0)
    >>> hmg.number_of_nodes
    7
    """

    def __init__(
        self,
        base_num_rows=0,
        base_num_cols=0,
        dx=1.0,
        xy_of_lower_left=(0.0, 0.0),
        orientation="horizontal",
        node_layout="hex",
        reorient_links=True,
        shape=None,
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
        xy_of_lower_left : tuple, optional
            Minimum x-of-node and y-of-node values. Depending on the grid
            no node may be present at this coordinate. Default is (0., 0.).
        xy_of_reference : tuple, optional
            Coordinate value in projected space of the reference point,
            `xy_of_lower_left`. Default is (0., 0.)
        orientation : string, optional
            One of the 3 cardinal directions in the grid, either 'horizontal'
            (default) or 'vertical'
        shape : string
            Either 'hex' (default) or 'rect'
        reorient_links : bool, optional
            Whether or not to re-orient all links to point between -45 deg
            and +135 deg clockwise from "north" (i.e., along y axis). default
            is True.

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

        # Pre-LL2.0: handle shape and node_layout keyword args
        # (note: in LL2.0, shape should refer ONLY to rows and columns)
        if type(shape) is str:  # "old" (LL < 2.0) usage
            node_layout = shape
        elif type(shape) is tuple:
            (base_num_rows, base_num_cols) = shape

        xy_of_lower_left = tuple(xy_of_lower_left)
        # Set number of nodes, and initialize if caller has given dimensions
        if base_num_rows * base_num_cols > 0:
            self._initialize(
                base_num_rows,
                base_num_cols,
                dx,
                xy_of_lower_left,
                orientation,
                node_layout,
                reorient_links,
                shape,
            )
        self._xy_of_lower_left = xy_of_lower_left
        super(HexModelGrid, self).__init__(**kwds)

    @property
    def xy_of_lower_left(self):
        """Return (x, y) of the reference point."""
        return self._xy_of_lower_left

    @xy_of_lower_left.setter
    def xy_of_lower_left(self, xy_of_lower_left):
        """Set a new value for the xy_of_lower_left."""
        dx = self.x_of_node[0] - xy_of_lower_left[0]
        dy = self.y_of_node[0] - xy_of_lower_left[1]
        self._xy_of_node -= (dx, dy)
        self._xy_of_lower_left = xy_of_lower_left

    def _initialize(
        self,
        base_num_rows,
        base_num_cols,
        dx,
        xy_of_lower_left,
        orientation,
        node_layout,
        reorient_links=True,
        shape=None,
    ):
        r"""Set up a hexagonal grid.

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
        xy_of_lower_left : tuple, optional
            (x, y) coordinates of the xy_of_lower_left.
            Default in ``__init__`` is (0., 0.)
        orientation : string
            Either 'horizontal' (default in ``__init__``) or 'vertical'
        shape : string
            Either 'hex' (default in ``__init__``) or 'rect'
        reorient_links : bool, optional
            Whether or not to re-orient all links to point between -45 deg
            and +135 deg clockwise from "north" (i.e., along y axis). default
            is True.

        Returns
        -------
        (none)

        Creates/modifies
        ----------------
        Creates and initializes and self._dx

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
        the number of rows (if 'vertical') will vary across the grid.
        Therefore, for hex-shaped grids we record only self._nrows for
        'horizontal' grids, and only self._ncols for 'vertical' grids.
        """
        if self._DEBUG_TRACK_METHODS:
            six.print_(
                "HexModelGrid._initialize("
                + str(base_num_rows)
                + ", "
                + str(base_num_cols)
                + ", "
                + str(dx)
                + ")"
            )

        # Pre-LL2.0: handle shape and node_layout keyword args
        # (note: in LL2.0, shape should refer ONLY to rows and columns)
        if type(shape) is str:  # "old" (LL < 2.0) usage
            node_layout = shape
            message = (
                "Use node_layout to specify 'rect' or 'hex' layout "
                + "(usage will be enforced in Landlab 2.0+)."
            )
            warn(message=message, category=DeprecationWarning)
        elif type(shape) is tuple:
            (base_num_rows, base_num_cols) = shape

        # Make sure the parameter *orientation* is correct
        assert (
            orientation[0].lower() == "h" or orientation[0].lower() == "v"
        ), 'orientation must be either "horizontal" (default) or "vertical"'

        # Make sure the parameter *node_layout* is correct
        assert (
            node_layout[0].lower() == "h" or node_layout[0].lower() == "r"
        ), 'node_layout must be either "hex" (default) or "rect"'

        # Create a set of hexagonally arranged points. These will be our nodes.
        if orientation[0].lower() == "h" and node_layout[0].lower() == "h":
            pts = HexModelGrid._hex_points_with_horizontal_hex(
                base_num_rows, base_num_cols, dx, xy_of_lower_left
            )
            self.orientation = "horizontal"
            self._nrows = base_num_rows
        elif orientation[0].lower() == "h" and node_layout[0].lower() == "r":
            pts = HexModelGrid._hex_points_with_horizontal_rect(
                base_num_rows, base_num_cols, dx, xy_of_lower_left
            )
            self.orientation = "horizontal"
            self._nrows = base_num_rows
            self._ncols = base_num_cols
            self._shape = (self._nrows, self._ncols)
            self._nodes = numpy.arange((self._nrows * self._ncols), dtype=int).reshape(
                self._shape
            )
        elif orientation[0].lower() == "v" and node_layout[0].lower() == "h":
            pts = HexModelGrid._hex_points_with_vertical_hex(
                base_num_rows, base_num_cols, dx, xy_of_lower_left
            )
            self.orientation = "vertical"
            self._ncols = base_num_cols
        else:
            pts = HexModelGrid._hex_points_with_vertical_rect(
                base_num_rows, base_num_cols, dx, xy_of_lower_left
            )
            self.orientation = "vertical"
            self._nrows = base_num_rows
            self._ncols = base_num_cols
            self._shape = (self._nrows, self._ncols)
            self._nodes = numpy.arange((self._nrows * self._ncols), dtype=int).reshape(
                self._shape
            )
            for col in range(self._ncols):
                base_node = (col // 2) + (col % 2) * ((self._ncols + 1) // 2)
                self._nodes[:, col] = numpy.arange(
                    base_node, self._nrows * self._ncols, self._ncols
                )

        # Call the VoronoiDelaunayGrid constructor to triangulate/Voronoi
        # the nodes into a grid.
        super(HexModelGrid, self)._initialize(pts[:, 0], pts[:, 1], reorient_links)

        # Handle special case of boundary nodes in rectangular node layout.
        # One pair of edges will have the nodes staggered. By default, only the
        # outer nodes will be assigned boundary status; we need the inner edge
        # nodes on these "ragged" edges also to be flagged as boundary nodes.
        if node_layout[0].lower() == "r":
            self._set_boundary_stat_at_rect_grid_ragged_edges(orientation, dx)

        # Remember grid spacing
        self._dx = dx

    def _set_boundary_stat_at_rect_grid_ragged_edges(self, orientation, dx):
        """Assign boundary status to all edge nodes along the 'ragged' edges.

        Handle special case of boundary nodes in rectangular node layout.
        One pair of edges will have the nodes staggered. By default, only the
        outer nodes will be assigned boundary status; we need the inner edge
        nodes on these "ragged" edges also to be flagged as boundary nodes.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> hg = HexModelGrid(3, 3, node_layout='rect', dx=2.0)
        >>> hg.status_at_node
        array([1, 1, 1, 1, 0, 1, 1, 1, 1], dtype=uint8)
        >>> hg = HexModelGrid(3, 3, node_layout='rect', orientation='vert')
        >>> hg.status_at_node
        array([1, 1, 1, 1, 1, 0, 1, 1, 1], dtype=uint8)
        >>> hg = HexModelGrid(4, 4, node_layout='rect', orientation='vert')
        >>> hg.status_at_node
        array([1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1], dtype=uint8)
        >>> hg.boundary_nodes
        array([ 0,  1,  2,  3,  4,  7,  8, 11, 12, 13, 14, 15])
        >>> hg = HexModelGrid(3, 4, node_layout='rect')
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
        self._boundary_nodes = numpy.where(self.status_at_node != CORE_NODE)[0]

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
    def _shift_to_lower_left(pts, xy_of_lower_left):
        xshift = xy_of_lower_left[0] - numpy.min(pts[:, 0])
        yshift = xy_of_lower_left[1] - numpy.min(pts[:, 1])
        pts[:, 0] += xshift
        pts[:, 1] += yshift
        return pts

    @staticmethod
    def _hex_points_with_horizontal_hex(num_rows, base_num_cols, dxh, xy_of_lower_left):
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
        xy_of_lower_left : tuple
            (x, y) coordinates of the xy_of_lower_left. Default is (0., 0.)
        Returns
        -------
        poinst : ndarray
            A 2D numpy array containing point (x,y) coordinates, and total
            number of points.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> points = HexModelGrid._hex_points_with_horizontal_hex(3, 2,
        ...                                                       1.0,
        ...                                                       (0., 0.))
        >>> len(points)
        7
        >>> points[1, :]
        array([ 1.5,  0. ])
        >>> points[:3, 0]
        array([ 0.5,  1.5,  0. ])
        """
        dxv = dxh * numpy.sqrt(3.0) / 2.0
        half_dxh = dxh / 2.0

        if numpy.mod(num_rows, 2) == 0:  # even number of rows
            npts = num_rows * base_num_cols + (num_rows * num_rows) // 4
        else:  # odd number of rows
            npts = num_rows * base_num_cols + ((num_rows - 1) // 2) * (
                (num_rows - 1) // 2
            )
        pts = numpy.zeros((npts, 2))
        middle_row = num_rows // 2
        extra_cols = 0

        xshift = 0
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

        # return pts
        return HexModelGrid._shift_to_lower_left(pts, xy_of_lower_left)

    @staticmethod
    def _hex_points_with_horizontal_rect(num_rows, num_cols, dxh, xy_of_lower_left):
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
        xy_of_lower_left : tuple
            (x, y) coordinates of the xy_of_lower_left. Default is (0., 0.)

        Returns
        -------
        points : ndarray of shape `(n_points, 2)`
            A 2D numpy array containing point (x, y) coordinates, and total
            number of points.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> points = HexModelGrid._hex_points_with_horizontal_rect(3, 3,
        ...                                                        1.0,
        ...                                                        (0., 0.))
        >>> len(points)
        9
        >>> points[1, :]
        array([ 1.,  0.])
        >>> points[:3, 0]
        array([ 0.,  1.,  2.])
        """
        dxv = dxh * numpy.sqrt(3.0) / 2.0
        half_dxh = dxh / 2.0

        npts = num_rows * num_cols
        pts = numpy.zeros((npts, 2))

        i = 0
        for r in range(num_rows):
            for c in range(num_cols):
                xshift = half_dxh * (r % 2)
                pts[i, 0] = c * dxh + xshift
                pts[i, 1] = r * dxv
                i += 1

        return HexModelGrid._shift_to_lower_left(pts, xy_of_lower_left)

    @staticmethod
    def _hex_points_with_vertical_hex(base_num_rows, num_cols, dxv, xy_of_lower_left):
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
        xy_of_lower_left : tuple
            (x, y) coordinates of the xy_of_lower_left. Default is (0., 0.)

        Returns
        -------
        points : ndarray of shape `(n_points, 2)`
            2D numpy array containing point (x,y) coordinates, and total
            number of points.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> points = HexModelGrid._hex_points_with_vertical_hex(2, 3,
        ...                                                     1.0,
        ...                                                     (0., 0.))
        >>> len(points)
        7
        >>> points[2, :]
        array([ 0.8660254,  0.       ])
        >>> points[:3, 1]
        array([ 0.5,  1.5,  0. ])
        """
        dxh = dxv * numpy.sqrt(3.0) / 2.0
        half_dxv = dxv / 2.0

        if numpy.mod(num_cols, 2) == 0:  # even number of columns
            npts = base_num_rows * num_cols + (num_cols * num_cols) // 4
        else:  # odd number of columns
            npts = base_num_rows * num_cols + ((num_cols - 1) // 2) * (
                (num_cols - 1) // 2
            )
        pts = numpy.zeros((npts, 2))
        middle_col = num_cols // 2
        extra_rows = 0

        yshift = 0
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

        return HexModelGrid._shift_to_lower_left(pts, xy_of_lower_left)

    @staticmethod
    def _hex_points_with_vertical_rect(num_rows, num_cols, dxv, xy_of_lower_left):
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
        xy_of_lower_left : tuple
            (x, y) coordinates of the xy_of_lower_left. Default is (0., 0.).

        Returns
        -------
        points : ndarray of shape `(n_points, 2)`
            2D numpy array containing point (x,y) coordinates, and total
            number of points.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> points = HexModelGrid._hex_points_with_vertical_rect(3, 3,
        ...                                                      1.0,
        ...                                                      (0., 0.))
        >>> len(points)
        9
        >>> points[1, :]
        array([ 0.,  1.])
        >>> points[:3, 1]
        array([ 0.,  1.,  2.])
        """
        dxh = dxv * numpy.sqrt(3.0) / 2.0
        half_dxv = dxv / 2.0

        npts = num_rows * num_cols
        pts = numpy.zeros((npts, 2))

        i = 0
        for c in range(num_cols):
            for r in range(num_rows):
                yshift = half_dxv * (c % 2)
                pts[i, 1] = r * dxv + yshift
                pts[i, 0] = c * dxh
                i += 1

        return HexModelGrid._shift_to_lower_left(pts, xy_of_lower_left)

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
        >>> grid = HexModelGrid(5, 5, node_layout='rect')
        >>> grid.number_of_node_columns
        5

        LLCATS: GINF NINF
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
        >>> grid = HexModelGrid(5, 5, node_layout='rect')
        >>> grid.number_of_node_rows
        5

        LLCATS: GINF NINF
        """
        return self._nrows

    @property
    def nodes_at_left_edge(self):
        """Get nodes along the left edge of a grid, if grid is rectangular.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import HexModelGrid
        >>> grid = HexModelGrid(3, 4, node_layout='rect')
        >>> grid.nodes_at_left_edge
        array([0, 4, 8])

        LLCATS: NINF BC SUBSET
        """
        try:
            return self._nodes[:, 0]
        except AttributeError:
            raise AttributeError("Only rectangular Hex grids have defined edges.")

    @property
    def nodes_at_right_edge(self):
        """Get nodes along the right edge of a grid, if grid is rectangular.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import HexModelGrid
        >>> grid = HexModelGrid(3, 4, node_layout='rect')
        >>> grid.nodes_at_right_edge
        array([ 3,  7, 11])

        LLCATS: NINF BC SUBSET
        """
        try:
            return self._nodes[:, -1]
        except AttributeError:
            raise AttributeError("Only rectangular Hex grids have defined edges.")

    @property
    def nodes_at_top_edge(self):
        """Get nodes along the top edge of a grid, if grid is rectangular.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import HexModelGrid
        >>> grid = HexModelGrid(3, 4, node_layout='rect')
        >>> grid.nodes_at_top_edge
        array([ 8,  9, 10, 11])

        LLCATS: NINF BC SUBSET
        """
        try:
            return self._nodes[-1, :]
        except AttributeError:
            raise AttributeError("Only rectangular Hex grids have defined edges.")

    @property
    def nodes_at_bottom_edge(self):
        """Get nodes along the bottom edge of a grid, if grid is rectangular.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import HexModelGrid
        >>> grid = HexModelGrid(3, 4, node_layout='rect')
        >>> grid.nodes_at_bottom_edge
        array([0, 1, 2, 3])

        LLCATS: NINF BC SUBSET
        """
        try:
            return self._nodes[0, :]
        except AttributeError:
            raise AttributeError("Only rectangular Hex grids have defined edges.")

    def node_row_and_column(self, node_id):
        """Row and column from node ID, FOR VERT RECT CONFIGURATION ONLY.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> grid = HexModelGrid(3, 4, node_layout='rect', orientation='vert')
        >>> grid.node_row_and_column(5)
        (1, 2)
        >>> grid = HexModelGrid(3, 5, node_layout='rect', orientation='vert')
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
        apothem = self._dx / 2.0
        # distance from node to each hexagon cell vertex
        radius = 2.0 * apothem / sqrt(3.0)

        # offsets from node x,y position
        offsets = zeros((6, 2))
        poly_verts = zeros((6, 2))

        # Figure out whether the orientation is horizontal or vertical
        if self.orientation[0] == "h":  # horizontal
            offsets[:, 0] = array([0.0, apothem, apothem, 0.0, -apothem, -apothem])
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
            offsets[:, 1] = array([apothem, 0.0, -apothem, -apothem, 0.0, apothem])

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
        plt.xlim([amin(self.node_x) - self._dx, amax(self.node_x) + self._dx])
        plt.ylim([amin(self.node_y) - self._dx, amax(self.node_y) + self._dx])

        return ax

    def set_watershed_boundary_condition_outlet_id(
        self, outlet_id, node_data, nodata_value=-9999.0
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
        node_data : field name or ndarray
            At-node field name or at-node data values to use for identifying
            watershed location.
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
        # get node_data if a field name
        node_data = self.return_array_or_field_values("node", node_data)

        # make ring of no data nodes
        self.status_at_node[self.boundary_nodes] = CLOSED_BOUNDARY

        # set no data nodes to inactive boundaries
        self.set_nodata_nodes_to_closed(node_data, nodata_value)

        # set the boundary condition (fixed value) at the outlet_node
        self.status_at_node[outlet_id] = FIXED_VALUE_BOUNDARY

    def set_watershed_boundary_condition(
        self, node_data, nodata_value=-9999.0, return_outlet_id=False
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
        node_data : field name or ndarray
            At-node field name or at-node data values to use for identifying
            watershed location.
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
        # get node_data if a field name
        node_data = self.return_array_or_field_values("node", node_data)

        # make ring of no data nodes
        self.status_at_node[self.boundary_nodes] = CLOSED_BOUNDARY

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
        self.status_at_node[outlet_loc] = FIXED_VALUE_BOUNDARY

        if return_outlet_id:
            return as_id_array(numpy.array([outlet_loc]))


def from_dict(param_dict):
    """
    Create a HexModelGrid from the dictionary-like object, *param_dict*.
    Required keys of the dictionary are NUM_ROWS, NUM_COLS. Raises a KeyError
    if either of these are missing.  If GRID_SPACING is given, use it as the
    HexModelGrid *dx* parameter, otherwise default to unit spacing.

    Deprecated in version 1.6.X. Will be removed in version 2.0.
    """
    msg = (
        "The non-class method version of 'from_dict' for RasterModelGrid "
        "was Deprecated in version 1.6.X. Will be removed in version 2.0."
    )
    warn(msg, DeprecationWarning)
    # Read and create a basic HexModelGrid
    try:
        n_rows = int(param_dict["NUM_ROWS"])
        n_cols = int(param_dict["NUM_COLS"])
        dx = float(param_dict.get("GRID_SPACING", 1.0))
    except KeyError:
        raise
    except ValueError:
        raise
    else:
        hg = HexModelGrid(n_rows, n_cols, dx)

    return hg
