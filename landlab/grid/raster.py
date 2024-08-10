#! /usr/env/python
"""A class used to create and manage regular raster grids for 2D numerical
models in Landlab.

Do NOT add new documentation here. Grid documentation is now built in a
semi- automated fashion. To modify the text seen on the web, edit the
files `docs/text_for_[gridfile].py.txt`.
"""
import contextlib

import numpy as np
import xarray as xr

from landlab.utils import structured_grid as sgrid
from landlab.utils.decorators import make_return_array_immutable

from ..core.utils import add_module_functions_to_class
from ..core.utils import as_id_array
from ..field import FieldError
from ..graph import DualUniformRectilinearGraph
from . import raster_funcs as rfuncs
from .base import ModelGrid
from .decorators import return_id_array
from .diagonals import DiagonalsMixIn
from .nodestatus import NodeStatus


def _node_has_boundary_neighbor(mg, id, method="d8"):
    """Test if a RasterModelGrid node is next to a boundary.

    Test if one of the neighbors of node *id* is a boundary node.

    Parameters
    ----------
    mg : RasterModelGrid
        Source grid
    node_id : int
        ID of node to test.
    method: string, optional
        default is d8 neighbor, other method is 'd4'

    Returns
    -------
    boolean
        ``True`` if node has a neighbor on the boundary, ``False`` otherwise.
    """
    for neighbor in mg.active_adjacent_nodes_at_node[id]:
        try:
            if mg.status_at_node[neighbor] != NodeStatus.CORE:
                return True
        except IndexError:
            return True
    if method == "d8":
        for neighbor in mg.diagonal_adjacent_nodes_at_node[id]:
            try:
                if mg.status_at_node[neighbor] != NodeStatus.CORE:
                    return True
            except IndexError:
                return True
    return False


_node_has_boundary_neighbor = np.vectorize(_node_has_boundary_neighbor, excluded=["mg"])


def grid_edge_is_closed_from_dict(boundary_conditions):
    """Get a list of closed-boundary status at grid edges.

    Get a list that indicates grid edges that are closed boundaries. The
    returned list provides a boolean that gives the boundary condition status
    for edges order as [*bottom*, *left*, *top*, *right*].

    *boundary_conditions* is a dict whose keys indicate edge location (as
    "bottom", "left", "top", "right") and values must be one of "open", or
    "closed". If an edge location key is missing, that edge is assumed to be
    *open*.

    Parameters
    ----------
    boundary_conditions : dict
        Boundary condition for grid edges.

    Returns
    -------
    list
        List of booleans indicating if an edge is a closed boundary.

    Examples
    --------
    >>> from landlab.grid.raster import grid_edge_is_closed_from_dict
    >>> grid_edge_is_closed_from_dict(dict(bottom="closed", top="open"))
    [False, False, False, True]
    >>> grid_edge_is_closed_from_dict({})
    [False, False, False, False]
    """
    for condition in boundary_conditions.values():
        if condition not in ["open", "closed"]:
            raise ValueError("%s: boundary condition type not understood", condition)

    return [
        boundary_conditions.get(loc, "open") == "closed"
        for loc in ["right", "top", "left", "bottom"]
    ]


class RasterModelGrid(DiagonalsMixIn, DualUniformRectilinearGraph, ModelGrid):
    """A 2D uniform rectilinear grid.

    Examples
    --------
    Create a uniform rectilinear grid that has 4 rows and 5 columns of nodes.
    Nodes along the edges will be *open*. That is, links connecting these
    nodes to core nodes are *active*.

    >>> from landlab import RasterModelGrid
    >>> rmg = RasterModelGrid((4, 5))
    >>> rmg.number_of_node_rows, rmg.number_of_node_columns
    (4, 5)
    >>> rmg.number_of_active_links
    17

    Set the nodes along the top edge of the grid to be *closed* boundaries.
    This means that any links touching these nodes will be *inactive*.

    >>> rmg = RasterModelGrid((4, 5), bc={"top": "closed"})
    >>> rmg.number_of_node_rows, rmg.number_of_node_columns
    (4, 5)
    >>> rmg.number_of_active_links
    14

    A `RasterModelGrid` can have different node spacings in the *x* and *y*
    directions.

    >>> grid = RasterModelGrid((4, 5), xy_spacing=(2, 1))
    >>> grid.dx, grid.dy
    (2.0, 1.0)
    >>> grid.node_y.reshape(grid.shape)
    array([[0.,  0.,  0.,  0.,  0.],
           [1.,  1.,  1.,  1.,  1.],
           [2.,  2.,  2.,  2.,  2.],
           [3.,  3.,  3.,  3.,  3.]])
    >>> grid.node_x.reshape(grid.shape)
    array([[0.,  2.,  4.,  6.,  8.],
           [0.,  2.,  4.,  6.,  8.],
           [0.,  2.,  4.,  6.,  8.],
           [0.,  2.,  4.,  6.,  8.]])
    """

    def __init__(
        self,
        shape,
        xy_spacing=1.0,
        xy_of_lower_left=(0.0, 0.0),
        xy_of_reference=(0.0, 0.0),
        xy_axis_name=("x", "y"),
        xy_axis_units="-",
        bc=None,
    ):
        """Create a 2D grid with equal spacing.

        Optionally takes numbers of rows and columns and cell size as
        inputs. If this are given, calls initialize() to set up the grid.
        At the moment, num_rows and num_cols MUST be specified. Both must be
        >=3 to allow correct automated setup of boundary conditions.

        Parameters
        ----------
        shape : tuple of int
            Shape of the grid in nodes as (nrows, ncols).
        xy_spacing : tuple or float, optional
            dx and dy spacing. Either provided as a float or a
            (dx, dy) tuple.
        xy_of_lower_left: tuple, optional
            (x, y) coordinates of the lower left corner.
        xy_of_reference : tuple, optional
            Coordinate value in projected space of the reference point,
            `xy_of_lower_left`. Default is (0., 0.)
        xy_axis_name: tuple of str
            Name to use for each axis.
        xy_axis_units: tuple of str, or str
            Units for coordinates of each axis.
        bc : dict, optional
            Edge boundary conditions.

        Returns
        -------
        RasterModelGrid
            A newly-created grid.

        Notes
        -----
        The option for NOT giving rows, cols, and dx no longer works,
        because the *field* init requires num_active_cells, etc., to be
        defined. Either we force users to give arguments on instantiation,
        or set it up such that one can create a zero-node grid.
        """
        shape = tuple(shape)
        xy_spacing = np.asfarray(np.broadcast_to(xy_spacing, 2))
        self._xy_of_lower_left = tuple(np.asfarray(xy_of_lower_left))

        if shape[0] <= 0 or shape[1] <= 0:
            raise ValueError("number of rows and columns must be positive")

        DualUniformRectilinearGraph.__init__(
            self, shape, spacing=xy_spacing[::-1], origin=self.xy_of_lower_left[::-1]
        )
        ModelGrid.__init__(
            self,
            xy_axis_name=xy_axis_name,
            xy_axis_units=xy_axis_units,
            xy_of_reference=xy_of_reference,
        )

        self._node_status = np.full(
            self.number_of_nodes, NodeStatus.CORE, dtype=np.uint8
        )
        self._node_status[self.perimeter_nodes] = NodeStatus.FIXED_VALUE

        if bc is None:
            bc = {"right": "open", "top": "open", "left": "open", "bottom": "open"}

        if "closed" in bc.values():
            self.set_closed_boundaries_at_grid_edges(*grid_edge_is_closed_from_dict(bc))

        self.looped_node_properties = {}

        # List of looped neighbor cells (all 8 neighbors) for
        # given *cell ids* can be created if requested by the user.
        self._looped_cell_neighbor_list = None

        # List of second ring looped neighbor cells (all 16 neighbors) for
        # given *cell ids* can be created if requested by the user.
        self._looped_second_ring_cell_neighbor_list_created = False

    def __repr__(self):
        return "RasterModelGrid({}, xy_spacing={}, xy_of_lower_left={})".format(
            repr(self.shape),
            repr((self.dx, self.dy)),
            repr((self.x_of_node.min(), self.y_of_node.min())),
        )

    def __setstate__(self, state_dict):
        """Set state for of RasterModelGrid from pickled state_dict."""
        if state_dict["type"] != "RasterModelGrid":
            assert TypeError("Saved model instance not of " "RasterModelGrid type.")

        xy_spacing = state_dict["xy_spacing"]
        shape = state_dict["shape"]
        xy_of_lower_left = state_dict["xy_of_lower_left"]
        xy_of_reference = state_dict["xy_of_reference"]
        xy_axis_name = state_dict["xy_axis_name"]
        xy_axis_units = state_dict["xy_axis_units"]

        status_at_node = state_dict["status_at_node"]

        RasterModelGrid.__init__(
            self,
            shape,
            xy_spacing=xy_spacing,
            xy_of_lower_left=xy_of_lower_left,
            xy_of_reference=xy_of_reference,
            xy_axis_name=xy_axis_name,
            xy_axis_units=xy_axis_units,
        )
        self.status_at_node = status_at_node

        # Add fields back to the grid
        fields = state_dict["fields"]
        for at in fields:
            for name in fields[at]:
                values = fields[at][name]["array"]
                units = fields[at][name]["units"]
                self.add_field(name, values, at=at, units=units)

    def __getstate__(self):
        """Get state for pickling."""
        state_dict = {}

        # save basic information about the shape and size of the grid
        state_dict = {
            "type": "RasterModelGrid",
            "xy_spacing": (self.dx, self.dy),
            "shape": self.shape,
            "xy_of_lower_left": self.xy_of_lower_left,
            "xy_of_reference": self.xy_of_reference,
            "xy_axis_name": self.axis_name,
            "xy_axis_units": self.axis_units,
            # save status information at nodes (status at link set based on status
            # at node
            "status_at_node": np.asarray(self._node_status),
        }

        groups = {}
        for at in ("node", "link", "patch", "corner", "face", "cell", "grid"):
            groups[at] = {}
            for name in self[at].keys():
                groups[at][name] = {
                    "array": self.field_values(at, name),
                    "units": self.field_units(at, name),
                }

        state_dict["fields"] = groups

        return state_dict

    @classmethod
    def from_dict(cls, params):
        """Create a RasterModelGrid from a dictionary.

        Parameters
        ----------
        params : dict_like
            Initialization parameters for a RasterModelGrid.

        Returns
        -------
        RasterModelGrid
            A newly-created grid.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid.from_dict({"shape": (3, 4), "bc": {"top": "closed"}})
        >>> grid.number_of_nodes
        12

        :meta landlab: info-grid
        """
        shape = params.pop("shape", None)
        return cls(shape, **params)

    @classmethod
    def from_dataset(cls, dataset):
        return cls(
            dataset["shape"],
            xy_spacing=dataset["xy_spacing"],
            xy_of_lower_left=dataset["xy_of_lower_left"],
        )

    def as_dataset(self, include="*", exclude=None, time=None):
        dataset = xr.Dataset(
            {
                "shape": (("dim",), list(self.shape)),
                "xy_spacing": (("dim",), [self.dx, self.dy]),
                "xy_of_lower_left": (("dim",), list(self.xy_of_lower_left)),
            },
            attrs={"grid_type": "uniform_rectilinear"},
        )
        return dataset.update(
            super().as_dataset(include=include, exclude=exclude, time=time)
        )

    @property
    def xy_of_lower_left(self):
        """Return (x, y) of the reference point."""
        return self._xy_of_lower_left

    @xy_of_lower_left.setter
    def xy_of_lower_left(self, xy_of_lower_left):
        """Set a new value for the xy_of_lower_left."""
        dx = self.xy_of_lower_left[0] - xy_of_lower_left[0]
        dy = self.xy_of_lower_left[1] - xy_of_lower_left[1]
        with self.thawed():
            self.x_of_node[:] -= dx
            self.y_of_node[:] -= dy
        self._xy_of_lower_left = tuple(np.asfarray(xy_of_lower_left))

    @property
    def cell_grid_shape(self):
        """Get the shape of the cellular grid (grid with only cells).

        Returns
        -------
        shape : tuple of ints
            The shape of the cellular grid as number of cell rows and cell
            columns.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4))
        >>> grid.cell_grid_shape
        (1, 2)

        :meta landlab: info-grid, info-cell
        """
        return (self.number_of_cell_rows, self.number_of_cell_columns)

    def _create_link_unit_vectors(self):
        """Make arrays to store the unit vectors associated with each link.

        Creates self.link_unit_vec_x and self.link_unit_vec_y. These contain,
        for each link, the x and y components of the link's unit vector (that
        is, the link's x and y dimensions if it were shrunk to unit length but
        retained its orientation). The length of these arrays is the number of
        links plus one. The last entry in each array is set to zero, and is
        used to handle references to "link -1" (meaning, a non-existent link,
        whose unit vector is (0,0)).

        Also builds arrays to store the unit-vector component sums for each
        node: node_unit_vector_sum_x and node_unit_vector_sum_y. These are
        designed to be used when mapping link vector values to nodes (one takes
        the average of the x- and y-components of all connected links).

        Notes
        -----
        .. note::

            Overrides ModelGrid._create_link_unit_vectors().

        Creates the following:

        *  `self.link_unit_vec_x`, `self.link_unit_vec_y` : `ndarray`
           x and y components of unit vectors at each link (extra 0
           entries at end)
        *  `self.node_vector_sum_x`, `self.node_vector_sum_y` : `ndarray`
           Sums of x & y unit vector components for each node. Sum is over all
           links connected to a given node.

        Examples
        --------
        In the example below, the first 8 links are vertical, and have unit
        vectors (0,1), whereas the remaining links are horizontal with (1,0).
        The middle columns have x-component vector sums equal to 2 (one
        horizontal inlink and one horizontal outlink), while the middle rows
        have y-component vector sums equal to 2 (one vertical inlink and one
        vertical outlink). The rest of the entries have ones, representing the
        left and right columns (only one horizontal link) and top and bottom
        rows (only one vertical link).

        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((3, 4), xy_spacing=(2.0, 2.0))

        >>> mg.unit_vector_at_link[:, 0]
        array([1.,  1.,  1.,  0.,  0.,  0.,  0.,
               1.,  1.,  1.,  0.,  0.,  0.,  0.,
               1.,  1.,  1.])
        >>> mg.unit_vector_at_link[:, 1]
        array([0.,  0.,  0.,  1.,  1.,  1.,  1.,
               0.,  0.,  0.,  1.,  1.,  1.,  1.,
               0.,  0.,  0.])

        >>> mg.unit_vector_at_node[:, 0]
        array([1.,  2.,  2.,  1.,  1.,  2.,  2.,  1.,  1.,  2.,  2.,  1.])
        >>> mg.unit_vector_at_node[:, 1]
        array([1.,  1.,  1.,  1.,  2.,  2.,  2.,  2.,  1.,  1.,  1.,  1.])
        """
        unit_vec_at_link = np.zeros((self.number_of_links + 1, 2), dtype=float)
        unit_vec_at_link[self.horizontal_links, 0] = 1.0
        unit_vec_at_link[self.vertical_links, 1] = 1.0

        self._unit_vec_at_node = unit_vec_at_link[self.links_at_node].sum(axis=1)
        self._unit_vec_at_link = unit_vec_at_link[:-1, :]

    @property
    def extent(self):
        """Extent of the grid in the y and x-dimensions.

        Return the y and x-dimension of the grid. Because boundary nodes
        don't have cells, the dimension of the grid is
        ``((num_rows - 1) * dy, (num_columns - 1) * dx)``, not
        ``(num_rows * dy, num_cols * dx)``.

        Returns
        -------
        (y_extent, x_extent) : tuple of float
            Length of the grid in the y and x-dimensions.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.extent
        (3.0, 4.0)

        >>> grid = RasterModelGrid((4, 5), xy_spacing=2.0)
        >>> grid.extent
        (6.0, 8.0)

        >>> grid = RasterModelGrid((4, 5), xy_spacing=(3, 2))
        >>> grid.extent
        (6.0, 12.0)

        :meta landlab: info-grid, quantity
        """
        # Method added 5/1/13 by DEJH, modified DEJH 4/3/14 to reflect fact
        # boundary nodes don't have defined
        return (
            (self.number_of_node_rows - 1) * self.dy,
            (self.number_of_node_columns - 1) * self.dx,
        )

    @property
    def number_of_interior_nodes(self):
        """Number of interior nodes.

        Returns the number of interior nodes on the grid, i.e., non-perimeter
        nodes. Compare self.number_of_core_nodes.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.number_of_interior_nodes
        6

        :meta landlab: info-node
        """
        return sgrid.interior_node_count(self.shape)

    @property
    def number_of_cell_columns(self):
        """Number of cell columns.

        Returns the number of columns, including boundaries.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.number_of_cell_columns
        3

        :meta landlab: info-grid, info-node
        """
        return self.shape[1] - 2

    @property
    def number_of_cell_rows(self):
        """Number of cell rows.

        Returns the number of rows, including boundaries.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.number_of_cell_rows
        2

        :meta landlab: info-grid, info-cell
        """
        return self.shape[0] - 2

    @property
    def cells_at_corners_of_grid(self):
        """Get array of cells in cellular grid (grid with only cells) corners.

        Return the IDs to the corner cells of the cellular grid, sorted by ID.

        Returns
        -------
        (4, ) ndarray
            Array of corner node IDs.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.cells_at_corners_of_grid
        array([0, 2, 3, 5])

        :meta landlab: info-grid, info-cell, subset
        """
        return sgrid.corners(self.cell_grid_shape)

    def is_point_on_grid(self, xcoord, ycoord):
        """Check if a point is on the grid.

        This method takes x, y coordinates and tests whether they lie within
        the grid. The limits of the grid are taken to be links connecting the
        boundary nodes. We perform a special test to detect looped boundaries.

        Coordinates can be ints or arrays of ints. If arrays, will return an
        array of the same length of boolean truth values.

        Parameters
        ----------
        xcoord : float or array_like
            The point's x-coordinate.
        ycoord : float or array_like
            The point's y-coordinate.

        Returns
        -------
        bool
            ``True`` if the point is on the grid. Otherwise, ``False``.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5), xy_spacing=(1, 2))
        >>> grid.is_point_on_grid(1, 1)
        True
        >>> grid.is_point_on_grid(
        ...     (1, 1, 1),
        ...     (1, 3.1, 6.1),
        ... )
        array([ True,  True, False])
        >>> grid.is_point_on_grid((-0.1, 0.1, 3.9, 4.1), (1, 1, 1, 1))
        array([False,  True,  True, False])

        :meta landlab: info-grid, quantity, subset
        """
        xcoord, ycoord = np.asarray(xcoord), np.asarray(ycoord)

        x_condition = (xcoord > 0.0) & (xcoord < (self.shape[1] - 1) * self.dx)
        y_condition = (ycoord > 0.0) & (ycoord < (self.shape[0] - 1) * self.dy)

        if np.all(self._node_status[self.nodes_at_left_edge] == 3) or np.all(
            self._node_status[self.nodes_at_right_edge] == 3
        ):
            try:
                x_condition[:] = 1
            except IndexError:
                x_condition = 1
        if np.all(self._node_status[self.nodes_at_top_edge] == 3) or np.all(
            self._node_status[self.nodes_at_bottom_edge] == 3
        ):
            try:
                y_condition[:] = 1
            except IndexError:
                y_condition = 1

        return x_condition & y_condition

    def nodes_around_point(self, xcoord, ycoord):
        """Get the nodes surrounding a point.

        Return IDs of the four nodes of the area around a point with
        coordinates *xcoord*, *ycoord*. Node IDs are returned
        counter-clockwise order starting from the southwest node.

        If either *xcoord* or *ycoord* are arrays the usual numpy broadcasting
        rules apply.

        Parameters
        ----------
        xcoord : float, array-like
            x-coordinate of point
        ycoord : float, array-like
            y-coordinate of point

        Returns
        -------
        (4, N) ndarray
            IDs of nodes around the point.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4))
        >>> grid.nodes_around_point(0.4, 1.2)
        array([4, 8, 9, 5])

        >>> grid.nodes_around_point([0.9, 1.1], 1.2)
        array([[ 4,  5],
               [ 8,  9],
               [ 9, 10],
               [ 5,  6]])

        >>> grid = RasterModelGrid((3, 4), xy_spacing=(1, 2))
        >>> grid.nodes_around_point(0.5, 1.5)
        array([0, 4, 5, 1])
        >>> grid = RasterModelGrid((3, 4))
        >>> grid.nodes_around_point(0.5, 1.5)
        array([4, 8, 9, 5])

        :meta landlab: info-node, subset
        """
        xcoord, ycoord = np.broadcast_arrays(xcoord, ycoord)

        # Method added 4/29/13 by DEJH, modified 9/24/13.
        id_ = ycoord // self.dy * self.number_of_node_columns + xcoord // self.dx
        try:
            id_ = int(id_)
        except TypeError:
            id_ = as_id_array(id_)
        return np.array(
            [
                id_,
                id_ + self.number_of_node_columns,
                id_ + self.number_of_node_columns + 1,
                id_ + 1,
            ]
        )

    def find_nearest_node(self, coords, mode="raise"):
        """Node nearest a point.

        Find the index to the node nearest the given x, y coordinates.
        Coordinates are provided as numpy arrays in the *coords* tuple.

        Use the *mode* keyword to specify what to do if the given coordinates
        are out-of-bounds. See :func:`np.ravel_multi_index` for a
        description of possible values for *mode*. Note that a coordinate is
        out-of-bounds if it is beyond one half the node spacing from the
        exterior nodes.

        Returns the indices of the nodes nearest the given coordinates.

        Parameters
        ----------
        coords : tuple of array-like
            Coordinates of points.
        mode : {'raise', 'wrap', 'clip'}, optional
            What to do if a point is off the grid.

        Returns
        -------
        array-like
            IDs of the nearest nodes.

        Notes
        -----
        For coordinates that are equidistant to two or more nodes, see
        the rounding rules for :func:`numpy.around`.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid((4, 5))
        >>> rmg.find_nearest_node([0.2, 0.2])
        0
        >>> rmg.find_nearest_node((np.array([1.6, 3.6]), np.array([2.3, 0.7])))
        array([12,  9])
        >>> rmg.find_nearest_node((-0.4999, 1.0))
        5

        :meta landlab: info-node, subset
        """
        return rfuncs.find_nearest_node(self, coords, mode=mode)

    def set_closed_boundaries_at_grid_edges(
        self, right_is_closed, top_is_closed, left_is_closed, bottom_is_closed
    ):
        """Set boundary not to be closed.

        Sets the status of nodes along the specified side(s) of a raster
        grid (bottom, right, top, and/or left) to ``BC_NODE_IS_CLOSED``.

        Arguments are booleans indicating whether the bottom, left, top, and
        right are closed (``True``) or not (``False``).

        For a closed boundary:

        *  the nodes are flagged ``BC_NODE_IS_CLOSED`` (status type 4)
        *  all links that connect to a ``BC_NODE_IS_CLOSED`` node are
           flagged as inactive (so they appear on link-based lists, but
           not active_link-based lists)

        This means that if you call the calc_grad_at_active_link
        method, links connecting to closed boundaries will be ignored: there
        can be no gradients or fluxes calculated, because the links that
        connect to that edge of the grid are not included in the calculation.
        So, setting a grid edge to BC_NODE_IS_CLOSED is a convenient way to
        impose a no-flux boundary condition. Note, however, that this applies
        to the grid as a whole, rather than a particular variable that you
        might use in your application. In other words, if you want a no-flux
        boundary in one variable but a different boundary condition for
        another, then use another method.

        Parameters
        ----------
        right_is_closed : boolean
            If ``True`` right-edge nodes are closed boundaries.
        top_is_closed : boolean
            If ``True`` top-edge nodes are closed boundaries.
        left_is_closed : boolean
            If ``True`` left-edge nodes are closed boundaries.
        bottom_is_closed : boolean
            If ``True`` bottom-edge nodes are closed boundaries.

        Notes
        -----
        Note that the four corners are treated as follows:

        *  bottom left = BOTTOM
        *  bottom right = BOTTOM
        *  top right = TOP
        *  top left = TOP

        This scheme is necessary for internal consistency with looped
        boundaries.

        Examples
        --------
        The following example sets the top and left boundaries as closed in a
        four-row by five-column grid that initially has all boundaries open
        and all boundary nodes coded as BC_NODE_IS_FIXED_VALUE (=1):

        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid((4, 5))  # rows, columns, spacing
        >>> rmg.number_of_active_links
        17
        >>> rmg.status_at_node.reshape(rmg.shape)
        array([[1, 1, 1, 1, 1],
               [1, 0, 0, 0, 1],
               [1, 0, 0, 0, 1],
               [1, 1, 1, 1, 1]], dtype=uint8)
        >>> rmg.set_closed_boundaries_at_grid_edges(True, True, False, False)
        >>> rmg.number_of_active_links
        12
        >>> rmg.status_at_node.reshape(rmg.shape)
        array([[1, 1, 1, 1, 1],
               [1, 0, 0, 0, 4],
               [1, 0, 0, 0, 4],
               [4, 4, 4, 4, 4]], dtype=uint8)

        :meta landlab: boundary-condition, subset
        """
        if bottom_is_closed:
            self._node_status[self.nodes_at_bottom_edge] = self.BC_NODE_IS_CLOSED

        if right_is_closed:
            self._node_status[self.nodes_at_right_edge[1:-1]] = self.BC_NODE_IS_CLOSED

        if top_is_closed:
            self._node_status[self.nodes_at_top_edge] = self.BC_NODE_IS_CLOSED

        if left_is_closed:
            self._node_status[self.nodes_at_left_edge[1:-1]] = self.BC_NODE_IS_CLOSED

        self.reset_status_at_node()

    def set_fixed_value_boundaries_at_grid_edges(
        self,
        right_is_fixed_val,
        top_is_fixed_val,
        left_is_fixed_val,
        bottom_is_fixed_val,
        value=None,
        value_of="topographic__elevation",
    ):
        """Create fixed values boundaries.

        Sets the status of nodes along the specified side(s) of a raster
        grid---bottom, right, top, and/or left---to NODE_IS_FIXED_VALUE

        Arguments are booleans indicating whether the bottom, right, top, and
        left sides are to be set (True) or not (False).

        *value* controls what values are held constant at these nodes. It can
        be either a float, an array of length number_of_fixed_nodes or
        number_of_nodes (total), or left blank. If left blank, the values will
        be set from the those already in the grid fields, according to
        'value_of'.

        *value_of* controls the name of the model field that contains the
        values. Remember, if you don't set value, the fixed values will be set
        from the field values ***at the time you call this method***. If no
        values are present in the field, the module will complain but accept
        this, warning that it will be unable to automatically update boundary
        conditions.

        The status of links (active or inactive) is automatically updated to
        reflect the changes.

        The following example sets the bottom and right boundaries as
        fixed-value in a four-row by five-column grid that initially has all
        boundaries closed (i.e., flagged as node_status=4):

        Parameters
        ----------
        bottom_is_fixed_val : boolean
            Set bottom edge as fixed boundary.
        left_is_fixed_val : boolean
            Set left edge as fixed boundary.
        top_is_fixed_val : boolean
            Set top edge as fixed boundary.
        right_is_fixed_val : boolean
            Set right edge as fixed boundary.
        value : float, array or None (default).
            Override value to be kept constant at nodes.
        value_of : string.
            The name of the grid field containing the values of interest.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid((4, 5), xy_spacing=(1, 1))
        >>> rmg.number_of_active_links
        17

        Put some arbitrary values in the grid fields:

        >>> import numpy as np
        >>> rmg.at_node["topographic__elevation"] = np.random.rand(20)
        >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> rmg.status_at_node.reshape(rmg.shape)
        array([[4, 4, 4, 4, 4],
               [4, 0, 0, 0, 4],
               [4, 0, 0, 0, 4],
               [4, 4, 4, 4, 4]], dtype=uint8)
        >>> rmg.set_fixed_value_boundaries_at_grid_edges(True, True, False, False)
        >>> rmg.number_of_active_links
        12
        >>> rmg.status_at_node.reshape(rmg.shape)
        array([[4, 4, 4, 4, 4],
               [4, 0, 0, 0, 1],
               [4, 0, 0, 0, 1],
               [1, 1, 1, 1, 1]], dtype=uint8)

        Note that the four corners are treated as follows:

        *  bottom left = BOTTOM
        *  bottom right = BOTTOM
        *  top right = TOP
        *  top left = TOP

        This scheme is necessary for internal consistency with looped
        boundaries.

        :meta landlab: boundary-condition, subset
        """
        bottom_edge = range(0, self.number_of_node_columns)
        right_edge = range(
            2 * self.number_of_node_columns - 1,
            self.number_of_nodes - 1,
            self.number_of_node_columns,
        )
        top_edge = range(
            (self.number_of_node_rows - 1) * self.number_of_node_columns,
            self.number_of_nodes,
        )
        left_edge = range(
            self.number_of_node_columns,
            self.number_of_nodes - self.number_of_node_columns,
            self.number_of_node_columns,
        )

        if bottom_is_fixed_val:
            self._node_status[bottom_edge] = NodeStatus.FIXED_VALUE

        if right_is_fixed_val:
            self._node_status[right_edge] = NodeStatus.FIXED_VALUE

        if top_is_fixed_val:
            self._node_status[top_edge] = NodeStatus.FIXED_VALUE

        if left_is_fixed_val:
            self._node_status[left_edge] = NodeStatus.FIXED_VALUE

        self.reset_status_at_node()

        # save some internal data to speed updating:
        self.fixed_value_node_properties = {}
        self.fixed_value_node_properties["boundary_node_IDs"] = as_id_array(
            np.where(self._node_status == NodeStatus.FIXED_VALUE)[0]
        )

        if value:
            if isinstance(value, (float, int)):
                values_to_use = float(value)
            elif isinstance(value, np.ndarray):
                if (
                    value.size
                    == self.fixed_value_node_properties["boundary_node_IDs"].size
                ):
                    values_to_use = value
                elif value.size == self.number_of_nodes:
                    values_to_use = value.take(
                        self.fixed_value_node_properties["boundary_node_IDs"]
                    )
                else:
                    raise TypeError(
                        "'value' must be of size nnodes or number of nodes " "to set!"
                    )
        else:
            try:
                values_to_use = self.at_node[value_of].take(
                    self.fixed_value_node_properties["boundary_node_IDs"]
                )
            except FieldError:
                pass  # we catch this case below
            else:
                # set a flag to indicate successful setting of internal values
                self.fixed_value_node_properties["internal_flag"] = True

        if not self.has_field("node", value_of):
            print(
                """
                *************************************************
                WARNING: set_fixed_value_boundaries_at_grid_edges
                has not been provided with a grid field name to
                allow internal boundary condition control. You
                will not be able to automate BC control with grid
                methods like update_boundary_nodes()!
                Not expecting this error? Try calling this method
                after loading the starting conditions into the
                grid fields.
                *************************************************
                """
            )

            # set a flag to indicate no internal values
            self.fixed_value_node_properties["internal_flag"] = False
        else:
            self.fixed_value_node_properties["internal_flag"] = True
            self.fixed_value_node_properties["fixed_value_of"] = value_of
        with contextlib.suppress(NameError):
            self.fixed_value_node_properties["values"] = values_to_use

    def set_looped_boundaries(self, top_bottom_are_looped, sides_are_looped):
        """Create wrap-around boundaries.

        Handles boundary conditions by setting corresponding parallel grid
        edges as looped "tracks_cell" (==3) status, linked to each other.
        If top_bottom_are_looped is True, the top and bottom edges will link
        to each other. If sides_are_looped is True, the left and right edges
        will link to each other.

        Looped boundaries are experimental, and not as yet well integrated into
        the Landlab framework. Many functions may not recognise them, or
        silently create unforeseen errors. Use at your own risk!

        Note that because of the symmetries this BC implies, the corner nodes
        are all paired with the bottom/top edges, not the sides.

        Parameters
        ----------
        top_bottom_are_looped : bool
            Top and bottom are wrap-around.
        sides_are_looped : bool
            Left and right sides are wrap-around.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid((4, 5))  # rows, columns, spacing
        >>> rmg.number_of_active_links
        17
        >>> rmg.status_at_node.reshape(rmg.shape)
        array([[1, 1, 1, 1, 1],
               [1, 0, 0, 0, 1],
               [1, 0, 0, 0, 1],
               [1, 1, 1, 1, 1]], dtype=uint8)
        >>> rmg.add_zeros("topographic__elevation", at="node")
        array([0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
               0.,  0.,  0.,  0.,  0.,  0.,  0.])
        >>> rmg.set_looped_boundaries(True, True)
        >>> rmg.looped_node_properties["boundary_node_IDs"]
        array([ 0,  1,  2,  3,  4,  5,  9, 10, 14, 15, 16, 17, 18, 19])
        >>> rmg.looped_node_properties["linked_node_IDs"]
        array([10, 11, 12, 13, 14,  8,  6, 13, 11,  5,  6,  7,  8,  9])

        :meta landlab: boundary-condition, subset
        """
        # Added DEJH Feb 2014
        # TODO: Assign BC_statuses also to *links*

        bottom_edge = np.array(range(0, self.number_of_node_columns))
        right_edge = np.array(
            range(
                2 * self.number_of_node_columns - 1,
                self.number_of_nodes - 1,
                self.number_of_node_columns,
            )
        )
        top_edge = np.array(
            range(
                (self.number_of_node_rows - 1) * self.number_of_node_columns,
                self.number_of_nodes,
            )
        )
        left_edge = np.array(
            range(
                self.number_of_node_columns,
                (self.number_of_nodes - self.number_of_node_columns),
                self.number_of_node_columns,
            )
        )
        these_boundary_IDs = np.array([])
        these_linked_nodes = np.array([])

        if top_bottom_are_looped:
            self._node_status[bottom_edge] = NodeStatus.LOOPED
            self._node_status[top_edge] = NodeStatus.LOOPED
            these_boundary_IDs = np.concatenate(
                (these_boundary_IDs, bottom_edge, top_edge)
            )
            these_linked_nodes = np.concatenate(
                (
                    these_linked_nodes,
                    top_edge - self.number_of_node_columns,
                    bottom_edge + self.number_of_node_columns,
                )
            )

        if sides_are_looped:
            self._node_status[right_edge] = NodeStatus.LOOPED
            self._node_status[left_edge] = NodeStatus.LOOPED
            these_boundary_IDs = np.concatenate(
                (these_boundary_IDs, left_edge, right_edge)
            )
            these_linked_nodes = np.concatenate(
                (these_linked_nodes, right_edge - 1, left_edge + 1)
            )

        self.reset_status_at_node()

        if not self.looped_node_properties:
            existing_IDs = np.array([])
            existing_links = np.array([])
        else:
            unrepeated_node_entries = np.logical_not(
                np.in1d(
                    self.looped_node_properties["boundary_node_IDs"], these_linked_nodes
                )
            )
            existing_IDs = self.looped_node_properties["boundary_node_IDs"][
                unrepeated_node_entries
            ]
            existing_links = self.looped_node_properties["linked_node_IDs"][
                unrepeated_node_entries
            ]

        self.looped_node_properties = {}
        all_the_IDs = np.concatenate((these_boundary_IDs, existing_IDs))
        ID_ordering = np.argsort(all_the_IDs)
        self.looped_node_properties["boundary_node_IDs"] = as_id_array(
            all_the_IDs[ID_ordering]
        )
        self.looped_node_properties["linked_node_IDs"] = as_id_array(
            np.concatenate((these_linked_nodes, existing_links))[ID_ordering]
        )

        if np.any(
            self._node_status[self.looped_node_properties["boundary_node_IDs"]] == 2
        ):
            raise AttributeError(
                "Switching a boundary between fixed gradient and looped will "
                "result in bad BC handling! Bailing out..."
            )

    def node_vector_to_raster(self, u, flip_vertically=False):
        """Unravel an array of node values.

        Converts node vector *u* to a 2D array and returns it, so that it
        can be plotted, output, etc.

        If the *flip_vertically* keyword is True, this function returns an
        array that has the rows in reverse order. This is useful for use in
        plot commands (such as the image display functions) that puts the
        first row at the top of the image. In the landlab coordinate system,
        the first row is thought to be at the bottom. Thus, a flipped matrix
        will plot in the landlab style with the first row at the bottom.

        The returned array is a view of *u*, not a copy.

        See also
        --------
        RasterModelGrid.nodes
            An equivalent property, but without the option to flip the grid.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid((4, 5))
        >>> u = rmg.zeros(centering="node")
        >>> u = u + range(0, len(u))
        >>> u
        array([ 0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.,
               11.,  12.,  13.,  14.,  15.,  16.,  17.,  18.,  19.])
        >>> ur = rmg.node_vector_to_raster(u)
        >>> ur
        array([[ 0.,   1.,   2.,   3.,   4.],
               [ 5.,   6.,   7.,   8.,   9.],
               [10.,  11.,  12.,  13.,  14.],
               [15.,  16.,  17.,  18.,  19.]])
        >>> ur = rmg.node_vector_to_raster(u, flip_vertically=True)
        >>> ur
        array([[15.,  16.,  17.,  18.,  19.],
               [10.,  11.,  12.,  13.,  14.],
               [ 5.,   6.,   7.,   8.,   9.],
               [ 0.,   1.,   2.,   3.,   4.]])

        :meta landlab: info-grid, info-node
        """
        return sgrid.reshape_array(self.shape, u, flip_vertically=flip_vertically)

    def cell_vector_to_raster(self, u, flip_vertically=False):
        """Unravel a 1D array.

        Converts cell vector u to a 2D array and returns it,
        so that it can be plotted, output, etc.

        If the optional argument flip_vertically=True, the function returns an
        array that has the rows in reverse order, for use in plot commands
        (such as the image display functions) that put the (0,0) axis at the
        top left instead of the bottom left.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid((4, 5))
        >>> u = rmg.zeros(centering="cell")
        >>> u = u + range(0, len(u))
        >>> u
        array([0.,  1.,  2.,  3.,  4.,  5.])
        >>> ur = rmg.cell_vector_to_raster(u)
        >>> ur
        array([[0.,  1.,  2.],
               [3.,  4.,  5.]])
        >>> ur = rmg.cell_vector_to_raster(u, flip_vertically=True)
        >>> ur
        array([[3.,  4.,  5.],
               [0.,  1.,  2.]])

        :meta landlab: info-grid, info-cell
        """
        return sgrid.reshape_array(
            (self.shape[0] - 2, self.shape[1] - 2), u, flip_vertically=flip_vertically
        )

    def roll_nodes_ud(self, data_name, shift, interior_only=False):
        """Roll (shift) specified data on nodes up or down in a raster grid.

        Similar to the Numpy roll() function, in that it shifts node values up
        by *shift* rows, wrapping the data in the top row(s) around to the
        bottom. If the *interior_only* is set, data along the left and right
        grid edges are not changed.

        Note that the contents of the *data_name* field are changed.

        Parameters
        ----------
        data_name : string
            Name of node-data item attached to grid.
        shift : int
            Number of rows to shift upward.
        interior_only : bool, optional
            If True, data along left and right edges are not shifted

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid((4, 3))
        >>> data = rmg.add_zeros("test_data", at="node")
        >>> data[:] = np.arange(12)
        >>> rmg.roll_nodes_ud("test_data", 1)
        >>> data.reshape(rmg.shape)
        array([[  9.,  10.,  11.],
               [  0.,   1.,   2.],
               [  3.,   4.,   5.],
               [  6.,   7.,   8.]])
        >>> rmg.roll_nodes_ud("test_data", 2)
        >>> data.reshape(rmg.shape)
        array([[  3.,   4.,   5.],
               [  6.,   7.,   8.],
               [  9.,  10.,  11.],
               [  0.,   1.,   2.]])
        >>> rmg.roll_nodes_ud("test_data", 1, interior_only=True)
        >>> data.reshape(rmg.shape)
        array([[  3.,   1.,   5.],
               [  6.,   4.,   8.],
               [  9.,   7.,  11.],
               [  0.,  10.,   2.]])

        :meta landlab: info-node
        """
        # Get the data
        data = self.at_node[data_name]

        # Get the IDs of the nodes in the top row, and number of rows and cols
        top_ids = self.nodes_at_top_edge
        ncols = self.number_of_node_columns
        nrows = self.number_of_node_rows

        # To handle "interior only" option, we use the variable *offset*,
        # which is zero if shifting everything, and 1 if shifting just the
        # interior -- we use this to go from column 1 to column N-2 (instead
        # of 0 to N-1) when interior_only is True.
        if interior_only:
            offset = 1
            top_ids = top_ids[1 : ncols - 1]
        else:
            offset = 0

        # Remember the top N rows
        top_rows_to_move = np.zeros((shift, ncols - 2 * offset))
        for i in range(0, shift):
            top_rows_to_move[shift - (i + 1), :] = data[top_ids - i * ncols]

        # Go row by row, starting from top
        for i in range(nrows - shift):
            to_row = nrows - (i + 1)
            from_row = to_row - shift
            data[ncols * to_row + offset : ncols * (to_row + 1) - offset] = data[
                ncols * from_row + offset : ncols * (from_row + 1) - offset
            ]

        # now replace the bottom *shift* rows
        for i in range(0, shift):
            data[ncols * i + offset : ncols * (i + 1) - offset] = top_rows_to_move[i, :]

    def node_has_boundary_neighbor(self, ids, method="d8"):
        """Check if nodes have neighbors that are boundary nodes.

        Checks to see if one of the eight neighbor nodes of node(s) with
        *id* has a boundary node.  Returns True if a node has a boundary node,
        False if all neighbors are interior.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((5, 5))
        >>> mg.node_has_boundary_neighbor(6)
        True
        >>> mg.node_has_boundary_neighbor(12)
        False
        >>> mg.node_has_boundary_neighbor([12, -1])
        array([False,  True])

        >>> mg.node_has_boundary_neighbor(25)
        Traceback (most recent call last):
        ...
        IndexError: index 25 is out of bounds for axis 0 with size 25

        :meta landlab: info-node, connectivity, boundary-condition
        """
        ans = _node_has_boundary_neighbor(self, ids, method=method)

        if ans.ndim == 0:
            return bool(ans)
        else:
            return ans

    @return_id_array
    def grid_coords_to_node_id(self, row, col, **kwds):
        """Convert node indices to node ID.

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

        Parameters
        ----------
        row : array-like
            Row of node.
        col : array-like
            Column of node.

        Returns
        -------
        ndarray
            Node IDs.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((4, 5))
        >>> mg.grid_coords_to_node_id(2, 3)
        13

        >>> mg.grid_coords_to_node_id([2, 0], [3, 4])
        array([13,  4])

        :meta landlab: info-node, subset, quantity
        """
        return np.ravel_multi_index((row, col), self.shape, **kwds)

    def calc_unit_normal_at_patch(self, elevs="topographic__elevation"):
        """Calculate and return the unit normal vector <a, b, c> to a patch.

        This method is not defined on a raster, as there is no unique unit
        normal for a square patch. Use
        `_calc_unit_normals_to_patch_subtriangles` instead.

        :meta landlab: info-patch, gradient
        """
        raise NotImplementedError(
            "This method is not defined on a raster, as there is no unique "
            "unit normal for a square patch. Use "
            "`_calc_unit_normals_to_patch_subtriangles` instead."
        )

    def calculate_slope_aspect_at_nodes_burrough(self, ids=None, vals="Elevation"):
        """Calculate topographic slope.

        Calculates the local topographic slope (i.e., the down-dip slope, and
        presented as positive), and the aspect (dip direction in degrees
        clockwise from north), at the given nodes, *ids*. All *ids* must be of
        core nodes.
        This method uses Burrough's 1998 Pg. 190 method similar to the methods
        used by ArcMap to calculate slope and aspect.

        If *ids* is not provided, the slope will be returned for nodes at all
        cells.

        *vals* is either the name of an existing grid field from which to draw
        topographic data, or an array of values to use. If an array of values
        is passed, it must be nnodes long.
        If *vals* is not provided, this method will default to trying to use
        the field 'Elevation'.

        Returns
        -------
        (slope, aspect) : tuple of float
            *slope*, a len(ids) array of slopes at each node provided.
            *aspect*, a len(ids) array of aspects at each node provided.

        Examples
        --------
        >>> import pytest
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4), xy_spacing=(4, 4))
        >>> z = np.array([0.0, 0.0, 0.0, 0.0, 3.0, 3.0, 3.0, 3, 6.0, 6.0, 6.0, 6.0])
        >>> slope, aspect = grid.calculate_slope_aspect_at_nodes_burrough(vals=z)
        >>> np.tan(slope)
        array([0.75,  0.75])
        >>> np.degrees(aspect)
        array([180.,  180.])

        We recommend using the following functions instead of this one:
        - :py:meth:`~landlab.grid.RasterModelGrid.calc_slope_at_node`
        - :py:meth:`~landlab.grid.RasterModelGrid.calc_aspect_at_node`
        Notice that :py:meth:`~landlab.grid.RasterModelGrid.calc_slope_at_node`
        and `:py:meth:`~landlab.grid.RasterModelGrid.calc_aspect_at_node` return
        values for all nodes, not just core nodes. In addition,
        `:py:meth:`~landlab.grid.RasterModelGrid.calc_aspect_at_node` returns
        compass-style angles in degrees.

        >>> np.tan(grid.calc_slope_at_node(elevs=z)[grid.core_nodes])
        array([0.75,  0.75])
        >>> grid.calc_aspect_at_node(elevs=z)[grid.core_nodes]
        array([180.,  180.])

        :meta landlab: info-node, surface, gradient
        """
        if ids is None:
            ids = self.node_at_cell
        if not isinstance(ids, np.ndarray):
            ids = np.array([ids])
        if isinstance(vals, str):
            vals = self.at_node[vals]
        else:
            if len(vals) != self.number_of_nodes:
                raise IndexError("*vals* was not of a compatible length!")

        neighbors = np.zeros([ids.shape[0], 4], dtype=int)
        diagonals = np.zeros([ids.shape[0], 4], dtype=int)
        # [right, top, left, bottom]
        neighbors[:] = self.active_adjacent_nodes_at_node[ids]
        # [topright, topleft, bottomleft, bottomright]
        diagonals[:] = self.diagonal_adjacent_nodes_at_node[ids]

        right = vals[neighbors[:, 0]]
        top = vals[neighbors[:, 1]]
        left = vals[neighbors[:, 2]]
        bottom = vals[neighbors[:, 3]]
        top_right = vals[diagonals[:, 0]]
        top_left = vals[diagonals[:, 1]]
        bottom_left = vals[diagonals[:, 2]]
        bottom_right = vals[diagonals[:, 3]]

        dz_dx = (
            (top_right + 2 * right + bottom_right) - (top_left + 2 * left + bottom_left)
        ) / (8.0 * self.dx)
        dz_dy = (
            (bottom_left + 2 * bottom + bottom_right) - (top_left + 2 * top + top_right)
        ) / (8.0 * self.dy)

        slope = np.zeros([ids.shape[0]], dtype=float)
        aspect = np.zeros([ids.shape[0]], dtype=float)
        slope = np.arctan(np.sqrt(dz_dx**2 + dz_dy**2))
        aspect = np.arctan2(dz_dy, -dz_dx)
        aspect = np.pi * 0.5 - aspect
        aspect[aspect < 0.0] = aspect[aspect < 0.0] + 2.0 * np.pi
        aspect[slope == 0.0] = -1.0

        return slope, aspect

    def save(self, path, names=None, format=None, at=None):
        """Save a grid and fields.

        If more than one field name is specified for names when saving to ARC
        ascii, multiple files will be produced, suffixed with the field names.

        When saving to netCDF (.nc), the fields are incorporated into the
        single named .nc file.

        Parameters
        ----------
        path : str
            Path to output file.
        names : iterable of strings, optional
            List of field names to save, defaults to all if not specified.
        format : {'netcdf', 'esri-ascii'}, optional
            Output file format. Guess from file extension if not given.
        at : str
            Grid element where values are defined.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> import os
        >>> from tempfile import mkdtemp

        >>> grid = RasterModelGrid((4, 5))
        >>> fname = os.path.join(mkdtemp(), "mysave.nc")
        >>> grid.save(fname)
        >>> os.path.isfile(fname)
        True
        >>> os.remove(fname)

        :meta landlab: info-grid
        """
        from ..io import write_esri_ascii
        from ..io.netcdf import write_netcdf

        format = format or _guess_format_from_name(path)
        path = _add_format_extension(path, format)

        if format == "netcdf":
            write_netcdf(path, self, format="NETCDF3_64BIT", names=names, at=at)
        elif format == "esri-ascii":
            write_esri_ascii(path, self, names=names)
        else:
            raise ValueError("format not understood")

    @property
    @make_return_array_immutable
    def looped_neighbors_at_cell(self):
        """For each cell in a raster, return the D8 neighboring cells, looping
        across grid boundaries as necessary.

        Returns lists of looped neighbor cell IDs of given *cell ids*.
        If *cell ids* are not given, it returns a 2D array of size
        (self.number_of_cells, 8).
        Order or neighbors is [ E, NE, N, NW, W, SW, S, SE ]

        Output is looped, regardless of boundary conditions! (see examples)

        Returns
        -------
        ndarray (num_cells, 8)
            The eight neighbors of each cell.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> neighbors = grid.looped_neighbors_at_cell
        >>> neighbors[1, :]
        array([2, 5, 4, 3, 0, 3, 4, 5])
        >>> neighbors[5, :]
        array([3, 0, 2, 1, 4, 1, 2, 0])
        >>> grid.looped_neighbors_at_cell[np.array([1, 5]), :]
        array([[2, 5, 4, 3, 0, 3, 4, 5],
               [3, 0, 2, 1, 4, 1, 2, 0]])

        :meta landlab: deprecated, info-cell, connectivity, boundary-condition
        """
        if self._looped_cell_neighbor_list is not None:
            return self._looped_cell_neighbor_list
        else:
            self._looped_cell_neighbor_list = self._create_looped_cell_neighbor_list()
            return self.looped_neighbors_at_cell

    def _create_looped_cell_neighbor_list(self):
        """Create a list of looped immediate cell neighbors (8 adjacent cells).

        Creates a list of looped immediate cell neighbors (*cell ids*) for each
        cell as a 2D array of size ( self.number_of_cells, 8 ).
        Order or neighbors is [ E, NE, N, NW, W, SW, S, SE ]

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> neighbors = grid._create_looped_cell_neighbor_list()
        >>> neighbors[1]
        array([2, 5, 4, 3, 0, 3, 4, 5])
        >>> neighbors[5]
        array([3, 0, 2, 1, 4, 1, 2, 0])
        """
        # CAUTION: Some terminology concerning cells in this module
        # is asynchronous to general understanding. This is intentionally
        # left as is until further discussion among dev group.
        # Any such instances are marked with (*TC - Terminoly Caution)
        nrows, ncols = self.cell_grid_shape
        interior_cells = sgrid.interior_nodes(self.cell_grid_shape)  # *TC
        cells_at_corners_of_grid = self.cells_at_corners_of_grid  # *TC

        # The cells along the edges minus the corner cells.
        top_edge_cells = self.cell_at_node[self.nodes[-2, :]][2:-2]
        bottom_edge_cells = self.cell_at_node[self.nodes[1, :]][2:-2]
        left_edge_cells = self.cell_at_node[self.nodes[:, 1]][2:-2]
        right_edge_cells = self.cell_at_node[self.nodes[:, -2]][2:-2]

        looped_cell_neighbors = np.empty([self.number_of_cells, 8], dtype=int)

        # order = [E,NE,N,NW,W,SW,S,SE]
        for cell in range(0, self.number_of_cells):
            if cell in interior_cells:
                neighbor_ = [
                    cell + 1,
                    cell + 1 + ncols,
                    cell + ncols,
                    cell + ncols - 1,
                    cell - 1,
                    cell - ncols - 1,
                    cell - ncols,
                    cell - ncols + 1,
                ]
            elif cell in bottom_edge_cells:
                neighbor_ = [
                    cell + 1,
                    cell + 1 + ncols,
                    cell + ncols,
                    cell + ncols - 1,
                    cell - 1,
                    cell + (nrows - 1) * ncols - 1,
                    cell + (nrows - 1) * ncols,
                    cell + (nrows - 1) * ncols + 1,
                ]
            elif cell in top_edge_cells:
                neighbor_ = [
                    cell + 1,
                    cell - (nrows - 1) * ncols + 1,
                    cell - (nrows - 1) * ncols,
                    cell - (nrows - 1) * ncols - 1,
                    cell - 1,
                    cell - ncols - 1,
                    cell - ncols,
                    cell - ncols + 1,
                ]
            elif cell in right_edge_cells:
                neighbor_ = [
                    cell - ncols + 1,
                    cell + 1,
                    cell + ncols,
                    cell + ncols - 1,
                    cell - 1,
                    cell - ncols - 1,
                    cell - ncols,
                    cell - 2 * ncols + 1,
                ]
            elif cell in left_edge_cells:
                neighbor_ = [
                    cell + 1,
                    cell + ncols + 1,
                    cell + ncols,
                    cell + 2 * ncols - 1,
                    cell + ncols - 1,
                    cell - 1,
                    cell - ncols,
                    cell - ncols + 1,
                ]
            elif cell == cells_at_corners_of_grid[0]:  # SW corner
                neighbor_ = [
                    cell + 1,
                    cell + ncols + 1,
                    cell + ncols,
                    cell + 2 * ncols - 1,
                    cell + ncols - 1,
                    cell + nrows * ncols - 1,
                    cell + (nrows - 1) * ncols,
                    cell + (nrows - 1) * ncols + 1,
                ]
            elif cell == cells_at_corners_of_grid[1]:  # SE corner
                neighbor_ = [
                    cell - ncols + 1,
                    cell + 1,
                    cell + ncols,
                    cell + ncols - 1,
                    cell - 1,
                    cell + (nrows - 1) * ncols - 1,
                    cell + (nrows - 1) * ncols,
                    cell + (nrows - 2) * ncols + 1,
                ]
            elif cell == cells_at_corners_of_grid[2]:  # NW corner
                neighbor_ = [
                    cell + 1,
                    cell - (nrows - 1) * ncols + 1,
                    cell - (nrows - 1) * ncols,
                    cell - (nrows - 2) * ncols - 1,
                    cell + ncols - 1,
                    cell - 1,
                    cell - ncols,
                    cell - ncols + 1,
                ]
            elif cell == cells_at_corners_of_grid[3]:  # NE corner
                neighbor_ = [
                    cell - ncols + 1,
                    cell - nrows * ncols + 1,
                    cell - (nrows - 1) * ncols,
                    cell - (nrows - 1) * ncols - 1,
                    cell - 1,
                    cell - ncols - 1,
                    cell - ncols,
                    cell - 2 * ncols + 1,
                ]
            looped_cell_neighbors[cell] = neighbor_

        return looped_cell_neighbors

    @property
    @make_return_array_immutable
    def second_ring_looped_neighbors_at_cell(self):
        """Get list of second ring looped neighbor cell IDs (all 16 neighbors).

        Returns lists of looped second ring neighbor cell IDs of
        given *cell ids*. If *cell ids* are not given, it returns
        a 2D array of size ( self.number_of_cells, 16 ).

        The cells are the 16 which encircle the nine true neighbor cells.
        Order of neighbors: Starts with E and goes counter clockwise

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((10, 10))
        >>> mg.second_ring_looped_neighbors_at_cell[36, :]
        array([38, 46, 54, 53, 52, 51, 50, 42, 34, 26, 18, 19, 20, 21, 22, 30])
        >>> mg.second_ring_looped_neighbors_at_cell[8, :]
        array([10, 18, 26, 25, 24, 31, 30, 22, 14,  6, 62, 63, 56, 57, 58,  2])

        ...take a look at the cell grid to understand why::

            [56, 57, 58, 59, 60, 61, 62, 63]
            [48, 49, 50, 51, 52, 53, 54, 55]
            [40, 41, 42, 43, 44, 45, 46, 47]
            [32, 33, 34, 35, 36, 37, 38, 39]
            [24, 25, 26, 27, 28, 29, 30, 31]
            [16, 17, 18, 19, 20, 21, 22, 23]
            [ 8,  9, 10, 11, 12, 13, 14, 15]
            [ 0,  1,  2,  3,  4,  5,  6,  7]

        :meta landlab: info-cell, connectivity, boundary-condition
        """
        if self._looped_second_ring_cell_neighbor_list_created:
            return self.second_ring_looped_cell_neighbor_list
        else:
            self.second_ring_looped_cell_neighbor_list = (
                self._create_second_ring_looped_cell_neighbor_list()
            )
            return self.second_ring_looped_neighbors_at_cell

    def _create_second_ring_looped_cell_neighbor_list(self):
        """Create list of looped second ring cell neighbors (16 cells).

        Creates a list of looped immediate cell neighbors for each cell
        as a 2D array of size ( self.number_of_cells, 16 ). Order or
        neighbors: Starts with E and goes counter clockwise
        """
        inf = self.looped_neighbors_at_cell
        second_ring = np.empty([self.number_of_cells, 16], dtype=int)
        order = np.arange(-1, 15)
        order[0] = 15
        for cell in range(0, self.number_of_cells):
            cell1, cell2, cell3, cell4 = (
                inf[cell][1],
                inf[cell][3],
                inf[cell][5],
                inf[cell][7],
            )
            ring_tw = np.concatenate(
                (
                    inf[cell1][0:4],
                    inf[cell2][2:6],
                    inf[cell3][4:8],
                    inf[cell4][6:8],
                    inf[cell4][0:2],
                )
            )[order]
            second_ring[cell] = ring_tw

        self._looped_second_ring_cell_neighbor_list_created = True
        return second_ring

    def set_watershed_boundary_condition(
        self,
        node_data,
        nodata_value=-9999.0,
        return_outlet_id=False,
        remove_disconnected=False,
        adjacency_method="D8",
    ):
        """Finds the node adjacent to a boundary node with the smallest value.
        This node is set as the outlet.  The outlet node must have a data
        value.  Can return the outlet id as a one element numpy array if
        return_outlet_id is set to True.

        All nodes with nodata_value are set to ``BC_NODE_IS_CLOSED``
        (grid.status_at_node == 4). All nodes with data values are set to
        ``BC_NODE_IS_CORE`` (grid.status_at_node == 0), with the exception that the
        outlet node is set to a ``BC_NODE_IS_FIXED_GRADIENT`` (grid.status_at_node == 1).

        Note that the outer ring (perimeter) of the raster is set to
        ``BC_NODE_IS_CLOSED``, even if there are nodes that have values. The only
        exception to this would be if the outlet node is on the perimeter, which
        is acceptable.

        This routine assumes that all of the nodata_values are on the outside of
        the data values. In other words, there are no islands of nodata_values
        surrounded by nodes with data.

        This also assumes that the grid has a single watershed (that is a single
        outlet node).

        If you are considering running flow routing algorithms on this model
        grid, it is recommended that you verify the absence of non-closed nodes
        surrounded only by closed nodes. The presence of these isolated non-
        closed nodes may result from clipping a raster to a polygon and will
        prevent the flow router from functioning as intended.

        This can be acomplished by setting the
        parameter: *remove_disconnected* to True (default is False).

        This will run the function:
        set_open_nodes_disconnected_from_watershed_to_closed
        which will find any isolated open nodes that have no neigbors in the
        main watershed and set them to closed. The adjacency method used
        to assess connectivity can be set to either 'D8'(default) or 'D4' using
        the flag *adjacency_method*.

        Finally, the developer has seen cases in which DEM data that has been
        filled results in a different outlet from DEM data which has not been
        filled.  Be aware that if you identify an outlet on a filled DEM, make
        sure that the filled DEM is what is being used for your modeling.
        Otherwise, this may find a different outlet.  To force the outlet
        location, use either set_watershed_boundary_condition_outlet_coords
        or set_watershed_boundary_condition_outlet_id.

        Parameters
        ----------
        node_data : field name or ndarray
            At-node field name or at-node data values to use for identifying
            watershed location.
        nodata_value : float, optional
            Value that indicates an invalid value.
        return_outlet_id : boolean, optional
            Indicates whether or not to return the id of the found outlet
        remove_disconnected : boolean, optional
            Indicates whether to search for and remove disconnected nodes.
        adjacency_method : string, optional. Default is 'D8'.
            Sets the connection method for use if remove_disconnected==True

        Returns
        --------
        outlet_loc : array
            Array of size 1 containing id of outlet location

        Examples
        --------

        The first example will use a 4,4 grid with node data values
        as illustrated::

            -9999. -9999. -9999. -9999.
            -9999.    67.     0. -9999.
            -9999.    67.    67. -9999.
            -9999. -9999. -9999. -9999.

        The second example will use a 4,4 grid with node data values
        as illustrated::

            -9999. -9999. -9999. -9999.
            -9999.    67.     0. -9999.
            -9999.    67.     67.   -2.
            -9999. -9999. -9999. -9999.

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid((4, 4))
        >>> node_data = np.array(
        ...     [
        ...         [-9999.0, -9999.0, -9999.0, -9999.0],
        ...         [-9999.0, 67.0, 67.0, -9999.0],
        ...         [-9999.0, 67.0, 0.0, -9999.0],
        ...         [-9999.0, -9999.0, -9999.0, -9999.0],
        ...     ]
        ... ).flatten()
        >>> out_id = rmg.set_watershed_boundary_condition(node_data, -9999.0, True)
        >>> out_id
        array([10])
        >>> rmg.status_at_node
        array([4, 4, 4, 4, 4, 0, 0, 4, 4, 0, 1, 4, 4, 4, 4, 4], dtype=uint8)
        >>> rmg2 = RasterModelGrid((4, 4))
        >>> node_data2 = np.array(
        ...     [
        ...         [-9999.0, -9999.0, -9999.0, -9999.0],
        ...         [-9999.0, 67.0, 67.0, -2.0],
        ...         [-9999.0, 67.0, 0.0, -9999.0],
        ...         [-9999.0, -9999.0, -9999.0, -9999.0],
        ...     ]
        ... ).flatten()
        >>> rmg2.set_watershed_boundary_condition(node_data2, -9999.0)
        >>> rmg2.status_at_node
        array([4, 4, 4, 4, 4, 0, 0, 1, 4, 0, 0, 4, 4, 4, 4, 4], dtype=uint8)

        The node data can also be provided as a model grid field.

        >>> rmg = RasterModelGrid((4, 4))
        >>> node_data = [
        ...     [-9999.0, -9999.0, -9999.0, -9999.0],
        ...     [-9999.0, 67.0, 67.0, -9999.0],
        ...     [-9999.0, 67.0, 0.0, -9999.0],
        ...     [-9999.0, -9999.0, -9999.0, -9999.0],
        ... ]
        >>> _ = rmg.add_field("topographic__elevation", node_data, at="node")
        >>> out_id = rmg.set_watershed_boundary_condition(
        ...     "topographic__elevation", -9999.0, True
        ... )
        >>> out_id
        array([10])
        >>> rmg.status_at_node
        array([4, 4, 4, 4, 4, 0, 0, 4, 4, 0, 1, 4, 4, 4, 4, 4], dtype=uint8)

        :meta landlab: boundary-condition
        """
        # get node_data if a field name
        node_data = self.return_array_or_field_values("node", node_data)

        # For this to be a watershed, need to make sure that there is a ring
        # of closed boundary nodes around the outside of the watershed,
        # barring the outlet location.  So enforce that all perimeter nodes
        # are inactive boundaries now, then set the outlet location later.
        # By enforcing the perimeter of closed values first, then fixing the
        # outlet later, it should be OK if the outlet is on the perimeter.
        self.set_closed_boundaries_at_grid_edges(True, True, True, True)

        # Set no data nodes to inactive boundaries.
        # This may be redundant, but must do in case there are no data
        # values that are not on the perimeter.
        self.set_nodata_nodes_to_closed(node_data, nodata_value)

        # need to find values that are not no_data

        # locs is a list that contains locations where
        # node data is not equal to the nodata value
        locs = np.where(node_data != nodata_value)
        if len(locs) < 1:
            raise ValueError("All data values are no_data values")

        # now find minimum of the data values
        min_val = np.min(node_data[locs])

        # now find where minimum values are
        min_locs = np.where(node_data == min_val)[0]

        # check all the locations with the minimum value to see if one
        # is adjacent to a boundary location.  If so, that will be the
        # watershed outlet.  If none of these points qualify, then
        # increase the minimum value and check again.  Keep checking
        # until a point next to the boundary is found.
        # NG I think the only way this would become an infinite loop
        # is if there are no interior nodes.  Should be checking for
        # this above.
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
                    potential_locs = min_locs[np.where(np.asarray(next_to_boundary))[0]]
                    raise ValueError(
                        "Grid has two potential outlet nodes."
                        "They have the following node IDs: \n"
                        + str(potential_locs)
                        + "\nUse the method set_watershed_boundary_condition_outlet_id "
                        "to explicitly select one of these "
                        "IDs as the outlet node."
                    )
                else:
                    outlet_loc = min_locs[np.where(next_to_boundary)[0][0]]

            # checked all of the min vals, (so done with inner while)
            # and none of the min values were outlet candidates
            if local_not_found:
                # need to find the next largest minimum value
                # first find the locations of all values greater
                # than the old minimum
                # not done with outer while
                locs = np.where((node_data > min_val) & (node_data != nodata_value))
                # now find new minimum of these values
                min_val = np.min(node_data[locs])
                min_locs = np.where(node_data == min_val)[0]
            else:
                # if locally found, it is also globally found
                # so done with outer while
                not_found = False

        # set outlet boundary condition
        self.status_at_node[outlet_loc] = NodeStatus.FIXED_VALUE

        if remove_disconnected:
            self.set_open_nodes_disconnected_from_watershed_to_closed(
                node_data=node_data,
                outlet_id=outlet_loc,
                nodata_value=nodata_value,
                adjacency_method=adjacency_method,
            )
        if return_outlet_id:
            return as_id_array(np.array([outlet_loc]))

    def set_open_nodes_disconnected_from_watershed_to_closed(
        self, node_data, outlet_id=None, nodata_value=-9999.0, adjacency_method="D8"
    ):
        """Identifys all non-closed nodes that are disconnected from the node
        given in.

        *outlet_id* and sets them as closed.

        If *outlet_id* is not given, the outlet will be identified as the node
        for which the status at the node is ``BC_NODE_IS_FIXED_VALUE``. If more than
        one node has this value, the algorithm will fail.

        If *outlet_id* is given, the algorithm will check that it is not a node
        with status of ``BC_NODE_IS_CLOSED``.

        The method supports both D4 and D8 (default) neighborhood evaluation in
        determining if a node is connected. This can be modified with the flag
        *adjacency_method*.

        This function can be run directly, or by setting the flag
        remove_disconnected to True in set_watershed_boundary_condition

        Parameters
        ----------
        node_data : field name or ndarray
            At-node field name or at-node data values to use for identifying
            watershed location.
        outlet_id : one element numpy array, optional.
            The node ID of the outlet that all open nodes must be connected to.
            If a node ID is provided, it does not need have the status
            ``BC_NODE_IS_FIXED_VALUE``. However, it must not have the status of
            ``BC_NODE_IS_CLOSED``.
        nodata_value : float, optional, default is -9999.
            Value that indicates an invalid value.
        adjacency_method : string, optional. Default is 'D8'.
            Sets the connection method.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> mg1 = RasterModelGrid((4, 6))
        >>> z1 = np.array(
        ...     [
        ...         [-9999.0, -9999.0, -9999.0, -9999.0, -9999.0, -9999.0],
        ...         [-9999.0, 67.0, 67.0, -9999.0, 50.0, -9999.0],
        ...         [-9999.0, 67.0, 0.0, -9999.0, -9999.0, -9999.0],
        ...         [-9999.0, -9999.0, -9999.0, -9999.0, -9999.0, -9999.0],
        ...     ]
        ... ).flatten()
        >>> mg2 = RasterModelGrid((4, 6))
        >>> z2 = np.array(
        ...     [
        ...         [-9999.0, -9999.0, -9999.0, -9999.0, -9999.0, -9999.0],
        ...         [-9999.0, 67.0, 67.0, -9999.0, 50.0, -9999.0],
        ...         [-9999.0, 67.0, 0.0, -9999.0, -9999.0, -9999.0],
        ...         [-9999.0, -9999.0, -9999.0, -9999.0, -9999.0, -9999.0],
        ...     ]
        ... ).flatten()
        >>> mg1.set_watershed_boundary_condition(z1, remove_disconnected=True)
        >>> mg2.set_watershed_boundary_condition(z2)
        >>> mg2.status_at_node.reshape(mg2.shape)
        array([[4, 4, 4, 4, 4, 4],
               [4, 0, 0, 4, 0, 4],
               [4, 0, 1, 4, 4, 4],
               [4, 4, 4, 4, 4, 4]], dtype=uint8)
        >>> mg2.set_open_nodes_disconnected_from_watershed_to_closed(z2)
        >>> np.allclose(mg1.status_at_node, mg2.status_at_node)
        True
        >>> np.allclose(z1, z2)
        True
        >>> mg2.status_at_node.reshape(mg2.shape)
        array([[4, 4, 4, 4, 4, 4],
               [4, 0, 0, 4, 4, 4],
               [4, 0, 1, 4, 4, 4],
               [4, 4, 4, 4, 4, 4]], dtype=uint8)
        >>> z1.reshape(mg1.shape)
        array([[-9999., -9999., -9999., -9999., -9999., -9999.],
               [-9999.,    67.,    67., -9999., -9999., -9999.],
               [-9999.,    67.,     0., -9999., -9999., -9999.],
               [-9999., -9999., -9999., -9999., -9999., -9999.]])

        :meta landlab: boundary-condition
        """
        # get node_data if a field name
        node_data = self.return_array_or_field_values("node", node_data)

        if outlet_id is None:
            # verify that there is one and only one node with the status
            # BC_NODE_IS_FIXED_VALUE.
            outlet_id = np.where(self.status_at_node == NodeStatus.FIXED_VALUE)[0]
            n_outlets = len(outlet_id)

            if n_outlets != 1:
                raise ValueError(
                    f"Incorrect number outlets ({n_outlets}). Grid must have one,"
                    " and only one, node with node status of BC_NODE_IS_FIXED_VALUE."
                )

            outlet_id = outlet_id[0]
        else:
            # check that the node status at the node given by outlet_id is not
            # BC_NODE_IS_CLOSED
            if self.status_at_node[outlet_id] == self.BC_NODE_IS_CLOSED:
                raise ValueError(
                    "The node given by outlet_id must not have the status: BC_NODE_IS_CLOSED"
                )

        if adjacency_method not in {"D4", "D8"}:
            raise ValueError(
                f"Method must be one of 'D4' or 'D8' (got {adjacency_method!r})"
            )

        # begin main code portion.
        # initialize list of core nodes and new nodes
        connected_nodes = [outlet_id]
        newNodes = connected_nodes

        # keep track of the number of nodes added in the previous itteration.
        numAdded = len(newNodes)

        # continue running until no new nodes are added.
        while numAdded > 0:
            # find all potential new nodes by filtering the nodes connected to
            # the most recent set of new nodes based on their status.
            connected_orthogonal_nodes = self.adjacent_nodes_at_node[newNodes]
            potentialNewNodes = list(
                connected_orthogonal_nodes[
                    self.status_at_node[connected_orthogonal_nodes]
                    != self.BC_NODE_IS_CLOSED
                ]
            )

            # if method is D8 (default), add the diagonal nodes.
            if adjacency_method == "D8":
                connected_diagonal_nodes = self.diagonal_adjacent_nodes_at_node[
                    newNodes
                ]
                potentialNewNodes.extend(
                    connected_diagonal_nodes[
                        self.status_at_node[connected_diagonal_nodes]
                        != self.BC_NODE_IS_CLOSED
                    ]
                )

            # filter new nodes further based on if they are already present in
            # the connected node list
            newNodes = list(set(potentialNewNodes) - set(connected_nodes))
            connected_nodes.extend(newNodes)

            # update number added, when this is zero, the loop will end
            numAdded = len(newNodes)

        # create an array that identifies the nodes that should be closed
        # of closed boundary nodes
        not_connected = np.array((0 * self.status_at_node) + 1)
        not_connected[np.array(connected_nodes)] = 0

        # identify those nodes that should be closed, but are not yet closed.
        is_not_connected_to_outlet = (self.status_at_node != self.BC_NODE_IS_CLOSED) & (
            not_connected == 1
        )

        # modify the node_data array to set those that are disconnected
        # to the no data value.
        node_data[is_not_connected_to_outlet] = nodata_value  #

        # finally update the status of the nodes based on the modified node_data.
        self.set_nodata_nodes_to_closed(node_data, nodata_value)

    def set_watershed_boundary_condition_outlet_coords(
        self, outlet_coords, node_data, nodata_value=-9999.0
    ):
        """Set the boundary conditions for a watershed. All nodes with
        nodata_value are set to ``BC_NODE_IS_CLOSED`` (grid.status_at_node == 4). All
        nodes with data values are set to CORE_NODES (grid.status_at_node ==
        0), with the exception that the outlet node is set to a
        BC_NODE_IS_FIXED_VALUE (grid.status_at_node == 1).

        Note that the outer ring of the raster is set to ``BC_NODE_IS_CLOSED``, even
        if there are nodes that have values.  The only exception to this would
        be if the outlet node is on the boundary, which is acceptable.

        Assumes that outlet is already known.

        This assumes that the grid has a single watershed.  If this is not
        the case this will not work.

        This must be passed the values of the outlet_row and outlet_column.
        Also takes node_data and optionally, nodata_value.

        Parameters
        ----------
        outlet_coords : list - two integer values
            row, column of outlet, NOT THE ABSOLUTE X AND Y LOCATIONS
        node_data : field name or ndarray
            At-node field name or at-node data values to use for identifying
            watershed location.
        nodata_value : float, optional
            Value that indicates an invalid value.

        Examples
        --------
        The example will use a 4,4 grid with node data values
        as illustrated::

            -9999. -9999. -9999. -9999.
            -9999.    67.     0. -9999.
            -9999.    67.    67. -9999.
            -9999. -9999. -9999. -9999.

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid((4, 4))
        >>> rmg.status_at_node
        array([1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1], dtype=uint8)
        >>> node_data = np.array(
        ...     [
        ...         [-9999.0, -9999.0, -9999.0, -9999.0],
        ...         [-9999.0, 67.0, 67.0, -9999.0],
        ...         [-9999.0, 67.0, 0.0, -9999.0],
        ...         [-9999.0, -9999.0, -9999.0, -9999.0],
        ...     ]
        ... ).flatten()
        >>> rmg.set_watershed_boundary_condition_outlet_coords(
        ...     (2, 2), node_data, -9999.0
        ... )
        >>> rmg.status_at_node.reshape(rmg.shape)
        array([[4, 4, 4, 4],
               [4, 0, 0, 4],
               [4, 0, 1, 4],
               [4, 4, 4, 4]], dtype=uint8)

        :meta landlab: boundary-condition
        """
        # get node_data if a field name
        node_data = self.return_array_or_field_values("node", node_data)

        # make ring of no data nodes
        self.set_closed_boundaries_at_grid_edges(True, True, True, True)

        # set no data nodes to inactive boundaries
        self.set_nodata_nodes_to_closed(node_data, nodata_value)

        # find the id of the outlet node
        outlet_node = self.grid_coords_to_node_id(outlet_coords[0], outlet_coords[1])
        # set the boundary condition (fixed value) at the outlet_node
        self.status_at_node[outlet_node] = NodeStatus.FIXED_VALUE

    def set_watershed_boundary_condition_outlet_id(
        self, outlet_id, node_data, nodata_value=-9999.0
    ):
        """Set the boundary conditions for a watershed. All nodes with
        nodata_value are set to ``BC_NODE_IS_CLOSED`` (4). All nodes with data values
        are set to ``BC_NODE_IS_CORE`` (0), with the exception that the outlet node is
        set to a ``BC_NODE_IS_FIXED_VALUE`` (1).

        Note that the outer ring of the raster is set to ``BC_NODE_IS_CLOSED``, even
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
        The example will use a 4,4 grid with node data values
        as illustrated:

            -9999. -9999. -9999. -9999.
            -9999.    67.     0. -9999.
            -9999.    67.    67. -9999.
            -9999. -9999. -9999. -9999.

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid((4, 4))
        >>> rmg.status_at_node.reshape(rmg.shape)
        array([[1, 1, 1, 1],
               [1, 0, 0, 1],
               [1, 0, 0, 1],
               [1, 1, 1, 1]], dtype=uint8)
        >>> node_data = np.array(
        ...     [
        ...         [-9999.0, -9999.0, -9999.0, -9999.0],
        ...         [-9999.0, 67.0, 67.0, -9999.0],
        ...         [-9999.0, 67.0, 0.0, -9999.0],
        ...         [-9999.0, -9999.0, -9999.0, -9999.0],
        ...     ]
        ... ).flatten()
        >>> outlet = rmg.set_watershed_boundary_condition_outlet_id(
        ...     10, node_data, -9999.0
        ... )
        >>> rmg.status_at_node.reshape(rmg.shape)
        array([[4, 4, 4, 4],
               [4, 0, 0, 4],
               [4, 0, 1, 4],
               [4, 4, 4, 4]], dtype=uint8)

        :meta landlab: boundary-condition
        """
        # get node_data if a field name
        node_data = self.return_array_or_field_values("node", node_data)

        # make ring of no data nodes
        self.set_closed_boundaries_at_grid_edges(True, True, True, True)

        # set no data nodes to inactive boundaries
        self.set_nodata_nodes_to_closed(node_data, nodata_value)

        # set the boundary condition (fixed value) at the outlet_node
        self.status_at_node[outlet_id] = NodeStatus.FIXED_VALUE


def _guess_format_from_name(path):
    """Get file format by name.

    Parameters
    ----------
    path : str
        Path to file.

    Returns
    -------
    str
        File format as a string.
    """
    import os

    fname = os.path.basename(path)

    if fname.endswith(".nc"):
        return "netcdf"
    elif fname.endswith(".asc"):
        return "esri-ascii"
    else:
        return None


def _add_format_extension(path, format):
    """Add format extension to a file name.

    Parameters
    ----------
    path : str
        File name.
    format : str
        File format.

    Returns
    -------
    str
        File name with the file-format extension added.
    """
    import os

    (base, ext) = os.path.splitext(path)
    if format == "netcdf":
        ext = ".nc"
    elif format == "esri-ascii":
        ext = ".asc"
    return base + ext


add_module_functions_to_class(RasterModelGrid, "raster_mappers.py", pattern="map_*")
add_module_functions_to_class(RasterModelGrid, "raster_gradients.py", pattern="calc_*")
add_module_functions_to_class(RasterModelGrid, "raster_divergence.py", pattern="calc_*")
add_module_functions_to_class(
    RasterModelGrid, "raster_set_status.py", pattern="set_status_at_node*"
)
