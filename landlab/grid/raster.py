#! /usr/env/python
"""
A class used to create and manage regular raster grids for 2D numerical models
in Landlab.

Do NOT add new documentation here. Grid documentation is now built in a semi-
automated fashion. To modify the text seen on the web, edit the files
`docs/text_for_[gridfile].py.txt`.
"""
from __future__ import absolute_import

from warnings import warn

import numpy as np
import six
from six.moves import range

from landlab.field.scalar_data_fields import FieldError
from landlab.grid.structured_quad import (
    cells as squad_cells,
    faces as squad_faces,
    links as squad_links,
)
from landlab.utils import structured_grid as sgrid
from landlab.utils.decorators import deprecated, make_return_array_immutable

from ..core.utils import add_module_functions_to_class, as_id_array
from ..io import write_esri_ascii
from ..io.netcdf import write_netcdf
from ..utils.decorators import cache_result_in_object
from . import raster_funcs as rfuncs
from .base import (
    BAD_INDEX_VALUE,
    CLOSED_BOUNDARY,
    CORE_NODE,
    FIXED_VALUE_BOUNDARY,
    LOOPED_BOUNDARY,
    ModelGrid,
)
from .decorators import return_id_array, return_readonly_id_array
from .diagonals import DiagonalsMixIn
from .warnings import (
    DeprecatedDxKeyword,
    DeprecatedOriginKeyword,
    DeprecatedRowsColsArguments,
    DeprecatedRowsColsKeywords,
    DeprecatedSpacingArgument,
    DeprecatedSpacingKeyword,
)


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
            if mg.status_at_node[neighbor] != CORE_NODE:
                return True
        except IndexError:
            return True
    if method == "d8":
        for neighbor in mg.diagonal_adjacent_nodes_at_node[id]:
            try:
                if mg.status_at_node[neighbor] != CORE_NODE:
                    return True
            except IndexError:
                return True
    return False


def _make_arg_into_array(arg):
    """Make an argument into an iterable.

    This function tests if the provided object is a Python list or a numpy
    array. If not, attempts to cast the object to a list. If it cannot, it will
    raise a TypeError.

    Parameters
    ----------
    arg : array_like
        Input array.

    Returns
    -------
    array_like
        The input array converted to an iterable.

    Examples
    --------
    >>> from landlab.grid.raster import _make_arg_into_array
    >>> _make_arg_into_array(1)
    [1]
    >>> _make_arg_into_array((1, ))
    [1]
    >>> _make_arg_into_array([1, 2])
    [1, 2]
    >>> import numpy as np
    >>> _make_arg_into_array(np.arange(3))
    array([0, 1, 2])
    """
    ids = arg
    if not isinstance(ids, list) and not isinstance(ids, np.ndarray):
        try:
            ids = list(ids)
        except TypeError:
            ids = [ids]
    return ids


_node_has_boundary_neighbor = np.vectorize(_node_has_boundary_neighbor, excluded=["mg"])


class RasterModelGridPlotter(object):

    """MixIn that provides plotting functionality.

    Inhert from this class to provide a ModelDataFields object with the
    method function, ``imshow``, that plots a data field.
    """

    def imshow(self, group, var_name, **kwds):
        """Plot a data field.

        This is a wrapper for `plot.imshow_grid`, and can take the same
        keywords. See that function for full documentation.

        Parameters
        ----------
        group : str
            Name of group.
        var_name : str
            Name of field

        See Also
        --------
        landlab.plot.imshow_grid

        LLCATS: GINF
        """
        from landlab.plot import imshow_grid

        kwds["values_at"] = group
        imshow_grid(self, var_name, **kwds)


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
    >>> grid_edge_is_closed_from_dict(dict(bottom='closed', top='open'))
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


def _old_style_args(args):
    """Test if arguments are the old-style RasterModelGrid __init__ method.

    The old way of initializing a :any:`RasterModelGrid` was like,

    .. code::
        grid = RasterModelGrid(n_rows, n_cols)

    The new way passes the grid shape as a tuple, like numpy functions,

    .. code::
        grid = RasterModelGrid((n_rows, n_cols))

    Parameters
    ----------
    args : iterable
        Arguments to a function.

    Examples
    --------
    >>> from landlab.grid.raster import _old_style_args
    >>> _old_style_args((4, 5))
    True
    >>> _old_style_args(((4, 5), ))
    False
    >>> _old_style_args(([4, 5], ))
    False
    """
    return len(args) in (2, 3) and isinstance(args[0], int)


def _parse_grid_shape_from_args(args):
    """Get grid shape from args.

    Parameters
    ----------
    args : iterable
        Arguments to a function.

    Examples
    --------
    >>> from landlab.grid.raster import _parse_grid_shape_from_args
    >>> _parse_grid_shape_from_args((3, 4))
    (3, 4)
    >>> _parse_grid_shape_from_args(((3, 4), ))
    (3, 4)
    """
    if _old_style_args(args):
        rows, cols = args[0], args[1]
    else:
        try:
            (rows, cols) = args[0]
        except ValueError:
            raise ValueError("grid shape must be tuple")
    return rows, cols


def _parse_grid_spacing_from_args(args):
    """Get grid spacing from args.

    Parameters
    ----------
    args : iterable
        Arguments to a function.

    Examples
    --------
    >>> from landlab.grid.raster import _parse_grid_spacing_from_args
    >>> _parse_grid_spacing_from_args((3, 4, 5))
    5
    >>> _parse_grid_spacing_from_args(((3, 4), 5))
    5
    """
    try:
        if _old_style_args(args):
            return args[2]
        else:
            return args[1]
    except IndexError:
        return None


class RasterModelGrid(DiagonalsMixIn, ModelGrid, RasterModelGridPlotter):

    """A 2D uniform rectilinear grid.

    Create a uniform rectilinear grid that has *num_rows* and *num_cols*
    of grid nodes, with a row and column spacing of *dx*.

    Use the *bc* keyword to specify boundary_conditions along the edge nodes
    of the grid. *bc* is a dict whose keys indicate edge location (as
    "bottom", "left", "top", "right") and values must be one of "open", or
    "closed". If an edge location key is missing, that edge is assumed to be
    *open*.

    *Deprecation Warning*: The keywords *spacing* and *origin* have been
    deprecated as of Landlab v1.5.5. They will be removed in v2.0.

    Parameters
    ----------
    shape : tuple of int
        Shape of the grid in nodes.
    xy_spacing : tuple or float, optional
        dx and dy spacing. Either provided as a float or a
        (dx, dy) tuple. Default is (1., 1.)
    xy_of_lower_left: tuple, optional
        (x, y) coordinates of the lower left corner. Default
        is (0., 0.)
    bc : dict, optional
        Edge boundary conditions.

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

    >>> rmg = RasterModelGrid((4, 5), bc={'top': 'closed'})
    >>> rmg.number_of_node_rows, rmg.number_of_node_columns
    (4, 5)
    >>> rmg.number_of_active_links
    14

    A `RasterModelGrid` can have different node spacings in the *x* and *y*
    directions.

    >>> grid = RasterModelGrid((4, 5), xy_spacing=(2, 1))
    >>> grid.dx, grid.dy
    (2.0, 1.0)
    >>> grid.node_y # doctest: +NORMALIZE_WHITESPACE
    array([ 0., 0., 0., 0., 0.,
            1., 1., 1., 1., 1.,
            2., 2., 2., 2., 2.,
            3., 3., 3., 3., 3.])
    >>> grid.node_x # doctest: +NORMALIZE_WHITESPACE
    array([ 0., 2., 4., 6., 8.,
            0., 2., 4., 6., 8.,
            0., 2., 4., 6., 8.,
            0., 2., 4., 6., 8.])

    Notes
    -----
    The option for NOT giving rows, cols, and dx no longer works,
    because the *field* init requires num_active_cells, etc., to be
    defined. Either we force users to give arguments on instantiation,
    or set it up such that one can create a zero-node grid.
    """

    @staticmethod
    def _parse_parameters(args, kwds):
        """Parse __init__ parameters into new-style parameters."""
        old = (args, kwds)
        deprecation_warning = None
        try:
            shape = kwds.pop("num_rows"), kwds.pop("num_cols")
        except KeyError:
            shape = _parse_grid_shape_from_args(args)
            if len(args) >= 2:
                deprecation_warning = DeprecatedRowsColsArguments
        else:
            deprecation_warning = DeprecatedRowsColsKeywords
            if args:
                raise ValueError(
                    "number of args must be 0 when using keywords for grid shape"
                )

        new_kwds = kwds.copy()
        new = ((shape,), new_kwds)

        xy_spacing = None
        if "dx" in kwds:
            deprecation_warning = DeprecatedDxKeyword
            xy_spacing = new_kwds.pop("dx")
        if "spacing" in kwds:
            deprecation_warning = DeprecatedSpacingKeyword
            xy_spacing = new_kwds.pop("spacing")
        if _parse_grid_spacing_from_args(args) is not None:
            deprecation_warning = DeprecatedSpacingArgument
            xy_spacing = _parse_grid_spacing_from_args(args)

        if "xy_spacing" not in kwds and xy_spacing is not None:
            try:
                new_kwds.setdefault("xy_spacing", xy_spacing[::-1])
            except TypeError:
                new_kwds.setdefault("xy_spacing", xy_spacing)

        if "origin" in kwds:
            deprecation_warning = DeprecatedOriginKeyword
            new_kwds.setdefault("xy_of_lower_left", new_kwds.pop("origin"))

        if deprecation_warning:
            warn(deprecation_warning("RasterModelGrid", old=old, new=new))

        return new

    def __init__(self, *args, **kwds):
        """Create a 2D grid with equal spacing.

        Optionally takes numbers of rows and columns and cell size as
        inputs. If this are given, calls initialize() to set up the grid.
        At the moment, num_rows and num_cols MUST be specified. Both must be
        >=3 to allow correct automated setup of boundary conditions.

        *Deprecation Warning*: The keywords *spacing* and *origin* have been
        deprecated as of Landlab v1.5.5. They will be removed in v2.0.
        Additionally, passing spacing as an argument has been deprecated.

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
        args, kwds = self._parse_parameters(args, kwds)

        shape = args[0]
        xy_spacing = np.broadcast_to(kwds.get("xy_spacing", 1.0), 2)
        xy_of_lower_left = kwds.get("xy_of_lower_left", (0.0, 0.0))

        if shape[0] <= 0 or shape[1] <= 0:
            raise ValueError("number of rows and columns must be positive")

        if len(xy_of_lower_left) != 2:
            raise ValueError("xy_of_lower_left must be size two")

        self._xy_of_lower_left = tuple(xy_of_lower_left)

        # Node Status
        self._node_status = np.empty(shape[0] * shape[1], dtype=np.uint8)

        # Set number of nodes, and initialize if caller has given dimensions
        self._initialize(
            shape[0], shape[1], xy_spacing, (xy_of_lower_left[0], xy_of_lower_left[1])
        )

        self.set_closed_boundaries_at_grid_edges(
            *grid_edge_is_closed_from_dict(kwds.pop("bc", {}))
        )

        super(RasterModelGrid, self).__init__(**kwds)

        self.looped_node_properties = {}

    def __repr__(self):
        return "RasterModelGrid({0}, xy_spacing={1}, xy_of_lower_left={2})".format(
            repr(self.shape), repr((self.dx, self.dy)), repr(self.xy_of_lower_left)
        )

    def __setstate__(self, state_dict):
        """Set state for of RasterModelGrid from pickled state_dict."""
        if state_dict["type"] != "RasterModelGrid":
            assert TypeError(("Saved model instance not of " "RasterModelGrid type."))

        dx = state_dict["dx"]
        dy = state_dict["dy"]
        shape = state_dict["shape"]
        origin = state_dict["origin"]
        num_rows = shape[0]
        num_cols = shape[1]

        self._node_status = np.empty(num_rows * num_cols, dtype=np.uint8)

        # Set number of nodes, and initialize if caller has given dimensions
        self._initialize(num_rows, num_cols, (dx, dy), origin)

        super(RasterModelGrid, self).__init__()

        self.looped_node_properties = {}

        # Recreate the state of the grid and the information it new
        # about itself

        # If angle of links existed, create them
        if state_dict["_angle_of_link_created"]:
            self._create_angle_of_link()

        # If patches existed, create them
        if state_dict["_patches_created"]:
            self.nodes_at_patch
            self.links_at_patch

        # If forced cell area existed.
        if state_dict["forced_cell_areas_created"]:
            self._create_cell_areas_array_force_inactive()

        # If neighbor list existed, create them
        if state_dict["neighbor_list_created"]:
            self._create_neighbor_list()

        # Set status at links and nodes
        self.status_at_node[:] = state_dict["status_at_node"]

        # Add fields back to the grid
        for group in state_dict["_groups"].keys():
            field_set = state_dict["_groups"][group]
            fields = field_set.keys()
            for field in fields:
                # add field to
                self.add_field(
                    group,
                    field,
                    field_set[field]["value_array"],
                    units=field_set[field]["units"],
                )
        self.bc_set_code = state_dict["bc_set_code"]

    def __getstate__(self):
        """Get state for pickling."""
        # initialize state_dict
        state_dict = {}

        # save basic information about the shape and size of the grid
        state_dict["type"] = "RasterModelGrid"
        state_dict["dx"] = self.dx
        state_dict["dy"] = self.dy
        state_dict["shape"] = self.shape
        state_dict["origin"] = (self.node_x[0], self.node_y[0])
        state_dict["_axis_name"] = self._axis_name
        state_dict["_axis_units"] = self._axis_units
        state_dict["_default_group"] = self._default_group

        # save information about things that might have been created
        state_dict["_angle_of_link_created"] = self._angle_of_link_created
        state_dict["_patches_created"] = self._patches_created
        state_dict["neighbor_list_created"] = self.neighbor_list_created
        try:
            if isinstance(self._forced_cell_areas, np.ndarray):
                state_dict["forced_cell_areas_created"] = True
            else:
                state_dict["forced_cell_areas_created"] = False
        except AttributeError:
            state_dict["forced_cell_areas_created"] = False

        # save status information at nodes (status at link set based on status
        # at node
        state_dict["status_at_node"] = np.asarray(self._node_status)

        # save all fields. This is the key part, since saving ScalarDataFields, breaks
        # pickle and/or dill

        groups = {}
        for group in self._groups.keys():
            field_set_dict = {}
            fields = self._groups[group].keys()
            for field in fields:
                field_dict = {}
                field_dict["value_array"] = self.field_values(group, field)
                field_dict["units"] = self.field_units(group, field)
                field_set_dict[field] = field_dict
            groups[group] = field_set_dict

        state_dict["_groups"] = groups

        # save BC set code
        state_dict["bc_set_code"] = self.bc_set_code

        # return state_dict
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
        >>> grid = RasterModelGrid.from_dict(
        ...     {'shape': (3, 4), 'bc': {'top': 'closed'}})
        >>> grid.number_of_nodes
        12

        LLCATS: GINF
        """
        shape = params.pop("shape", None)
        return cls(shape, **params)

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

    def _initialize(self, num_rows, num_cols, xy_spacing, xy_of_lower_left):
        """Set up a raster grid.

        Sets up a *num_rows* by *num_cols* grid with cell *spacing*  and
        (by default) regular boundaries (that is, all perimeter cells are
        boundaries and all interior cells are active).

        Spacing may be provided as a value for square cells or as a tuple
        (dx, dy) for rectangular cells.

        The lower left corner is set through *xy_of_lower_left*, which is a
        (lower left corner x, lower left corner y) tuple.

        To be consistent with unstructured grids, the raster grid is
        managed not as a 2D array but rather as a set of vectors that
        describe connectivity information between nodes, links, active links,
        cells, active cells, faces, patches, junctions, and corners.

        By default, all interior nodes are set to active, and all perimeter
        nodes are set as fixed value, open boundaries (type 1, see supporting
        documentation).

        Note that by default, a RasterModelGrid ONLY has links to
        orthogonal neighboring nodes. However, if you wish to work with the
        diagonal links (e.g., D8 flow routing), these functions are available
        as methods, and the diagonal links can readily be created after
        initialization.

        Examples
        --------
        >>> from landlab import RasterModelGrid, BAD_INDEX_VALUE
        >>> numrows = 20          # number of rows in the grid
        >>> numcols = 30          # number of columns in the grid
        >>> dx = 10.0             # grid cell spacing
        >>> rmg = RasterModelGrid((numrows, numcols), xy_spacing=dx)
        >>> (rmg.number_of_nodes, rmg.number_of_cells, rmg.number_of_links,
        ...  rmg.number_of_active_links)
        (600, 504, 1150, 1054)
        >>> rmg = RasterModelGrid((4, 5))
        >>> (rmg.number_of_nodes, rmg.number_of_cells, rmg.number_of_links,
        ...  rmg.number_of_active_links)
        (20, 6, 31, 17)
        >>> rmg.status_at_node # doctest: +NORMALIZE_WHITESPACE
        array([1, 1, 1, 1, 1,
               1, 0, 0, 0, 1,
               1, 0, 0, 0, 1,
               1, 1, 1, 1, 1], dtype=uint8)
        >>> rmg.node_at_cell # doctest: +NORMALIZE_WHITESPACE
        array([ 6,  7,  8,
               11, 12, 13])
        >>> rmg.node_at_link_tail # doctest: +NORMALIZE_WHITESPACE
        array([ 0,  1,  2,  3,  0,  1,  2,  3,  4,  5,  6,  7,  8,  5,  6,  7,
                8,  9, 10, 11, 12, 13, 10, 11, 12, 13, 14, 15, 16, 17, 18])
        >>> rmg.node_at_link_head # doctest: +NORMALIZE_WHITESPACE
        array([ 1,  2,  3,  4,  5,  6,  7,  8,  9,  6,  7,  8,  9, 10, 11, 12,
               13, 14, 11, 12, 13, 14, 15, 16, 17, 18, 19, 16, 17, 18, 19])
        >>> rmg.face_at_link[20]
        12
        >>> rmg.active_links # doctest: +NORMALIZE_WHITESPACE
        array([ 5,  6,  7,  9, 10, 11, 12, 14, 15, 16, 18, 19, 20, 21, 23, 24,
               25])
        """
        if isinstance(xy_spacing, (float, int)):
            xy_spacing = (xy_spacing, xy_spacing)

        # Basic info about raster size and shape
        self._nrows = num_rows
        self._ncols = num_cols

        self._dx, self._dy = float(xy_spacing[0]), float(xy_spacing[1])
        self.cellarea = self._dy * self._dx

        self._node_at_cell = sgrid.node_at_cell(self.shape)
        self._cell_at_node = squad_cells.cell_id_at_nodes(self.shape).reshape((-1,))

        # We need at least one row or column of boundary cells on each
        # side, so the grid has to be at least 3x3
        assert np.min((num_rows, num_cols)) >= 3

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
        (x_of_node, y_of_node) = sgrid.node_coords(
            (num_rows, num_cols),
            (self._dy, self._dx),
            (xy_of_lower_left[1], xy_of_lower_left[0]),
        )

        self._xy_of_node = np.hstack(
            (x_of_node.reshape((-1, 1)), y_of_node.reshape((-1, 1)))
        )

        # Node boundary/active status:
        # Next, we set up an array of "node status" values, which indicate
        # whether a given node is an active, non-boundary node, or some type of
        # boundary. Here we default to having all perimeter nodes be active
        # fixed-value boundaries.
        self._node_status[:] = sgrid.status_at_node(
            self.shape, boundary_status=FIXED_VALUE_BOUNDARY
        )

        # Cell lists:
        # For all cells, we create a list of the corresponding node ID for
        # each cell.
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

        self._links_at_node = squad_links.links_at_node(self.shape)

        # Link lists:
        # For all links, we encode the "tail" and "head" nodes, and the face
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
        # 22      23      24      25      26
        #  |       |       |       |       |
        #  *--18-->*--19-->*--20-->*--21-->*
        #  ^       ^       ^       ^       ^
        # 13      14      15      16      17
        #  |       |       |       |       |
        #  *---9-->*--10-->*--11-->*--12-->*
        #  ^       ^       ^       ^       ^
        #  4       5       6       7       8
        #  |       |       |       |       |
        #  *---0-->*---1-->*---2-->*---3-->*
        #
        #   create the tail-node and head-node lists
        (node_at_link_tail, node_at_link_head) = sgrid.node_index_at_link_ends(
            self.shape
        )
        self._nodes_at_link = np.hstack(
            (node_at_link_tail.reshape((-1, 1)), node_at_link_head.reshape((-1, 1)))
        )

        # Sort them by midpoint coordinates
        self._sort_links_by_midpoint()

        # Flag indicating whether we have created patches
        self._patches_created = False

        # Flag indicating whether we have created link angles
        self._angle_of_link_created = False

        # set up the list of active links
        self._reset_link_status_list()

        # Create 2D array containing, for each node, direction of connected
        # link (1=incoming, -1=outgoing, 0=no link present at this position)
        # needs to come after BC setting
        self._create_link_dirs_at_node()

        #   set up link unit vectors and node unit-vector sums
        self._create_link_unit_vectors()

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
        self._face_at_link = sgrid.face_at_link(self.shape, actives=self.active_links)
        self._create_cell_areas_array()

        # List of neighbors for each cell: we will start off with no
        # list. If a caller requests it via active_adjacent_nodes_at_node or
        # _create_neighbor_list, we'll create it if necessary.
        self.neighbor_list_created = False
        self._neighbor_node_dict = {}

        # List of looped neighbor cells (all 8 neighbors) for
        # given *cell ids* can be created if requested by the user.
        self._looped_cell_neighbor_list = None

        # List of second ring looped neighbor cells (all 16 neighbors) for
        # given *cell ids* can be created if requested by the user.
        self._looped_second_ring_cell_neighbor_list_created = False

    def _setup_nodes(self):
        self._nodes = np.arange(self.number_of_nodes, dtype=int).reshape(self.shape)
        return self._nodes

    @property
    @make_return_array_immutable
    def nodes(self):
        """Get a shaped array of nodes.

        Returns
        -------
        ndarray
            Node IDs in an array shaped as *number_of_node_rows* by
            *number_of_node_columns*.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4))
        >>> grid.nodes
        array([[ 0,  1,  2,  3],
               [ 4,  5,  6,  7],
               [ 8,  9, 10, 11]])

        You can't change node ids.

        >>> grid.nodes[0] = 99 # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        ValueError: assignment destination is read-only

        LLCATS: NINF
        """
        return super(RasterModelGrid, self).nodes

    @property
    def nodes_at_right_edge(self):
        """Get nodes along the right edge of a grid.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4))
        >>> vals = np.array([ 0,  1,  2,  3,
        ...                   4,  5,  6,  7,
        ...                   8,  9, 10, 11])
        >>> vals[grid.nodes_at_right_edge]
        array([ 3,  7, 11])

        LLCATS: NINF BC SUBSET
        """
        return self.nodes[:, -1]

    @property
    def nodes_at_top_edge(self):
        """Get nodes along the top edge of a grid.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4))
        >>> vals = np.array([ 0,  1,  2,  3,
        ...                   4,  5,  6,  7,
        ...                   8,  9, 10, 11])
        >>> vals[grid.nodes_at_top_edge]
        array([ 8,  9, 10, 11])

        LLCATS: NINF BC SUBSET
        """
        return self.nodes[-1, :]

    @property
    def nodes_at_left_edge(self):
        """Get nodes along the left edge of a grid.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4))
        >>> vals = np.array([ 0,  1,  2,  3,
        ...                   4,  5,  6,  7,
        ...                   8,  9, 10, 11])
        >>> vals[grid.nodes_at_left_edge]
        array([0, 4, 8])

        LLCATS: NINF BC SUBSET
        """
        return self.nodes[:, 0]

    @property
    def nodes_at_bottom_edge(self):
        """Get nodes along the bottom edge of a grid.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4))
        >>> vals = np.array([ 0,  1,  2,  3,
        ...                   4,  5,  6,  7,
        ...                   8,  9, 10, 11])
        >>> vals[grid.nodes_at_bottom_edge]
        array([0, 1, 2, 3])

        LLCATS: NINF BC SUBSET
        """
        return self.nodes[0, :]

    def nodes_at_edge(self, edge):
        """Get edge nodes by edge name.

        Parameters
        ----------
        edge : {'right', 'top', 'left', 'bottom'}
            Edge location.

        Returns
        -------
        slice
            Slice of the nodes on an edge.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4))
        >>> vals = np.array([ 0,  1,  2,  3,
        ...                   4,  5,  6,  7,
        ...                   8,  9, 10, 11])
        >>> vals[grid.nodes_at_edge('left')]
        array([0, 4, 8])

        LLCATS: NINF BC SUBSET
        """
        if edge not in ("right", "top", "left", "bottom"):
            raise ValueError("value for edge not understood")
        return getattr(self, "nodes_at_{edge}_edge".format(edge=edge))

    def _create_cell_areas_array(self):
        """Set up array of cell areas.

        This method supports the creation of the array that stores cell areas.
        It is not meant to be called manually.
        """
        self._area_of_cell = np.full(
            self.number_of_cells, self.dx * self.dy, dtype=float
        )
        return self._area_of_cell

    def _create_cell_areas_array_force_inactive(self):
        """Set up array cell areas including extra cells for perimeter nodes.

        This method supports the creation of the array that stores cell areas.
        It differs from _create_cell_areas_array in that it forces ALL nodes to
        have a surrounding cell, which is not actually the case for the generic
        perimeter node (these are unbounded). This is only possible because the
        grid is a raster.
        It is not meant to be called manually.
        """
        self._forced_cell_areas = np.full(self.shape, self.dx * self.dy, dtype=float)
        self._forced_cell_areas[(0, -1), :] = 0.0
        self._forced_cell_areas[:, (0, -1)] = 0.0
        self._forced_cell_areas.shape = (-1,)
        return self._forced_cell_areas

    @property
    def shape(self):
        """Get the shape of the grid.

        Returns
        -------
        shape : tuple of ints
            The shape of the grid as number of node rows and node columns.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4))
        >>> grid.shape
        (3, 4)

        LLCATS: GINF NINF
        """
        return (self.number_of_node_rows, self.number_of_node_columns)

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

        LLCATS: GINF CINF
        """
        return (self.number_of_cell_rows, self.number_of_cell_columns)

    @property
    def dx(self):
        """Get node spacing in the column direction.

        Returns
        -------
        float
            Spacing of node columns.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.dx
        1.0
        >>> grid = RasterModelGrid((4, 5), xy_spacing=2.0)
        >>> grid.dx
        2.0

        LLCATS: GINF MEAS
        """
        return self._dx

    @property
    def dy(self):
        """Get node spacing in the row direction.

        Note in a RasterModelGrid, dy==dx.

        Returns
        -------
        float
            Spacing of node rows.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.dy
        1.0
        >>> grid = RasterModelGrid((4, 5), xy_spacing=(4, 2))
        >>> grid.dy
        2.0

        LLCATS: GINF MEAS
        """
        return self._dy

    @deprecated(use="vals[links_at_node]*active_link_dirs_at_node", version=1.0)
    def _active_links_at_node(self, *args):
        """_active_links_at_node([node_ids])
        Active links of a node.

        Parameters
        ----------
        node_ids : int or list of ints
                   ID(s) of node(s) for which to find connected active links

        Returns
        -------
        (4, N) ndarray
            The ids of active links attached to grid nodes with
            *node_ids*. If *node_ids* is not given, return links for all of the
            nodes in the grid. Link ids are listed in clockwise order starting
            with the south link. Diagonal links are never returned.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid((3, 4))
        >>> rmg.links_at_node[5]
        array([ 8, 11,  7,  4])
        >>> rmg._active_links_at_node((5, 6))
        array([[ 4,  5],
               [ 7,  8],
               [11, 12],
               [ 8,  9]])
        >>> rmg._active_links_at_node()
        array([[-1, -1, -1, -1, -1,  4,  5, -1, -1, 11, 12, -1],
               [-1, -1, -1, -1, -1,  7,  8,  9, -1, -1, -1, -1],
               [-1,  4,  5, -1, -1, 11, 12, -1, -1, -1, -1, -1],
               [-1, -1, -1, -1,  7,  8,  9, -1, -1, -1, -1, -1]])

        array([[-1, -1, -1, -1, -1,  0,  1, -1, -1,  2,  3, -1],
               [-1, -1, -1, -1, -1,  4,  5,  6, -1, -1, -1, -1],
               [-1,  0,  1, -1, -1,  2,  3, -1, -1, -1, -1, -1],
               [-1, -1, -1, -1,  4,  5,  6, -1, -1, -1, -1, -1]])

        LLCATS: DEPR LINF NINF
        """
        active_links_at_node = self.links_at_node.copy()
        active_links_at_node[self.active_link_dirs_at_node == 0] = -1
        active_links_at_node = active_links_at_node[:, (3, 2, 1, 0)]

        if len(args) == 0:
            return active_links_at_node.T
        elif len(args) == 1:
            node_ids = np.broadcast_arrays(args[0])[0]
            return active_links_at_node[node_ids, :].T
        else:
            raise ValueError("only zero or one arguments accepted")

    @property
    @make_return_array_immutable
    def horizontal_links(self):
        """
        LLCATS: LINF
        """
        try:
            return self._horizontal_links
        except AttributeError:
            self._horizontal_links = squad_links.horizontal_link_ids(self.shape)
            return self._horizontal_links

    @property
    @make_return_array_immutable
    def vertical_links(self):
        """
        LLCATS: LINF
        """
        try:
            return self._vertical_links
        except AttributeError:
            self._vertical_links = squad_links.vertical_link_ids(self.shape)
            return self._vertical_links

    @property
    @return_readonly_id_array
    def patches_at_node(self):
        """Get array of patches attached to nodes.

        Returns a (N, 4) array of the patches associated with each node in the
        grid.
        The four possible patches are returned in order CCW from east, i.e.,
        NE, NW, SW, SE.

        Missing patches are indexed -1.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((3, 3))
        >>> mg.patches_at_node
        array([[ 0, -1, -1, -1],
               [ 1,  0, -1, -1],
               [-1,  1, -1, -1],
               [ 2, -1, -1,  0],
               [ 3,  2,  0,  1],
               [-1,  3,  1, -1],
               [-1, -1, -1,  2],
               [-1, -1,  2,  3],
               [-1, -1,  3, -1]])

        LLCATS: PINF NINF CONN
        """
        try:
            return self.node_patch_matrix
        except AttributeError:
            self.node_patch_matrix = np.full((self.number_of_nodes, 4), -1, dtype=int)
            self.node_patch_matrix[:, 2][
                np.setdiff1d(
                    np.arange(self.number_of_nodes),
                    np.union1d(self.nodes_at_left_edge, self.nodes_at_bottom_edge),
                )
            ] = np.arange(self.number_of_patches)
            self.node_patch_matrix[:, 3][
                np.setdiff1d(
                    np.arange(self.number_of_nodes),
                    np.union1d(self.nodes_at_right_edge, self.nodes_at_bottom_edge),
                )
            ] = np.arange(self.number_of_patches)
            self.node_patch_matrix[:, 1][
                np.setdiff1d(
                    np.arange(self.number_of_nodes),
                    np.union1d(self.nodes_at_left_edge, self.nodes_at_top_edge),
                )
            ] = np.arange(self.number_of_patches)
            self.node_patch_matrix[:, 0][
                np.setdiff1d(
                    np.arange(self.number_of_nodes),
                    np.union1d(self.nodes_at_right_edge, self.nodes_at_top_edge),
                )
            ] = np.arange(self.number_of_patches)
            # we no longer blank out any patches that have a closed node as any
            # vertex, per modern LL style. Instead, we will make a closed/open
            # mask
            self._patches_created = True
            return self.node_patch_matrix

    @property
    @return_readonly_id_array
    def nodes_at_patch(self):
        """Get array of nodes of a patch.

        Returns the four nodes at the corners of each patch in a regular grid.
        Shape of the returned array is (nnodes, 4). Returns in order CCW from
        east, i.e., [NE, NW, SW, SE].

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((3, 3))
        >>> mg.nodes_at_patch
        array([[4, 3, 0, 1],
               [5, 4, 1, 2],
               [7, 6, 3, 4],
               [8, 7, 4, 5]])

        LLCATS: NINF PINF CONN
        """
        self._patches_created = True
        base = np.arange(self.number_of_patches)
        bottom_left_corner = base + base // (self._ncols - 1)
        return np.column_stack(
            (
                bottom_left_corner + self._ncols + 1,
                bottom_left_corner + self._ncols,
                bottom_left_corner,
                bottom_left_corner + 1,
            )
        )

    @property
    @return_readonly_id_array
    def links_at_patch(self):
        """Get array of links defining each patch.

        Examples
        --------
        >>> mg = RasterModelGrid((3, 4))
        >>> mg.links_at_patch
        array([[ 4,  7,  3,  0],
               [ 5,  8,  4,  1],
               [ 6,  9,  5,  2],
               [11, 14, 10,  7],
               [12, 15, 11,  8],
               [13, 16, 12,  9]])

        LLCATS: PINF LINF CONN
        """
        self._patches_created = True
        base = np.arange(self.number_of_patches)
        bottom_edge = base + (base // (self._ncols - 1)) * self._ncols
        return np.column_stack(
            (
                bottom_edge + self._ncols,
                bottom_edge + 2 * self._ncols - 1,
                bottom_edge + self._ncols - 1,
                bottom_edge,
            )
        )

    @property
    @return_readonly_id_array
    def patches_at_link(self):
        """Get array of patches adjoined to each link.

        Missing paches are indexed as -1.

        Examples
        --------
        >>> mg = RasterModelGrid((3, 4))
        >>> mg.patches_at_link
        array([[ 0, -1],
               [ 1, -1],
               [ 2, -1],
               [ 0, -1],
               [ 0,  1],
               [ 1,  2],
               [ 2, -1],
               [ 0,  3],
               [ 1,  4],
               [ 2,  5],
               [ 3, -1],
               [ 3,  4],
               [ 4,  5],
               [ 5, -1],
               [ 3, -1],
               [ 4, -1],
               [ 5, -1]])

        LLCATS: PINF LINF CONN
        """
        try:
            return self._patches_at_link
        except AttributeError:
            self._create_patches_at_link()
            return self._patches_at_link

    def _create_patches_at_link(self):
        from .cfuncs import create_patches_at_element

        self._patches_created = True
        self._patches_at_link = np.empty((self.number_of_links, 2), dtype=int)
        self._patches_at_link.fill(-1)
        create_patches_at_element(
            self.links_at_patch, self.number_of_links, self._patches_at_link
        )

    # a sort of the links will be performed here once we have corners

    def _create_link_dirs_at_node(self):
        """Make array with link directions at each node

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid((3, 4))
        >>> rmg._links_at_node
        array([[ 0,  3, -1, -1],
               [ 1,  4,  0, -1],
               [ 2,  5,  1, -1],
               [-1,  6,  2, -1],
               [ 7, 10, -1,  3],
               [ 8, 11,  7,  4],
               [ 9, 12,  8,  5],
               [-1, 13,  9,  6],
               [14, -1, -1, 10],
               [15, -1, 14, 11],
               [16, -1, 15, 12],
               [-1, -1, 16, 13]])
        >>> rmg._link_dirs_at_node
        array([[-1, -1,  0,  0],
               [-1, -1,  1,  0],
               [-1, -1,  1,  0],
               [ 0, -1,  1,  0],
               [-1, -1,  0,  1],
               [-1, -1,  1,  1],
               [-1, -1,  1,  1],
               [ 0, -1,  1,  1],
               [-1,  0,  0,  1],
               [-1,  0,  1,  1],
               [-1,  0,  1,  1],
               [ 0,  0,  1,  1]], dtype=int8)
        """
        # Create arrays for link-at-node information
        self._link_dirs_at_node = np.zeros((self.number_of_nodes, 4), dtype=np.int8)
        num_links_per_row = (self.number_of_node_columns * 2) - 1
        # Sweep over all links
        node_at_link_tail = self.node_at_link_tail
        node_at_link_head = self.node_at_link_head
        for lk in range(self.number_of_links):
            # Find the orientation
            is_horiz = (lk % num_links_per_row) < (self.number_of_node_columns - 1)
            # Find the IDs of the tail and head nodes
            t = node_at_link_tail[lk]
            h = node_at_link_head[lk]

            # If the link is horizontal, the index (row) in the links_at_node
            # array should be 0 (east) for the tail node, and 2 (west) for the
            # head node.
            # If vertical, the index should be 1 (north) for the tail node and
            # 3 (south) for the head node.
            if is_horiz:
                tail_index = 0
                head_index = 2
            else:
                tail_index = 1
                head_index = 3

            # Add this link to the list for this node, set the direction
            # (outgoing, indicated by -1), and increment the number found so
            # far
            self._link_dirs_at_node[t][tail_index] = -1
            self._link_dirs_at_node[h][head_index] = 1

        # setup the active link equivalent
        # self._active_link_dirs_at_node = self._link_dirs_at_node.copy()
        # inactive_links = (self.status_at_link[self.links_at_node] ==
        #                   INACTIVE_LINK)
        # inactive_links[self.link_dirs_at_node == 0] = False
        # self._active_link_dirs_at_node[inactive_links] = 0

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
        >>> mg = RasterModelGrid((3, 4), xy_spacing=(2., 2.))

        >>> mg.unit_vector_at_link[:, 0] # doctest: +NORMALIZE_WHITESPACE
        array([ 1.,  1.,  1.,  0.,  0.,  0.,  0.,
                1.,  1.,  1.,  0.,  0.,  0.,  0.,
                1.,  1.,  1.])
        >>> mg.unit_vector_at_link[:, 1] # doctest: +NORMALIZE_WHITESPACE
        array([ 0.,  0.,  0.,  1.,  1.,  1.,  1.,
                0.,  0.,  0.,  1.,  1.,  1.,  1.,
                0.,  0.,  0.])

        >>> mg.unit_vector_at_node[:, 0]
        array([ 1.,  2.,  2.,  1.,  1.,  2.,  2.,  1.,  1.,  2.,  2.,  1.])
        >>> mg.unit_vector_at_node[:, 1]
        array([ 1.,  1.,  1.,  1.,  2.,  2.,  2.,  2.,  1.,  1.,  1.,  1.])
        """
        unit_vec_at_link = np.zeros((self.number_of_links + 1, 2), dtype=float)
        unit_vec_at_link[self.horizontal_links, 0] = 1.0
        unit_vec_at_link[self.vertical_links, 1] = 1.0

        self._unit_vec_at_node = unit_vec_at_link[self.links_at_node].sum(axis=1)
        self._unit_vec_at_link = unit_vec_at_link[:-1, :]

    def _setup_link_at_face(self):
        """Set up links associated with faces.

        Returns an array of the link IDs for the links which intersect the
        faces specificed by *face_id*. *face_id* can be either a scalar or an
        array.

        Parameters
        ----------
        face_id : int
            Face of a cell.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((4, 5))
        >>> mg.link_at_face[0]
        5
        >>> mg.link_at_face[(0, 4, 13), ]
        array([ 5, 10, 21])
        """
        self._link_at_face = squad_faces.link_at_face(self.shape)
        return self._link_at_face

    def _create_face_at_link(self):
        """Set up array of faces associated with links.

        Return an array of the face IDs for the faces that intersect the links
        specified by *link_id*. *link_id* can be either a scalar or array. If
        *link_id* is not given, return the faces of all links.

        If a link does not have an associated face (e.g., some inactive links),
        that entry in the returned array is set to `BAD_INDEX_VALUE`.

        Parameters
        ----------
        link_id : array-like, optional
            Grid links.

        Examples
        --------
        >>> from landlab import RasterModelGrid, BAD_INDEX_VALUE
        >>> rmg = RasterModelGrid((4, 5))
        >>> rmg.face_at_link[5]
        0
        >>> faces = rmg.face_at_link[(0, 1, 15, 19, 12, 26), ]
        >>> faces[faces == BAD_INDEX_VALUE] = -1
        >>> faces
        array([-1, -1,  8, 11,  6, -1])
        """
        self._face_at_link = squad_faces.face_at_link(self.shape)
        return self._face_at_link

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

        >>> grid = RasterModelGrid((4, 5), xy_spacing=2.)
        >>> grid.extent
        (6.0, 8.0)

        >>> grid = RasterModelGrid((4, 5), xy_spacing=(3, 2))
        >>> grid.extent
        (6.0, 12.0)

        LLCATS: GINF MEAS
        """
        # Method added 5/1/13 by DEJH, modified DEJH 4/3/14 to reflect fact
        # boundary nodes don't have defined
        return (
            (self.number_of_node_rows - 1) * self._dy,
            (self.number_of_node_columns - 1) * self._dx,
        )

    @property
    def grid_xdimension(self):
        """Length of the grid in the x-dimension.

        Return the x-dimension of the grid. Because boundary nodes don't have
        cells, the dimension of the grid is num_cols-1, not num_cols.

        Returns
        -------
        float
            Length of the grid in the x-dimension.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.grid_xdimension
        4.0

        >>> grid = RasterModelGrid((4, 5), xy_spacing=0.5)
        >>> grid.grid_xdimension
        2.0

        >>> grid = RasterModelGrid((4, 5), xy_spacing=(3, 2))
        >>> grid.grid_xdimension
        12.0

        LLCATS: GINF MEAS
        """
        return (self.number_of_node_columns - 1) * self._dx

    @property
    def grid_ydimension(self):
        """Length of the grid in the y-dimension.

        Return the y-dimension of the grid. Because boundary nodes don't have
        cells, the dimension of the grid is num_rows-1, not num_rows.

        Returns
        -------
        float
            Length of the grid in the y-dimension.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.grid_ydimension
        3.0

        >>> grid = RasterModelGrid((4, 5), xy_spacing=0.5)
        >>> grid.grid_ydimension
        1.5

        >>> grid = RasterModelGrid((4, 5), xy_spacing=(3, 2))
        >>> grid.grid_ydimension
        6.0

        LLCATS: GINF MEAS
        """
        # Method added 5/1/13 by DEJH, modified DEJH 4/3/14, as above.
        return (self.number_of_node_rows - 1) * self._dy

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

        LLCATS: NINF
        """
        return sgrid.interior_node_count(self.shape)

    @property
    def number_of_node_columns(self):
        """Number of node columns.

        Returns the number of columns, including boundaries.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.number_of_node_columns
        5

        LLCATS: GINF NINF
        """
        return self._ncols

    @property
    def number_of_node_rows(self):
        """Number of node rows.

        Returns the number of rows, including boundaries.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.number_of_node_rows
        4

        LLCATS: GINF NINF
        """
        return self._nrows

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

        LLCATS: GINF NINF
        """
        return self._ncols - 2

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

        LLCATS: GINF CINF
        """
        return self._nrows - 2

    @property
    def number_of_patches(self):
        """Number of patches.

        Returns the number of patches over the grid.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.number_of_patches
        12

        LLCATS: PINF
        """
        return (self._nrows - 1) * (self._ncols - 1)

    @property
    @deprecated(use="nodes_at_corners_of_grid", version=1.0)
    def corner_nodes(self):
        """
        LLCATS: DEPR GINF NINF SUBSET
        """
        return self.nodes_at_corners_of_grid

    @property
    def nodes_at_corners_of_grid(self):
        """Get array of the nodes in grid corners.

        Return the IDs to the corner nodes of the grid, sorted by ID.

        Returns
        -------
        (4, ) ndarray
            Array of corner node IDs.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.nodes_at_corners_of_grid
        array([ 0,  4, 15, 19])

        LLCATS: GINF NINF SUBSET
        """
        return sgrid.corners((self._nrows, self._ncols))

    @property
    @deprecated(use="cells_at_corners_of_grid", version=1.0)
    def corner_cells(self):
        """
        LLCATS: DEPR GINF CINF SUBSET
        """
        return self.cells_at_corners_of_grid

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

        LLCATS: GINF CINF SUBSET
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
        >>> grid.is_point_on_grid((1, 1, 1,), (1, 3.1, 6.1))
        array([ True,  True, False], dtype=bool)
        >>> grid.is_point_on_grid((-.1, .1, 3.9, 4.1), (1, 1, 1, 1))
        array([False,  True,  True, False], dtype=bool)

        LLCATS: GINF MEAS SUBSET
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
        >>> grid.nodes_around_point(.4, 1.2)
        array([4, 8, 9, 5])

        >>> grid.nodes_around_point([.9, 1.1], 1.2)
        array([[ 4,  5],
               [ 8,  9],
               [ 9, 10],
               [ 5,  6]])

        >>> grid = RasterModelGrid((3, 4), xy_spacing=(1, 2))
        >>> grid.nodes_around_point(.5, 1.5)
        array([0, 4, 5, 1])
        >>> grid = RasterModelGrid((3, 4))
        >>> grid.nodes_around_point(.5, 1.5)
        array([4, 8, 9, 5])

        LLCATS: NINF SUBSET
        """
        xcoord, ycoord = np.broadcast_arrays(xcoord, ycoord)

        # Method added 4/29/13 by DEJH, modified 9/24/13.
        id_ = ycoord // self._dy * self.number_of_node_columns + xcoord // self._dx
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
        >>> rmg.find_nearest_node((np.array([1.6, 3.6]), np.array([2.3, .7])))
        array([12,  9])
        >>> rmg.find_nearest_node((-.4999, 1.))
        5

        LLCATS: NINF SUBSET
        """
        return rfuncs.find_nearest_node(self, coords, mode=mode)

    @property
    @cache_result_in_object()
    def length_of_link(self):
        """Get lengths of links.

        Return the link lengths in the grid, as a nlinks-long array.

        Returns
        -------
        (4, N) ndarray
            Link lengths.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 3), xy_spacing=(4, 3))
        >>> grid.length_of_link
        array([ 4.,  4.,  3.,  3.,  3.,  4.,  4.,  3.,  3.,  3.,  4.,  4.])

        LLCATS: LINF MEAS
        """
        length_of_link = np.empty(self.number_of_links, dtype=float)
        length_of_link[self.vertical_links] = self.dy
        length_of_link[self.horizontal_links] = self.dx

        return length_of_link

    def _create_length_of_link(self):
        """Calculate link lengths for a raster grid.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4), xy_spacing=(3, 2))
        >>> grid._create_length_of_link() # doctest: +NORMALIZE_WHITESPACE
        array([ 3., 3., 3.,
                2., 2., 2., 2.,
                3., 3., 3.,
                2., 2., 2., 2.,
                3., 3., 3.])

        >>> grid = RasterModelGrid((3, 3), xy_spacing=(2, 1))
        >>> grid._create_length_of_link() # doctest: +NORMALIZE_WHITESPACE
        array([ 2., 2.,
                1., 1., 1.,
                2., 2.,
                1., 1., 1.,
                2., 2.])

        Notes
        -----
        We initially set all lengths to dy. Then we loop over each row, setting
        the horizontal links in that row to dx.
        """
        if self._link_length is None:
            self._link_length = np.empty(self.number_of_links, dtype=float)
            self._link_length[self.vertical_links] = self.dy
            self._link_length[self.horizontal_links] = self.dx

        return self._link_length

    def set_closed_boundaries_at_grid_edges(
        self, right_is_closed, top_is_closed, left_is_closed, bottom_is_closed
    ):
        """Set boundary not to be closed.

        Sets the status of nodes along the specified side(s) of a raster
        grid (bottom, right, top, and/or left) to ``CLOSED_BOUNDARY``.

        Arguments are booleans indicating whether the bottom, left, top, and
        right are closed (``True``) or not (``False``).

        For a closed boundary:

        *  the nodes are flagged ``CLOSED_BOUNDARY`` (status type 4)
        *  all links that connect to a ``CLOSED_BOUNDARY`` node are
           flagged as inactive (so they appear on link-based lists, but
           not active_link-based lists)

        This means that if you call the calc_grad_at_active_link
        method, links connecting to closed boundaries will be ignored: there
        can be no gradients or fluxes calculated, because the links that
        connect to that edge of the grid are not included in the calculation.
        So, setting a grid edge to CLOSED_BOUNDARY is a convenient way to
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
        and all boundary nodes coded as FIXED_VALUE_BOUNDARY (=1):

        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid((4, 5)) # rows, columns, spacing
        >>> rmg.number_of_active_links
        17
        >>> rmg.status_at_node # doctest: +NORMALIZE_WHITESPACE
        array([1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1],
              dtype=uint8)
        >>> rmg.set_closed_boundaries_at_grid_edges(True, True, False, False)
        >>> rmg.number_of_active_links
        12
        >>> rmg.status_at_node # doctest: +NORMALIZE_WHITESPACE
        array([1, 1, 1, 1, 1, 1, 0, 0, 0, 4, 1, 0, 0, 0, 4, 4, 4, 4, 4, 4],
              dtype=uint8)

        LLCATS: BC SUBSET
        """
        if self._DEBUG_TRACK_METHODS:
            six.print_("ModelGrid.set_closed_boundaries_at_grid_edges")

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

        if bottom_is_closed:
            self._node_status[bottom_edge] = CLOSED_BOUNDARY

        if right_is_closed:
            self._node_status[right_edge] = CLOSED_BOUNDARY

        if top_is_closed:
            self._node_status[top_edge] = CLOSED_BOUNDARY

        if left_is_closed:
            self._node_status[left_edge] = CLOSED_BOUNDARY

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
        grid---bottom, right, top, and/or left---to FIXED_VALUE_BOUNDARY.

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
        conditions (and such methods, e.g.,
        ``RasterModelGrid.update_boundary_nodes()``, will raise exceptions
        if you try).

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
        >>> rmg.at_node['topographic__elevation'] = np.random.rand(20)
        >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> rmg.status_at_node # doctest: +NORMALIZE_WHITESPACE
        array([4, 4, 4, 4, 4,
               4, 0, 0, 0, 4,
               4, 0, 0, 0, 4,
               4, 4, 4, 4, 4], dtype=uint8)
        >>> rmg.set_fixed_value_boundaries_at_grid_edges(
        ...     True, True, False, False)
        >>> rmg.number_of_active_links
        12
        >>> rmg.status_at_node # doctest: +NORMALIZE_WHITESPACE
        array([4, 4, 4, 4, 4,
               4, 0, 0, 0, 1,
               4, 0, 0, 0, 1,
               1, 1, 1, 1, 1], dtype=uint8)

        Note that the four corners are treated as follows:

        *  bottom left = BOTTOM
        *  bottom right = BOTTOM
        *  top right = TOP
        *  top left = TOP

        This scheme is necessary for internal consistency with looped
        boundaries.

        LLCATS: BC SUBSET
        """
        if self._DEBUG_TRACK_METHODS:
            six.print_("ModelGrid.set_closed_boundaries_at_grid_edges")

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
            self._node_status[bottom_edge] = FIXED_VALUE_BOUNDARY

        if right_is_fixed_val:
            self._node_status[right_edge] = FIXED_VALUE_BOUNDARY

        if top_is_fixed_val:
            self._node_status[top_edge] = FIXED_VALUE_BOUNDARY

        if left_is_fixed_val:
            self._node_status[left_edge] = FIXED_VALUE_BOUNDARY

        self.reset_status_at_node()

        # save some internal data to speed updating:
        self.fixed_value_node_properties = {}
        self.fixed_value_node_properties["boundary_node_IDs"] = as_id_array(
            np.where(self._node_status == FIXED_VALUE_BOUNDARY)[0]
        )

        if value:
            if type(value) == float or type(value) == int:
                values_to_use = float(value)
            elif type(value) == np.ndarray:
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
            six.print_(
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
        try:
            self.fixed_value_node_properties["values"] = values_to_use
        except NameError:
            pass  # the flag will catch this case

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
        >>> rmg = RasterModelGrid((4, 5)) # rows, columns, spacing
        >>> rmg.number_of_active_links
        17
        >>> rmg.status_at_node # doctest: +NORMALIZE_WHITESPACE
        array([1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1],
              dtype=uint8)
        >>> rmg.add_zeros('topographic__elevation', at='node')
        array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
                0.,  0.,  0.,  0.,  0.,  0.,  0.])
        >>> rmg.set_looped_boundaries(True, True)
        >>> rmg.looped_node_properties['boundary_node_IDs']
        array([ 0,  1,  2,  3,  4,  5,  9, 10, 14, 15, 16, 17, 18, 19])
        >>> rmg.looped_node_properties['linked_node_IDs']
        array([10, 11, 12, 13, 14,  8,  6, 13, 11,  5,  6,  7,  8,  9])

        LLCATS: BC SUBSET
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
            self._node_status[bottom_edge] = LOOPED_BOUNDARY
            self._node_status[top_edge] = LOOPED_BOUNDARY
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
            self._node_status[right_edge] = LOOPED_BOUNDARY
            self._node_status[left_edge] = LOOPED_BOUNDARY
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

    @deprecated(use="_update_links_nodes_cells_to_new_BCs", version=1.0)
    def update_boundary_nodes(self):
        """Update the boundary nodes.

        This method updates all the boundary nodes in the grid field on which
        they are set (i.e., it updates the field
        rmg.at_node[rmg.fixed_gradient_node_properties['fixed_gradient_of']]).
        It currently works only with fixed value (type 1) and fixed gradient
        (type 2) conditions. Looping must be handled internally to a component,
        and is not dealt with here.

        LLCATS: DEPR NINF BC
        """
        try:
            self.fixed_value_node_properties["boundary_node_IDs"]
        except AttributeError:
            # no fixed value boundaries have been set
            pass
        else:
            assert self.fixed_value_node_properties["internal_flag"], (
                "Values were not supplied to the method that set the "
                "boundary conditions! You cant update automatically!"
            )
            values_val = self.at_node[
                self.fixed_value_node_properties["fixed_value_of"]
            ]
            values_val[
                self.fixed_value_node_properties["boundary_node_IDs"]
            ] = self.fixed_value_node_properties["values"]

        try:
            values_grad = self.at_node[
                self.fixed_gradient_node_properties["fixed_gradient_of"]
            ]
            values_grad[self.fixed_gradient_node_properties["boundary_node_IDs"]] = (
                values_grad[self.fixed_gradient_node_properties["anchor_node_IDs"]]
                + self.fixed_gradient_node_properties["values_to_add"]
            )
        except AttributeError:
            # no fixed grad boundaries have been set
            pass

    # DEJH believes this needs deprecating, but it's pretty hard wired into
    # the flow router. So I've restored it for now.
    def _calculate_gradients_at_d8_active_links(self, node_values):
        """Calculate gradients over D8 active links.

        MAY 16: Landlab's handling of diagonal links may soon be enhanced;
        methods like this may be soon superceded.

        Parameters
        ----------
        node_values : ndarray
            Values at nodes.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> import numpy as np
        >>> grid = RasterModelGrid((3, 4), xy_spacing=(4, 3))
        >>> z = np.array([3., 3., 3., 3.,
        ...               3., 3., 0., 0.,
        ...               3., 0., 0., 0.])
        >>> grid._calculate_gradients_at_d8_active_links(z)
        ...     # doctest: +NORMALIZE_WHITESPACE
        array([ 0.  , -1.  ,  0.  , -0.75,  0.  , -1.  ,  0.  ,  0.  , -0.6 ,
                0.  , -0.6 ,  0.  , -0.6 ,  0.  , 0. ])

        LLCATS: LINF GRAD
        """
        active_links = self.active_d8
        diagonal_links = squad_links.is_diagonal_link(self.shape, active_links)
        active_links = active_links[~diagonal_links]

        vertical_links = squad_links.is_vertical_link(self.shape, active_links)
        horizontal_links = squad_links.is_horizontal_link(self.shape, active_links)

        diffs = np.diff(node_values[self.nodes_at_link[self.active_links]], axis=1)

        diffs[vertical_links] /= self.dy
        diffs[horizontal_links] /= self.dx

        diag_dist = np.sqrt(self.dy ** 2.0 + self.dx ** 2.0)
        diagonal_link_slopes = (
            np.diff(node_values[self.nodes_at_diagonal[self.active_diagonals]], axis=1)
            / diag_dist
        )

        return np.concatenate((diffs.flatten(), diagonal_link_slopes.flatten()))

    def _calculate_gradients_at_d8_links(self, node_values):
        """Calculate gradients over all D8 links.

        MAY 16: Landlab's handling of diagonal links may soon be enhanced;
        methods like this may be soon superceded.

        Parameters
        ----------
        node_values : ndarray
            Values at nodes.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> import numpy as np
        >>> grid = RasterModelGrid((3, 4), xy_spacing=(4, 3))
        >>> z = np.array([3., 3., 3., 3.,
        ...               3., 3., 0., 0.,
        ...               3., 0., 0., 0.])
        >>> grid._calculate_gradients_at_d8_links(z)
        ...     # doctest: +NORMALIZE_WHITESPACE
        array([ 0. ,  0. , -0.6,  0. , -0.6, -0.6, -0.6,  0. , -0.6,  0. ,  0. ,
                0. ])

        LLCATS: LINF GRAD
        """
        diag_dist = np.sqrt(self.dy ** 2.0 + self.dx ** 2.0)
        diagonal_link_slopes = (
            node_values[self.nodes_at_diagonal[:, 1]]
            - node_values[self.nodes_at_diagonal[:, 0]]
        ) / diag_dist

        return diagonal_link_slopes

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
        >>> u = rmg.zeros(centering='node')
        >>> u = u + range(0, len(u))
        >>> u # doctest: +NORMALIZE_WHITESPACE
        array([ 0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.,
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

        LLCATS: GINF NINF
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

        LLCATS: GINF CINF
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
        >>> data = rmg.add_zeros('test_data', at='node')
        >>> data[:] = np.arange(12)
        >>> rmg.roll_nodes_ud('test_data', 1)
        >>> data # doctest: +NORMALIZE_WHITESPACE
        array([ 9.,  10.,  11.,   0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,
                8.])
        >>> rmg.roll_nodes_ud('test_data', 2)
        >>> data # doctest: +NORMALIZE_WHITESPACE
        array([ 3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.,  11.,   0.,   1.,
                2.])
        >>> rmg.roll_nodes_ud('test_data', 1, interior_only=True)
        >>> data # doctest: +NORMALIZE_WHITESPACE
        array([ 3.,   1.,   5.,   6.,   4.,   8.,   9.,   7.,  11.,   0.,  10.,
                2.])

        LLCATS: NINF
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

    def REMOVE_create_neighbor_list(self, bad_index=BAD_INDEX_VALUE):
        """Create list of neighbor node IDs.

        Creates a list of IDs of neighbor nodes for each node, as a
        2D array. Only record neighbor nodes that are on the other end of an
        *active* link. Nodes attached to *inactive* links or neighbor nodes
        that would be outside of the grid are given an ID of
        :const:`~landlab.grid.base.BAD_INDEX_VALUE`.

        Neighbors are ordered as [*right*, *top*, *left*, *bottom*].
        """
        # assert(self.neighbor_list_created == False)
        # this method can now be called to create multiple neighbor lists with
        # different BAD_INDEX_VALUES
        # note self.nieghbor_nodes is no longer created... but nobody should be
        # calling it direct anyway.

        neighbor_nodes = sgrid.neighbor_node_array(
            self.shape,
            closed_boundary_nodes=self.closed_boundary_nodes,
            open_boundary_nodes=self.open_boundary_nodes,
            inactive=bad_index,
        ).T

        self.neighbor_list_created = True
        return neighbor_nodes

    @deprecated(use="node_has_boundary_neighbor", version=1.0)
    def has_boundary_neighbor(self, ids, method="d8"):
        """
        LLCATS: DEPR NINF CONN BC
        """
        return self.node_has_boundary_neighbor(ids, method=method)

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
        array([False,  True], dtype=bool)

        >>> mg.node_has_boundary_neighbor(25)
        Traceback (most recent call last):
            ...
        IndexError: index 25 is out of bounds for axis 0 with size 25

        LLCATS: NINF CONN BC
        """
        ans = _node_has_boundary_neighbor(self, ids, method=method)

        if ans.ndim == 0:
            return bool(ans)
        else:
            return ans

    @deprecated(use="node_is_core", version=1.0)
    def is_core(self, *args):
        """
        LLCATS: DEPR NINF BC
        """
        return self.node_is_core(*args)

    def node_is_core(self, *args):
        """node_is_core([ids])
        Check if a node is a core node.

        Returns an boolean array of truth values for each node ID provided;
        True if the node is a core node, False otherwise.
        If no IDs are provided, method returns a boolean array for every node.

        (Core status is typically indicated by a value of 0 in node_status.)

        LLCATS: NINF BC
        """
        # NG changed this.
        # Modified DEJH May 2014 to accept simulaneous tests of multiple nodes;
        # should still be back-conmpatible.
        try:
            node_ids = args[0]
        except IndexError:  # return all nodes
            return np.equal(self._node_status, CORE_NODE)
        else:
            return np.equal(self._node_status[node_ids], CORE_NODE)

    @deprecated(use="nodes_are_all_core", version=1.0)
    def are_all_interior(self, IDs):
        """Check if nodes are interior.

        Returns a single boolean truth value, True if all nodes with *IDs* are
        interior nodes, False if not.

        LLCATS: DEPR NINF BC
        """
        return np.all(np.equal(self._node_status[IDs], CORE_NODE))

    @deprecated(use="nodes_are_all_core", version=1.0)
    def are_all_core(self, ids):
        """
        LLCATS: DEPR NINF BC
        """
        return self.nodes_are_all_core(ids)

    def nodes_are_all_core(self, ids):
        """Check if nodes are all core.

        Returns a single boolean truth value, True if all nodes with *IDs* are
        core nodes, False if not.

        Parameters
        ----------
        ids : array-like
            Grid nodes.

        Returns
        -------
        boolean
            ``True`` if all the given nodes are *core* nodes.

        LLCATS: NINF BC
        """
        return np.all(np.equal(self._node_status[ids], CORE_NODE))

    @deprecated(use="no replacement", version=1.0)
    def face_connecting_cell_pair(self, cell_a, cell_b):
        """Get the face that connects two cells.

        Returns an array of face indices that *cell_a* and *cell_b* share.
        If the cells do not share any faces, returns an empty array.

        Examples
        --------
        >>> import pytest
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((4, 5))
        >>> with pytest.deprecated_call():
        ...     mg.face_connecting_cell_pair(0, 1)
        array([4])
        >>> with pytest.deprecated_call():
        ...     mg.face_connecting_cell_pair(0, 2).size  # empty array returned
        0

        LLCATS: DEPR FINF CINF CONN
        """
        cell_faces = self.faces_at_cell[[cell_a, cell_b]]
        return as_id_array(
            np.intersect1d(cell_faces[0], cell_faces[1], assume_unique=True)
        )

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

        LLCATS: NINF SUBSET MEAS
        """
        return np.ravel_multi_index((row, col), self.shape, **kwds)

    def _create_face_width(self):
        """Set up array of face widths.

        Produces an array of length nfaces containing the face width.

        Returns
        -------
        ndarray of float
            Width of faces (listed as horizontal, then vertical).

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 3))
        >>> grid.width_of_face
        array([ 1.,  1.,  1.,  1.])
        """
        n_horizontal_faces = (self.shape[0] - 2) * (self.shape[1] - 1)

        self._face_width = np.empty(squad_faces.number_of_faces(self.shape))
        self._face_width[:n_horizontal_faces] = self.dx
        self._face_width[n_horizontal_faces:] = self.dy
        return self._face_width

    def _unit_test(self):
        """Stub for adding unit tests to RasterModelGrid."""
        pass

    def calc_unit_normal_at_patch(self, elevs="topographic__elevation"):
        """Calculate and return the unit normal vector <a, b, c> to a patch.

        This method is not defined on a raster, as there is no unique unit
        normal for a square patch. Use
        `_calc_unit_normals_to_patch_subtriangles` instead.

        LLCATS: PINF GRAD
        """
        raise NotImplementedError(
            "This method is not defined on a raster, as there is no unique "
            "unit normal for a square patch. Use "
            "`_calc_unit_normals_to_patch_subtriangles` instead."
        )

    @deprecated(use="calc_aspect_at_node", version=1.0)
    def calculate_aspect_at_nodes_bestFitPlane(self, id, val):
        """Aspect at nodes.

        .. codeauthor:: Katy Barnhart <katherine.barnhart@colorado.edu>

        Calculates the aspect at each node based on the elevation of
        the node and its neighbors using a best fit plane calculated
        using single value decomposition.

        Parameters
        ----------
        id : array-like
            ID of nodes at which to calculate the aspect.
        val : ndarray
            Elevation at all nodes

        Returns
        -------
        ndarray
            Aspect at the nodes given by id

        LLCATS: DEPR NINF SURF
        """
        # additional note, KRB has written three codes in raster.py
        # one to calculate slope, one to calculate aspect, and one
        # to calculate both

        # get the list of neighboring nodes for the nodes given by id
        n = self.active_adjacent_nodes_at_node[id]
        a = []

        # for each node in id make a list with the node id and the ids of
        # its neighbors.

        # determine the values for the x, y, and z coordinates of each node,
        # pass these to rfuncs.calculate_slope_aspect_bfp to calculate the
        # slope and aspect.

        indBool = n != BAD_INDEX_VALUE

        for i in range(len(id)):
            # make a list of the neighbor nodes and
            # check that none of the nodes are bad

            ns = list(n[0][indBool[0]])
            ns.append(id[i])

            x = self.node_x[ns]
            y = self.node_y[ns]
            z = val[ns]
            slope, aspect = rfuncs.calculate_slope_aspect_bfp(x, y, z)
            a.append(aspect)
            del ns
        # return aspect alone
        return a

    @deprecated(use="calc_slope_at_node", version=1.0)
    def calculate_slope_at_nodes_bestFitPlane(self, id, val):
        """Slope of best-fit plane at nodes.

        .. codeauthor:: Katy Barnhart <katherine.barnhart@colorado.edu>

        Calculates the slope at each node based on the elevation of
        the node and its neighbors using a best fit plane calculated
        using single value decomposition.

        Parameters
        ----------
        id : array-like
            ID of nodes at which to calculate the aspect
        val : ndarray
            Elevation at all nodes

        Returns
        -------
        ndarray
            Slope at the nodes given by id

        LLCATS: DEPR NINF GRAD SURF
        """
        #
        # additional note, KRB has written three codes in raster.py
        # one to calculate slope, one to calculate aspect, and one
        # to calculate both

        # get the list of neighboring nodes for the nodes given by id
        n = self.active_adjacent_nodes_at_node[id]
        s = []

        # for each node in id make a list with the node id and the ids of
        # its neighbors.

        # determine the values for the x, y, and z coordinates of each node,
        # pass these to rfuncs.calculate_slope_aspect_bfp to calculate the
        # slope and aspect.

        indBool = n != BAD_INDEX_VALUE

        for i in range(len(id)):
            # make a list of the neighbor nodes and
            # check that none of the nodes are bad

            ns = list(n[0][indBool[0]])
            ns.append(id[i])

            x = self.node_x[ns]
            y = self.node_y[ns]
            z = val[ns]

            slope, _ = rfuncs.calculate_slope_aspect_bfp(x, y, z)
            s.append(slope)
            del ns
        # return slope alone
        return s

    @deprecated(use="calc_slope_at_node, calc_aspect_at_node", version=1.0)
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
        >>> z = np.array([0., 0., 0., 0.,
        ...               3., 3., 3., 3,
        ...               6., 6., 6., 6.])
        >>> with pytest.deprecated_call():
        ...     (slope, aspect) = (
        ...              grid.calculate_slope_aspect_at_nodes_burrough(vals=z))
        >>> np.tan(slope)
        array([ 0.75,  0.75])
        >>> np.degrees(aspect)
        array([ 180.,  180.])

        This method is *deprecated*. Use ``calc_slope_at_node`` and
        ``calc_aspect_at_node`` instead. Notice that ``calc_slope_at_node``
        and ``calc_aspect_at_node`` return values for all nodes, not just
        core nodes. In addition, ``calc_aspect_at_node`` returns compass-style
        angles in degrees.

        >>> np.tan(grid.calc_slope_at_node(elevs=z)[grid.core_nodes])
        array([ 0.75,  0.75])
        >>> grid.calc_aspect_at_node(elevs=z)[grid.core_nodes]
        array([ 180.,  180.])

        LLCATS: DEPR NINF SURF GRAD
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
        ) / (8.0 * self._dx)
        dz_dy = (
            (bottom_left + 2 * bottom + bottom_right) - (top_left + 2 * top + top_right)
        ) / (8.0 * self._dy)

        slope = np.zeros([ids.shape[0]], dtype=float)
        aspect = np.zeros([ids.shape[0]], dtype=float)
        slope = np.arctan(np.sqrt(dz_dx ** 2 + dz_dy ** 2))
        aspect = np.arctan2(dz_dy, -dz_dx)
        aspect = np.pi * 0.5 - aspect
        aspect[aspect < 0.0] = aspect[aspect < 0.0] + 2.0 * np.pi
        aspect[slope == 0.0] = -1.0

        return slope, aspect

    @deprecated(use="calc_slope_at_node, calc_aspect_at_node", version=1.0)
    def calculate_slope_aspect_at_nodes_best_fit_plane(self, nodes, val):
        r"""Calculate slope aspect.

        Slope aspect of best-fit plane at nodes.

        .. codeauthor:: Katy Barnhart <katherine.barnhart@colorado.edu>

        .. note::

            THIS CODE HAS ISSUES (SN 25-Sept-14): This code didn't perform
            well on a NS facing elevation profile. Please check
            slope_aspect_routines_comparison.py under landlab\examples before
            using this.  Suggested alternative:
            calculate_slope_aspect_at_nodes_burrough

        Calculates both the slope and aspect at each node based on the
        elevation of the node and its neighbors using a best fit plane
        calculated using single value decomposition.

        Parameters
        ----------
        nodes : array-like
            ID of nodes at which to calculate the aspect
        val : ndarray
            Elevation at all nodes

        Returns
        -------
        (slope, aspect) : tuple of floats
            Tuple containing (*slope*, *aspect*)

        LLCATS: DEPR NINF GRAD SURF
        """
        # additional note, KRB has written three codes in raster.py
        # one to calculate slope, one to calculate aspect, and one
        # to calculate both

        # get the list of neighboring nodes for the nodes given by id
        node_neighbors = self.active_adjacent_nodes_at_node[nodes]
        aspects = []
        slopes = []

        # for each node in id make a list with the node id and the ids of
        # its neighbors.

        # determine the values for the x, y, and z coordinates of each node,
        # pass these to rfuncs.calculate_slope_aspect_bfp to calculate the
        # slope and aspect.

        indBool = node_neighbors != BAD_INDEX_VALUE

        for id_ in range(len(nodes)):
            # make a list of the neighbor nodes and
            # check that none of the nodes are bad
            neighbors = list(node_neighbors[0, indBool[0]])
            neighbors.append(nodes[id_])

            node_x = self.node_x[neighbors]
            node_y = self.node_y[neighbors]
            node_z = val[neighbors]
            slope, aspect = rfuncs.calculate_slope_aspect_bfp(node_x, node_y, node_z)
            aspects.append(aspect)
            slopes.append(slope)

            del neighbors
        return slopes, aspects

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

        LLCATS: GINF
        """
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
        """
        For each cell in a raster, return the D8 neighboring cells, looping
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

        LLCATS: DEPR CINF CONN BC
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

        ...take a look at the cell grid to understand why:
        [56, 57, 58, 59, 60, 61, 62, 63]
        [48, 49, 50, 51, 52, 53, 54, 55]
        [40, 41, 42, 43, 44, 45, 46, 47]
        [32, 33, 34, 35, 36, 37, 38, 39]
        [24, 25, 26, 27, 28, 29, 30, 31]
        [16, 17, 18, 19, 20, 21, 22, 23]
        [ 8,  9, 10, 11, 12, 13, 14, 15]
        [ 0,  1,  2,  3,  4,  5,  6,  7]

        LLCATS: CINF CONN BC
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

        Creates a list of looped immediate cell neighbors for each cell as a
        2D array of size ( self.number_of_cells, 16 ).
        Order or neighbors: Starts with E and goes counter clockwise
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

    # This is deprecated as link status should only be changed by
    # changing node status.
    def set_fixed_link_boundaries_at_grid_edges(
        self,
        right_is_fixed,
        top_is_fixed,
        left_is_fixed,
        bottom_is_fixed,
        link_value=None,
        node_value=None,
        fixed_node_value_of="topographic__elevation",
        fixed_link_value_of="topographic__slope",
    ):
        """Create fixed link boundaries at the grid edges.

        Sets the status of links along the specified side(s) of a raster
        grid--- bottom vertical links, right horizontal, top vertical links,
        and/or left horizontal links ---to FIXED_LINK.

        By definition, fixed links exist between fixed gradient nodes
        (status_at_node == 2) and core nodes (status_at_node == 0). Because the
        outer ring of nodes are fixed gradient (status_at_node == 2), the links
        between them are inactive (status_at_link == 4) and are not set using
        this function (the inactive links are the top and bottom horizontal
        edge links, and left and right edge vertical edge links.)

        Arguments are booleans indicating whether the bottom, right, top, and
        left sides are to be set (True) or not (False).

        *node_value* controls what values are held constant at the fixed
        gradient nodes (status_at_node == 2). It can be either a float, an
        array of length number_of_fixed_nodes or number_of_nodes (total), or
        left blank. If left blank, the values will be set from the those
        already in the grid fields, according to 'fixed_node_value_of'.

        *link_value* controls what values are held constant at the fixed
        links (status_at_link == 2). It can be either a float, an array of
        length number_of_fixed_links or number_of_links (total), or
        left blank. If left blank, the values will be set from the those
        already in the grid fields, according to 'fixed_link_value_of'.

        *fixed_node_value_of* controls the name of the model field that
        contains the node values. Remember, if you don't set value, the fixed
        gradient node values will be set from the field values ***at the time
        you call this method***. If no values are present in the field, the
        module will complain but accept this, warning that it will be unable to
        automatically update boundary conditions (and such methods, e.g.,
        ``RasterModelGrid.update_boundary_nodes()``, will raise exceptions
        if you try).

        *fixed_link_value_of* controls the name of the model field that
        contains the fixed link values. Remember, if you don't set value, the
        fixed link values will be set from the field values ***at the time you
        call this method***. If no values are present in the field, the module
        will complain but accept this, warning that it will be unable to
        automatically update boundary conditions (and such methods, e.g.,
        ``RasterModelGrid.update_boundary_nodes()``, will raise exceptions
        if you try).

        The following example sets the bottom and right link boundaries as
        fixed-value in a four-row by nine-column grid that initially has all
        boundaries set to fixed_gradient (nodes, i.e. flagged at
        (status_at_node == 2) and fixed_link (links, i.e., flagged as
        (status_at_link == 2).

        Parameters
        ----------
        right_is_fixed : boolean
            Set right edge  horizontal links as fixed boundary.
        top_is_fixed : boolean
            Set top edge vertical links as fixed boundary.
        left_is_fixed : boolean
            Set left edge horizontal links as fixed boundary.
        bottom_is_fixed : boolean
            Set bottom edge vertical links as fixed boundary.
        link_value : float, array or None (default).
            Override value to be kept constant at links.
        node_value : float, array or None (default).
            Override value to be kept constant at nodes.
        fixed_node_value_of : string.
            The name of the grid field containing the values of interest at
            nodes.
        fixed_link_value_of : string.
            The name of the grid field containing the values of interest at
            links.

        Examples
        --------

        The following grid is used in the example::

            *--I--->*--I--->*--I--->*--I--->*--I--->*--I--->*--I--->*--I--->*
            ^       ^       ^       ^       ^       ^       ^       ^       ^
            I       X       X       X       X       X       X       X       I
            |       |       |       |       |       |       |       |       |
            *--X--->o       o       o       o       o       o       o--X--->*
            ^       ^       ^       ^       ^       ^       ^       ^       ^
            I       |       |       |       |       |       |       |       I
            |       |       |       |       |       |       |       |       |
            *--X--->o       o       o       o       o       o       o--X--->*
            ^       ^       ^       ^       ^       ^       ^       ^       ^
            I       X       X       X       X       X       X       X       I
            |       |       |       |       |       |       |       |       |
            *--I--->*--I--->*--I--->*--I--->*--I--->*--I--->*--I--->*--I--->*

        .. note::

          Links set to :any:`ACTIVE_LINK` are not indicated in this diagram.

        ``*`` indicates the nodes that are set to
        :any:`FIXED_GRADIENT BOUNDARY`

        ``o`` indicates the nodes that are set to :any:`CORE_NODE`

        ``I`` indicates the links that are set to :any:`INACTIVE_LINK`

        ``X`` indicates the links that are set to :any:`FIXED_LINK`

        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid((4, 9), xy_spacing=1.0)
        >>> import numpy as np
        >>> z = np.arange(0, rmg.number_of_nodes)
        >>> s = np.arange(0, rmg.number_of_links)
        >>> rmg['node']['topographic__elevation'] = z
        >>> rmg['link']['topographic__slope'] = s
        >>> rmg.set_fixed_link_boundaries_at_grid_edges(True, True, True, True)
        >>> rmg.status_at_node # doctest: +NORMALIZE_WHITESPACE
        array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0,
               0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2], dtype=uint8)
        >>> rmg.status_at_link # doctest: +NORMALIZE_WHITESPACE
        array([4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 4, 2, 0, 0, 0,
               0, 0, 0, 2, 4, 0, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 2,
               4, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4], dtype=uint8)
        >>> rmg.fixed_link_properties['fixed_gradient_of']
        'topographic__slope'
        >>> rmg.fixed_gradient_node_properties['fixed_gradient_of']
        'topographic__elevation'

        LLCATS: BC SUBSET
        """
        # Fixed link boundaries are found between core nodes (node_status==0)
        # and fixed gradient nodes (node_status==2). To assure these conditions
        # are met, we store link and node boundary IDs in arrays...
        fixed_nodes = np.array([], dtype=int)

        # Based on the inputs, we then assign boundary status. Starting
        # from the right edge (east edge) we look to see if the boolean input
        # is True or False. If true, we find the appropriate links and nodes
        # and set them to the boundary condition of FIXED_GRADIENT_BOUNDARY
        # for nodes and FIXED_LINK for links.
        if right_is_fixed:
            fixed_nodes = np.append(fixed_nodes, self.nodes_at_right_edge)

        if top_is_fixed:
            fixed_nodes = np.append(fixed_nodes, self.nodes_at_top_edge)

        if left_is_fixed:
            fixed_nodes = np.append(fixed_nodes, self.nodes_at_left_edge)

        if bottom_is_fixed:
            fixed_nodes = np.append(fixed_nodes, self.nodes_at_bottom_edge)

        fixed_nodes = np.unique(fixed_nodes)
        fixed_nodes = np.sort(fixed_nodes)

        self.status_at_node[fixed_nodes] = self.BC_NODE_IS_FIXED_GRADIENT

        # Now we are testing to see what values will be assigned to these
        # boundaries

        fixed_links = self.fixed_links

        # For links, the default is topographic slope ('topographic__slope')
        # First, see if there is a scalar value...
        if link_value is None:

            # if not, we assign the link values from the field of
            # 'topographic_slope'. If it does not exists, an error will be
            # kicked out.
            assigned_link_values = self.at_link[fixed_link_value_of][fixed_links]

        else:

            # If there IS a scalar link value, it is instead set here.
            assigned_link_values = np.ones(fixed_links.size) * link_value

        # For nodes, the default value is 'topographic__elevation'.
        if node_value is None:

            # If no scalar is found, the field values for
            # 'topographic__elevation' are used. If this field does not
            # exist, an error will be returned.
            assigned_node_values = self.at_node[fixed_node_value_of][fixed_nodes]
        else:

            # If there is a scalar, it is instead set here.
            assigned_node_values = np.ones(fixed_nodes.size) * node_value

        # Now we will set the attributes using a Python dictionary.
        # First, nodes.
        try:
            # Simply testing to make sure no boundary conditions exist...
            self.fixed_gradient_node_properties["boundary_node_IDs"]
        except AttributeError:

            # If they don't exist, we set them here.
            self.fixed_gradient_node_properties = {}

            # Setting the node ids in the dictionary.
            self.fixed_gradient_node_properties["boundary_node_IDs"] = fixed_nodes

            # What gradient was assigned to the nodes? That is set here.
            self.fixed_gradient_node_properties[
                "fixed_gradient_of"
            ] = fixed_node_value_of

            # Assigned gradient values at the nodes set in the dictionary.
            self.fixed_gradient_node_properties[
                "boundary_node_gradients"
            ] = assigned_node_values

        # Then, links
        try:
            # First, test to make sure no boundary conditions exist.
            self.fixed_link_properties["boundary_link_IDs"]

        except AttributeError:

            # If they don't exist, we set them here
            self.fixed_link_properties = {}

            # Setting the link IDs in the dictionary
            self.fixed_link_properties["boundary_link_IDs"] = fixed_links

            # What gradient are we assigning to the links? That is set here.
            self.fixed_link_properties["fixed_gradient_of"] = fixed_link_value_of

            # Assigned gradient values at the links set in the dictionary
            self.fixed_link_properties["boundary_link_gradients"] = assigned_link_values

    def set_watershed_boundary_condition(
        self,
        node_data,
        nodata_value=-9999.0,
        return_outlet_id=False,
        remove_disconnected=False,
        adjacency_method="D8",
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

        Note that the outer ring (perimeter) of the raster is set to
        CLOSED_BOUNDARY, even if there are nodes that have values. The only
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
        as illustrated:

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
        >>> rmg = RasterModelGrid((4,4))
        >>> node_data = np.array([-9999., -9999., -9999., -9999.,
        ...                      -9999.,    67.,    67., -9999.,
        ...                      -9999.,    67.,     0., -9999.,
        ...                      -9999., -9999., -9999., -9999.])
        >>> out_id = rmg.set_watershed_boundary_condition(node_data, -9999.,
        ...                                              True)
        >>> out_id
        array([10])
        >>> rmg.status_at_node
        array([4, 4, 4, 4, 4, 0, 0, 4, 4, 0, 1, 4, 4, 4, 4, 4], dtype=uint8)
        >>> rmg2 = RasterModelGrid((4,4))
        >>> node_data2 = np.array([-9999., -9999., -9999., -9999.,
        ...                      -9999.,    67.,    67.,    -2.,
        ...                      -9999.,    67.,     0., -9999.,
        ...                      -9999., -9999., -9999., -9999.])
        >>> rmg2.set_watershed_boundary_condition(node_data2, -9999.)
        >>> rmg2.status_at_node
        array([4, 4, 4, 4, 4, 0, 0, 1, 4, 0, 0, 4, 4, 4, 4, 4], dtype=uint8)

        The node data can also be provided as a model grid field.

        >>> rmg = RasterModelGrid((4,4))
        >>> node_data = np.array([-9999., -9999., -9999., -9999.,
        ...                      -9999.,    67.,    67., -9999.,
        ...                      -9999.,    67.,     0., -9999.,
        ...                      -9999., -9999., -9999., -9999.])
        >>> _ = rmg.add_field('node', 'topographic__elevation', node_data)
        >>> out_id = rmg.set_watershed_boundary_condition('topographic__elevation',
        ...                                               -9999.,
        ...                                               True)
        >>> out_id
        array([10])
        >>> rmg.status_at_node
        array([4, 4, 4, 4, 4, 0, 0, 4, 4, 0, 1, 4, 4, 4, 4, 4], dtype=uint8)

        LLCATS: BC
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
        self.status_at_node[outlet_loc] = FIXED_VALUE_BOUNDARY

        if remove_disconnected:
            self.set_open_nodes_disconnected_from_watershed_to_closed(
                node_data=node_data,
                outlet_id=as_id_array(np.array([outlet_loc])),
                nodata_value=nodata_value,
                adjacency_method=adjacency_method,
            )
        if return_outlet_id:
            return as_id_array(np.array([outlet_loc]))

    def set_open_nodes_disconnected_from_watershed_to_closed(
        self, node_data, outlet_id=None, nodata_value=-9999.0, adjacency_method="D8"
    ):
        """
        Identifys all non-closed nodes that are disconnected from the node given in
        *outlet_id* and sets them as closed.

        If *outlet_id* is not given, the outlet will be identified as the node
        for which the status at the node is FIXED_VALUE_BOUNDARY. If more than
        one node has this value, the algorithm will fail.

        If *outlet_id* is given, the algorithm will check that it is not a node
        with status of CLOSED_BOUNDARY.

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
            FIXED_VALUE_BOUNDARY. However, it must not have the status of
            CLOSED_BOUNDARY.
        nodata_value : float, optional, default is -9999.
            Value that indicates an invalid value.
        adjacency_method : string, optional. Default is 'D8'.
            Sets the connection method.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> mg1 = RasterModelGrid((4,6))
        >>> z1 = np.array([-9999., -9999., -9999.,  -9999., -9999., -9999.,
        ...                -9999.,    67.,    67.,  -9999.,    50., -9999.,
        ...                -9999.,    67.,     0.,  -9999., -9999., -9999.,
        ...                -9999., -9999., -9999.,  -9999., -9999., -9999.])
        >>> mg2 = RasterModelGrid((4,6))
        >>> z2 = np.array([-9999., -9999., -9999.,  -9999., -9999., -9999.,
        ...                -9999.,    67.,    67.,  -9999.,    50., -9999.,
        ...                -9999.,    67.,     0.,  -9999., -9999., -9999.,
        ...                -9999., -9999., -9999.,  -9999., -9999., -9999.])
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

        LLCATS: BC
        """
        # get node_data if a field name
        node_data = self.return_array_or_field_values("node", node_data)

        if outlet_id is None:
            # verify that there is one and only one node with the status
            # FIXED_VALUE_BOUNDARY.
            possible_outlets = np.where(self.status_at_node == FIXED_VALUE_BOUNDARY)[0]

            if len(possible_outlets) > 1:
                raise ValueError(
                    "Model grid must only have one node with node status of FIXED_VALUE_BOUNDARY. This grid has %r"
                    % len(possible_outlets)
                )
            if len(possible_outlets) < 1:
                raise ValueError(
                    "Model grid must only have one node with node status of FIXED_VALUE_BOUNDARY. This grid has none"
                )

            outlet_id = possible_outlets

        elif outlet_id.size != 1 or (isinstance(outlet_id, np.ndarray) is False):
            # check that the value given by outlet_id is an integer
            raise ValueError("outlet_id must be a length 1 numpy array")
        else:
            # check that the node status at the node given by outlet_id is not
            # CLOSED_BOUNDARY
            if self.status_at_node[outlet_id] == CLOSED_BOUNDARY:
                raise ValueError(
                    "The node given by outlet_id must not have the status: CLOSED_BOUNDARY"
                )

        # now test that the method given is either 'D8' or 'D4'
        if adjacency_method != "D8":
            assert (
                adjacency_method == "D4"
            ), "Method must be either 'D8'(default) or 'D4'"

        # begin main code portion.
        # initialize list of core nodes and new nodes
        connected_nodes = list(outlet_id)
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
                    self.status_at_node[connected_orthogonal_nodes] != CLOSED_BOUNDARY
                ]
            )

            # if method is D8 (default), add the diagonal nodes.
            if adjacency_method == "D8":
                connected_diagonal_nodes = self.diagonal_adjacent_nodes_at_node[
                    newNodes
                ]
                potentialNewNodes.extend(
                    connected_diagonal_nodes[
                        self.status_at_node[connected_diagonal_nodes] != CLOSED_BOUNDARY
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
        is_not_connected_to_outlet = (self.status_at_node != CLOSED_BOUNDARY) & (
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
        """
        Set the boundary conditions for a watershed.
        All nodes with nodata_value are set to CLOSED_BOUNDARY
        (grid.status_at_node == 4). All nodes with data values
        are set to CORE_NODES (grid.status_at_node == 0), with
        the exception that the outlet node is set to a
        FIXED_VALUE_BOUNDARY (grid.status_at_node == 1).

        Note that the outer ring of the raster is set to CLOSED_BOUNDARY, even
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
        as illustrated:

            -9999. -9999. -9999. -9999.
            -9999.    67.     0. -9999.
            -9999.    67.    67. -9999.
            -9999. -9999. -9999. -9999.

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid((4,4))
        >>> rmg.status_at_node
        array([1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1], dtype=uint8)
        >>> node_data = np.array([-9999., -9999., -9999., -9999.,
        ...                      -9999.,    67.,    67., -9999.,
        ...                      -9999.,    67.,     0., -9999.,
        ...                      -9999., -9999., -9999., -9999.])
        >>> rmg.set_watershed_boundary_condition_outlet_coords(
        ...     (2, 2), node_data, -9999.)
        >>> rmg.status_at_node
        array([4, 4, 4, 4, 4, 0, 0, 4, 4, 0, 1, 4, 4, 4, 4, 4], dtype=uint8)

        LLCATS: BC
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
        self.status_at_node[outlet_node] = FIXED_VALUE_BOUNDARY

    def set_watershed_boundary_condition_outlet_id(
        self, outlet_id, node_data, nodata_value=-9999.0
    ):
        """
        Set the boundary conditions for a watershed.
        All nodes with nodata_value are set to CLOSED_BOUNDARY (4).
        All nodes with data values are set to CORE_NODES (0), with the
        exception that the outlet node is set to a FIXED_VALUE_BOUNDARY (1).

        Note that the outer ring of the raster is set to CLOSED_BOUNDARY, even
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
        >>> rmg = RasterModelGrid((4,4))
        >>> rmg.status_at_node
        array([1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1], dtype=uint8)
        >>> node_data = np.array([-9999., -9999., -9999., -9999.,
        ...                      -9999.,    67.,    67., -9999.,
        ...                      -9999.,    67.,     0., -9999.,
        ...                      -9999., -9999., -9999., -9999.])
        >>> outlet = rmg.set_watershed_boundary_condition_outlet_id(
        ...     10, node_data, -9999.)
        >>> rmg.status_at_node
        array([4, 4, 4, 4, 4, 0, 0, 4, 4, 0, 1, 4, 4, 4, 4, 4], dtype=uint8)

        LLCATS: BC
        """
        # get node_data if a field name
        node_data = self.return_array_or_field_values("node", node_data)

        # make ring of no data nodes
        self.set_closed_boundaries_at_grid_edges(True, True, True, True)

        # set no data nodes to inactive boundaries
        self.set_nodata_nodes_to_closed(node_data, nodata_value)

        # set the boundary condition (fixed value) at the outlet_node
        self.status_at_node[outlet_id] = FIXED_VALUE_BOUNDARY


def _is_closed_boundary(boundary_string):
    """Check if boundary string indicates a closed boundary.

    Helper function, probably depreciated due to changes in BC handling
    procedures (DEJH, May 14).
    """
    return boundary_string.lower() == "closed"


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


def from_dict(param_dict):
    """Create a RasterModelGrid from a dict-like object.

    Create a RasterModelGrid from the dictionary-like object, *param_dict*.
    Required keys of the dictionary are NUM_ROWS, NUM_COLS. Raises a KeyError
    if either of these are missing is given, use it as the
    HexModelGrid *dx* parameter, otherwise default to unit spacing.

    Deprecated in version 1.6.X. Will be removed in version 2.0.
    """
    msg = (
        "The non-class method version of 'from_dict' for RasterModelGrid "
        "was Deprecated in version 1.6.X. Will be removed in version 2.0."
    )
    warn(msg, DeprecationWarning)

    # Read and create basic raster grid
    try:
        nrows = int(param_dict["NUM_ROWS"])
        ncols = int(param_dict["NUM_COLS"])
        spacing = float(param_dict.get("GRID_SPACING", 1.0))
    except KeyError:
        raise
    except ValueError:
        raise
    else:
        grid = RasterModelGrid((nrows, ncols), xy_spacing=spacing)

    # Set boundaries
    left_boundary_type = param_dict.get("LEFT_BOUNDARY", "open")
    right_boundary_type = param_dict.get("RIGHT_BOUNDARY", "open")
    top_boundary_type = param_dict.get("TOP_BOUNDARY", "open")
    bottom_boundary_type = param_dict.get("BOTTOM_BOUNDARY", "open")
    grid.set_closed_boundaries_at_grid_edges(
        _is_closed_boundary(right_boundary_type),
        _is_closed_boundary(top_boundary_type),
        _is_closed_boundary(left_boundary_type),
        _is_closed_boundary(bottom_boundary_type),
    )

    # Return the created and initialized grid
    return grid


add_module_functions_to_class(RasterModelGrid, "raster_mappers.py", pattern="map_*")
add_module_functions_to_class(RasterModelGrid, "raster_gradients.py", pattern="calc_*")
add_module_functions_to_class(
    RasterModelGrid, "raster_set_status.py", pattern="set_status_at_node*"
)
