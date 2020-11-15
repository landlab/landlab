#! /usr/env/python
"""Python implementation of ModelGrid, a base class used to create and manage
grids for 2D numerical models.

Do NOT add new documentation here. Grid documentation is now built in a
semi-automated fashion. To modify the text seen on the web, edit the
files `docs/text_for_[gridfile].py.txt`.
"""
import fnmatch
from functools import lru_cache

import numpy as np
import xarray as xr

from landlab.utils.decorators import make_return_array_immutable

from ..core import load_params
from ..core.utils import add_module_functions_to_class
from ..field.graph_field import GraphFields
from ..layers.eventlayers import EventLayersMixIn
from ..layers.materiallayers import MaterialLayersMixIn
from ..utils.decorators import cache_result_in_object
from . import grid_funcs as gfuncs
from .decorators import (
    override_array_setitem_and_reset,
    return_id_array,
    return_readonly_id_array,
)
from .linkstatus import LinkStatus, set_status_at_link
from .nodestatus import NodeStatus

#: Indicates an index is, in some way, *bad*.
BAD_INDEX_VALUE = -1
# DEJH thinks the user should be able to override this value if they want

# Map names grid elements to the ModelGrid attribute that contains the count
# of that element in the grid.
_ARRAY_LENGTH_ATTRIBUTES = {
    "node": "number_of_nodes",
    "patch": "number_of_patches",
    "link": "number_of_links",
    "corner": "number_of_corners",
    "face": "number_of_faces",
    "cell": "number_of_cells",
    "active_link": "number_of_active_links",
    "active_face": "number_of_active_faces",
    "core_node": "number_of_core_nodes",
    "core_cell": "number_of_core_cells",
}

# Fields whose sizes can not change.
_SIZED_FIELDS = {"node", "link", "patch", "corner", "face", "cell"}


def _sort_points_into_quadrants(x, y, nodes):
    """Divide x, y points into quadrants.

    Divide points with locations given in the *x*, and *y* arrays into north,
    south, east, and west quadrants. Returns nodes contained in quadrants
    (west, east, north, south).

    Parameters
    ----------
    x : array_like
        X-coordinates of points.
    y : array_like
        Y-coordinates of points.
    nodes : array_like
        Nodes associated with points.

    Returns
    -------
    tuple of array_like
        Tuple of nodes in each coordinate. Nodes are grouped as
        (*east*, *north*, *west*, *south*).

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.base import _sort_points_into_quadrants
    >>> x = np.array([0, 1, 0, -1])
    >>> y = np.array([1, 0, -1, 0])
    >>> nodes = np.array([1, 2, 3, 4])
    >>> _sort_points_into_quadrants(x, y, nodes)
    (array([2]), array([1]), array([4]), array([3]))
    """
    above_x_axis = y > 0
    right_of_y_axis = x > 0
    closer_to_y_axis = np.abs(y) >= np.abs(x)

    north_nodes = nodes[above_x_axis & closer_to_y_axis]
    south_nodes = nodes[(~above_x_axis) & closer_to_y_axis]
    east_nodes = nodes[right_of_y_axis & (~closer_to_y_axis)]
    west_nodes = nodes[(~right_of_y_axis) & (~closer_to_y_axis)]

    return (east_nodes, north_nodes, west_nodes, south_nodes)


def _default_axis_names(n_dims):
    """Name of each axis.

    Parameters
    ----------
    n_dims : int
        Number of spatial dimensions.

    Returns
    -------
    tuple of str
        Name of each axis.

    Examples
    --------
    >>> from landlab.grid.base import _default_axis_names
    >>> _default_axis_names(1)
    ('x',)
    >>> _default_axis_names(2)
    ('x', 'y')
    >>> _default_axis_names(3)
    ('x', 'y', 'z')
    """
    _DEFAULT_NAMES = ("x", "y", "z")
    return _DEFAULT_NAMES[:n_dims]


def _default_axis_units(n_dims):
    """Unit names for each axis.

    Parameters
    ----------
    n_dims : int
        Number of spatial dimensions.

    Returns
    -------
    tuple of str
        Units of each axis.

    Examples
    --------
    >>> from landlab.grid.base import _default_axis_units
    >>> _default_axis_units(1)
    ('-',)
    >>> _default_axis_units(2)
    ('-', '-')
    >>> _default_axis_units(3)
    ('-', '-', '-')
    """
    return ("-",) * n_dims


def find_true_vector_from_link_vector_pair(L1, L2, b1x, b1y, b2x, b2y):
    r"""Separate a pair of links with vector values into x and y components.

    The concept here is that a pair of adjacent links attached to a node are
    projections of a 'true' but unknown vector. This function finds and returns
    the x and y components of this true vector. The trivial case is the
    situation in which the two links are orthogonal and aligned with the grid
    axes, in which case the vectors of these two links *are* the x and y
    components.

    Parameters
    ----------
    L1, L2 : float
        Values (magnitudes) associated with the two links
    b1x, b1y, b2x, b2y : float
        Unit vectors of the two links

    Returns
    -------
    ax, ay : float
        x and y components of the 'true' vector

    Notes
    -----
    The function does an inverse vector projection. Suppose we have a given
    'true' vector :math:`a`, and we want to project it onto two other lines
    with unit vectors (b1x,b1y) and (b2x,b2y). In the context of Landlab,
    the 'true' vector is some unknown vector quantity, which might for
    example represent the local water flow velocity. The lines represent two
    adjacent links in the grid.

    Let :math:`\mathbf{a}` be the true vector, :math:`\mathbf{B}` be a
    different vector with unit vector :math:`\mathbf{b}`, and :math:`L`
    be the scalar projection of *a* onto *B*. Then,

    ..math::

        L = \mathbf{a} \dot \mathbf{b} = a_x b_x + a_y b_y,

    where :math:`(a_x,a_y)` are the components of **a** and :math:`(b_x,b_y)`
    are the components of the unit vector **b**.

    In this case, we know *b* (the link unit vector), and we want to know the
    *x* and *y* components of **a**. The problem is that we have one equation
    and two unknowns (:math:`a_x` and :math:`a_y`). But we can solve this if
    we have *two* vectors, both of which are projections of **a**. Using the
    subscripts 1 and 2 to denote the two vectors, we can obtain equations for
    both :math:`a_x` and :math:`a_y`:

    ..math::

        a_x = L_1 / b_{1x} - a_y b_{1y} / b_{1x}

        a_y = L_2 / b_{2y} - a_x b_{2x} / b_{2y}

    Substituting the second into the first,

    ..math::

        a_x = [L_1/b_{1x}-L_2 b_{1y}/(b_{1x} b_{2y})] / [1-b_{1y} b_{2x}/(b_{1x} b_{2y})]

    Hence, we find the original vector :math:`(a_x,a_y)` from two links with
    unit vectors :math:`(b_{1x},b_{1y})` and :math:`(b_{2x},b_{2y})` and
    associated values :math:`L_1` and :math:`L_2`.

    Note that the above equations require that :math:`b_{1x}>0` and
    :math:`b_{2y}>0`. If this isn't the case, we invert the order of the two
    links, which requires :math:`b_{2x}>0` and :math:`b_{1y}>0`. If none of
    these conditions is met, then we have a degenerate case.

    Examples
    --------
    The following example represents the active links in a 7-node hexagonal
    grid, with just one core node. The 'true' vector has a magnitude of 5 units
    and an orientation of 30 degrees, pointing up and to the right (i.e., the
    postive-x and postive-y quadrant), so that its vector components are 4 (x)
    and 3 (y) (in other words, it is a 3-4-5 triangle). The values assigned to
    L below are the projection of that true vector onto the six link
    vectors. The algorithm should recover the correct vector component
    values of 4 and 3. The FOR loop examines each pair of links in turn.

    >>> import numpy as np
    >>> from landlab.grid.base import find_true_vector_from_link_vector_pair
    >>> bx = np.array([0.5, -0.5, -1., -0.5, 1., 0.5])
    >>> by = np.array([0.866, 0.866, 0., -0.866, 0., -0.866])
    >>> L = np.array([4.6, 0.6, -4., -4.6, 4., -0.6])
    >>> for i in range(5):
    ...     ax, ay = find_true_vector_from_link_vector_pair(
    ...         L[i], L[i+1], bx[i], by[i], bx[i+1], by[i+1])
    ...     round(ax,1), round(ay,1)
    (4.0, 3.0)
    (4.0, 3.0)
    (4.0, 3.0)
    (4.0, 3.0)
    (4.0, 3.0)
    """
    assert (b1x != 0 and b2y != 0) or (b2x != 0 and b1y != 0), "Improper unit vectors"

    if b1x != 0.0 and b2y != 0.0:
        ax = (L1 / b1x - L2 * (b1y / (b1x * b2y))) / (1.0 - (b1y * b2x) / (b1x * b2y))
        ay = L2 / b2y - ax * (b2x / b2y)
    elif b2x != 0.0 and b1y != 0.0:
        ax = (L2 / b2x - L1 * (b2y / (b2x * b1y))) / (1.0 - (b2y * b1x) / (b2x * b1y))
        ay = L1 / b1y - ax * (b1x / b1y)

    return ax, ay


class ModelGrid(GraphFields, EventLayersMixIn, MaterialLayersMixIn):

    """Base class for 2D structured or unstructured grids for numerical models.

    The idea is to have at least two inherited
    classes, RasterModelGrid and DelaunayModelGrid, that can create and
    manage grids. To this might be added a GenericModelGrid, which would
    be an unstructured polygonal grid that doesn't necessarily obey or
    understand the Delaunay triangulation, but rather simply accepts
    an input grid from the user. Also a :class:`~.HexModelGrid` for hexagonal.

    Attributes
    ----------
    at_node : dict-like
        Values at nodes.
    at_cell : dict-like
        Values at cells.
    at_link : dict-like
        Values at links.
    at_face : dict-like
        Values at faces.
    at_grid: dict-like
        Global values
    BAD_INDEX : int
        Indicates a grid element is undefined.
    BC_NODE_IS_CORE : int
        Indicates a node is *core*.
    BC_NODE_IS_FIXED_VALUE : int
        Indicates a boundary node has a fixed value.
    BC_NODE_IS_FIXED_GRADIENT : int
        Indicates a boundary node has a fixed gradient.
    BC_NODE_IS_LOOPED : int
        Indicates a boundary node is wrap-around.
    BC_NODE_IS_CLOSED : int
        Indicates a boundary node is closed
    BC_LINK_IS_ACTIVE : int
        Indicates a link is *active*, and can carry flux.
    BC_LINK_IS_FIXED : int
        Indicates a link has a fixed gradient value, and behaves as a boundary
    BC_LINK_IS_INACTIVE : int
        Indicates a link is *inactive*, and cannot carry flux.

    Other Parameters
    ----------------
    axis_name : tuple, optional
        Name of axes
    axis_units : tuple, optional
        Units of coordinates
    """

    #: Indicates a node is *bad index*.
    BAD_INDEX = BAD_INDEX_VALUE

    #: Indicates a node is *core*.
    BC_NODE_IS_CORE = NodeStatus.CORE
    #: Indicates a boundary node has a fixed value.
    BC_NODE_IS_FIXED_VALUE = NodeStatus.FIXED_VALUE
    #: Indicates a boundary node has a fixed gradient.
    BC_NODE_IS_FIXED_GRADIENT = NodeStatus.FIXED_GRADIENT
    #: Indicates a boundary node is wrap-around.
    BC_NODE_IS_LOOPED = NodeStatus.LOOPED
    #: Indicates a boundary node is closed
    BC_NODE_IS_CLOSED = NodeStatus.CLOSED

    #: Indicates a link is *active*, and can carry flux
    BC_LINK_IS_ACTIVE = LinkStatus.ACTIVE
    #: Indicates a link has a fixed gradient value, and behaves as a boundary
    BC_LINK_IS_FIXED = LinkStatus.FIXED
    #: Indicates a link is *inactive*, and cannot carry flux
    BC_LINK_IS_INACTIVE = LinkStatus.INACTIVE

    #: Grid elements on which fields can be placed.
    VALID_LOCATIONS = ("node", "link", "patch", "corner", "face", "cell", "grid")

    at_node = {}  # : Values defined at nodes
    at_link = {}  # : Values defined at links
    at_patch = {}  # : Values defined at patches
    at_corner = {}  # : Values defined at corners
    at_face = {}  # : Values defined at faces
    at_cell = {}  # : Values defined at cells
    at_grid = {}  # : Values defined at grid

    @classmethod
    def from_file(cls, file_like):
        """Create grid from a file-like object.

        File to load either as a file-like object, path to an existing file, or
        the contents of a file as a string.

        Parameters
        ----------
        file_like :
            File-like object, filepath, or string.

        Examples
        --------
        >>> from io import StringIO
        >>> from landlab import RasterModelGrid
        >>> filelike = StringIO('''
        ... shape:
        ...     - 3
        ...     - 4
        ... xy_spacing: 2
        ... ''')
        >>> grid = RasterModelGrid.from_file(filelike)
        >>> grid.x_of_node
        array([ 0.,  2.,  4.,  6.,  0.,  2.,  4.,  6.,  0.,  2.,  4.,  6.])
        >>> grid.y_of_node
        array([ 0.,  0.,  0.,  0.,  2.,  2.,  2.,  2.,  4.,  4.,  4.,  4.])
        """
        params = load_params(file_like)
        return cls.from_dict(params)

    @classmethod
    def from_dict(cls, params):
        """Create grid from dictionary.

        Parameters
        ----------
        params : dictionary
            Dictionary of required parameters to create a model grid.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> params = {"shape": (3,4), "xy_spacing": 2}
        >>> grid = RasterModelGrid.from_dict(params)
        >>> grid.x_of_node
        array([ 0.,  2.,  4.,  6.,  0.,  2.,  4.,  6.,  0.,  2.,  4.,  6.])
        >>> grid.y_of_node
        array([ 0.,  0.,  0.,  0.,  2.,  2.,  2.,  2.,  4.,  4.,  4.,  4.])
        """
        return cls(**params)

    def __init__(self, **kwds):
        axis_units = kwds.pop("xy_axis_units", "-")
        axis_name = kwds.pop("xy_axis_name", ("x", "y"))

        super().__init__()

        self.new_field_location("node", self.number_of_nodes)
        self.new_field_location("link", self.number_of_links)
        self.new_field_location("patch", self.number_of_patches)
        self.new_field_location("corner", self.number_of_corners)
        self.new_field_location("face", self.number_of_faces)
        self.new_field_location("cell", self.number_of_cells)
        self.new_field_location("grid", None)
        self.default_group = "node"

        # self.axis_name = kwds.get("axis_name", _default_axis_names(self.ndim))
        # self.axis_units = kwds.get("axis_units", _default_axis_units(self.ndim))

        self._ref_coord = tuple(kwds.get("xy_of_reference", (0.0, 0.0)))
        self._link_length = None
        self._all_node_distances_map = None
        self._all_node_azimuths_map = None
        self.bc_set_code = 0

        self._axis_units = tuple(np.broadcast_to(axis_units, self.ndim))
        self._axis_name = tuple(np.broadcast_to(axis_name, self.ndim))

    def fields(self, include="*", exclude=None):
        """List of fields held by the grid.

        Parameters
        ----------
        include : str, or iterable of str, optional
            Glob-style pattern for field names to include.
        exclude : str, or iterable of str, optional
            Glob-style pattern for field names to exclude.

        Returns
        -------
        set
            Filtered set of canonical field names held by the grid

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4))
        >>> _ = grid.add_full("elevation", 3.0, at="node")
        >>> _ = grid.add_full("elevation", 4.0, at="link")
        >>> _ = grid.add_full("temperature", 5.0, at="node")

        >>> sorted(grid.fields())
        ['at_link:elevation', 'at_node:elevation', 'at_node:temperature']
        >>> sorted(grid.fields(include="at_node*"))
        ['at_node:elevation', 'at_node:temperature']
        >>> sorted(grid.fields(include="at_node*", exclude="*temp*"))
        ['at_node:elevation']
        """
        if isinstance(include, str):
            include = [include]
        if isinstance(exclude, str):
            exclude = [exclude]

        canonical_names = set()
        for at in self.groups | set(["layer"]):
            canonical_names.update(["at_{0}:{1}".format(at, name) for name in self[at]])

        names = set()
        for pattern in include:
            names.update(fnmatch.filter(canonical_names, pattern))
        for pattern in exclude or []:
            names.difference_update(fnmatch.filter(canonical_names, pattern))

        return names

    def as_dataset(self, include="*", exclude=None):
        """Create an xarray Dataset representation of a grid.

        Parameters
        ----------
        include : str or iterable or str
            Glob-style patterns of fields to include in the dataset.
        exclude : str or iterable or str
            Glob-style patterns of fields to exclude from the dataset.

        Returns
        -------
        Dataset
            An xarray Dataset representation of a *ModelGrid*.
        """
        names = self.fields(include=include, exclude=exclude)

        layer_names = set([name for name in names if name.startswith("at_layer")])
        names.difference_update(layer_names)

        data = {}
        for name in names:
            dim, field_name = name[len("at_") :].split(":")
            data[name] = ((dim,), getattr(self, "at_" + dim)[field_name])

        for name in layer_names:
            dim, field_name = name[len("at_") :].split(":")
            data[name] = (("layer", "cell"), self.at_layer[field_name])

        data["status_at_node"] = (("node",), self.status_at_node)

        return xr.Dataset(data)

    @property
    def xy_of_reference(self):
        """Return the coordinates (x, y) of the reference point.

        For RasterModelGrid and HexModelGrid the reference point is the
        minimum of x_of_node and of y_of_node. By default it is (0, 0). For
        VoronoiDelaunayGrid the reference point is (0, 0). For RadialModelGrid
        it is the (x, y) of the center point.

        The intention of these coordinates is to provide a method to store
        the large float values of projected coordinates.

        Example
        -------

        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid((4, 5),
        ...       xy_of_reference = (12345, 678910))
        >>> rmg.xy_of_reference
        (12345, 678910)
        >>> rmg.xy_of_reference = (98765, 43210)
        >>> rmg.xy_of_reference
        (98765, 43210)
        """
        return self._ref_coord

    @xy_of_reference.setter
    def xy_of_reference(self, new_xy_of_reference):
        """Set a new value for the model grid xy_of_reference."""
        self._ref_coord = (new_xy_of_reference[0], new_xy_of_reference[1])

    @property
    def ndim(self):
        """Number of spatial dimensions of the grid.

        LLCATS: GINF
        """
        return 2

    @property
    @override_array_setitem_and_reset("reset_status_at_node")
    def status_at_node(self):
        """Get array of the boundary status for each node.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import LinkStatus, NodeStatus, RasterModelGrid
        >>> mg = RasterModelGrid((4, 5))
        >>> mg.status_at_node.reshape((4, 5))
        array([[1, 1, 1, 1, 1],
               [1, 0, 0, 0, 1],
               [1, 0, 0, 0, 1],
               [1, 1, 1, 1, 1]], dtype=uint8)
        >>> np.any(mg.status_at_link == LinkStatus.FIXED)
        False

        >>> mg.status_at_node[mg.nodes_at_left_edge] = NodeStatus.FIXED_GRADIENT
        >>> mg.status_at_node.reshape((4, 5))
        array([[2, 1, 1, 1, 1],
               [2, 0, 0, 0, 1],
               [2, 0, 0, 0, 1],
               [2, 1, 1, 1, 1]], dtype=uint8)
        >>> np.any(mg.status_at_link == LinkStatus.FIXED)  # links auto-update
        True

        LLCATS: NINF BC
        """
        return self._node_status

    @status_at_node.setter
    def status_at_node(self, new_status):
        """Set the array of node boundary statuses."""
        self._node_status[:] = new_status[:]
        self.reset_status_at_node()

    @property
    @cache_result_in_object()
    @return_readonly_id_array
    def active_adjacent_nodes_at_node(self):
        """Adjacent nodes for each grid node.

        For each grid node, get the adjacent nodes ordered
        counterclockwise starting from the positive x axis.

        Examples
        --------
        >>> from landlab import RasterModelGrid, HexModelGrid
        >>> grid = RasterModelGrid((4, 5))

        >>> grid.active_adjacent_nodes_at_node[(-1, 6, 2), ]
        array([[-1, -1, -1, -1],
               [ 7, 11,  5,  1],
               [-1,  7, -1, -1]])

        Setting a node to closed causes all links touching it to
        be inactive.

        >>> grid.status_at_node[6] = grid.BC_NODE_IS_CLOSED
        >>> grid.active_adjacent_nodes_at_node[(-1, 6, 2), ]
        array([[-1, -1, -1, -1],
               [-1, -1, -1, -1],
               [-1,  7, -1, -1]])

        >>> grid.active_adjacent_nodes_at_node[7]
        array([ 8, 12, -1,  2])
        >>> grid.active_adjacent_nodes_at_node[2]
        array([-1,  7, -1, -1])

        >>> grid = HexModelGrid((3, 2))
        >>> grid.status_at_node[0] = grid.BC_NODE_IS_CLOSED
        >>> grid.active_adjacent_nodes_at_node
        array([[-1, -1, -1, -1, -1, -1],
               [-1,  3, -1, -1, -1, -1],
               [ 3, -1, -1, -1, -1, -1],
               [ 4,  6,  5,  2, -1,  1],
               [-1,  3, -1, -1, -1, -1],
               [-1, -1,  3, -1, -1, -1],
               [-1,  3, -1, -1, -1, -1]])

        LLCATS: NINF CONN BC
        """
        return np.choose(
            self.status_at_link[self.links_at_node] == LinkStatus.ACTIVE,
            (-1, self.adjacent_nodes_at_node),
        )

    @property
    @make_return_array_immutable
    @cache_result_in_object()
    def active_link_dirs_at_node(self):
        """Link flux directions at each node: 1=incoming flux, -1=outgoing
        flux, 0=no flux. Note that inactive links receive zero, but active and
        fixed links are both reported normally.

        Returns
        -------
        (NODES, LINKS) ndarray of int
            Link directions relative to the nodes of a grid. The shape of the
            matrix will be number of nodes rows by max number of links per
            node. A zero indicates no link at this position.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 3))
        >>> grid.status_at_node[grid.nodes_at_left_edge] = grid.BC_NODE_IS_CLOSED
        >>> grid.active_link_dirs_at_node # doctest: +NORMALIZE_WHITESPACE
        array([[ 0,  0,  0,  0], [ 0, -1,  0,  0], [ 0,  0,  0,  0],
               [ 0,  0,  0,  0], [-1, -1,  0,  1], [ 0,  0,  1,  0],
               [ 0,  0,  0,  0], [-1, -1,  0,  1], [ 0,  0,  1,  0],
               [ 0,  0,  0,  0], [ 0,  0,  0,  1], [ 0,  0,  0,  0]],
               dtype=int8)

        LLCATS: NINF LINF CONN
        """
        return np.choose(
            self.link_status_at_node == LinkStatus.ACTIVE, (0, self.link_dirs_at_node)
        )

    @property
    @make_return_array_immutable
    @cache_result_in_object()
    def link_status_at_node(self):
        return self.status_at_link[self.links_at_node]

    @property
    @return_readonly_id_array
    @cache_result_in_object()
    def core_nodes(self):
        """Get array of core nodes.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((4, 5))
        >>> mg.core_nodes
        array([ 6,  7,  8, 11, 12, 13])

        LLCATS: NINF BC
        """
        return np.where(self.status_at_node == NodeStatus.CORE)[0]

    @property
    @return_readonly_id_array
    def boundary_nodes(self):
        """Get array of boundary nodes.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((4, 5))
        >>> mg.boundary_nodes
        array([ 0,  1,  2,  3,  4,  5,  9, 10, 14, 15, 16, 17, 18, 19])

        LLCATS: NINF BC
        """
        try:
            return self._boundary_nodes
        except AttributeError:
            (boundary_node_ids,) = np.where(self._node_status != NodeStatus.CORE)
            return boundary_node_ids

    @property
    @return_readonly_id_array
    def open_boundary_nodes(self):
        """Get array of open boundary nodes.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((4, 5))
        >>> for edge in (mg.nodes_at_left_edge, mg.nodes_at_right_edge,
        ...              mg.nodes_at_bottom_edge):
        ...     mg.status_at_node[edge] = mg.BC_NODE_IS_CLOSED
        >>> mg.open_boundary_nodes
        array([16, 17, 18])

        LLCATS: NINF BC
        """
        (open_boundary_node_ids,) = np.where(
            (self._node_status != self.BC_NODE_IS_CLOSED)
            & (self._node_status != self.BC_NODE_IS_CORE)
        )
        return open_boundary_node_ids

    @property
    @return_readonly_id_array
    def closed_boundary_nodes(self):
        """Get array of closed boundary nodes.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((4, 5))
        >>> mg.status_at_node[mg.nodes_at_top_edge] = mg.BC_NODE_IS_CLOSED
        >>> mg.closed_boundary_nodes
        array([15, 16, 17, 18, 19])

        LLCATS: NINF BC
        """
        (closed_boundary_node_ids,) = np.where(
            self._node_status == self.BC_NODE_IS_CLOSED
        )
        return closed_boundary_node_ids

    @property
    @return_readonly_id_array
    def fixed_gradient_boundary_nodes(self):
        """Get array of fixed gradient boundary nodes.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((4, 5))
        >>> mg.status_at_node[mg.nodes_at_top_edge] = mg.BC_NODE_IS_FIXED_GRADIENT
        >>> mg.fixed_gradient_boundary_nodes
        array([15, 16, 17, 18, 19])

        LLCATS: NINF BC
        """
        (fixed_gradient_boundary_node_ids,) = np.where(
            self._node_status == self.BC_NODE_IS_FIXED_GRADIENT
        )
        return fixed_gradient_boundary_node_ids

    @property
    @return_readonly_id_array
    def fixed_gradient_boundary_node_fixed_link(self):
        """An array of the fixed_links connected to fixed gradient boundary
        nodes.

        Note that on a raster, some nodes (notably the corners) can be
        `NodeStatus.FIXED_GRADIENT`, but not have a true `LinkStatus.FIXED` neighboring
        link. In such cases, the link returned will be a closed link joining
        the corner node to a neighboring `NodeStatus.FIXED_GRADIENT` node (see
        example).

        An AssertionError will be raised if for some reason a
        `NodeStatus.FIXED_GRADIENT` node exists which has neither a
        `NodeStatus.FIXED_GRADIENT` neighbor, or a `LinkStatus.FIXED`.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4))
        >>> leftedge = grid.nodes_at_left_edge
        >>> grid.status_at_node[leftedge] = grid.BC_NODE_IS_FIXED_GRADIENT
        >>> grid.fixed_gradient_boundary_nodes
        array([0, 4, 8])
        >>> grid.fixed_gradient_boundary_node_fixed_link
        array([ 3,  7, 10])
        """
        try:
            return self._fixed_gradient_boundary_node_links
        except AttributeError:
            self._create_fixed_gradient_boundary_node_links()
            return self._fixed_gradient_boundary_node_links

    @property
    @return_readonly_id_array
    def fixed_gradient_boundary_node_anchor_node(self):
        """Returns the node at the other end of the fixed link for a fixed
        gradient boundary node.

        Degenerate `NodeStatus.FIXED_GRADIENT` nodes (e.g., corners) are handled as
        in :func:`fixed_gradient_boundary_node_fixed_link`, by pointing to a
        neighboring `NodeStatus.FIXED_GRADIENT` node.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4))
        >>> leftedge = grid.nodes_at_left_edge
        >>> grid.status_at_node[leftedge] = grid.BC_NODE_IS_FIXED_GRADIENT
        >>> grid.fixed_gradient_boundary_nodes
        array([0, 4, 8])
        >>> grid.fixed_gradient_boundary_node_fixed_link
        array([ 3,  7, 10])
        >>> grid.fixed_gradient_boundary_node_anchor_node
        array([4, 5, 4])
        """
        try:
            return self._fixed_gradient_boundary_node_anchor_node
        except AttributeError:
            self._create_fixed_gradient_boundary_node_anchor_node()
            return self._fixed_gradient_boundary_node_anchor_node

    def _create_fixed_gradient_boundary_node_links(self):
        """Builds a data structure to hold the fixed_links which control the
        values of any `NodeStatus.FIXED_GRADIENT` nodes in the grid.

        An AssertionError will be raised if for some reason a
        `NodeStatus.FIXED_GRADIENT` node exists which has neither a
        `NodeStatus.FIXED_GRADIENT` neighbor, or a `LinksStatus.FIXED`.
        """
        self._fixed_grad_links_created = True
        self._fixed_gradient_boundary_node_links = np.empty_like(
            self.fixed_gradient_boundary_nodes, dtype=int
        )
        fix_nodes = self.fixed_gradient_boundary_nodes
        neighbor_links = self.links_at_node[fix_nodes]  # -1s
        boundary_exists = self.link_dirs_at_node[fix_nodes]
        # next line retains -1 indexes
        link_stat_badind = self.status_at_link[neighbor_links] == LinkStatus.FIXED
        true_connection = np.logical_and(link_stat_badind, boundary_exists)
        true_fix_nodes = true_connection.sum(axis=1).astype(bool)
        self._fixed_gradient_boundary_node_links[true_fix_nodes] = neighbor_links[
            true_connection
        ]
        # resolve any corner nodes
        neighbor_nodes = self.adjacent_nodes_at_node[fix_nodes]  # BAD_INDEX_VALUEs
        neighbor_nodes[neighbor_nodes == self.BAD_INDEX] = -1
        fixed_grad_neighbor = np.logical_and(
            (self.status_at_node[neighbor_nodes] == NodeStatus.FIXED_GRADIENT),
            boundary_exists,
        )
        # ^True when NodeStatus.FIXED_GRADIENT for real
        # winnow it down to only one possibility for fixed_grad neighbor:
        which_neighbor = np.argmax(fixed_grad_neighbor, axis=1)
        indexing_range = np.arange(fixed_grad_neighbor.shape[0])
        a_link_to_fixed_grad = neighbor_links[indexing_range, which_neighbor]
        corners = np.logical_not(true_fix_nodes)
        assert np.all(fixed_grad_neighbor[indexing_range, which_neighbor][corners])
        self._fixed_gradient_boundary_node_links[corners] = a_link_to_fixed_grad[
            corners
        ]

    def _create_fixed_gradient_boundary_node_anchor_node(self):
        """Builds a data structure to hold the nodes which anchor the values of
        any `NodeStatus.FIXED_GRADIENT` nodes in the grid, i.e., those at the other
        ends of the `LinkStatus.FIXED`.

        An AssertionError will be raised if for some reason a
        `NodeStatus.FIXED_GRADIENT` node exists which has neither a
        `NodeStatus.FIXED_GRADIENT` neighbor, or a `LinkStatus.FIXED`.
        """
        self._fixed_grad_links_created = True
        fix_grad_nodes = self.fixed_gradient_boundary_nodes
        self._fixed_gradient_boundary_node_anchor_node = np.empty_like(fix_grad_nodes)
        heads_and_tails = np.empty((fix_grad_nodes.size, 2))
        which_one = np.empty_like(heads_and_tails, dtype=bool)
        heads_and_tails[:, 0] = self.node_at_link_head[
            self.fixed_gradient_boundary_node_fixed_link
        ]
        heads_and_tails[:, 1] = self.node_at_link_tail[
            self.fixed_gradient_boundary_node_fixed_link
        ]
        which_one[:, 0] = heads_and_tails[:, 0] == fix_grad_nodes
        which_one[:, 1] = heads_and_tails[:, 1] == fix_grad_nodes
        assert np.all(which_one.sum(axis=1) == 1)
        self._fixed_gradient_boundary_node_anchor_node = heads_and_tails[
            np.logical_not(which_one)
        ]

    @property
    @return_readonly_id_array
    @cache_result_in_object()
    def fixed_value_boundary_nodes(self):
        """Get array of fixed value boundary nodes.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))

        Initially all the perimeter nodes are fixed value boundary.

        >>> grid.fixed_value_boundary_nodes
        array([ 0,  1,  2,  3,  4, 5,  9, 10, 14, 15, 16, 17, 18, 19])

        Set left, right, and bottom edges to closed.

        >>> for edge in (grid.nodes_at_left_edge, grid.nodes_at_right_edge,
        ...              grid.nodes_at_bottom_edge):
        ...     grid.status_at_node[edge] = grid.BC_NODE_IS_CLOSED

        Now nodes on just the top edge are fixed.

        >>> grid.fixed_value_boundary_nodes
        array([16, 17, 18])

        LLCATS: NINF BC
        """
        return np.where(self._node_status == NodeStatus.FIXED_VALUE)[0]

    @property
    @return_readonly_id_array
    @cache_result_in_object()
    def active_faces(self):
        """Get array of active faces.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4))
        >>> grid.active_faces
        array([0, 1, 2, 3, 4, 5, 6])

        >>> grid.status_at_node[6] = grid.BC_NODE_IS_CLOSED
        >>> grid.active_faces
        array([0, 2, 5])

        LLCATS: FINF BC
        """
        return self.face_at_link[self.active_links]

    @property
    @return_readonly_id_array
    @cache_result_in_object()
    def active_links(self):
        """Get array of active links.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4))
        >>> grid.active_links
        array([ 4,  5,  7,  8,  9, 11, 12])

        LLCATS: LINF BC
        """
        return np.where(self.status_at_link == LinkStatus.ACTIVE)[0]

    @property
    @return_readonly_id_array
    @cache_result_in_object()
    def fixed_links(self):
        """Get array of fixed links.

        Examples
        --------
        >>> from landlab import NodeStatus, RasterModelGrid
        >>> grid = RasterModelGrid((3, 4))
        >>> grid.status_at_node # doctest: +NORMALIZE_WHITESPACE
        array([1, 1, 1, 1,
               1, 0, 0, 1,
               1, 1, 1, 1], dtype=uint8)
        >>> grid.fixed_links.size
        0

        >>> grid.status_at_node[:4] = NodeStatus.FIXED_GRADIENT
        >>> grid.status_at_node # doctest: +NORMALIZE_WHITESPACE
        array([2, 2, 2, 2,
               1, 0, 0, 1,
               1, 1, 1, 1], dtype=uint8)
        >>> grid.fixed_links
        array([4, 5])

        LLCATS: LINF BC
        """
        return np.where(self.status_at_link == LinkStatus.FIXED)[0]

    @return_id_array
    def link_with_node_status(self, status_at_tail=None, status_at_head=None):
        """Links with a given node status.

        Parameters
        ----------
        status_at_tail : NodeStatus, optional
            Status of the link tail node.
        status_at_head : NodeStatus, optional
            Status of the link head node.

        Returns
        -------
        array of int
            Links with the given tail and head node statuses.

        Examples
        --------
        >>> from landlab import RasterModelGrid, NodeStatus
        >>> grid = RasterModelGrid((4, 5))

        >>> grid.status_at_node[13] = NodeStatus.FIXED_VALUE
        >>> grid.status_at_node[2] = NodeStatus.CLOSED
        >>> grid.link_with_node_status(
        ...     status_at_tail=NodeStatus.CORE, status_at_head=NodeStatus.CORE
        ... )
        array([10, 11, 14, 15, 19])
        >>> grid.link_with_node_status(
        ...     status_at_tail=NodeStatus.CORE, status_at_head=NodeStatus.FIXED_VALUE
        ... )
        array([12, 16, 20, 23, 24])
        >>> grid.link_with_node_status(
        ...     status_at_tail=NodeStatus.FIXED_VALUE, status_at_head=NodeStatus.CORE
        ... )
        array([ 5,  7,  9, 18])

        >>> grid.link_with_node_status(status_at_head=NodeStatus.CORE)
        array([ 5,  6,  7,  9, 10, 11, 14, 15, 18, 19])
        >>> grid.link_with_node_status(status_at_tail=NodeStatus.CORE)
        array([10, 11, 12, 14, 15, 16, 19, 20, 23, 24])
        >>> grid.link_with_node_status()
        array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
               17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30])
        """
        masks = []
        if status_at_tail is not None:
            masks.append(self.status_at_node[self.node_at_link_tail] == status_at_tail)
        if status_at_head is not None:
            masks.append(self.status_at_node[self.node_at_link_head] == status_at_head)

        if len(masks) == 0:
            return np.arange(self.number_of_links, dtype=int)
        if len(masks) == 1:
            return np.where(masks[0])[0]
        else:
            return np.where(masks[0] & masks[1])[0]

    @return_id_array
    def link_with_angle(self, angle, in_degrees=False):
        """Return array of IDs of links with given angle.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> grid = HexModelGrid((3, 3))
        >>> grid.link_with_angle(0.0)
        array([  0,  1,  8,  9, 10, 17, 18])
        >>> grid.link_with_angle(60.0, in_degrees=True)
        array([  3,  5,  7, 11, 13, 15])
        >>> grid.link_with_angle(2.0944)  # 120 degrees
        array([  2,  4,  6, 12, 14, 16])
        >>> len(grid.link_with_angle(0.5236))  # no links at 30 deg
        0
        >>> grid = HexModelGrid((3, 3), orientation='vertical')
        >>> grid.link_with_angle(30.0, in_degrees=True)
        array([  1,  3,  8, 10, 15, 17])
        >>> grid.link_with_angle(1.5708)  # 90 degrees
        array([ 2,  5,  6,  9, 12, 13, 16])
        >>> grid.link_with_angle(330.0, in_degrees=True)
        array([ 0,  4,  7, 11, 14, 18])
        >>> len(grid.link_with_angle(60.0, in_degrees=True))  # none at 60 deg
        0
        """
        if in_degrees:
            angle = np.deg2rad(angle % 360.0)
        return np.where(np.isclose(self.angle_of_link, angle))[0]

    @property
    @cache_result_in_object()
    @return_readonly_id_array
    def node_at_core_cell(self):
        """Get array of nodes associated with core cells.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))

        Initially each cell's node is core.

        >>> grid.node_at_core_cell
        array([ 6,  7,  8,
               11, 12, 13])

        Setting a node to closed causes means its cell is also
        "closed".

        >>> grid.status_at_node[8] = grid.BC_NODE_IS_CLOSED
        >>> grid.node_at_core_cell
        array([ 6,  7, 11, 12, 13])

        LLCATS: NINF CINF BC CONN
        """
        return np.where(self.status_at_node == NodeStatus.CORE)[0]

    @property
    @make_return_array_immutable
    @cache_result_in_object()
    def core_cells(self):
        """Get array of core cells.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))

        Initially all of the cells are "core".

        >>> grid.core_cells
        array([0, 1, 2,
               3, 4, 5])

        Setting a node to closed causes its cell to no longer be core.

        >>> grid.status_at_node[8] = grid.BC_NODE_IS_CLOSED
        >>> grid.core_cells
        array([0, 1, 3, 4, 5])

        LLCATS: CINF BC
        """
        return self.cell_at_node[self.core_nodes]

    @property
    def number_of_active_faces(self):
        """Total number of active faces.

        Returns
        -------
        int
            Total number of active faces in the grid.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4))
        >>> grid.number_of_active_faces
        7

        The number of active faces is updated when a node status changes.

        >>> grid.status_at_node[6] = grid.BC_NODE_IS_CLOSED
        >>> grid.number_of_active_faces
        3

        LLCATS: FINF BC
        """
        return self.active_faces.size

    @property
    def number_of_core_nodes(self):
        """Number of core nodes.

        The number of core nodes on the grid (i.e., excluding all boundary
        nodes).

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.number_of_core_nodes
        6

        >>> grid.status_at_node[7] = grid.BC_NODE_IS_CLOSED
        >>> grid.number_of_core_nodes
        5

        LLCATS: NINF BC
        """
        return self.core_nodes.size

    @property
    def number_of_core_cells(self):
        """Number of core cells.

        A core cell excludes all boundary cells.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.number_of_core_cells
        6

        >>> grid.status_at_node[7] = grid.BC_NODE_IS_CLOSED
        >>> grid.number_of_core_cells
        5

        LLCATS: CINF BC
        """
        return self.core_cells.size

    @property
    def number_of_active_links(self):
        """Number of active links.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((4, 5))
        >>> mg.number_of_active_links
        17
        >>> for edge in (mg.nodes_at_left_edge, mg.nodes_at_right_edge,
        ...              mg.nodes_at_bottom_edge):
        ...     mg.status_at_node[edge] = mg.BC_NODE_IS_CLOSED
        >>> mg.number_of_active_links
        10

        LLCATS: LINF BC
        """
        return self.active_links.size

    @property
    def number_of_fixed_links(self):
        """Number of fixed links.

        Examples
        --------
        >>> from landlab import NodeStatus, RasterModelGrid
        >>> mg = RasterModelGrid((4, 5))
        >>> mg.number_of_fixed_links
        0
        >>> mg.status_at_node[mg.nodes_at_top_edge] = NodeStatus.FIXED_GRADIENT
        >>> mg.number_of_fixed_links
        3

        LLCATS: LINF BC
        """
        return self.fixed_links.size

    def number_of_elements(self, name):
        """Number of instances of an element.

        Get the number of instances of a grid element in a grid.

        Parameters
        ----------
        name : {'node', 'cell', 'link', 'face', 'core_node', 'core_cell',
                'active_link', 'active_face'}
            Name of the grid element.

        Returns
        -------
        int
            Number of elements in the grid.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((4, 5))
        >>> mg.number_of_elements('node')
        20
        >>> mg.number_of_elements('core_cell')
        6
        >>> mg.number_of_elements('link')
        31
        >>> mg.number_of_elements('active_link')
        17
        >>> mg.status_at_node[8] = mg.BC_NODE_IS_CLOSED
        >>> mg.number_of_elements('link')
        31
        >>> mg.number_of_elements('active_link')
        13

        LLCATS: GINF
        """
        try:
            return getattr(self, _ARRAY_LENGTH_ATTRIBUTES[name])
        except KeyError:
            raise TypeError("{name}: element name not understood".format(name=name))

    @make_return_array_immutable
    def node_axis_coordinates(self, axis=0):
        """Get the coordinates of nodes along a particular axis.

        Return node coordinates from a given *axis* (defaulting to 0). Axis
        numbering is the same as that for numpy arrays. That is, the zeroth
        axis is along the rows, and the first along the columns.

        Parameters
        ----------
        axis : int, optional
            Coordinate axis.

        Returns
        -------
        ndarray
            Coordinates of nodes for a given axis.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.node_axis_coordinates(0) # doctest: +NORMALIZE_WHITESPACE
        array([ 0., 0., 0., 0., 0.,
                1., 1., 1., 1., 1.,
                2., 2., 2., 2., 2.,
                3., 3., 3., 3., 3.])
        >>> grid.node_axis_coordinates(1) # doctest: +NORMALIZE_WHITESPACE
        array([ 0., 1., 2., 3., 4.,
                0., 1., 2., 3., 4.,
                0., 1., 2., 3., 4.,
                0., 1., 2., 3., 4.])

        LLCATS: GINF NINF MEAS
        """
        AXES = ("node_y", "node_x")
        try:
            return getattr(self, AXES[axis])
        except IndexError:
            raise ValueError("'axis' entry is out of bounds")

    @property
    def axis_units(self):
        """Get units for each axis.

        Returns
        -------
        tuple of str
            The units (as a string) for each of a grid's coordinates.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((4, 5), xy_spacing=(3., 2.))
        >>> mg.axis_units
        ('-', '-')
        >>> mg.axis_units = ("degrees_north", "degrees_east")
        >>> mg.axis_units
        ('degrees_north', 'degrees_east')
        >>> mg.axis_units = "m"
        >>> mg.axis_units
        ('m', 'm')

        LLCATS: GINF
        """
        return self._axis_units

    @axis_units.setter
    def axis_units(self, new_units):
        """Set the units for each coordinate axis."""
        new_units = np.broadcast_to(new_units, self.ndim)
        if len(new_units) != self.ndim:
            raise ValueError("length of units does not match grid dimension")
        self._axis_units = tuple(new_units)

    @property
    def axis_name(self):
        """Get the name of each coordinate axis.

        Returns
        -------
        tuple of str
            The names of each axis.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.axis_name
        ('x', 'y')
        >>> grid.axis_name = ('lon', 'lat')
        >>> grid.axis_name
        ('lon', 'lat')

        LLCATS: GINF
        """
        return self._axis_name

    @axis_name.setter
    def axis_name(self, new_names):
        """Set the names of a grid's coordinate axes.

        Raises
        ------
        ValueError
            If the number of dimension do not match.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.axis_name = ("y", "x")
        >>> grid.axis_name = ("lon", "lat")
        >>> grid.axis_name
        ('lon', 'lat')
        """
        new_names = np.broadcast_to(new_names, self.ndim)
        if len(new_names) != self.ndim:
            raise ValueError("length of names does not match grid dimension")
        self._axis_name = tuple(new_names)

    @property
    @make_return_array_immutable
    @cache_result_in_object()
    def status_at_link(self):
        """Get array of the status of all links.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((4, 5))
        >>> mg.status_at_node[mg.nodes_at_left_edge] = mg.BC_NODE_IS_CLOSED
        >>> mg.status_at_node[mg.nodes_at_right_edge] = mg.BC_NODE_IS_FIXED_GRADIENT
        >>> mg.status_at_link # doctest: +NORMALIZE_WHITESPACE
        array([4, 4, 4, 4, 4, 0, 0, 0, 4, 4, 0, 0, 2, 4, 0, 0, 0, 4, 4, 0, 0,
               2, 4, 0, 0, 0, 4, 4, 4, 4, 4], dtype=uint8)

        LLCATS: BC LINF
        """
        return set_status_at_link(self.status_at_node[self.nodes_at_link])

    @property
    @make_return_array_immutable
    @cache_result_in_object()
    def angle_of_link_about_head(self):
        """Find and return the angle of a link about the node at the link head.

        Because links have direction, their angle can be specified as an angle
        about either the node at the link head, or the node at the link tail.
        The default behaviour of `angle_of_link` is to return the angle about
        the link tail, but this method gives the angle about the link head.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> import numpy as np

        >>> grid = HexModelGrid((3, 2), node_layout="hex")
        >>> np.round(grid.angle_of_link[:3] / np.pi * 3.0)
        array([ 0., 2.,  1.])
        >>> np.round(grid.angle_of_link_about_head[:3] / np.pi * 3.0)  # 60 deg segments
        array([ 3.,  5.,  4.])

        LLCATS: LINF MEAS
        """
        angles = np.arctan2(-np.sin(self.angle_of_link), -np.cos(self.angle_of_link))
        return np.mod(angles, 2.0 * np.pi, out=angles)

    def resolve_values_on_links(self, link_values, out=None):
        """Resolve the xy-components of links.

        Resolves values provided defined on links into the x and y directions.
        Returns values_along_x, values_along_y

        LLCATS: LINF
        """
        return gfuncs.resolve_values_on_links(self, link_values, out=out)

    def link_at_node_is_upwind(self, values, out=None):
        """Return a boolean the same shape as :func:`links_at_node` which flags
        links which are upwind of the node as True.

        link_at_node_is_upwind iterates across the grid and identifies the link
        values at each link connected to a node. It then uses the
        link_dirs_at_node data structure to identify links bringing flux into
        the node. It then return a boolean array the same shape as
        links_at_node flagging these links. e.g., for a raster, the returned
        array will be shape (nnodes, 4).

        Parameters
        ----------
        values : str or array
            Name of variable field defined at links, or array of values at
            links.
        out : ndarray, optional
            Buffer to place mapped values into or `None` to create a new array.
            Must be correct shape and boolean dtype.

        Returns
        -------
        ndarray
            Boolean of which links are upwind at nodes.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid

        >>> rmg = RasterModelGrid((3, 4))
        >>> rmg.at_link['grad'] = np.array([-1., -2., -1.,
        ...                                 -2., -3., -4., -5.,
        ...                                 -1., -2., -1.,
        ...                                 -1., -2., -3., -4.,
        ...                                 -1., -2., -1.])
        >>> rmg.link_at_node_is_upwind('grad')
        array([[False, False, False, False],
               [False, False,  True, False],
               [False, False,  True, False],
               [False, False,  True, False],
               [False, False, False,  True],
               [False, False,  True,  True],
               [False, False,  True,  True],
               [False, False,  True,  True],
               [False, False, False,  True],
               [False, False,  True,  True],
               [False, False,  True,  True],
               [False, False,  True,  True]], dtype=bool)

        LLCATS: LINF NINF CONN
        """
        if out is None:
            out = np.empty_like(self.links_at_node, dtype=bool)
        else:
            assert out.shape is self.links_at_node.shape
            assert out.dtype is bool
        if type(values) is str:
            vals = self.at_link[values]
        else:
            assert len(values) == self.number_of_links
            vals = values
        values_at_links = vals[self.links_at_node] * self.link_dirs_at_node
        # this procedure makes incoming links NEGATIVE
        np.less(values_at_links, 0.0, out=out)

        return out

    def link_at_node_is_downwind(self, values, out=None):
        """Return a boolean the same shape as :func:`links_at_node` which flags
        links which are downwind of the node as True.

        link_at_node_is_downwind iterates across the grid and identifies the
        link values at each link connected to a node. It then uses the
        link_dirs_at_node data structure to identify links carrying flux out of
        the node. It then return a boolean array the same shape as
        links_at_node flagging these links. e.g., for a raster, the returned
        array will be shape (nnodes, 4).

        Parameters
        ----------
        values : str or array
            Name of variable field defined at links, or array of values at
            links.
        out : ndarray, optional
            Buffer to place mapped values into or `None` to create a new array.
            Must be correct shape and boolean dtype.

        Returns
        -------
        ndarray
            Boolean of which links are downwind at nodes.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid

        >>> rmg = RasterModelGrid((3, 4))
        >>> rmg.at_link['grad'] = np.array([-1., -2., -1.,
        ...                                 -2., -3., -4., -5.,
        ...                                 -1., -2., -1.,
        ...                                 -1., -2., -3., -4.,
        ...                                 -1., -2., -1.])
        >>> rmg.link_at_node_is_downwind('grad')
        array([[ True,  True, False, False],
               [ True,  True, False, False],
               [ True,  True, False, False],
               [False,  True, False, False],
               [ True,  True, False, False],
               [ True,  True, False, False],
               [ True,  True, False, False],
               [False,  True, False, False],
               [ True, False, False, False],
               [ True, False, False, False],
               [ True, False, False, False],
               [False, False, False, False]], dtype=bool)

        LLCATS: LINF NINF CONN
        """
        if out is None:
            out = np.empty_like(self.links_at_node, dtype=bool)
        else:
            assert out.shape is self.links_at_node.shape
            assert out.dtype is bool
        if type(values) is str:
            vals = self.at_link[values]
        else:
            assert len(values) == self.number_of_links
            vals = values
        values_at_links = vals[self.links_at_node] * self.link_dirs_at_node
        # this procedure makes incoming links NEGATIVE
        np.greater(values_at_links, 0.0, out=out)

        return out

    def upwind_links_at_node(self, values, bad_index=-1):
        """Return an (nnodes, X) shape array of link IDs of which links are
        upwind of each node, according to *values* (field or array).

        X is the maximum upwind links at any node. Nodes with fewer upwind
        links than this have additional slots filled with *bad_index*. Links
        are ordered anticlockwise from east.

        Parameters
        ----------
        values : str or array
            Name of variable field defined at links, or array of values at
            links.
        bad_index : int
            Index to place in array indicating no link.

        Returns
        -------
        ndarray
            Array of upwind link IDs

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid

        >>> rmg = RasterModelGrid((3, 4))
        >>> rmg.at_link['grad'] = np.array([-1., -2., -1.,
        ...                                 -2., -3., -4., -5.,
        ...                                 -1., -2., -1.,
        ...                                 -1., -2., -3., -4.,
        ...                                 -1., -2., -1.])
        >>> rmg.upwind_links_at_node('grad', bad_index=-1)
        array([[-1, -1],
               [ 0, -1],
               [ 1, -1],
               [ 2, -1],
               [ 3, -1],
               [ 7,  4],
               [ 8,  5],
               [ 9,  6],
               [10, -1],
               [14, 11],
               [15, 12],
               [16, 13]])

        LLCATS: LINF NINF CONN
        """
        if type(values) is str:
            vals = self.at_link[values]
        else:
            assert len(values) == self.number_of_links
            vals = values
        values_at_links = vals[self.links_at_node] * self.link_dirs_at_node
        # this procedure makes incoming links NEGATIVE
        unordered_IDs = np.where(values_at_links < 0.0, self.links_at_node, bad_index)
        bad_IDs = unordered_IDs == bad_index
        nnodes = self.number_of_nodes
        flat_sorter = np.argsort(bad_IDs, axis=1) + self.links_at_node.shape[
            1
        ] * np.arange(nnodes).reshape((nnodes, 1))
        big_ordered_array = unordered_IDs.ravel()[flat_sorter].reshape(
            self.links_at_node.shape
        )
        cols_to_cut = int(bad_IDs.sum(axis=1).min())

        if cols_to_cut > 0:
            return big_ordered_array[:, :-cols_to_cut]
        else:
            return big_ordered_array

    def downwind_links_at_node(self, values, bad_index=-1):
        """Return an (nnodes, X) shape array of link IDs of which links are
        downwind of each node, according to *values* (array or field).

        X is the maximum downwind links at any node. Nodes with fewer downwind
        links than this have additional slots filled with *bad_index*. Links
        are ordered anticlockwise from east.

        Parameters
        ----------
        values : str or array
            Name of variable field defined at links, or array of values at
            links.
        bad_index : int
            Index to place in array indicating no link.

        Returns
        -------
        ndarray
            Array of upwind link IDs

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid

        >>> rmg = RasterModelGrid((3, 4))
        >>> rmg.at_link['grad'] = np.array([-1., -2., -1.,
        ...                                 -2., -3., -4., -5.,
        ...                                 -1., -2., -1.,
        ...                                 -1., -2., -3., -4.,
        ...                                 -1., -2., -1.])
        >>> rmg.downwind_links_at_node('grad', bad_index=rmg.BAD_INDEX)
        array([[ 0,  3],
               [ 1,  4],
               [ 2,  5],
               [ 6, -1],
               [ 7, 10],
               [ 8, 11],
               [ 9, 12],
               [13, -1],
               [14, -1],
               [15, -1],
               [16, -1],
               [-1, -1]])

        LLCATS: LINF NINF CONN
        """
        if type(values) is str:
            vals = self.at_link[values]
        else:
            assert len(values) == self.number_of_links
            vals = values
        values_at_links = vals[self.links_at_node] * self.link_dirs_at_node
        # this procedure makes incoming links NEGATIVE
        unordered_IDs = np.where(values_at_links > 0.0, self.links_at_node, bad_index)
        bad_IDs = unordered_IDs == bad_index
        nnodes = self.number_of_nodes
        flat_sorter = np.argsort(bad_IDs, axis=1) + self.links_at_node.shape[
            1
        ] * np.arange(nnodes).reshape((nnodes, 1))
        big_ordered_array = unordered_IDs.ravel()[flat_sorter].reshape(
            self.links_at_node.shape
        )
        cols_to_cut = int(bad_IDs.sum(axis=1).min())

        if cols_to_cut > 0:
            return big_ordered_array[:, :-cols_to_cut]
        else:
            return big_ordered_array

    @property
    @make_return_array_immutable
    def patches_present_at_node(self):
        """A boolean array, False where a patch has a closed node or is
        missing.

        The array is the same shape as :func:`patches_at_node`, and is designed
        to mask it.

        Note that in cases where patches may have more than 3 nodes (e.g.,
        rasters), a patch is considered still present as long as at least 3
        open nodes are present.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((3, 3))
        >>> mg.status_at_node[mg.nodes_at_top_edge] = mg.BC_NODE_IS_CLOSED
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
        >>> mg.patches_present_at_node
        array([[ True, False, False, False],
               [ True,  True, False, False],
               [False,  True, False, False],
               [False, False, False,  True],
               [False, False,  True,  True],
               [False, False,  True, False],
               [False, False, False, False],
               [False, False, False, False],
               [False, False, False, False]], dtype=bool)
        >>> 1 in mg.patches_at_node * mg.patches_present_at_node
        True
        >>> 2 in mg.patches_at_node * mg.patches_present_at_node
        False

        LLCATS: PINF NINF
        """
        try:
            return self._patches_present_mask
        except AttributeError:
            self.patches_at_node
            self._reset_patch_status()
            return self._patches_present_mask

    @property
    @make_return_array_immutable
    def patches_present_at_link(self):
        """A boolean array, False where a patch has a closed node or is
        missing.

        The array is the same shape as :func:`patches_at_link`, and is designed
        to mask it.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((3, 3))
        >>> mg.status_at_node[mg.nodes_at_top_edge] = mg.BC_NODE_IS_CLOSED
        >>> mg.patches_at_link
        array([[-1,  0],
               [-1,  1],
               [ 0, -1],
               [ 1,  0],
               [-1,  1],
               [ 0,  2],
               [ 1,  3],
               [ 2, -1],
               [ 3,  2],
               [-1,  3],
               [ 2, -1],
               [ 3, -1]])
        >>> mg.patches_present_at_link
        array([[False,  True],
               [False,  True],
               [ True, False],
               [ True,  True],
               [False,  True],
               [ True, False],
               [ True, False],
               [False, False],
               [False, False],
               [False, False],
               [False, False],
               [False, False]], dtype=bool)
        >>> 1 in mg.patches_at_link * mg.patches_present_at_link
        True
        >>> 2 in mg.patches_at_link * mg.patches_present_at_link
        False

        LLCATS: PINF LINF
        """
        try:
            return self._patches_present_link_mask
        except AttributeError:
            self.patches_at_node
            self._reset_patch_status()
            return self._patches_present_link_mask

    @property
    @make_return_array_immutable
    def number_of_patches_present_at_node(self):
        """Return the number of patches at a node without a closed node.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((3, 3))
        >>> mg.status_at_node[mg.nodes_at_top_edge] = mg.BC_NODE_IS_CLOSED
        >>> mg.patches_present_at_node
        array([[ True, False, False, False],
               [ True,  True, False, False],
               [False,  True, False, False],
               [False, False, False,  True],
               [False, False,  True,  True],
               [False, False,  True, False],
               [False, False, False, False],
               [False, False, False, False],
               [False, False, False, False]], dtype=bool)
        >>> mg.number_of_patches_present_at_node
        array([1, 2, 1, 1, 2, 1, 0, 0, 0])

        LLCATS: PINF NINF BC
        """
        try:
            return self._number_of_patches_present_at_node
        except AttributeError:
            self.patches_at_node
            self._reset_patch_status()
            return self._number_of_patches_present_at_node

    @property
    @make_return_array_immutable
    def number_of_patches_present_at_link(self):
        """Return the number of patches at a link without a closed node.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((3, 3))
        >>> mg.status_at_node[mg.nodes_at_top_edge] = mg.BC_NODE_IS_CLOSED
        >>> mg.patches_present_at_link
        array([[False,  True],
               [False,  True],
               [ True, False],
               [ True,  True],
               [False,  True],
               [ True, False],
               [ True, False],
               [False, False],
               [False, False],
               [False, False],
               [False, False],
               [False, False]], dtype=bool)
        >>> mg.number_of_patches_present_at_link
        array([1, 1, 1, 2, 1, 1, 1, 0, 0, 0, 0, 0])

        LLCATS: PINF LINF BC
        """
        try:
            return self._number_of_patches_present_at_link
        except AttributeError:
            self.patches_at_node
            self._reset_patch_status()
            return self._number_of_patches_present_at_link

    def _reset_patch_status(self):
        """Creates the array which stores patches_present_at_node.

        Call whenever boundary conditions are updated on the grid.
        """
        from landlab import RasterModelGrid, VoronoiDelaunayGrid

        node_status_at_patch = self.status_at_node[self.nodes_at_patch]
        if isinstance(self, RasterModelGrid):
            max_nodes_at_patch = 4
        elif isinstance(self, VoronoiDelaunayGrid):
            max_nodes_at_patch = 3
        else:
            max_nodes_at_patch = (self.nodes_at_patch > -1).sum(axis=1)
        any_node_at_patch_closed = (node_status_at_patch == self.BC_NODE_IS_CLOSED).sum(
            axis=1
        ) > (max_nodes_at_patch - 3)
        absent_patches = any_node_at_patch_closed[self.patches_at_node]
        bad_patches = np.logical_or(absent_patches, self.patches_at_node == -1)
        self._patches_present_mask = np.logical_not(bad_patches)
        self._number_of_patches_present_at_node = np.sum(
            self._patches_present_mask, axis=1
        )
        absent_patches = any_node_at_patch_closed[self.patches_at_link]
        bad_patches = np.logical_or(absent_patches, self.patches_at_link == -1)
        self._patches_present_link_mask = np.logical_not(bad_patches)
        self._number_of_patches_present_at_link = np.sum(
            self._patches_present_link_mask, axis=1
        )

    def calc_hillshade_at_node(
        self,
        alt=45.0,
        az=315.0,
        slp=None,
        asp=None,
        unit="degrees",
        elevs="topographic__elevation",
    ):
        """Get array of hillshade.

        .. codeauthor:: Katy Barnhart <katherine.barnhart@colorado.edu>

        Parameters
        ----------
        alt : float
            Sun altitude (from horizon) - defaults to 45 degrees
        az : float
            Sun azimuth (CW from north) - defaults to 315 degrees
        slp : float
            slope of cells at surface - optional
        asp : float
            aspect of cells at surface (from north) - optional (with slp)
        unit : string
            'degrees' (default) or 'radians' - only needed if slp and asp
                                                are not provided

        If slp and asp are both not specified, 'elevs' must be provided as
        a grid field name (defaults to 'topographic__elevation') or an
        nnodes-long array of elevation values. In this case, the method will
        calculate local slopes and aspects internally as part of the hillshade
        production.

        Returns
        -------
        ndarray of float
            Hillshade at each cell.

        Notes
        -----
        code taken from GeospatialPython.com example from December 14th, 2014
        DEJH found what looked like minor sign problems, and adjusted to follow
        the `ArcGIS algorithm
        <http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#/How_Hillshade_works/009z000000z2000000/>`.

        Remember when plotting that bright areas have high values. cmap='Greys'
        will give an apparently inverted color scheme. *cmap='gray'* has white
        associated with the high values, so is recommended for plotting.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid

        >>> mg = RasterModelGrid((5, 5), xy_spacing=1.)
        >>> z = mg.x_of_node * np.tan(60. * np.pi / 180.)
        >>> mg.calc_hillshade_at_node(elevs=z, alt=30., az=210.)
        array([ 0.625,  0.625,  0.625,  0.625,  0.625,  0.625,  0.625,  0.625,
                0.625,  0.625,  0.625,  0.625,  0.625,  0.625,  0.625,  0.625,
                0.625,  0.625,  0.625,  0.625,  0.625,  0.625,  0.625,  0.625,
                0.625])

        LLCATS: NINF SURF
        """
        if slp is not None and asp is not None:
            if unit == "degrees":
                (alt, az, slp, asp) = (
                    np.radians(alt),
                    np.radians(az),
                    np.radians(slp),
                    np.radians(asp),
                )
            elif unit == "radians":
                if alt > np.pi / 2.0 or az > 2.0 * np.pi:
                    print(
                        "Assuming your solar properties are in degrees, "
                        "but your slopes and aspects are in radians..."
                    )
                    (alt, az) = (np.radians(alt), np.radians(az))
                    # ...because it would be super easy to specify radians,
                    # but leave the default params alone...
            else:
                raise TypeError("unit must be 'degrees' or 'radians'")
        elif slp is None and asp is None:
            if unit == "degrees":
                (alt, az) = (np.radians(alt), np.radians(az))
            elif unit == "radians":
                pass
            else:
                raise TypeError("unit must be 'degrees' or 'radians'")
            slp, slp_comps = self.calc_slope_at_node(elevs, return_components=True)

            asp = self.calc_aspect_at_node(
                slope_component_tuple=slp_comps, unit="radians"
            )
        else:
            raise TypeError("Either both slp and asp must be set, or neither!")

        shaded = np.sin(alt) * np.cos(slp) + np.cos(alt) * np.sin(slp) * np.cos(
            az - asp
        )

        return shaded.clip(0.0)

    @property
    @lru_cache()
    @make_return_array_immutable
    def cell_area_at_node(self):
        """Cell areas in a nnodes-long array.

        Zeros are entered at all perimeter nodes, which lack cells.

        Returns
        -------
        ndarray
            Cell areas as an n_nodes-long array.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5), xy_spacing=(3, 4))
        >>> grid.status_at_node[7] = grid.BC_NODE_IS_CLOSED
        >>> grid.cell_area_at_node
        array([  0.,   0.,   0.,   0.,   0.,
                 0.,  12.,  12.,  12.,   0.,
                 0.,  12.,  12.,  12.,   0.,
                 0.,   0.,   0.,   0.,   0.])

        LLCATS: CINF NINF CONN
        """
        cell_area_at_node = np.zeros(self.number_of_nodes, dtype=float)
        cell_area_at_node[self.node_at_cell] = self.area_of_cell
        return cell_area_at_node

    def reset_status_at_node(self):
        attrs = [
            "_active_link_dirs_at_node",
            "_status_at_link",
            "_active_links",
            "_fixed_links",
            "_activelink_fromnode",
            "_activelink_tonode",
            "_active_faces",
            "_core_nodes",
            "_core_cells",
            "_fixed_links",
            "_active_adjacent_nodes_at_node",
            "_fixed_value_boundary_nodes",
            "_node_at_core_cell",
            "_link_status_at_node",
        ]

        for attr in attrs:
            try:
                del self.__dict__[attr]
            except KeyError:
                pass
        try:
            self.bc_set_code += 1
        except AttributeError:
            self.bc_set_code = 0
        try:
            del self.__dict__["__node_active_inlink_matrix"]
        except KeyError:
            pass
        try:
            del self.__dict__["__node_active_outlink_matrix"]
        except KeyError:
            pass

    def set_nodata_nodes_to_closed(self, node_data, nodata_value):
        """Make no-data nodes closed boundaries.

        Sets node status to `BC_NODE_IS_CLOSED` for all nodes whose value
        of *node_data* is equal to the *nodata_value*.

        Any links connected to `BC_NODE_IS_CLOSED` nodes are automatically
        set to `LinkStatus.INACTIVE` boundary.

        Parameters
        ----------
        node_data : ndarray
            Data values.
        nodata_value : float
            Value that indicates an invalid value.

        Examples
        --------

        The following example uses the following grid::

          *--I--->o------>o------>o
          ^       ^       ^       ^
          I       I       |       |
          |       |       |       |
          *--I--->*--I--->o------>o
          ^       ^       ^       ^
          I       I       I       I
          |       |       |       |
          *--I--->*--I--->*--I--->*

        .. note::

          Links set to `LinkStatus.ACTIVE` are not shown in this diagram.

        ``*`` indicates the nodes that are set to `NodeStatus.CLOSED`

        ``o`` indicates the nodes that are set to `NodeStatus.CORE`

        ``I`` indicates the links that are set to `LinkStatus.INACTIVE`

        >>> import numpy as np
        >>> import landlab as ll
        >>> mg = ll.RasterModelGrid((3, 4))
        >>> mg.status_at_node
        array([1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1], dtype=uint8)
        >>> h = np.array([-9999, -9999, -9999, -9999, -9999, -9999, 12345.,
        ...     0., -9999, 0., 0., 0.])
        >>> mg.set_nodata_nodes_to_closed(h, -9999)
        >>> mg.status_at_node
        array([4, 4, 4, 4, 4, 4, 0, 1, 4, 1, 1, 1], dtype=uint8)

        LLCATS: BC NINF
        """
        # Find locations where value equals the NODATA code and set these nodes
        # as inactive boundaries.
        nodata_locations = np.nonzero(node_data == nodata_value)
        self.status_at_node[nodata_locations] = self.BC_NODE_IS_CLOSED

    def set_nodata_nodes_to_fixed_gradient(self, node_data, nodata_value):
        """Make no-data nodes fixed gradient boundaries.

        Set node status to `BC_NODE_IS_FIXED_VALUE` for all nodes
        whose value of *node_data* is equal to *nodata_value*.

        Any links between `BC_NODE_IS_FIXED_GRADIENT` nodes and
        `BC_NODE_IS_CORE` are automatically set to `LinkStatus.FIXED` boundary
        status.

        Parameters
        ----------
        node_data : ndarray
            Data values.
        nodata_value : float
            Value that indicates an invalid value.

        Examples
        --------

        The following examples use this grid::

          *--I--->*--I--->*--I--->*--I--->*--I--->*--I--->*--I--->*--I--->*
          ^       ^       ^       ^       ^       ^       ^       ^       ^
          I       I       I       X       X       X       X       X       I
          |       |       |       |       |       |       |       |       |
          *--I--->*--I--->*--X--->o       o       o       o       o--X--->*
          ^       ^       ^       ^       ^       ^       ^       ^       ^
          I       I       I       |       |       |       |       |       I
          |       |       |       |       |       |       |       |       |
          *--I--->*--I--->*--X--->o       o       o       o       o--X--->*
          ^       ^       ^       ^       ^       ^       ^       ^       ^
          I       I       I       X       X       X       X       X       I
          |       |       |       |       |       |       |       |       |
          *--I--->*--I--->*--I--->*--I--->*--I--->*--I--->*--I--->*--I--->*

        .. note::

            Links set to `LinkStatus.ACTIVE` are not shown in this diagram.

        ``X`` indicates the links that are set to `LinkStatus.FIXED`

        ``I`` indicates the links that are set to `LinkStatus.INACTIVE`

        ``o`` indicates the nodes that are set to `NodeStatus.CORE`

        ``*`` indicates the nodes that are set to
              `BC_NODE_IS_FIXED_GRADIENT`

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid((4, 9))
        >>> rmg.status_at_node # doctest: +NORMALIZE_WHITESPACE
        array([1, 1, 1, 1, 1, 1, 1, 1, 1,
               1, 0, 0, 0, 0, 0, 0, 0, 1,
               1, 0, 0, 0, 0, 0, 0, 0, 1,
               1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=uint8)

        >>> z = rmg.zeros(at='node')
        >>> z = np.array([
        ...     -99., -99., -99., -99., -99., -99., -99., -99., -99.,
        ...     -99., -99., -99.,   0.,   0.,   0.,   0.,   0., -99.,
        ...     -99., -99., -99.,   0.,   0.,   0.,   0.,   0., -99.,
        ...     -99., -99., -99., -99., -99., -99., -99., -99., -99.])

        >>> rmg.set_nodata_nodes_to_fixed_gradient(z, -99)
        >>> rmg.status_at_node # doctest: +NORMALIZE_WHITESPACE
        array([2, 2, 2, 2, 2, 2, 2, 2, 2,
               2, 2, 2, 0, 0, 0, 0, 0, 2,
               2, 2, 2, 0, 0, 0, 0, 0, 2,
               2, 2, 2, 2, 2, 2, 2, 2, 2], dtype=uint8)

        >>> rmg.status_at_link # doctest: +NORMALIZE_WHITESPACE
        array([4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 4,
               4, 4, 2, 0, 0, 0, 0, 2, 4, 4, 4, 0, 0, 0, 0, 0, 4,
               4, 4, 2, 0, 0, 0, 0, 2, 4, 4, 4, 2, 2, 2, 2, 2, 4,
               4, 4, 4, 4, 4, 4, 4, 4], dtype=uint8)

        LLCATS: BC NINF
        """
        # Find locations where value equals the NODATA code and set these nodes
        # as inactive boundaries.
        nodata_locations = np.nonzero(node_data == nodata_value)
        self.status_at_node[nodata_locations] = NodeStatus.FIXED_GRADIENT

    @property
    def unit_vector_sum_xcomponent_at_node(self):
        """Get array of x-component of unit vector sums at each node.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 3))
        >>> len(grid.unit_vector_sum_xcomponent_at_node) == grid.number_of_nodes
        True
        >>> grid.unit_vector_sum_xcomponent_at_node
        array([ 1.,  2.,  1.,  1.,  2.,  1.,  1.,  2.,  1.])

        LLCATS: NINF MEAS
        """
        return self.unit_vector_at_node[:, 0]

    @property
    def unit_vector_sum_ycomponent_at_node(self):
        """Get array of y-component of unit vector sums at each node.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 3))
        >>> len(grid.unit_vector_sum_ycomponent_at_node) == grid.number_of_nodes
        True
        >>> grid.unit_vector_sum_ycomponent_at_node
        array([ 1.,  1.,  1.,  2.,  2.,  2.,  1.,  1.,  1.])

        LLCATS: NINF MEAS
        """
        return self.unit_vector_at_node[:, 1]

    def map_link_vector_to_nodes(self, q):
        r"""Map data defined on links to nodes.

        Given a variable defined on links, breaks it into x and y components
        and assigns values to nodes by averaging each node's attached links.

        Parameters
        ----------
        q : ndarray of floats (1D, length = number of links in grid)
            Variable defined on links

        Returns
        -------
        ndarray, ndarray
            x and y components of variable mapped to nodes (1D,
            length = number of nodes)

        See Also
        --------
        _create_link_unit_vectors : sets up unit vectors at links and unit-vector
                                  sums at nodes

        Notes
        -----
        THIS ALGORITHM IS NOT CORRECT AND NEEDS TO BE CHANGED!

        The concept here is that q contains a vector variable that is defined
        at each link. The magnitude is given by the value of q, and the
        direction is given by the orientation of the link, as described by
        its unit vector.

        To map the link-vector values to the nodes, we break the values into
        x- and y-components according to each link's unit vector. The
        x-component of q at a node is a weighted sum of the x-components of the
        links that are attached to that node. A good way to appreciate this
        is by example. Consider a 3x4 raster grid::

            8--14---9--15--10--16--11
            |       |       |       |
            4       5       6       7
            |       |       |       |
            4--11---5---12--6---13--7
            |       |       |       |
            0       1       2       3
            |       |       |       |
            0---8---1---9---2--10---3

        Imagine that for each node, we were to add up the unit vector
        components for each connected link; in other words, add up all the x
        components of the unit vectors associated with each link, and add up
        all the y components. Here's what that would look like for the above
        grid ("vsx" and "vsy" stand for "vector sum x" and "vector sum y"):

        *  Corner nodes (0, 3, 8, 11): vsx = 1, vsy = 1
        *  Bottom and top nodes (1-2, 9-10): vsx = 2, vsy = 1
        *  Left and right nodes (4, 7): vsx = 1, vsy = 2
        *  All others: vsx = 2, vsy = 2

        The process of creating unit-vector sums at nodes is handled by
        ModelGrid._create_link_unit_vectors() (and, for raster grids, by the
        overriding method RasterModelGrid._create_link_unit_vectors()). The node
        unit-vector sums are then stored in self.node_unit_vector_sum_x and
        self.node_unit_vector_sum_y.

        How would you use this? Suppose you have a vector variable q defined at
        links. What's the average at the nodes? We'll define the average as
        follows.  The terminology here is: :math:`q = (u,v)` represents the
        vector quantity defined at links, :math:`Q = (U,V)` represents its
        definition at nodes, :math:`(m,n)` represents the unit vector
        components at a link, and :math:`(S_x,S_y)` represents the unit-vector
        sum at a given node.

        .. math::

            U_i = \sum_{j=1}^{L_i} q_j m_j / S_{xi}
            V_i = \sum_{j=1}^{L_i} q_j n_j / S_{yi}

        Suppose that the vector q is uniform and equal to one.
        Then, at node 0 in the above grid, this works out to::

            U_0 = (q_0 m_0) / 1 + (q_8 m_8) / 1 = (1 0)/ 1 + (1 1)/1 = 1
            V_0 = (q_0 n_0) / 1 + (q_8 n_8) / 1 = (1 1) / 1 + (1 0) / 1 = 1

        At node 1, in the bottom row but not a corner, we add up the values
        of **q** associated with THREE links. The x-vector sum of these links
        is 2 because there are two horizontal links, each with an x- unit
        vector value of unity.  The y-vector sum is 1 because only one of the
        three (link #1) has a non-zero y component (equal to one). Here is
        how the numbers work out::

            U_1 = (q_1 m_1) / 2 + (q_8 m_8) / 2 + (q_9 m_9) / 2
                = (1 0) / 2 + (1 1) / 2 + (1 1) / 2 = 1
            V_1 = (q_1 n_1) / 1 + (q_8 n_8) / 1 + (q_9 n_9) / 1
                = (1 1) / 1 + (1 0) / 1 + (1 0) / 1 = 1

        At node 5, in the interior, there are four connected links (two
        in-links and two out-links; two horizontal and two vertical). So, we
        add up the q values associated with all four::

            U_5 = (q_1 m_1) / 2 + (q_5 m_5) / 2 + (q_11 m_11) / 2 + (q_12 m_12) / 2
                = (1 0) / 2 + (1 0) / 2 + (1 1) / 2 + (1 1) / 2 = 1

            V_5 = (q_1 n_1) / 2 + (q_5 n_5) / 2 + (q_11 n_11) / 2 + (q_12 n_12) / 2
                = (1 1) / 2 + (1 1) / 2 + (1 0) / 2 + (1 0) / 2 = 1

        To do this calculation efficiently, we use the following algorithm::

            FOR each row in _node_inlink_matrix (representing one inlink @ each
            node)
                Multiply the link's q value by its unit x component ...
                ... divide by node's unit vector sum in x ...
                ... and add it to the node's total q_x
                Multiply the link's q value by its unit y component ...
                ... divide by node's unit vector sum in y ...
                ... and add it to the node's total q_y

        Examples
        --------

        **Example 1**

        q[:] = 1. Vector magnitude is :math:`\sqrt{2}`, direction is
        :math:`(1,1)`.

        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4), xy_spacing=(2., 2.))
        >>> grid.unit_vector_at_node
        array([[ 1.,  1.],
               [ 2.,  1.],
               [ 2.,  1.],
               [ 1.,  1.],
               [ 1.,  2.],
               [ 2.,  2.],
               [ 2.,  2.],
               [ 1.,  2.],
               [ 1.,  1.],
               [ 2.,  1.],
               [ 2.,  1.],
               [ 1.,  1.]])
        >>> q = grid.ones(at='link')
        >>> grid.map_link_vector_to_nodes(q)
        array([[ 1.,  1.],
               [ 1.,  1.],
               [ 1.,  1.],
               [ 1.,  1.],
               [ 1.,  1.],
               [ 1.,  1.],
               [ 1.,  1.],
               [ 1.,  1.],
               [ 1.,  1.],
               [ 1.,  1.],
               [ 1.,  1.],
               [ 1.,  1.]])

        **Example 2**

        Vector magnitude is 5, angle is 30 degrees from horizontal,
        forming a 3-4-5 triangle.

        >>> import numpy as np
        >>> q = np.array([4., 4., 4., 3., 3., 3., 3.,
        ...               4., 4., 4., 3., 3., 3., 3.,
        ...               4., 4., 4])
        >>> grid.map_link_vector_to_nodes(q)
        array([[ 4.,  3.],
               [ 4.,  3.],
               [ 4.,  3.],
               [ 4.,  3.],
               [ 4.,  3.],
               [ 4.,  3.],
               [ 4.,  3.],
               [ 4.,  3.],
               [ 4.,  3.],
               [ 4.,  3.],
               [ 4.,  3.],
               [ 4.,  3.]])

        ..todo::

            Fix and finish example 3 below.

        Example 3: Hexagonal grid with vector as above. Here, q is
        pre-calculated to have the right values to represent a uniform
        vector with magnitude 5 and orientation 30 degrees counter-clockwise
        from horizontal.

        LLCATS: NINF LINF CONN MAP
        """
        # Break the link-based vector input variable, q, into x- and
        # y-components.
        # Notes:
        #   1) We make the arrays 1 element longer than the number of links,
        #       so that references to -1 in the node-link matrices will refer
        #       to the last element of these two arrays, which will contain
        #       zeros. (Same trick as in the flux divergence functions)
        #   2) This requires memory allocation. Because this function might be
        #       called repeatedly, it would be good to find a way to
        #       pre-allocate to improve speed.
        qx = q * self.unit_vector_at_link[:, 0]
        qy = q * self.unit_vector_at_link[:, 1]

        active_links_at_node = self.link_dirs_at_node != 0
        unit_vec_at_node = np.empty((self.number_of_nodes, 2), dtype=float)
        unit_vec_at_node[:, 0] = (qx[self.links_at_node] * active_links_at_node).sum(
            axis=1
        )
        unit_vec_at_node[:, 1] = (qy[self.links_at_node] * active_links_at_node).sum(
            axis=1
        )

        return np.divide(
            unit_vec_at_node, self.unit_vector_at_node, out=unit_vec_at_node
        )

    def node_is_boundary(self, ids, boundary_flag=None):
        """Check if nodes are boundary nodes.

        Check if nodes at given *ids* are boundary nodes. Use the
        *boundary_flag* to specify a particular boundary type status flag.

        Parameters
        ----------
        ids : ndarray
            Node IDs to check.
        boundary_flag : int, optional
            A boundary type to check for.

        Returns
        -------
        ndarray
            Array of booleans indicating if nodes are boundary nodes.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((4, 5))
        >>> mg.node_is_boundary([0, 6])
        array([ True, False], dtype=bool)
        >>> mg.node_is_boundary([0, 6], boundary_flag=mg.BC_NODE_IS_CLOSED)
        array([False, False], dtype=bool)

        LLCATS: NINF BC
        """
        if boundary_flag is None:
            return ~(self._node_status[ids] == NodeStatus.CORE)
        else:
            return self._node_status[ids] == boundary_flag

    def calc_distances_of_nodes_to_point(
        self, coord, get_az=None, node_subset=None, out_distance=None, out_azimuth=None
    ):
        """Get distances for nodes to a given point.

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

        .. note::

            Angles are returned in radians but measured clockwise from
            north.

        Parameters
        ----------
        coord : tuple of float
            Coodinates of point as (x, y).
        get_az: {None, 'angles', 'displacements'}, optional
            Optionally calculate azimuths as either angles or displacements.
            The calculated values will be returned along with the distances
            as the second item of a tuple.
        node_subset : array_like, optional
            Calculate distances on a subset of grid nodes. The default is to
            calculate distances from the provided points to all nodes.
        out_distance : array_like, optional
            If provided, put the calculated distances here. Otherwise,
            create a new array.
        out_azimuth : array_like, optional
            If provided, put the calculated distances here. Otherwise,
            create a new array.

        Returns
        -------
        ndarray or tuple of ndarray
            If *get_az* is ``None`` return the array of distances. Otherwise,
            return a tuple of distances and azimuths.

        Notes
        -----
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

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))

        Calculate distances from point at (2., 1.) to a subset of nodes on
        the grid.

        >>> grid.calc_distances_of_nodes_to_point((2, 1),
        ...     node_subset=(2, 6, 7, 8, 12))
        array([ 1.,  1.,  0.,  1.,  1.])

        Calculate distances from a point to all nodes on the grid.

        >>> dist = grid.calc_distances_of_nodes_to_point((2, 1))
        >>> dist.shape == (grid.number_of_nodes, )
        True
        >>> dist.take((2, 6, 7, 8, 12))
        array([ 1.,  1.,  0.,  1.,  1.])

        Put the distances into a buffer.

        >>> out = np.empty(grid.number_of_nodes, dtype=float)
        >>> dist = grid.calc_distances_of_nodes_to_point((2, 1),
        ...     out_distance=out)
        >>> out is dist
        True
        >>> out.take((2, 6, 7, 8, 12))
        array([ 1.,  1.,  0.,  1.,  1.])

        Calculate azimuths along with distances. The azimuths are calculated
        in radians but measured clockwise from north.

        >>> (_, azim) = grid.calc_distances_of_nodes_to_point((2, 1),
        ...     get_az='angles')
        >>> azim.take((2, 6, 7, 8, 12)) * 180. / np.pi
        array([ 180.,  270.,    0.,   90.,    0.])
        >>> (_, azim) = grid.calc_distances_of_nodes_to_point((2, 1),
        ...     get_az='angles', node_subset=(1, 3, 11, 13))
        >>> azim * 180. / np.pi
        array([ 225.,  135.,  315.,   45.])

        When calculating displacements, the first row contains displacements
        in x and the second displacements in y.

        >>> (_, azim) = grid.calc_distances_of_nodes_to_point((2, 1),
        ...     get_az='displacements', node_subset=(2, 6, 7, 8, 12))
        >>> azim
        array([[ 0., -1.,  0.,  1.,  0.],
               [-1.,  0.,  0.,  0.,  1.]])

        LLCATS: NINF MEAS
        """
        if len(coord) != 2:
            raise ValueError("coordinate must iterable of length 2")

        if get_az not in (None, "displacements", "angles"):
            raise ValueError("get_az not understood")

        if node_subset is not None and np.any(np.isnan(node_subset)):
            node_subset = None

        if node_subset is not None:
            if not isinstance(node_subset, np.ndarray):
                node_subset = np.array(node_subset)
            node_subset = node_subset.reshape((-1,))
            len_subset = node_subset.size
        else:
            len_subset = self.number_of_nodes

        if out_distance is None:
            out_distance = np.empty(len_subset, dtype=np.float)
        if out_distance.size != len_subset:
            raise ValueError("output array size mismatch for distances")

        if get_az is not None:
            if get_az == "displacements":
                az_shape = (2, len_subset)
            else:
                az_shape = (len_subset,)
            if out_azimuth is None:
                out_azimuth = np.empty(az_shape, dtype=np.float)
            if out_azimuth.shape != az_shape:
                raise ValueError("output array mismatch for azimuths")

        azimuths_as_displacements = np.empty((2, self.number_of_nodes))
        dummy_nodes_1 = np.empty(self.number_of_nodes)
        dummy_nodes_2 = np.empty(self.number_of_nodes)
        dummy_nodes_3 = np.empty(self.number_of_nodes)

        if node_subset is None:
            azimuths_as_displacements[0] = self.node_x - coord[0]
            azimuths_as_displacements[1] = self.node_y - coord[1]
        else:
            azimuths_as_displacements[0, :len_subset] = (
                self.node_x[node_subset] - coord[0]
            )
            azimuths_as_displacements[1, :len_subset] = (
                self.node_y[node_subset] - coord[1]
            )

        np.square(
            azimuths_as_displacements[0, :len_subset], out=dummy_nodes_1[:len_subset]
        )
        np.square(
            azimuths_as_displacements[1, :len_subset], out=dummy_nodes_2[:len_subset]
        )
        np.add(
            dummy_nodes_1[:len_subset],
            dummy_nodes_2[:len_subset],
            out=dummy_nodes_3[:len_subset],
        )
        np.sqrt(dummy_nodes_3[:len_subset], out=out_distance)

        if get_az:
            if get_az == "displacements":
                out_azimuth[:] = azimuths_as_displacements[:, :len_subset]
            elif get_az == "angles":
                np.arctan2(
                    azimuths_as_displacements[0, :len_subset],
                    azimuths_as_displacements[1, :len_subset],
                    out=out_azimuth[:len_subset],
                )

                less_than_zero = np.empty(self.number_of_nodes, dtype=bool)
                np.less(out_azimuth, 0.0, out=less_than_zero[:len_subset])
                out_azimuth[less_than_zero[:len_subset]] += 2.0 * np.pi

            return out_distance, out_azimuth
        else:
            return out_distance

    @property
    def all_node_distances_map(self):
        """Get distances from every node to every other node.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4))
        >>> distances = grid.all_node_distances_map

        The shape of the array is ``number_of_nodes`` by ``number_of_nodes``
        and distance from a node to itself is zero.

        >>> distances.shape == (grid.number_of_nodes, grid.number_of_nodes)
        True
        >>> distances.diagonal()
        array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])

        The distances from the first node to all nodes in its row and all the
        nodes in its column.

        >>> distances[0, :4]
        array([ 0.,  1.,  2.,  3.])
        >>> distances[0, ::4]
        array([ 0.,  1.,  2.])

        LLCATS: NINF MEAS
        """
        if self._all_node_distances_map is None:
            self._create_all_node_distances_azimuths_maps()
        return self._all_node_distances_map

    @property
    def all_node_azimuths_map(self):
        """Get azimuths from every node to every other node.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 4))
        >>> angles = grid.all_node_azimuths_map

        The shape of the array is ``number_of_nodes`` by ``number_of_nodes``
        and azimuth from a node to itself is zero.

        >>> angles.shape == (grid.number_of_nodes, grid.number_of_nodes)
        True
        >>> angles.diagonal()
        array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])

        Angles are measured in radians and increase clockwise starting at
        north.

        >>> angles *= 180. / np.pi
        >>> angles[0, :4]
        array([  0.,  90.,  90.,  90.])
        >>> angles[0, ::4]
        array([ 0.,  0.,  0.])
        >>> angles[0, ::5]
        array([  0.,  45.,  45.])

        LLCATS: NINF MEAS
        """
        if self._all_node_azimuths_map is None:
            self._create_all_node_distances_azimuths_maps()
        return self._all_node_azimuths_map

    def _create_all_node_distances_azimuths_maps(self):
        """Build distance-azimuth maps.

        This function creates and stores in the grid field two ``nnodes`` by
        ``nnodes`` arrays that map the distances and azimuths of all nodes
        in the grid to all nodes in the grid.

        This is useful if your module needs to make repeated lookups of
        distances between the same nodes, but does potentially use up a lot
        of memory so should be used with caution.

        The map is symmetrical, so it does not matter whether rows are
        "from" or "to".

        The arrays are called:
        - ``self.all_node_distances_map``
        - ``self.all_node_azimuths_map``

        Returns
        -------
        tuple of ndarrays
            Tuple of (distances, azimuths)
        """
        self._all_node_distances_map = np.empty(
            (self.number_of_nodes, self.number_of_nodes)
        )
        self._all_node_azimuths_map = np.empty(
            (self.number_of_nodes, self.number_of_nodes)
        )

        node_coords = np.empty((self.number_of_nodes, 2))
        node_coords[:, 0] = self.node_x
        node_coords[:, 1] = self.node_y

        for i in range(self.number_of_nodes):
            (
                self._all_node_distances_map[i, :],
                self._all_node_azimuths_map[i, :],
            ) = self.calc_distances_of_nodes_to_point(
                (node_coords[i, 0], node_coords[i, 1]), get_az="angles"
            )

        assert np.all(self._all_node_distances_map >= 0.0)

        return self._all_node_distances_map, self._all_node_azimuths_map

    # def node_has_boundary_neighbor(self, ids):
    def node_has_boundary_neighbor(self):
        """Check if ModelGrid nodes have neighbors that are boundary nodes.

        Checks to see if one of the eight neighbor nodes of node(s) with
        *id* has a boundary node.  Returns True if a node has a boundary node,
        False if all neighbors are interior.

        Parameters
        ----------
        ids : int, or iterable of int
            ID of node to test.

        Returns
        -------
        boolean
            ``True`` if node has a neighbor with a boundary ID,
            ``False`` otherwise.

        Examples
        --------

        ::

                0,  1,  2,  3,
              4,  5,  6,  7,  8,
            9, 10,  11, 12, 13, 14,
              15, 16, 17, 18, 19,
                20, 21, 22, 23

        >>> from landlab import HexModelGrid
        >>> grid = HexModelGrid((5, 4))
        >>> grid.node_has_boundary_neighbor()
        array([ True,  True,  True,  True,  True,  True,  True,  True,  True,
                True,  True, False, False,  True,  True,  True,  True,  True,
                True,  True,  True,  True,  True,  True], dtype=bool)

        >>> grid.node_has_boundary_neighbor()[6]
        True
        >>> grid.node_has_boundary_neighbor()[12]
        False
        >>> grid.node_has_boundary_neighbor()[((12, 0),)]
        array([False,  True], dtype=bool)

        LLCATS: NINF CONN BC
        """
        status_of_neighbor = self._node_status[self.adjacent_nodes_at_node]
        neighbor_not_core = status_of_neighbor != NodeStatus.CORE
        bad_neighbor = self.adjacent_nodes_at_node == self.BAD_INDEX
        neighbor_not_core[bad_neighbor] = False
        return np.any(neighbor_not_core, axis=1)

        # node_has_boundary_neighbor = np.any(neighbor_not_core, axis=1)
        # return node_has_boundary_neighbor[ids]


add_module_functions_to_class(ModelGrid, "mappers.py", pattern="map_*")
# add_module_functions_to_class(ModelGrid, 'gradients.py',
#                               pattern='calculate_*')
add_module_functions_to_class(ModelGrid, "gradients.py", pattern="calc_*")
add_module_functions_to_class(ModelGrid, "divergence.py", pattern="calc_*")
