#! /usr/env/python
"""
Python implementation of ModelGrid, a base class used to create and manage
grids for 2D numerical models.

Data Fields in ModelGrid
------------------------
:class:`~.ModelGrid` inherits from the :class:`~.ModelDataFields` class. This
provides `~.ModelGrid`, and its subclasses, with the ability to, optionally,
store data values that are associated with the different types grid elements
(nodes, cells, etc.). In particular, as part of ``ModelGrid.__init__()``,
data field *groups* are added to the `ModelGrid` that provide containers to
put data fields into. There is one group for each of the eight grid elements
(node, cell, link, face, core_node, core_cell, active_link, and active_face).

To access these groups, use the same methods as accessing groups with
`~.ModelDataFields`. ``ModelGrid.__init__()`` adds the following attributes to
itself that provide access to the values groups:

.. autosummary::
    :toctree: generated/
    :nosignatures:

    ~landlab.grid.base.ModelGrid.at_node
    ~landlab.grid.base.ModelGrid.at_cell
    ~landlab.grid.base.ModelGrid.at_link
    ~landlab.grid.base.ModelGrid.at_face
    ~landlab.grid.base.ModelGrid.at_core_node
    ~landlab.grid.base.ModelGrid.at_core_cell
    ~landlab.grid.base.ModelGrid.at_active_link
    ~landlab.grid.base.ModelGrid.at_active_face

Each of these attributes returns a ``dict``-like object whose keys are value
names as strings and values are numpy arrays that gives quantities at
grid elements.


Create Field Arrays
+++++++++++++++++++
:class:`~.ModelGrid` inherits several useful methods for creating new data
fields and adding new data fields to a ModelGrid instance. Methods to add or
create a new data array follow the ``numpy`` syntax for creating arrays. The
folowing methods create and, optionally, initialize new arrays. These arrays
are of the correct size but a new field will not be added to the field:

.. autosummary::
    :toctree: generated/
    :nosignatures:

    ~landlab.field.grouped.ModelDataFields.empty
    ~landlab.field.grouped.ModelDataFields.ones
    ~landlab.field.grouped.ModelDataFields.zeros

Add Fields to a ModelGrid
+++++++++++++++++++++++++
Unlike with the equivalent numpy functions, these do not take a size argument
as the size of the returned arrays is determined from the size of the
ModelGrid. However, the keyword arguments are the same as those of the numpy
equivalents.

The following methods will create a new array and add a reference to that
array to the ModelGrid:

.. autosummary::
    :toctree: generated/
    :nosignatures:

    ~landlab.field.grouped.ModelDataFields.add_empty
    ~landlab.field.grouped.ModelDataFields.add_ones
    ~landlab.field.grouped.ModelDataFields.add_zeros
    ~landlab.field.grouped.ModelDataFields.add_field

These methods operate in the same way as the previous set except that, in
addition to creating a new array, the newly-created array is added to the
ModelGrid. The calling signature is the same but with the addition of an
argument that gives the name of the new field as a string. The additional
method, :meth:`~.ModelDataFields.add_field`, adds a previously allocation
array to the ModelGrid. If the array is of the incorrect size it will raise
``ValueError``.

Query Fields
++++++++++++
Use the following methods/attributes get information about the stored data
fields:

.. autosummary::
    :toctree: generated/
    :nosignatures:

    ~landlab.field.grouped.ModelDataFields.size
    ~landlab.field.grouped.ModelDataFields.keys
    ~landlab.field.grouped.ModelDataFields.has_group
    ~landlab.field.grouped.ModelDataFields.has_field
    ~landlab.field.grouped.ModelDataFields.groups

    i.e., call, e.g. mg.has_field('node', 'my_field_name')

Notes
-----
It is important that when creating a new grid class that inherits from
``ModelGrid``, to call ``ModelGrid.__init__()`` in the new grid's
``__init__()``. For example, the new class's __init__ should contain the
following code,

.. code-block:: python

    class NewGrid(ModelGrid):
        def __init__(self, *args, **kwds):
            ModelGrid.__init__(self, **kwds)
            # Code that initializes the NewGrid

Without this, the new grid class will not have the ``at_*`` attributes.

Examples
--------
Although the following examples use a :class:`~.RasterModelGrid`, they apply
equally to any grid that inherits from :class:`~.ModelGrid`.  The new grid
comes with a set of pre-defined value groups. One group for each grid element.
Use the groups attribute to see the group names.

>>> from landlab import RasterModelGrid
>>> grid = RasterModelGrid(3, 3)
>>> groups = list(grid.groups)
>>> groups.sort()
>>> groups # doctest: +NORMALIZE_WHITESPACE
['active_face', 'active_link', 'cell', 'core_cell', 'core_node', 'face',
 'link', 'node']

Create Field Arrays
+++++++++++++++++++
If you just want to create an array but not add it to the grid, you can use
the :meth:`~.ModelGrid.ones` method.

>>> grid.ones(centering='node')
array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])
>>> list(grid.at_node.keys()) # Nothing has been added to the grid
[]

Add Field Arrays
++++++++++++++++
Use the ``add_*`` methods to add value arrays attached to grid elements. Each
of these methods accepts two arguments. The first is name of the grid element
where values are associated and the second the name of the quantity. The
quantity name must be unique within a group but the same quantity can appear
in multiple goups.

>>> list(grid.at_node.keys()) # There a no values defined at grid nodes
[]
>>> z = grid.add_ones('node', 'topographic__elevation')

We now see that the array has been added to the grid as a reference to the
array returned by ``add_ones``.

>>> list(grid.at_node.keys())
['topographic__elevation']
>>> grid.at_node['topographic__elevation']
array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])
>>> z is grid.at_node['topographic__elevation']
True

To add a previously created array to the grid, use the
:meth:`~.ModelGrid.add_field` method but be aware that it must be of the
correct size (if it's not a ``ValueError`` will be raised).

>>> grid.has_field('node', 'air__temperature')
False
>>> import numpy as np
>>> t = np.zeros(9)
>>> t is grid.add_field('node', 'air__temperature', t)
True
>>> grid.has_field('node', 'air__temperature')
True
>>> grid.has_field('cell', 'air__temperature')
False
>>> t is grid.at_node['air__temperature']
True
"""

import numpy
import numpy as np
import warnings

import six
from six.moves import range

from landlab.testing.decorators import track_this_method
from landlab.utils import count_repeated_values
from landlab.utils.decorators import make_return_array_immutable, deprecated
from landlab.field import ModelDataFields, ScalarDataFields
from landlab.field.scalar_data_fields import FieldError
from landlab.core.model_parameter_dictionary import MissingKeyError
from . import grid_funcs as gfuncs
from ..core.utils import as_id_array
from ..core.utils import add_module_functions_to_class
from .decorators import (override_array_setitem_and_reset, return_id_array,
                         return_readonly_id_array)

#: Indicates an index is, in some way, *bad*.
BAD_INDEX_VALUE = numpy.iinfo(numpy.int32).max
# DEJH thinks the user should be able to override this value if they want

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

#: Indicates a node is *core*.
CORE_NODE = 0

#: Indicates a boundary node is has a fixed values.
FIXED_VALUE_BOUNDARY = 1

#: Indicates a boundary node is has a fixed gradient.
FIXED_GRADIENT_BOUNDARY = 2

#: Indicates a boundary node is wrap-around.
TRACKS_CELL_BOUNDARY = 3

#: Indicates a boundary node is closed
CLOSED_BOUNDARY = 4

# Define the link types

#: Indicates a link is *active*, and can carry flux
ACTIVE_LINK = 0

#: Indicates a link has a fixed (gradient) value, & behaves as a boundary
FIXED_LINK = 2

#: Indicates a link is *inactive*, and cannot carry flux
INACTIVE_LINK = 4

BOUNDARY_STATUS_FLAGS_LIST = [
    FIXED_VALUE_BOUNDARY,
    FIXED_GRADIENT_BOUNDARY,
    TRACKS_CELL_BOUNDARY,
    CLOSED_BOUNDARY,
]
BOUNDARY_STATUS_FLAGS = set(BOUNDARY_STATUS_FLAGS_LIST)

LINK_STATUS_FLAGS_LIST = [
    ACTIVE_LINK,
    FIXED_LINK,
    INACTIVE_LINK,
]
LINK_STATUS_FLAGS = set(LINK_STATUS_FLAGS_LIST)


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
        Tuple of nodes in each coordinate. Nodes are grouped as (*west*,
        *east*, *north*, *south*).

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.base import _sort_points_into_quadrants
    >>> x = np.array([0, 1, 0, -1])
    >>> y = np.array([1, 0, -1, 0])
    >>> nodes = np.array([1, 2, 3, 4])
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
    ('y', 'x')
    >>> _default_axis_names(3)
    ('z', 'y', 'x')
    """
    _DEFAULT_NAMES = ('z', 'y', 'x')
    return _DEFAULT_NAMES[- n_dims:]


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
    return ('-', ) * n_dims


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
    assert ((b1x != 0 and b2y != 0) or (b2x != 0 and b1y != 0)), \
        'Improper unit vectors'

    if b1x != 0. and b2y != 0.:
        ax = (L1 / b1x - L2 * (b1y / (b1x * b2y))) / \
            (1. - (b1y * b2x) / (b1x * b2y))
        ay = L2 / b2y - ax * (b2x / b2y)
    elif b2x != 0. and b1y != 0.:
        ax = (L2 / b2x - L1 * (b2y / (b2x * b1y))) / \
            (1. - (b2y * b1x) / (b2x * b1y))
        ay = L1 / b1y - ax * (b1x / b1y)

    return ax, ay


class ModelGrid(ModelDataFields):
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
    at_core_node : dict-like
        Values at core nodes.
    at_core_cell : dict-like
        Values at core cells.
    at_active_link : dict-like
        Values at active links.
    at_active_face : dict-like
        Values at active faces.

    Other Parameters
    ----------------
    axis_name : tuple, optional
        Name of axes
    axis_units : tuple, optional
        Units of coordinates
    """
    # Debugging flags (if True, activates some output statements)
    _DEBUG_VERBOSE = False
    _DEBUG_TRACK_METHODS = False

    at_node = {}  # : Values defined at nodes
    at_cell = {}  # : Values defined at cells
    at_link = {}  # : Values defined at links
    at_face = {}  # : Values defined at faces
    at_core_node = {}  # : Values defined at core nodes
    at_core_cell = {}  # : Values defined at core cells
    at_active_link = {}  # : Values defined at active links
    at_active_face = {}  # : Values defined at active faces

    # : Nodes on the other end of links pointing into a node.
    node_inlink_matrix = numpy.array([], dtype=numpy.int32)
    # : Nodes on the other end of links pointing out of a node.
    node_outlink_matrix = numpy.array([], dtype=numpy.int32)

    def __init__(self, **kwds):
        super(ModelGrid, self).__init__()
        for element_name in _ARRAY_LENGTH_ATTRIBUTES:
            array_length = self.number_of_elements(element_name)
            try:
                self.new_field_location(element_name, array_length)
            except AttributeError:
                pass

        self.axis_name = kwds.get('axis_name', _default_axis_names(self.ndim))
        self.axis_units = kwds.get(
            'axis_units', _default_axis_units(self.ndim))

        self._link_length = None
        self._all_node_distances_map = None
        self._all_node_azimuths_map = None
        self._node_unit_vector_sum_x = None
        self._node_unit_vector_sum_y = None
        self._link_unit_vec_x = None
        self._link_unit_vec_y = None

    def _initialize(self):
        raise NotImplementedError('_initialize')

    @property
    def ndim(self):
        """Number of spatial dimensions of the grid"""
        return 2

    @property
    @override_array_setitem_and_reset('update_links_nodes_cells_to_new_BCs')
    def status_at_node(self):
        """Get array of the boundary status for each node."""
        return self._node_status

    @status_at_node.setter
    def status_at_node(self, new_status):
        """Set the array of node boundary statuses."""
        self._node_status[:] = new_status[:]
        self.update_links_nodes_cells_to_new_BCs()

    @property
    @make_return_array_immutable
    def neighbors_at_node(self):
        """Get neighboring nodes.

        Examples
        --------
        >>> from landlab import RasterModelGrid, BAD_INDEX_VALUE
        >>> grid = RasterModelGrid((4, 3))
        >>> neighbors = grid.neighbors_at_node.copy()
        >>> neighbors[neighbors == BAD_INDEX_VALUE] = -1
        >>> neighbors # doctest: +NORMALIZE_WHITESPACE
        array([[ 1,  3, -1, -1], [ 2,  4,  0, -1], [-1,  5,  1, -1],
               [ 4,  6, -1,  0], [ 5,  7,  3,  1], [-1,  8,  4,  2],
               [ 7,  9, -1,  3], [ 8, 10,  6,  4], [-1, 11,  7,  5],
               [10, -1, -1,  6], [11, -1,  9,  7], [-1, -1, 10,  8]])
        """
        return self._neighbors_at_node

    @property
    def node_at_cell(self):
        """Node ID associated with grid cells.

        Examples
        --------
        >>> from landlab import RasterModelGrid, BAD_INDEX_VALUE
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.node_at_cell # doctest: +NORMALIZE_WHITESPACE
        array([ 6,  7,  8,
               11, 12, 13])
        """
        return self._node_at_cell

    @property
    def cell_at_node(self):
        """Node ID associated with grid cells.

        Examples
        --------
        >>> from landlab import RasterModelGrid, BAD_INDEX_VALUE
        >>> grid = RasterModelGrid((4, 5))
        >>> ids = grid.cell_at_node
        >>> ids[ids == BAD_INDEX_VALUE] = -1
        >>> ids # doctest: +NORMALIZE_WHITESPACE
        array([-1, -1, -1, -1, -1,
               -1,  0,  1,  2, -1,
               -1,  3,  4,  5, -1,
               -1, -1, -1, -1, -1])
        """
        return self._cell_at_node

    @property
    @return_readonly_id_array
    def active_nodes(self):
        """Get array of active nodes.

        Node IDs of all active (core & open boundary) nodes.
        core_nodes will return just core nodes.
        """
        (active_node_ids, ) = numpy.where(self._node_status != CLOSED_BOUNDARY)
        return active_node_ids

    @property
    @return_readonly_id_array
    def core_nodes(self):
        """Get array of core nodes."""
        try:
            return self._core_nodes
        except AttributeError:
            (core_node_ids, ) = numpy.where(self._node_status == CORE_NODE)
            return core_node_ids

    @property
    @return_readonly_id_array
    def boundary_nodes(self):
        """Get array of boundary nodes."""
        try:
            return self._boundary_nodes
        except:
            (boundary_node_ids, ) = numpy.where(self._node_status != CORE_NODE)
            return boundary_node_ids

    @property
    @return_readonly_id_array
    def open_boundary_nodes(self):
        """Get array of open boundary nodes."""
        (open_boundary_node_ids, ) = numpy.where(
            (self._node_status != CLOSED_BOUNDARY) &
            (self._node_status != CORE_NODE))
        return open_boundary_node_ids

    @property
    @return_readonly_id_array
    def closed_boundary_nodes(self):
        """Get array of closed boundary nodes."""
        (closed_boundary_node_ids, ) = numpy.where(
            self._node_status == CLOSED_BOUNDARY)
        return closed_boundary_node_ids

    @property
    @return_readonly_id_array
    def fixed_gradient_boundary_nodes(self):
        """Get array of fixed gradient boundary nodes."""
        (fixed_gradient_boundary_node_ids, ) = numpy.where(
            self._node_status == FIXED_GRADIENT_BOUNDARY)
        return fixed_gradient_boundary_node_ids

    @property
    @return_readonly_id_array
    def fixed_value_boundary_nodes(self):
        """Get array of fixed value boundary nodes."""
        (fixed_value_boundary_node_ids, ) = numpy.where(
            self._node_status == FIXED_VALUE_BOUNDARY)
        return fixed_value_boundary_node_ids

    @property
    def active_links(self):
        """Get array of active links."""
        try:
            return self.active_link_ids
        except AttributeError:
            self._reset_link_status_list()
            return self.active_link_ids

    @property
    def fixed_links(self):
        """Get array of fixed links."""
        try:
            return self.fixed_link_ids
        except AttributeError:
            self._reset_link_status_list()
            return self.fixed_link_ids

    @property
    @return_readonly_id_array
    def node_at_core_cell(self):
        """Get array of nodes associated with core cells."""
        (core_cell_ids, ) = numpy.where(self._node_status == CORE_NODE)
        return core_cell_ids

    @property
    def core_cells(self):
        """Get array of core cells."""
        return self._core_cells

    @property
    def node_at_link_head(self):
        """Get array of the node at each link head (*to-node*)."""
        return self._node_at_link_head

    @property
    def node_at_link_tail(self):
        """Get array of the node at each link tail (*from-node*)."""
        return self._node_at_link_tail

    @property
    def face_at_link(self):
        """Get array of faces associated with links."""
        return self.link_face

    @property
    def number_of_nodes(self):
        """Total number of nodes."""
        return self._num_nodes

    @property
    def number_of_cells(self):
        """Total number of cells."""
        return self._num_cells

    @property
    def number_of_links(self):
        """Total number of links."""
        return self._num_links

    @property
    def number_of_faces(self):
        """Total number of faces."""
        return self._num_faces

    @property
    def number_of_active_nodes(self):
        """Number of active nodes.

        An active node is either a core node or an open-boundary node.
        """
        return self._num_active_nodes

    @property
    def number_of_core_nodes(self):
        """Number of core nodes.

        A core node is a non-boundary node
        """
        return self._num_core_nodes

    @property
    def number_of_active_cells(self):
        """Number of active cells.

        An active cell includes any possible boundary cells.
        """
        return self._num_active_cells

    @property
    def number_of_core_cells(self):
        """Number of core cells.

        A core cell excludes all boundary cells.
        """
        return self._num_core_cells

    @property
    def number_of_active_links(self):
        """Number of active links."""
        try:
            return self._num_active_links
        except AttributeError:
            self._reset_link_status_list()
            return self._num_active_links

    @property
    def number_of_active_faces(self):
        """Number of active faces."""
        try:
            return self._num_active_faces
        except AttributeError:
            self._reset_link_status_list()
            return self._num_active_faces

    @property
    def number_of_fixed_links(self):
        """Number of fixed links."""
        try:
            return self._num_fixed_links
        except AttributeError:
            self._reset_link_status_list()
            return self._num_fixed_links

    @property
    def number_of_fixed_faces(self):
        """Number of fixed faces."""
        try:
            return self._num_fixed_faces
        except AttributeError:
            self._reset_link_status_list()
            return self._num_fixed_faces

    def number_of_elements(self, element_name):
        """Number of instances of an element.

        Get the number of instances of a grid element in a grid.

        Parameters
        ----------
        element_name : {'node', 'cell', 'link', 'face', 'core_node',
            'core_cell', 'active_link', 'active_face'}
            Name of the grid element.

        Returns
        -------
        int
            Number of elements in the grid.
        """
        try:
            return getattr(self, _ARRAY_LENGTH_ATTRIBUTES[element_name])
        except KeyError:
            raise TypeError('element name not understood')

    @property
    @make_return_array_immutable
    def node_x(self):
        """Get array of the x-coordinates of nodes."""
        return self._node_x

    @property
    @make_return_array_immutable
    def node_y(self):
        """Get array of the y-coordinates of nodes."""
        return self._node_y

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
        """
        AXES = ('node_y', 'node_x')
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
        """
        return self._axis_units

    @axis_units.setter
    def axis_units(self, new_units):
        """Set the units for each coordinate axis."""
        if len(new_units) != self.ndim:
            raise ValueError('length of units does not match grid dimension')
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
        ('y', 'x')
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
        >>> grid.axis_name = ('lon', 'lat')
        >>> grid.axis_name
        ('lon', 'lat')
        """
        if len(new_names) != self.ndim:
            raise ValueError('length of names does not match grid dimension')
        self._axis_name = tuple(new_names)

    @property
    @make_return_array_immutable
    def status_at_link(self):
        """Get array of the status of all links."""
        return self._link_status

    @status_at_node.setter
    def status_at_node(self, new_status_array):
        self._node_status[:] = new_status_array[:]
        self.update_links_nodes_cells_to_new_BCs()

    def active_links_at_node(self, *args):
        """active_links_at_node([node_ids])
        Active links of a node.

        Parameters
        ----------
        node_ids : int or list of ints
                   ID(s) of node(s) for which to find connected active links

        Returns
        -------
        (M, N) ndarray
            The ids of active links attached to grid nodes with
            *node_ids*. If *node_ids* is not given, return links for all of
            the nodes in the grid. M is the number of rows in the grid's
            node_active_inlink_matrix, which can vary depending on the type
            and structure of the grid; in a hex grid, for example, it is 6.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> hmg = HexModelGrid(3, 2)
        >>> hmg.active_links_at_node(3)
        array([[-1],
               [-1],
               [-1],
               [-1],
               [-1],
               [-1],
               [ 0],
               [ 1],
               [ 2],
               [ 3],
               [ 4],
               [ 5]])
        >>> hmg.active_links_at_node()
        array([[ 3,  5,  2, -1,  4,  1,  0],
               [-1, -1, -1, -1, -1, -1, -1],
               [-1, -1, -1, -1, -1, -1, -1],
               [-1, -1, -1, -1, -1, -1, -1],
               [-1, -1, -1, -1, -1, -1, -1],
               [-1, -1, -1, -1, -1, -1, -1],
               [-1, -1, -1,  0, -1, -1, -1],
               [-1, -1, -1,  1, -1, -1, -1],
               [-1, -1, -1,  2, -1, -1, -1],
               [-1, -1, -1,  3, -1, -1, -1],
               [-1, -1, -1,  4, -1, -1, -1],
               [-1, -1, -1,  5, -1, -1, -1]])
        """
        if len(args) == 0:
            return numpy.vstack((self.node_active_inlink_matrix,
                                 self.node_active_outlink_matrix))
        elif len(args) == 1:
            node_ids = numpy.broadcast_arrays(args[0])[0]
            return numpy.vstack(
                (self.node_active_inlink_matrix[:, node_ids],
                 self.node_active_outlink_matrix[:, node_ids])
            ).reshape(2 * numpy.size(self.node_active_inlink_matrix, 0), -1)
        else:
            raise ValueError('only zero or one arguments accepted')

    def active_links_at_node2(self, *args):
        """active_links_at_node2([node_ids])
        Get active links attached to nodes.

        Parameters
        ----------
        node_ids : int or list of ints (optional)
                   ID(s) of node(s) for which to find connected active links.
                   (Default: all nodes)

        Returns
        -------
        (M, N) ndarray
            The link IDs of active links attached to grid nodes with
            *node_ids*. If *node_ids* is not given, return links for all of
            the nodes in the grid. M is the number of rows in the grid's
            node_active_inlink_matrix, which can vary depending on the type
            and structure of the grid; in a hex grid, for example, it is 6.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> hmg = HexModelGrid(3, 2)
        >>> hmg.active_links_at_node2(3)
        array([[-1],
               [-1],
               [-1],
               [-1],
               [-1],
               [-1],
               [ 0],
               [ 1],
               [ 2],
               [ 3],
               [ 4],
               [ 5]])
        >>> hmg.active_links_at_node2()
        array([[ 3,  5,  2, -1,  4,  1,  0],
               [-1, -1, -1, -1, -1, -1, -1],
               [-1, -1, -1, -1, -1, -1, -1],
               [-1, -1, -1, -1, -1, -1, -1],
               [-1, -1, -1, -1, -1, -1, -1],
               [-1, -1, -1, -1, -1, -1, -1],
               [-1, -1, -1,  0, -1, -1, -1],
               [-1, -1, -1,  1, -1, -1, -1],
               [-1, -1, -1,  2, -1, -1, -1],
               [-1, -1, -1,  3, -1, -1, -1],
               [-1, -1, -1,  4, -1, -1, -1],
               [-1, -1, -1,  5, -1, -1, -1]])

        """
        if len(args) == 0:
            return numpy.vstack((self.node_active_inlink_matrix2,
                                 self.node_active_outlink_matrix2))
        elif len(args) == 1:
            node_ids = numpy.broadcast_arrays(args[0])[0]
            return numpy.vstack(
                (self.node_active_inlink_matrix2[:, node_ids],
                 self.node_active_outlink_matrix2[:, node_ids])
            ).reshape(2 * numpy.size(self.node_active_inlink_matrix2, 0), -1)
        else:
            raise ValueError('only zero or one arguments accepted')

    def create_node_array_zeros(self, name=None, **kwds):
        """Create an array of the given type, filled with zeros.

        Returns a 1D numpy array the same length as the number of nodes. If
        user gives optional argument *name*, we add this data to the grid with
        the specified name and return a reference to it; otherwise, we just
        create and return a 1D numpy array.

        Parameters
        ----------
        name : str, optional
            Name of the quantity.

        Returns
        -------
        ndarray
            The newly created array.

        See Also
        --------
        zeros

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid(4,5)
        >>> mydata = np.arange(20, dtype=float)
        >>> rmg.create_node_array_zeros('topographic__elevation')
        array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
                0.,  0.,  0.,  0.,  0.,  0.,  0.])
        >>> rmg.at_node['topographic__elevation'] = mydata
        >>> rmg.at_node['topographic__elevation']
        ...     # doctest: +NORMALIZE_WHITESPACE
        array([  0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,
                10., 11.,  12.,  13.,  14.,  15.,  16.,  17.,  18.,  19.])
        """
        if name is None:
            return numpy.zeros(self.number_of_nodes, **kwds)
        else:
            self.add_zeros('node', name, **kwds)
            return self.at_node[name]

    def create_active_link_array_zeros(self, name=None):
        """Array, filled with zeros, for values at active links.

        Returns a 1D numpy array the same length as the number of nodes. If
        user gives optional argument 'name', we add this data to the grid with
        the specified name and return a reference to it; otherwise, we just
        create and return a 1D numpy array.

        Parameters
        ----------
        name : str
            Name of the quantity.

        Returns
        -------
        ndarray
            The newly created array.

        See Also
        --------
        zeros

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid(3, 3)
        >>> grid.create_active_link_array_zeros()
        array([ 0.,  0.,  0.,  0.])
        """
        if name is None:
            return numpy.zeros(self.number_of_active_links)
        else:
            self.add_zeros('link', name)
            return self.at_link[name]

    def zeros(self, **kwds):
        """Array, filled with zeros, for a given element.

        Returns a numpy array of zeros that is the same length as the number
        of nodes in the grid. Use the *centering* keyword to return an
        array for other elements of the grid. *centering* is a string that is
        one of *node*, *cell*, *link*, or *face*.

        All other keywords are the same as for the numpy zeros function.

        Parameters
        ----------
        centering : str, optional
            Grid element on which the values are defined.

        Returns
        -------
        ndarray
            A newly-allocated array.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.zeros()
        array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
                0.,  0.,  0.,  0.,  0.,  0.,  0.])
        """
        centering = kwds.pop('centering', 'node')
        try:
            return numpy.zeros(self.number_of_elements(centering), **kwds)
        except KeyError:
            raise TypeError(centering)

    def empty(self, **kwds):
        """Array, filled with unititialized values, for a given element.

        Returns a numpy array of uninitialized values that is the same length
        as the number of nodes in the grid. Use the *centering* keyword to
        return an array for other elements of the grid. *centering* is a
        string that is one of *node*, *cell*, *link*, or *face*.

        All other keywords are the same as for the numpy zeros function.

        Parameters
        ----------
        centering : str, optional
            Grid element on which the values are defined.

        Returns
        -------
        ndarray
            A newly-allocated array.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> len(grid.empty())
        20
        """
        centering = kwds.pop('centering', 'node')
        try:
            return numpy.empty(self.number_of_elements(centering), **kwds)
        except KeyError:
            raise TypeError(centering)

    def ones(self, **kwds):
        """Array, filled with ones, for a given element.

        Returns a numpy array of ones that is the same length as the number
        of nodes in the grid. Use the *centering* keyword to return an
        array for other elements of the grid. *centering* is a string that is
        one of *node*, *cell*, *link*, or *face*.

        All other keywords are the same as for the numpy zeros function.

        Parameters
        ----------
        centering : str, optional
            Grid element on which the values are defined.

        Returns
        -------
        ndarray
            A newly-allocated array.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.zeros(dtype=int)
        array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        """
        centering = kwds.pop('centering', 'node')
        try:
            return numpy.ones(self.number_of_elements(centering), **kwds)
        except KeyError:
            raise TypeError(centering)

    def resolve_values_on_links(self, link_values, out=None):
        """Resolve the xy-components of links.

        Resolves values provided defined on links into the x and y directions.
        Returns values_along_x, values_along_y
        """
        return gfuncs.resolve_values_on_links(self, link_values, out=out)

    def resolve_values_on_active_links(self, link_values, out=None):
        """Resolve the xy-components of active links.

        Resolves values provided defined on active links into the x and y
        directions.
        Returns values_along_x, values_along_y
        """
        return gfuncs.resolve_values_on_active_links(self, link_values,
                                                     out=out)

    def node_slopes_using_patches(self, elevs='topographic__elevation',
                                  unit='degrees', return_components=False):
        """
        trial run to extract average local slopes at nodes by the average slope
        of its surrounding patches. DEJH 10/1/14
        elevs either a field name or an nnodes-array.
        unit is 'degrees' or 'radians'.
        If return_components=False (the default), returns the slope magnitude.
        If return_components=True, returns the slope magnitude, then the vector
        (a tuple) of the slope components in the x, y directions. Note the
        slope components will always be returned as rise/run.
        If closed nodes were present in the original array, their values will
        be masked.
        """
        dummy_patch_nodes = numpy.empty(
            (self.nodes_at_patch.shape[0] + 1, self.nodes_at_patch.shape[1]),
            dtype=int)
        dummy_patch_nodes[:-1, :] = self.nodes_at_patch[:]
        dummy_patch_nodes[-1, :] = -1

        # Now any ref to a null node will be -1 in this new
        # (N, patch_max_dim, 4or3) array.
        nodes_on_patches = dummy_patch_nodes[self.patches_at_node()][:, :, :3]
        # Note: we truncate the array to be [N, patch_max_dim,3]; we only
        # need 3 pts per patch, if we're working on a raster

        # Using the wrong values in -1 won't matter, as we'll mask with
        # nodes_on_patches at the end
        node_elevs = numpy.ones((nodes_on_patches.shape[0],
                                 nodes_on_patches.shape[1], 3, 3),
                                dtype=float)

        mask_from_nop = nodes_on_patches[:, :, 0] == -1
        node_elevs[:, :, :, 0] = self.node_x[nodes_on_patches]
        node_elevs[:, :, :, 1] = self.node_y[nodes_on_patches]
        c = numpy.ma.array(numpy.linalg.det(node_elevs), mask=mask_from_nop)
        try:
            node_elevs[:, :, :, 2] = self.at_node[elevs][nodes_on_patches]
        except TypeError:
            node_elevs[:, :, :, 2] = elevs[nodes_on_patches]
        node_elevs[:, :, :, 1] = 1.
        b = numpy.linalg.det(node_elevs)
        node_elevs[:, :, :, 1] = self.node_y[nodes_on_patches]
        node_elevs[:, :, :, 0] = 1.
        a = numpy.linalg.det(node_elevs)

        mask_from_nop = nodes_on_patches[:, :, 0] == -1
        grad_x = -a / c
        grad_y = -b / c  # ...still for each patch
        mean_grad_x = numpy.mean(grad_x, axis=1)
        mean_grad_y = numpy.mean(grad_y, axis=1)

        slope_mag = numpy.arctan(numpy.sqrt(mean_grad_x**2 + mean_grad_y**2))

        if unit == 'radians':
            if not return_components:
                return slope_mag
            else:
                return slope_mag, (mean_grad_x, mean_grad_y)
        if unit == 'degrees':
            if not return_components:
                return 180. / numpy.pi * slope_mag
            else:
                return 180. / numpy.pi * slope_mag, (mean_grad_x, mean_grad_y)
        else:
            raise TypeError("unit must be 'degrees' or 'radians'")

    def node_slopes(self, **kwargs):
        """Array of slopes at nodes.

        This method is an alias for :any:`node_slopes_using_patches`.

        Parameters
        ----------
        elevs : str or ndarray, optional
            Field name or array of node values.
        unit : {'degrees', 'radians'}
            Units for slopes.
        """
        return self.node_slopes_using_patches(**kwargs)

    def aspect(self, slope_component_tuple=None,
               elevs='topographic__elevation', unit='degrees'):
        """Get array of aspect of a surface.

        Calculates at returns the aspect of a surface. Aspect is returned as
        radians clockwise of north, unless input parameter units is set to
        'degrees'.

        If slope_component_tuple is provided, i.e., (slope_x, slope_y), the
        aspect will be calculated from these data.

        If it is not, it will be derived from elevation data at the nodes,
        which can either be a string referring to a grid field (default:
        'topographic__elevation'), or an nnodes-long numpy array of the
        values themselves.
        """
        if slope_component_tuple:
            assert type(slope_component_tuple) == tuple
            assert len(slope_component_tuple) == 2
        else:
            try:
                elev_array = self.at_node[elevs]
            except MissingKeyError:
                assert elevs.size == self.number_of_nodes
                elev_array = elevs
            _, slope_component_tuple = self.node_slopes_using_patches(
                elevs=elev_array, return_components=True)
        angle_from_x_ccw = numpy.arctan2(
            slope_component_tuple[1], slope_component_tuple[0])
        angle_from_N_cw = -(angle_from_x_ccw + numpy.pi / 2.) % (2 * numpy.pi)
        if unit == 'degrees':
            return 180. / numpy.pi * angle_from_N_cw
        elif unit == 'radians':
            return angle_from_N_cw
        else:
            raise TypeError("unit must be 'degrees' or 'radians'")

    def hillshade(self, alt=45., az=315., slp=None, asp=None, unit='degrees',
                  elevs='topographic__elevation'):
        """Get array of hillshade.

        .. codeauthor:: Katy Barnhart <katherine.barnhart@colorado.edu>

        Parameters
        ----------
        alt : float
            Sun altitude (from horizon) - defaults to 45 degrees
        az : float
            Sun azimuth (from north) - defaults to 315 degrees
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
        the ArcGIS algorithm: http://help.arcgis.com/en/arcgisdesktop/10.0/
        help/index.html#/How_Hillshade_works/009z000000z2000000/ .

        NB: grid.node_slopes_using_patches() returns slopes as RADIANS. Use
        caution.

        Remember when plotting that bright areas have high values. cmap='Greys'
        will give an apparently inverted color scheme. *cmap='gray'* has white
        associated with the high values, so is recommended for plotting.

        """
        if slp is not None and asp is not None:
            if unit == 'degrees':
                (alt, az, slp, asp) = (numpy.radians(alt), numpy.radians(az),
                                       numpy.radians(slp), numpy.radians(asp))
            elif unit == 'radians':
                if alt > numpy.pi / 2. or az > 2. * numpy.pi:
                    six.print_(
                        'Assuming your solar properties are in degrees, '
                        'but your slopes and aspects are in radians...')
                    (alt, az) = (numpy.radians(alt), numpy.radians(az))
                    # ...because it would be super easy to specify radians,
                    # but leave the default params alone...
            else:
                raise TypeError("unit must be 'degrees' or 'radians'")
        elif slp is None and asp is None:
            if unit == 'degrees':
                (alt, az) = (numpy.radians(alt), numpy.radians(az))
            elif unit == 'radians':
                pass
            else:
                raise TypeError("unit must be 'degrees' or 'radians'")
            slp, slp_comps = self.node_slopes_using_patches(
                elevs, unit='radians', return_components=True)
            asp = self.aspect(slope_component_tuple=slp_comps, unit='radians')
        else:
            raise TypeError('Either both slp and asp must be set, or neither!')

        shaded = (
            numpy.sin(alt) * numpy.cos(slp) +
            numpy.cos(alt) * numpy.sin(slp) * numpy.cos(az - asp)
        )

        return shaded

    def calculate_flux_divergence_at_core_nodes(self, active_link_flux,
                                                net_unit_flux=None):
        r"""Get array of flux divergence for core nodes.

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

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid(4, 5, 1.0)
        >>> u = [0., 1., 2., 3., 0.,
        ...      1., 2., 3., 2., 3.,
        ...      0., 1., 2., 1., 2.,
        ...      0., 0., 2., 2., 0.]
        >>> u = np.array(u)
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

        >>> divflux = rmg.zeros(centering='core_cell') # outside loop

        Then do this inside the loop:

        >>> divflux = rmg.calculate_flux_divergence_at_core_nodes(
        ...     flux, divflux)

        In this case, the function will not have to create the divflux array.

        Note this method is untested with looped boundary conditions.
        """

        if self._DEBUG_TRACK_METHODS:
            six.print_('ModelGrid.calculate_flux_divergence_at_core_nodes')

        assert (len(active_link_flux) == self.number_of_active_links), \
            "incorrect length of active_link_flux array"

        # If needed, create net_unit_flux array
        if net_unit_flux is None:
            net_unit_flux = numpy.zeros(self.number_of_core_nodes)
        else:
            net_unit_flux[:] = 0.

        assert (len(net_unit_flux)) == self.number_of_core_nodes

        node_net_unit_flux = self.calculate_flux_divergence_at_nodes(
            active_link_flux)

        node_at_core_cell = self.node_at_cell[self.core_cells]
        net_unit_flux = node_net_unit_flux[node_at_core_cell]

        return net_unit_flux

    @track_this_method
    def calculate_flux_divergence_at_nodes(self, active_link_flux, out=None):
        """Flux divergence at nodes.

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
        return gfuncs.calculate_flux_divergence_at_nodes(self,
                                                         active_link_flux,
                                                         out=out)

    @property
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
        >>> grid = RasterModelGrid((4, 5), spacing=(3, 4))
        >>> grid.status_at_node[7] = CLOSED_BOUNDARY
        >>> grid.cell_area_at_node
        array([  0.,   0.,   0.,   0.,   0.,
                 0.,  12.,  12.,  12.,   0.,
                 0.,  12.,  12.,  12.,   0.,
                 0.,   0.,   0.,   0.,   0.])
        """
        try:
            return self._cell_area_at_node
        except AttributeError:
            return self._setup_cell_areas_array_force_inactive()

    @property
    def face_widths(self):
        """Width of grid faces."""
        try:
            return self._face_widths
        except AttributeError:
            return self._setup_face_widths()

    def _setup_cell_areas_array_force_inactive(self):
        """Set up an array of cell areas that is n_nodes long.

        Sets up an array of cell areas that is nnodes long. Nodes that have
        cells receive the area of that cell. Nodes which do not, receive
        zeros.
        """
        _cell_area_at_node_zero = numpy.zeros(self.number_of_nodes,
                                              dtype=float)
        _cell_area_at_node_zero[self.node_at_cell] = self.area_of_cell
        self._cell_area_at_node = _cell_area_at_node_zero
        return self._cell_area_at_node

    def get_active_link_connecting_node_pair(self, node1, node2):
        """Get the active link that connects a pair of nodes.

        Returns the ID number of the active link that connects the given pair
        of nodes, or BAD_INDEX_VALUE if not found.
        This method is slow, and can only take single ints as *node1* and
        *node2*. It should ideally be overridden for optimal functionality in
        more specialized grid modules (e.g., raster).

        Examples
        --------
        >>> import landlab as ll
        >>> rmg = ll.RasterModelGrid(4, 5)
        >>> rmg.get_active_link_connecting_node_pair(8, 3)
        array([2])
        """
        active_link = BAD_INDEX_VALUE
        for alink in range(0, self.number_of_active_links):
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
    @make_return_array_immutable
    def area_of_cell(self):
        """Get areas of grid cells.

        Returns
        -------
        ndarray
            Areas of cells, in ID order.
        """
        # return self._area_of_cell
        try:
            return self._area_of_cell
        except AttributeError:
            return self._setup_cell_areas_array()

    @property
    def link_length(self):
        """Get lengths of links.

        Returns
        -------
        ndarray
            Lengths of all links, in ID order.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.link_length
        array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,
                1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,
                1.,  1.,  1.,  1.,  1.])
        >>> len(grid.link_length) == grid.number_of_links
        True
        """
        if self._link_length is None:
            return self._calculate_link_length()
        else:
            return self._link_length

    def min_active_link_length(self):
        """Get length of the shortest active link.

        Returns
        -------
        float
            Length of the shortest active link in the grid.
        """
        return numpy.amin(self.link_length[self.active_link_ids])

    def max_active_link_length(self):
        """Get the length of the longest active link.

        Returns
        -------
        float
            Length of the longest active link in the grid.
        """
        return numpy.amax(self.link_length[self.active_link_ids])

    def _calculate_link_length(self):
        """Get array of the lengths of all links.

        Calculates, returns, and stores as a property of the grid the lengths
        of all the links in the grid.
        """
        if self._link_length is None:
            self._link_length = self.empty(centering='link', dtype=float)

        diff_x = (self.node_x[self.node_at_link_tail] -
                  self.node_x[self.node_at_link_head])
        diff_y = (self.node_y[self.node_at_link_tail] -
                  self.node_y[self.node_at_link_head])
        numpy.sqrt(diff_x ** 2 + diff_y ** 2, out=self._link_length)

        return self._link_length

    def assign_upslope_vals_to_active_links(self, u, v=None):
        """Assign upslope node value to link.

        Assigns to each active link the value of *u* at whichever of its
        neighbors has a higher value of *v*. If *v* is omitted, uses *u* for
        both. The order of the link values is by link ID.

        Parameters
        ----------
        u : array-like
            Node values to assign to links.
        v : array-like, optional
            Node values to test for upslope-ness.

        Returns
        -------
        ndarray
            Values at active links.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> import numpy as np
        >>> grid = RasterModelGrid(3, 3)
        >>> u = np.arange(9.)
        >>> grid.assign_upslope_vals_to_active_links(u)
        array([ 4.,  7.,  4.,  5.])
        """
        if v is None:
            v = numpy.array((0., ))

        fv = numpy.zeros(self.number_of_active_links)
        if len(v) < len(u):
            for i in range(0, self.number_of_active_links):
                fv[i] = max(u[self.activelink_fromnode[i]],
                            u[self.activelink_tonode[i]])
        else:
            for i in range(0, self.number_of_active_links):
                if (v[self.activelink_fromnode[i]] >
                        v[self.activelink_tonode[i]]):
                    fv[i] = u[self.activelink_fromnode[i]]
                else:
                    fv[i] = u[self.activelink_tonode[i]]
        return fv

    def _reset_link_status_list(self):
        """Create of reset a list of links statuses.

        Creates or resets a list of link statuses. We do this by sweeping
        through the given lists of from and to nodes, and checking the status
        of these as given in the node_status list. A link is active if both its
        nodes are core, or if one is core and the other is fixed value.
        A link is inactive if either node is closed.
        A link is fixed if either node is fixed gradient.

        Note that by default, any link which has been previously set as fixed
        will remain so, and if a closed-core node pair is found at each of its
        ends, the closed node will be converted to a fixed gradient node. If
        you want to close a node which has a fixed link already connected to
        it, first change the link status to inactive.

        A further test is performed to ensure that the final maps of node and
        link status are internally consistent.
        """
        if self._DEBUG_TRACK_METHODS:
            six.print_('ModelGrid._reset_link_status_list')

        try:
            already_fixed = self._link_status == FIXED_LINK
        except AttributeError:
            already_fixed = numpy.zeros(self.number_of_links, dtype=bool)

        fromnode_status = self._node_status[self.node_at_link_tail]
        tonode_status = self._node_status[self.node_at_link_head]

        if not numpy.all((fromnode_status[already_fixed] ==
                          FIXED_GRADIENT_BOUNDARY) |
                         (tonode_status[already_fixed] ==
                          FIXED_GRADIENT_BOUNDARY)):
            assert numpy.all(np.logical_not((fromnode_status[already_fixed] ==
                                             CLOSED_BOUNDARY) &
                                            (tonode_status[already_fixed] ==
                                             CLOSED_BOUNDARY)))
            fromnode_status[already_fixed] = numpy.where(
                (fromnode_status[already_fixed] == CLOSED_BOUNDARY) &
                (tonode_status[already_fixed] == CORE_NODE),
                FIXED_GRADIENT_BOUNDARY,
                fromnode_status[already_fixed])
            tonode_status[already_fixed] = numpy.where(
                (tonode_status[already_fixed] == CLOSED_BOUNDARY) &
                (fromnode_status[already_fixed] == CORE_NODE),
                FIXED_GRADIENT_BOUNDARY,
                tonode_status[already_fixed])
            print("""
                  Remember, fixed_links are dominant over node statuses.
                  Your grid may have had an incompatibility between
                  fixed_links and closed nodes, which has been resolved by
                  converting the closed nodes to fixed gradient nodes. If
                  you were trying to deliberately close a node which had
                  once been set to fixed gradient, you need to open the
                  links before changing the node statuses. If you were
                  setting a node to fixed_value, you can ignore this
                  message.
                  """)

        active_links = (((fromnode_status == CORE_NODE) & ~
                         (tonode_status == CLOSED_BOUNDARY)) |
                        ((tonode_status == CORE_NODE) & ~
                         (fromnode_status == CLOSED_BOUNDARY)))
        # ...this still includes things that will become fixed_link

        fixed_links = ((((fromnode_status == FIXED_GRADIENT_BOUNDARY) &
                         (tonode_status == CORE_NODE)) |
                        ((tonode_status == FIXED_GRADIENT_BOUNDARY) &
                         (fromnode_status == CORE_NODE))) |
                       already_fixed)
        
        fixed_link_fixed_val = (((fromnode_status == FIXED_VALUE_BOUNDARY) |
                                 (tonode_status == FIXED_VALUE_BOUNDARY)) &
                                already_fixed)                        
        # these are the "special cases", where the user is probably trying to
        # adjust an individual fixed_link back to fixed value. We'll allow it:
        fixed_links[fixed_link_fixed_val] = False

        try:
            self._link_status.fill(4)
        except AttributeError:
            self._link_status = numpy.empty(self.number_of_links, dtype=int)
            self._link_status.fill(4)

        self._link_status[active_links] = 0

        self._link_status[fixed_links] = 2

        active_links = self._link_status == 0  # now it's correct
        (self.active_link_ids, ) = numpy.where(active_links)
        (self.fixed_link_ids, ) = numpy.where(fixed_links)
        self.active_link_ids = as_id_array(self.active_link_ids)
        self.fixed_link_ids = as_id_array(self.fixed_link_ids)

        self._num_active_links = (active_links).sum()
        self._num_active_faces = self._num_active_links
        self._num_fixed_links = fixed_links.sum()
        self._num_fixed_faces = self._num_fixed_links
        self.activelink_fromnode = self.node_at_link_tail[active_links]
        self.activelink_tonode = self.node_at_link_head[active_links]

        # Set up active inlink and outlink matrices
        self._setup_active_inlink_and_outlink_matrices()

    def _reset_lists_of_nodes_cells(self):
        """Create of reset lists of nodes and cells based on their status.

        Creates or resets various lists of nodes and cells based on their
        statuses. Call this function whenever you make changes to the
        boundary conditions in the grid.
        The updated attributes and arrays are:
        * _num_active_nodes
        * _num_active_cells
        * _num_core_nodes
        * _num_core_cells
        * activecell_node *
        * corecell_node *
        * active_cells
        * core_cells
        * _boundary_nodes

        Examples
        --------
        >>> import landlab
        >>> grid = landlab.RasterModelGrid((4, 5))
        >>> grid.status_at_node[7] = landlab.CLOSED_BOUNDARY
        >>> grid.core_cells
        array([0, 2, 3, 4, 5])
        """
        (self._core_nodes, ) = numpy.where(self._node_status == CORE_NODE)
        self._num_core_nodes = self._core_nodes.size

        self._core_cells = self.cell_at_node[self._core_nodes]
        self._num_core_cells = self._core_cells.size

        self._boundary_nodes = as_id_array(
            numpy.where(self._node_status != CORE_NODE)[0])

    def update_links_nodes_cells_to_new_BCs(self):
        """Update grid element connectivity, status.

        This method updates all of the various lists and attributes governed
        by node status (e.g., core nodes, active links, etc) when you change
        node statuses. Call it if your method or driver makes changes to the
        boundary conditions of nodes in the grid.
        """
        self._reset_link_status_list()
        self._reset_lists_of_nodes_cells()
        try:
            if self.diagonal_list_created:
                self.diagonal_list_created = False
        except AttributeError:
            pass
        try:
            if self.neighbor_list_created:
                self.neighbor_list_created = False
        except AttributeError:
            pass

    @deprecated
    def set_nodata_nodes_to_inactive(self, node_data, nodata_value):
        """Make no-data nodes inactive.

        Set the status to CLOSED_BOUNDARY for all nodes whose value
        of node_data is equal to the nodata_value.

        .. note:: Deprecated since version 0.6.
            Deprecated due to out of date terminology;
            use :func:`set_nodata_nodes_to_closed` instead.

        Parameters
        ----------
        node_data : ndarray
            Data values.
        nodata_value : float
            Value that indicates an invalid value.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid(3, 4, 1.0)
        >>> mg.status_at_node
        array([1, 1, 1, 1,
               1, 0, 0, 1,
               1, 1, 1, 1], dtype=int8)
        >>> h = np.array([-9999, -9999, -9999, -9999,
        ...               -9999, -9999, 12345.,   0.,
        ...               -9999,    0.,     0.,   0.])
        >>> mg.set_nodata_nodes_to_inactive(h, -9999)
        >>> mg.status_at_node
        array([4, 4, 4, 4,
               4, 4, 0, 1,
               4, 1, 1, 1], dtype=int8)
        """
        self.set_nodata_nodes_to_closed(node_data, nodata_value)

    def set_nodata_nodes_to_closed(self, node_data, nodata_value):
        """Make no-data nodes closed boundaries.

        Sets node status to :any:`CLOSED_BOUNDARY` for all nodes whose value
        of *node_data* is equal to the *nodata_value*.

        Any links connected to :any:`CLOSED_BOUNDARY` nodes are automatically
        set to :any:`INACTIVE_LINK` boundary.

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

          Links set to :any:`ACTIVE_LINK` are not shown in this diagram.

        ``*`` indicates the nodes that are set to :any:`CLOSED_BOUNDARY`

        ``o`` indicates the nodes that are set to :any:`CORE_NODE`

        ``I`` indicates the links that are set to :any:`INACTIVE_LINK`

        >>> import numpy as np
        >>> import landlab as ll
        >>> mg = ll.RasterModelGrid(3, 4, 1.0)
        >>> mg.status_at_node
        array([1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1], dtype=int8)
        >>> h = np.array([-9999, -9999, -9999, -9999, -9999, -9999, 12345.,
        ...     0., -9999, 0., 0., 0.])
        >>> mg.set_nodata_nodes_to_closed(h, -9999)
        >>> mg.status_at_node
        array([4, 4, 4, 4, 4, 4, 0, 1, 4, 1, 1, 1], dtype=int8)
        """
        # Find locations where value equals the NODATA code and set these nodes
        # as inactive boundaries.
        nodata_locations = numpy.nonzero(node_data == nodata_value)
        self._node_status[nodata_locations] = CLOSED_BOUNDARY

        # Recreate the list of active cell IDs
        self.update_links_nodes_cells_to_new_BCs()

    def set_nodata_nodes_to_fixed_gradient(self, node_data, nodata_value):
        """Make no-data nodes fixed gradient boundaries.

        Set node status to :any:`FIXED_GRADIENT_BOUNDARY` for all nodes
        whose value of *node_data* is equal to *nodata_value*.

        Any links between :any:`FIXED_GRADIENT_BOUNDARY` nodes and
        :any:`CORE_NODES` are automatically set to :any:`FIXED_LINK` boundary
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

            Links set to :any:`ACTIVE_LINK` are not shown in this diagram.

        ``X`` indicates the links that are set to :any:`FIXED_LINK`

        ``I`` indicates the links that are set to :any:`INACTIVE_LINK`

        ``o`` indicates the nodes that are set to :any:`CORE_NODE`

        ``*`` indicates the nodes that are set to
              :any:`FIXED_GRADIENT_BOUNDARY`

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> rmg = RasterModelGrid(4, 9)
        >>> rmg.status_at_node # doctest: +NORMALIZE_WHITESPACE
        array([1, 1, 1, 1, 1, 1, 1, 1, 1,
               1, 0, 0, 0, 0, 0, 0, 0, 1,
               1, 0, 0, 0, 0, 0, 0, 0, 1,
               1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=int8)

        >>> z = rmg.create_node_array_zeros()
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
               2, 2, 2, 2, 2, 2, 2, 2, 2], dtype=int8)

        >>> rmg.status_at_link # doctest: +NORMALIZE_WHITESPACE
        array([4, 4, 4, 2, 2, 2, 2, 2, 4,
               4, 4, 4, 0, 0, 0, 0, 0, 4,
               4, 4, 4, 2, 2, 2, 2, 2, 4,
               4, 4, 4, 4, 4, 4, 4, 4,
               4, 4, 2, 0, 0, 0, 0, 2,
               4, 4, 2, 0, 0, 0, 0, 2,
               4, 4, 4, 4, 4, 4, 4, 4])
        """
        # Find locations where value equals the NODATA code and set these nodes
        # as inactive boundaries.
        nodata_locations = numpy.nonzero(node_data == nodata_value)
        self._node_status[nodata_locations] = FIXED_GRADIENT_BOUNDARY

        # Recreate the list of active cell IDs
        self.update_links_nodes_cells_to_new_BCs()

    def max_of_link_end_node_values(self, node_data):
        """Maximum value at the end of links.

        For each active link, finds and returns the maximum value of node_data
        at either of the two ends. Use this, for example, if you want to find
        the maximum value of water depth at linked pairs of nodes (by passing
        in an array of water depth values at nodes).

        Parameters
        ----------
        node_data : ndarray
            Values at grid nodes.

        Returns
        -------
        ndarray :
            Maximum values whose length is the number of active links.

        Examples
        --------
        >>> import numpy as np
        >>> import landlab as ll
        >>> mg = ll.RasterModelGrid(3, 4, 1.0)
        >>> h = np.array([2.,2.,8.,0.,8.,0.,3.,0.,5.,6.,8.,3.])
        >>> mg.max_of_link_end_node_values(h)
        array([ 2.,  8.,  6.,  8.,  8.,  3.,  3.])
        """
        return numpy.maximum(node_data[self.activelink_fromnode],
                             node_data[self.activelink_tonode])

    def calculate_numbers_of_node_neighbors(self):
        """Number of neighbor nodes.

        Calculates the number of neighboring nodes for each node, and returns
        the result as a 1D numpy array. Used to find the maximum number of
        neighbors, so that inlink and outlink matrices can be dimensioned
        accordingly. Assumes that self.number_of_nodes, self.node_at_link_tail,
        and self.node_at_link_head have already been set up.

        Algorithm works by simply looping through all links; for each, the
        endpoints are neighbors of one another, so we increment the number of
        neighbors for both the endpoint nodes.
        """
        num_nbrs = numpy.zeros(self.number_of_nodes, dtype=int)
        for link in range(self.number_of_links):
            num_nbrs[self.node_at_link_tail[link]] += 1
            num_nbrs[self.node_at_link_head[link]] += 1
        return num_nbrs

    def _setup_inlink_and_outlink_matrices(self):
        """Create data structured for number of inlinks and outlinks.

        Creates data structures to record the numbers of inlinks and outlinks
        for each node. An inlink of a node is simply a link that has the
        node as its "to" node, and an outlink is a link that has the node
        as its "from".

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
        link number -1, Python takes the value of the last element in the
        array, which is zero. By doing it this way, we maintain the efficiency
        that comes with the use of numpy. Again, more info can be found in the
        description of the flux divergence functions.
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
        tonodes = self.node_at_link_head
        self.node_numinlink = numpy.bincount(tonodes,
                                             minlength=self.number_of_nodes)

        counts = count_repeated_values(self.node_at_link_head)
        for (count, (tonodes, link_ids)) in enumerate(counts):
            self.node_inlink_matrix[count][tonodes] = link_ids

        # Set up the outlink arrays
        fromnodes = self.node_at_link_tail
        self.node_numoutlink = numpy.bincount(fromnodes,
                                              minlength=self.number_of_nodes)
        counts = count_repeated_values(self.node_at_link_tail)
        for (count, (fromnodes, link_ids)) in enumerate(counts):
            self.node_outlink_matrix[count][fromnodes] = link_ids

    def _setup_active_inlink_and_outlink_matrices(self):
        """Create data structures for number of active inlinks and outlinks.

        Creates data structures to record the numbers of active inlinks and
        active outlinks for each node. These data structures are equivalent to
        the "regular" inlink and outlink matrices, except that it uses the IDs
        of active links (only).

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> hg = HexModelGrid(3, 2)
        >>> hg.node_numactiveinlink
        array([1, 1, 1, 0, 1, 1, 1])
        >>> hg.node_active_inlink_matrix2
        array([[ 3,  5,  2, -1,  4,  1,  0],
               [-1, -1, -1, -1, -1, -1, -1],
               [-1, -1, -1, -1, -1, -1, -1],
               [-1, -1, -1, -1, -1, -1, -1],
               [-1, -1, -1, -1, -1, -1, -1],
               [-1, -1, -1, -1, -1, -1, -1]])
        >>> hg.node_numactiveoutlink
        array([0, 0, 0, 6, 0, 0, 0])
        >>> hg.node_active_outlink_matrix2
        array([[-1, -1, -1,  0, -1, -1, -1],
               [-1, -1, -1,  1, -1, -1, -1],
               [-1, -1, -1,  2, -1, -1, -1],
               [-1, -1, -1,  3, -1, -1, -1],
               [-1, -1, -1,  4, -1, -1, -1],
               [-1, -1, -1,  5, -1, -1, -1]])
        """
        # Create active in-link and out-link matrices.
        self.node_active_inlink_matrix = - numpy.ones(
            (self.max_num_nbrs, self.number_of_nodes), dtype=numpy.int)
        self.node_active_outlink_matrix = - numpy.ones(
            (self.max_num_nbrs, self.number_of_nodes), dtype=numpy.int)

        # Set up the inlink arrays
        tonodes = self.activelink_tonode
        self.node_numactiveinlink = as_id_array(numpy.bincount(
            tonodes, minlength=self.number_of_nodes))

        counts = count_repeated_values(self.activelink_tonode)
        for (count, (tonodes, active_link_ids)) in enumerate(counts):
            self.node_active_inlink_matrix[count][tonodes] = active_link_ids

        # Set up the outlink arrays
        fromnodes = self.activelink_fromnode
        self.node_numactiveoutlink = as_id_array(numpy.bincount(
            fromnodes, minlength=self.number_of_nodes))
        counts = count_repeated_values(self.activelink_fromnode)
        for (count, (fromnodes, active_link_ids)) in enumerate(counts):
            self.node_active_outlink_matrix[count][fromnodes] = active_link_ids

        # THE FOLLOWING IS MEANT TO REPLACE THE ABOVE CODE, USING LINK IDS
        # FOR ACTIVE LINKS (ONLY), INSTEAD OF "ACTIVE LINK IDS". THE POINT IS
        # TO HAVE JUST ONE ID/NUMBERING SYSTEM FOR LINKS, RATHER THAN A
        # SEPARATE NUMBERING SYSTEM FOR ACTIVE LINKS
        # GT JUNE 2015
        # TODO: CLEAN THIS UP

        # Create AN ALTERNATIVE VERSION OF active in-link and out-link
        # matrices, WHICH WILL EVENTUALLY REPLACE THE ONE ABOVE (AND BE
        # RENAMED TO GET RID OF THE "2")
        # TODO: MAKE THIS CHANGE ONCE CODE THAT USES IT HAS BEEN PREPPED
        self.node_active_inlink_matrix2 = - numpy.ones(
            (self.max_num_nbrs, self.number_of_nodes), dtype=numpy.int)
        self.node_active_outlink_matrix2 = - numpy.ones(
            (self.max_num_nbrs, self.number_of_nodes), dtype=numpy.int)

        # Set up the inlink arrays
        tonodes = self.node_at_link_head[self.active_links]
        self.node_numactiveinlink = as_id_array(numpy.bincount(
            tonodes, minlength=self.number_of_nodes))

        # OK, HERE WE HAVE TO MAKE A CHANGE, BECAUSE THE INDICES RETURNED BY
        # count_repeated_values ARE "ACTIVE LINK INDICES", WHICH WE ARE NO
        # LONGER USING. HAVE TO TURN THESE BACK INTO LINK IDS. I THINK WE CAN
        # DO THIS BY CHANGING active_link_ids TO
        # self.active_links[active_link_ids] BUT HAVEN'T MADE THIS CHANGE YET.
        # NEED TO WORK THROUGH EXAMPLE 3,2 HMG
        counts = count_repeated_values(
            self.node_at_link_head[self.active_links])
        for (count, (tonodes, active_link_ids)) in enumerate(counts):
            self.node_active_inlink_matrix2[count][
                tonodes] = self.active_links[active_link_ids]

        # Set up the outlink arrays
        fromnodes = self.node_at_link_tail[self.active_links]
        self.node_numactiveoutlink = as_id_array(numpy.bincount(
            fromnodes, minlength=self.number_of_nodes))
        counts = count_repeated_values(self.activelink_fromnode)
        for (count, (fromnodes, active_link_ids)) in enumerate(counts):
            self.node_active_outlink_matrix2[count][
                fromnodes] = self.active_links[active_link_ids]

    def _make_link_unit_vectors(self):
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

        Creates the following:
        * ``self.link_unit_vec_x``, ``self.link_unit_vec_y`` : ndarray
            x and y components of unit vectors at each link (extra 0
            entries @ end)
        * ``self.node_vector_sum_x``, ``self.node_vector_sum_y`` : ndarray
            Sums of x & y unit vector components for each node. Sum is over all
            links connected to a given node.

        Examples
        --------
        The example below is a seven-node hexagonal grid, with six nodes around
        the perimeter and one node (#3) in the interior. There are four
        horizontal links with unit vector (1,0), and 8 diagonal links with
        unit vector (+/-0.5, +/-sqrt(3)/2) (note: sqrt(3)/2 ~ 0.866).

        .. note::

            This example assumes that the triangulation places links in a
            certain order. Because the order is arbitrary, this might break on
            different platforms. If that happens, the example needs to be
            made generic somehow ...

        >>> import landlab as ll
        >>> hmg = ll.HexModelGrid(3, 2, 2.0)
        >>> hmg.link_unit_vec_x # doctest: +NORMALIZE_WHITESPACE
        array([ 0.5, -0.5, -1. , -0.5,  1. ,  0.5,  0.5, -1. , -0.5,  0.5,
                1. , -0.5,  0. ])
        >>> hmg.link_unit_vec_y
        array([ 0.8660254,  0.8660254,  0.       , -0.8660254,  0.       ,
               -0.8660254, -0.8660254,  0.       , -0.8660254, -0.8660254,
                0.       , -0.8660254,  0.       ])
        >>> hmg.node_unit_vector_sum_x
        array([ 2.,  2.,  2.,  4.,  2.,  2.,  2.])
        >>> hmg.node_unit_vector_sum_y
        array([ 1.73205081,  1.73205081,  1.73205081,  3.46410162,  1.73205081,
                1.73205081,  1.73205081])
        """
        # Create the arrays for unit vectors for each link. These each get an
        # additional array element at the end with the value zero. This allows
        # any references to "link ID -1" in the node_inlink_matrix and
        # node_outlink_matrix to refer to the zero value in this extra element,
        # so that when we're summing up link unit vectors, or multiplying by a
        # nonexistent unit vector, we end up just treating these as zero.
        self._link_unit_vec_x = numpy.zeros(self.number_of_links + 1)
        self._link_unit_vec_y = numpy.zeros(self.number_of_links + 1)

        # Calculate the unit vectors using triangle similarity and the
        # Pythagorean Theorem.
        dx = self.node_x[self.node_at_link_head] - \
            self.node_x[self.node_at_link_tail]
        dy = self.node_y[self.node_at_link_head] - \
            self.node_y[self.node_at_link_tail]
        self._link_unit_vec_x[:self.number_of_links] = dx / self.link_length
        self._link_unit_vec_y[:self.number_of_links] = dy / self.link_length

        # While we're at it, calculate the unit vector sums for each node.
        # These will be useful in averaging link-based vectors at the nodes.
        self._node_unit_vector_sum_x = numpy.zeros(self.number_of_nodes)
        self._node_unit_vector_sum_y = numpy.zeros(self.number_of_nodes)
        max_num_inlinks_per_node = numpy.size(self.node_inlink_matrix, 0)
        for i in range(max_num_inlinks_per_node):
            self._node_unit_vector_sum_x += abs(
                self._link_unit_vec_x[self.node_inlink_matrix[i, :]])
            self._node_unit_vector_sum_y += abs(
                self._link_unit_vec_y[self.node_inlink_matrix[i, :]])
            self._node_unit_vector_sum_x += abs(
                self._link_unit_vec_x[self.node_outlink_matrix[i, :]])
            self._node_unit_vector_sum_y += abs(
                self._link_unit_vec_y[self.node_outlink_matrix[i, :]])

    @property
    def link_unit_vec_x(self):
        """Get array of x-component of unit vector for links.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 3))
        >>> len(grid.link_unit_vec_x) == grid.number_of_links + 1
        True
        >>> grid.link_unit_vec_x # doctest: +NORMALIZE_WHITESPACE
        array([ 0.,  0.,  0.,  0.,  0.,  0.,
                1.,  1.,  1.,  1.,  1.,  1.,  0.])
        """
        if self._link_unit_vec_x is None:
            self._make_link_unit_vectors()
        return self._link_unit_vec_x

    @property
    def link_unit_vec_y(self):
        """Get array of y-component of unit vector for links.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 3))
        >>> len(grid.link_unit_vec_y) == grid.number_of_links + 1
        True
        >>> grid.link_unit_vec_y # doctest: +NORMALIZE_WHITESPACE
        array([ 1.,  1.,  1.,  1.,  1.,  1.,
                0.,  0.,  0.,  0.,  0.,  0.,  0.])
        """
        if self._link_unit_vec_y is None:
            self._make_link_unit_vectors()
        return self._link_unit_vec_y

    @property
    def node_unit_vector_sum_x(self):
        """Get array of x-component of unit vector sums at each node.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 3))
        >>> len(grid.node_unit_vector_sum_x) == grid.number_of_nodes
        True
        >>> grid.node_unit_vector_sum_x
        array([ 1.,  2.,  1.,  1.,  2.,  1.,  1.,  2.,  1.])
        """
        if self._node_unit_vector_sum_x is None:
            self._make_link_unit_vectors()
        return self._node_unit_vector_sum_x

    @property
    def node_unit_vector_sum_y(self):
        """Get array of y-component of unit vector sums at each node.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 3))
        >>> len(grid.node_unit_vector_sum_y) == grid.number_of_nodes
        True
        >>> grid.node_unit_vector_sum_y
        array([ 1.,  1.,  1.,  2.,  2.,  2.,  1.,  1.,  1.])
        """
        if self._node_unit_vector_sum_y is None:
            self._make_link_unit_vectors()
        return self._node_unit_vector_sum_y

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
        _make_link_unit_vectors : sets up unit vectors at links and unit-vector
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
        is by example. Consider a 3x4 raster grid:

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

            Corner nodes (0, 3, 8, 11): vsx = 1, vsy = 1
            Bottom and top nodes (1-2, 9-10): vsx = 2, vsy = 1
            Left and right nodes (4, 7): vsx = 1, vsy = 2
            All others: vsx = 2, vsy = 2

        The process of creating unit-vector sums at nodes is handled by
        ModelGrid._make_link_unit_vectors() (and, for raster grids, by the
        overriding method RasterModelGrid._make_link_unit_vectors()). The node
        unit-vector sums are then stored in self.node_unit_vector_sum_x and
        self.node_unit_vector_sum_y.

        How would you use this? Suppose you have a vector variable q defined at
        links. What's the average at the nodes? We'll define the average as
        follows.  The terminology here is: :math:`q = (u,v)` represents the
        vector quantity defined at links, :math:`Q = (U,V)` represents its
        definition at nodes, :math:`(m,n)` represents the unit vector
        components at a link, and :math:`(S_x,S_y)` represents the unit-vector
        sum at a given node.

        ..math::

            U_i = \sum_{j=1}^{L_i} q_j m_j / S_{xi}
            V_i = \sum_{j=1}^{L_i} q_j n_j / S_{yi}

        Suppose that the vector q is uniform and equal to one.
        Then, at node 0 in the above grid, this works out to:

        ..math::

            U_0 = (q_0 m_0) / 1 + (q_8 m_8) / 1 = (1 0)/ 1 + (1 1)/1 = 1
            V_0 = (q_0 n_0) / 1 + (q_8 n_8) / 1 = (1 1) / 1 + (1 0) / 1 = 1

        At node 1, in the bottom row but not a corner, we add up the values
        of **q** associated with THREE links. The x-vector sum of these links
        is 2 because there are two horizontal links, each with an x- unit
        vector value of unity.  The y-vector sum is 1 because only one of the
        three (link #1) has a non-zero y component (equal to one). Here is
        how the numbers work out:

        ..math::

            U_1 = (q_1 m_1) / 2 + (q_8 m_8) / 2 + (q_9 m_9) / 2
                = (1 0) / 2 + (1 1) / 2 + (1 1) / 2 = 1
            V_1 = (q_1 n_1) / 1 + (q_8 n_8) / 1 + (q_9 n_9) / 1
                = (1 1) / 1 + (1 0) / 1 + (1 0) / 1 = 1

        At node 5, in the interior, there are four connected links (two
        in-links and two out-links; two horizontal and two vertical). So, we
        add up the q values associated with all four:

        U_5 = (q_1 m_1) / 2 + (q_5 m_5) / 2 + (q_11 m_11) / 2 + (q_12 m_12) / 2
            = (1 0) / 2 + (1 0) / 2 + (1 1) / 2 + (1 1) / 2 = 1

        V_5 = (q_1 n_1) / 2 + (q_5 n_5) / 2 + (q_11 n_11) / 2 + (q_12 n_12) / 2
            = (1 1) / 2 + (1 1) / 2 + (1 0) / 2 + (1 0) / 2 = 1

        To do this calculation efficiently, we use the following algorithm:

        FOR each row in node_inlink_matrix (representing one inlink @ each
        node)
        * Multiply the link's q value by its unit x component ...
          ... divide by node's unit vector sum in x ...
          ... and add it to the node's total q_x
        * Multiply the link's q value by its unit y component ...
          ... divide by node's unit vector sum in y ...
          ... and add it to the node's total q_y

        Examples
        --------

        Example 1: q[:] = 1. Vector magnitude is :math:`\sqrt{2}`, direction is
        :math:`(1,1)`.

        >>> import numpy as np
        >>> import landlab as ll
        >>> rmg = ll.RasterModelGrid(3, 4, 2.0)
        >>> rmg.node_unit_vector_sum_x
        array([ 1.,  2.,  2.,  1.,  1.,  2.,  2.,  1.,  1.,  2.,  2.,  1.])
        >>> rmg.node_unit_vector_sum_y
        array([ 1.,  1.,  1.,  1.,  2.,  2.,  2.,  2.,  1.,  1.,  1.,  1.])
        >>> q = np.ones(rmg.number_of_links)
        >>> nvx, nvy = rmg.map_link_vector_to_nodes(q)
        >>> nvx
        array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])
        >>> nvy
        array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])

        Example 2: Vector magnitude is 5, angle is 30 degrees from horizontal,
        forming a 3-4-5 triangle.
        >>> q[:8] = 3.
        >>> q[8:] = 4.
        >>> nvx, nvy = rmg.map_link_vector_to_nodes(q)
        >>> nvx
        array([ 4.,  4.,  4.,  4.,  4.,  4.,  4.,  4.,  4.,  4.,  4.,  4.])
        >>> nvy
        array([ 3.,  3.,  3.,  3.,  3.,  3.,  3.,  3.,  3.,  3.,  3.,  3.])

        ..todo::

            Fix and finish example 3 below.

        Example 3: Hexagonal grid with vector as above. Here, q is
        pre-calculated to have the right values to represent a uniform
        vector with magnitude 5 and orientation 30 degrees counter-clockwise
        from horizontal.
        # >>> hmg = ll.HexModelGrid(3, 2, 2.0)
        # >>> q = np.array([4.598, 0.598, -4., -4.598, 4., -0.598, -0.598,
        # ...     -4., -4.598, -0.598, 4., -4.598])
        """

        # Create the arrays to hold the node-based values of the x and y
        # components of the vector (q)
        node_vec_x = numpy.zeros(self.number_of_nodes)
        node_vec_y = numpy.zeros(self.number_of_nodes)

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
        qx = numpy.zeros(self.number_of_links + 1)
        qy = numpy.zeros(self.number_of_links + 1)
        qx[:self.number_of_links] = q * \
            self.link_unit_vec_x[:self.number_of_links]
        qy[:self.number_of_links] = q * \
            self.link_unit_vec_y[:self.number_of_links]

        # Loop over each row in the node_inlink_matrix and node_outlink_matrix.
        # This isn't a big loop! In a raster grid, these have only two rows
        # each; in an unstructured grid, it depends on the grid geometry;
        # for a hex grid, there are up to 6 rows.
        n_matrix_rows = numpy.size(self.node_inlink_matrix, 0)
        for i in range(n_matrix_rows):
            node_vec_x += qx[self.node_inlink_matrix[i, :]]
            node_vec_x += qx[self.node_outlink_matrix[i, :]]
            node_vec_y += qy[self.node_inlink_matrix[i, :]]
            node_vec_y += qy[self.node_outlink_matrix[i, :]]
        node_vec_x /= self.node_unit_vector_sum_x
        node_vec_y /= self.node_unit_vector_sum_y

        return node_vec_x, node_vec_y

    def display_grid(self, draw_voronoi=False):
        """Display the grid."""
        import matplotlib.pyplot as plt

        # Plot nodes, colored by boundary vs interior
        plt.plot(self.node_x[self.core_nodes],
                 self.node_y[self.core_nodes], 'go')
        plt.plot(self.node_x[self.boundary_nodes],
                 self.node_y[self.boundary_nodes], 'ro')

        # Draw links
        for i in range(self.number_of_links):
            plt.plot([self.node_x[self.node_at_link_tail[i]],
                      self.node_x[self.node_at_link_head[i]]],
                     [self.node_y[self.node_at_link_tail[i]],
                      self.node_y[self.node_at_link_head[i]]], 'k-')

        # Draw active links
        for link in self.active_link_ids:
            plt.plot([self.node_x[self.node_at_link_tail[link]],
                      self.node_x[self.node_at_link_head[link]]],
                     [self.node_y[self.node_at_link_tail[link]],
                      self.node_y[self.node_at_link_head[link]]], 'g-')

        # If caller asked for a voronoi diagram, draw that too
        if draw_voronoi:
            from scipy.spatial import Voronoi, voronoi_plot_2d
            pts = numpy.zeros((self.number_of_nodes, 2))
            pts[:, 0] = self.node_x
            pts[:, 1] = self.node_y
            vor = Voronoi(pts)
            voronoi_plot_2d(vor)

        plt.show()

    def is_boundary(self, ids, boundary_flag=None):
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
        """
        if boundary_flag is None:
            return ~ (self._node_status[ids] == CORE_NODE)
        else:
            return self._node_status[ids] == boundary_flag

    def _assign_boundary_nodes_to_grid_sides(self):
        """Assign boundary nodes to a quadrant.

        For each boundary node, determines whether it belongs to the left,
        right, top or bottom of the grid, based on its distance from the grid's
        centerpoint (mean (x,y) position). Returns lists of nodes on each of
        the four grid sides. Assumes self.status_at_node, self.number_of_nodes,
        self.boundary_nodes, self._node_x, and self._node_y have been
        initialized.

        Examples
        --------
        >>> import landlab as ll
        >>> m = ll.HexModelGrid(5, 3, 1.0)
        >>> [l,r,t,b] = m._assign_boundary_nodes_to_grid_sides()
        >>> l
        array([ 7, 12,  3])
        >>> r
        array([11, 15,  6])
        >>> t
        array([16, 18, 17])
        >>> b
        array([0, 2, 1])
        """
        # Calculate x and y distance from centerpoint
        diff_x = self.node_x[self.boundary_nodes] - numpy.mean(self.node_x)
        diff_y = self.node_y[self.boundary_nodes] - numpy.mean(self.node_y)

        return _sort_points_into_quadrants(diff_x, diff_y, self.boundary_nodes)

    @deprecated
    def set_inactive_boundaries(self, bottom_is_inactive, right_is_inactive,
                                top_is_inactive, left_is_inactive):
        """Set boundaries to inactive.

        .. note:: Deprecated since version 0.6.
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

        Examples
        --------
        >>> import landlab as ll
        >>> rmg = ll.HexModelGrid(5, 3, 1.0) # rows, columns, spacing
        >>> rmg.number_of_active_links
        30
        >>> rmg.status_at_node # doctest: +NORMALIZE_WHITESPACE
        array([1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1],
              dtype=int8)
        >>> rmg.set_inactive_boundaries(False, False, True, True)
        >>> rmg.number_of_active_links
        21
        >>> rmg.status_at_node # doctest: +NORMALIZE_WHITESPACE
        array([1, 1, 1, 4, 0, 0, 1, 4, 0, 0, 0, 1, 4, 0, 0, 1, 4, 4, 4],
               dtype=int8)
        """
        if self._DEBUG_TRACK_METHODS:
            six.print_('ModelGrid.set_inactive_boundaries')

        [left_edge, right_edge, top_edge, bottom_edge] = \
            self._assign_boundary_nodes_to_grid_sides()

        if bottom_is_inactive:
            self._node_status[bottom_edge] = CLOSED_BOUNDARY
        else:
            self._node_status[bottom_edge] = FIXED_VALUE_BOUNDARY

        if right_is_inactive:
            self._node_status[right_edge] = CLOSED_BOUNDARY
        else:
            self._node_status[right_edge] = FIXED_VALUE_BOUNDARY

        if top_is_inactive:
            self._node_status[top_edge] = CLOSED_BOUNDARY
        else:
            self._node_status[top_edge] = FIXED_VALUE_BOUNDARY

        if left_is_inactive:
            self._node_status[left_edge] = CLOSED_BOUNDARY
        else:
            self._node_status[left_edge] = FIXED_VALUE_BOUNDARY

        self.update_links_nodes_cells_to_new_BCs()

    def set_closed_nodes(self, nodes):
        """Make nodes closed boundaries.

        Sets the given nodes' boundary condition statuses to CLOSED (==4),
        and resets the list of active links to reflect any changes.
        """
        self._node_status[nodes] = CLOSED_BOUNDARY
        self.update_links_nodes_cells_to_new_BCs()

    def get_distances_of_nodes_to_point(self, coord, get_az=None,
                                        node_subset=None,
                                        out_distance=None, out_azimuth=None):
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

        >>> grid.get_distances_of_nodes_to_point((2, 1),
        ...     node_subset=(2, 6, 7, 8, 12))
        array([ 1.,  1.,  0.,  1.,  1.])

        Calculate distances from a point to all nodes on the grid.

        >>> dist = grid.get_distances_of_nodes_to_point((2, 1))
        >>> dist.shape == (grid.number_of_nodes, )
        True
        >>> dist.take((2, 6, 7, 8, 12))
        array([ 1.,  1.,  0.,  1.,  1.])

        Put the distances into a buffer.

        >>> out = np.empty(grid.number_of_nodes, dtype=float)
        >>> dist = grid.get_distances_of_nodes_to_point((2, 1),
        ...     out_distance=out)
        >>> out is dist
        True
        >>> out.take((2, 6, 7, 8, 12))
        array([ 1.,  1.,  0.,  1.,  1.])

        Calculate azimuths along with distances. The azimuths are calculated
        in radians but measured clockwise from north.

        >>> (_, azim) = grid.get_distances_of_nodes_to_point((2, 1),
        ...     get_az='angles')
        >>> azim.take((2, 6, 7, 8, 12)) * 180. / np.pi
        array([ 180.,  270.,    0.,   90.,    0.])
        >>> (_, azim) = grid.get_distances_of_nodes_to_point((2, 1),
        ...     get_az='angles', node_subset=(1, 3, 11, 13))
        >>> azim * 180. / np.pi
        array([ 225.,  135.,  315.,   45.])

        When calculating displacements, the first row contains displacements
        in x and the second displacements in y.

        >>> (_, azim) = grid.get_distances_of_nodes_to_point((2, 1),
        ...     get_az='displacements', node_subset=(2, 6, 7, 8, 12))
        >>> azim
        array([[ 0., -1.,  0.,  1.,  0.],
               [-1.,  0.,  0.,  0.,  1.]])
        """
        if len(coord) != 2:
            raise ValueError('coordinate must iterable of length 2')

        if get_az not in (None, 'displacements', 'angles'):
            raise ValueError('get_az not understood')

        if node_subset is not None and numpy.any(numpy.isnan(node_subset)):
            node_subset = None

        if node_subset is not None:
            if not isinstance(node_subset, numpy.ndarray):
                node_subset = numpy.array(node_subset)
            node_subset = node_subset.reshape((-1, ))
            len_subset = node_subset.size
        else:
            len_subset = self.number_of_nodes

        if out_distance is None:
            out_distance = numpy.empty(len_subset, dtype=numpy.float)
        if out_distance.size != len_subset:
            raise ValueError('output array size mismatch for distances')

        if get_az is not None:
            if get_az == 'displacements':
                az_shape = (2, len_subset)
            else:
                az_shape = (len_subset, )
            if out_azimuth is None:
                out_azimuth = numpy.empty(az_shape, dtype=numpy.float)
            if out_azimuth.shape != az_shape:
                raise ValueError('output array mismatch for azimuths')

        azimuths_as_displacements = numpy.empty((2, self.number_of_nodes))
        dummy_nodes_1 = numpy.empty(self.number_of_nodes)
        dummy_nodes_2 = numpy.empty(self.number_of_nodes)
        dummy_nodes_3 = numpy.empty(self.number_of_nodes)

        if node_subset is None:
            azimuths_as_displacements[0] = (self.node_x - coord[0])
            azimuths_as_displacements[1] = (self.node_y - coord[1])
        else:
            azimuths_as_displacements[0, :len_subset] = (
                self.node_x[node_subset] - coord[0])
            azimuths_as_displacements[1, :len_subset] = (
                self.node_y[node_subset] - coord[1])

        numpy.square(azimuths_as_displacements[0, :len_subset],
                     out=dummy_nodes_1[:len_subset])
        numpy.square(azimuths_as_displacements[1, :len_subset],
                     out=dummy_nodes_2[:len_subset])
        numpy.add(dummy_nodes_1[:len_subset], dummy_nodes_2[:len_subset],
                  out=dummy_nodes_3[:len_subset])
        numpy.sqrt(dummy_nodes_3[:len_subset], out=out_distance)

        if get_az:
            if get_az == 'displacements':
                out_azimuth[:] = azimuths_as_displacements[:, :len_subset]
            elif get_az == 'angles':
                numpy.arctan2(
                    azimuths_as_displacements[0, :len_subset],
                    azimuths_as_displacements[1, :len_subset],
                    out=out_azimuth[:len_subset])

                less_than_zero = numpy.empty(self.number_of_nodes, dtype=bool)
                numpy.less(out_azimuth, 0., out=less_than_zero[:len_subset])
                out_azimuth[less_than_zero[:len_subset]] += 2. * numpy.pi

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
        """
        if self._all_node_distances_map is None:
            self.build_all_node_distances_azimuths_maps()
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
        """
        if self._all_node_azimuths_map is None:
            self.build_all_node_distances_azimuths_maps()
        return self._all_node_azimuths_map

    def build_all_node_distances_azimuths_maps(self):
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
        self._all_node_distances_map = numpy.empty((self.number_of_nodes,
                                                    self.number_of_nodes))
        self._all_node_azimuths_map = numpy.empty((self.number_of_nodes,
                                                   self.number_of_nodes))

        node_coords = numpy.empty((self.number_of_nodes, 2))
        node_coords[:, 0] = self.node_x
        node_coords[:, 1] = self.node_y

        for i in range(self.number_of_nodes):
            (self._all_node_distances_map[i, :],
             self._all_node_azimuths_map[i, :]) = (
                 self.get_distances_of_nodes_to_point(
                     (node_coords[i, 0], node_coords[i, 1]), get_az='angles'))

        assert numpy.all(self._all_node_distances_map >= 0.)

        return self._all_node_distances_map, self._all_node_azimuths_map


add_module_functions_to_class(ModelGrid, 'mappers.py', pattern='map_*')
add_module_functions_to_class(ModelGrid, 'gradients.py', pattern='calculate_*')


if __name__ == '__main__':
    import doctest
    doctest.testmod()
