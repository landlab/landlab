#! /usr/env/python
"""
Python implementation of ModelGrid, a base class used to create and manage
grids for 2D numerical models.

Getting Information about a Grid
--------------------------------
The following attributes, properties, and methods provide data about the grid,
its geometry, and the connectivity among the various elements. Each grid
element has an ID number, which is also its position in an array that
contains information about that type of element. For example, the *x*
coordinate of node 5 would be found at `grid.node_x[5]`.

The naming of grid-element arrays is *attribute*`_at_`*element*, where
*attribute* is the name of the data in question, and *element* is the element
to which the attribute applies. For example, the property `node_at_cell`
contains the ID of the node associated with each cell. For example,
`node_at_cell[3]` contains the *node ID* of the node associated with cell 3.
The *attribute* is singular if there is only one value per element; for
example, there is only one node associated with each cell. It is plural when
there are multiple values per element; for example, the `faces_at_cell` array
contains multiple faces for each cell. Exceptions to these general rules are
functions that return indices of a subset of all elements of a particular type.
For example, you can obtain an array with IDs of only the core nodes using
`core_nodes`, while `active_links` provides an array of IDs of active links
(only). Finally, attributes that represent a measurement of something, such as
the length of a link or the surface area of a cell, are described using `_of_`,
as in the example `area_of_cell`.

Information about the grid as a whole
+++++++++++++++++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.hex.HexModelGrid.axis_name
    ~landlab.grid.hex.HexModelGrid.axis_units
    ~landlab.grid.hex.HexModelGrid.from_dict
    ~landlab.grid.hex.HexModelGrid.hexplot
    ~landlab.grid.hex.HexModelGrid.move_origin
    ~landlab.grid.hex.HexModelGrid.node_axis_coordinates
    ~landlab.grid.hex.HexModelGrid.number_of_active_faces
    ~landlab.grid.hex.HexModelGrid.number_of_active_links
    ~landlab.grid.hex.HexModelGrid.number_of_cells
    ~landlab.grid.hex.HexModelGrid.number_of_core_cells
    ~landlab.grid.hex.HexModelGrid.number_of_core_nodes
    ~landlab.grid.hex.HexModelGrid.number_of_elements
    ~landlab.grid.hex.HexModelGrid.number_of_faces
    ~landlab.grid.hex.HexModelGrid.number_of_fixed_links
    ~landlab.grid.hex.HexModelGrid.number_of_links
    ~landlab.grid.hex.HexModelGrid.number_of_node_columns
    ~landlab.grid.hex.HexModelGrid.number_of_node_rows
    ~landlab.grid.hex.HexModelGrid.number_of_nodes
    ~landlab.grid.hex.HexModelGrid.number_of_patches
    ~landlab.grid.hex.HexModelGrid.save

Information about nodes
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.hex.HexModelGrid.active_link_dirs_at_node
    ~landlab.grid.hex.HexModelGrid.all_node_azimuths_map
    ~landlab.grid.hex.HexModelGrid.all_node_distances_map
    ~landlab.grid.hex.HexModelGrid.boundary_nodes
    ~landlab.grid.hex.HexModelGrid.cell_area_at_node
    ~landlab.grid.hex.HexModelGrid.cell_at_node
    ~landlab.grid.hex.HexModelGrid.closed_boundary_nodes
    ~landlab.grid.hex.HexModelGrid.core_nodes
    ~landlab.grid.hex.HexModelGrid.downwind_links_at_node
    ~landlab.grid.hex.HexModelGrid.fixed_gradient_boundary_nodes
    ~landlab.grid.hex.HexModelGrid.fixed_value_boundary_nodes
    ~landlab.grid.hex.HexModelGrid.link_at_node_is_downwind
    ~landlab.grid.hex.HexModelGrid.link_at_node_is_upwind
    ~landlab.grid.hex.HexModelGrid.link_dirs_at_node
    ~landlab.grid.hex.HexModelGrid.links_at_node
    ~landlab.grid.hex.HexModelGrid.neighbors_at_node
    ~landlab.grid.hex.HexModelGrid.node_at_cell
    ~landlab.grid.hex.HexModelGrid.node_at_core_cell
    ~landlab.grid.hex.HexModelGrid.node_at_link_head
    ~landlab.grid.hex.HexModelGrid.node_at_link_tail
    ~landlab.grid.hex.HexModelGrid.node_axis_coordinates
    ~landlab.grid.hex.HexModelGrid.node_is_boundary
    ~landlab.grid.hex.HexModelGrid.node_x
    ~landlab.grid.hex.HexModelGrid.node_y
    ~landlab.grid.hex.HexModelGrid.nodes
    ~landlab.grid.hex.HexModelGrid.nodes_at_patch
    ~landlab.grid.hex.HexModelGrid.number_of_core_nodes
    ~landlab.grid.hex.HexModelGrid.number_of_links_at_node
    ~landlab.grid.hex.HexModelGrid.number_of_node_columns
    ~landlab.grid.hex.HexModelGrid.number_of_node_rows
    ~landlab.grid.hex.HexModelGrid.number_of_nodes
    ~landlab.grid.hex.HexModelGrid.open_boundary_nodes
    ~landlab.grid.hex.HexModelGrid.patches_at_node
    ~landlab.grid.hex.HexModelGrid.set_nodata_nodes_to_closed
    ~landlab.grid.hex.HexModelGrid.set_nodata_nodes_to_fixed_gradient
    ~landlab.grid.hex.HexModelGrid.status_at_node
    ~landlab.grid.hex.HexModelGrid.unit_vector_sum_xcomponent_at_node
    ~landlab.grid.hex.HexModelGrid.unit_vector_sum_ycomponent_at_node
    ~landlab.grid.hex.HexModelGrid.upwind_links_at_node
    ~landlab.grid.hex.HexModelGrid.x_of_node
    ~landlab.grid.hex.HexModelGrid.y_of_node

Information about links
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.hex.HexModelGrid.active_link_dirs_at_node
    ~landlab.grid.hex.HexModelGrid.active_links
    ~landlab.grid.hex.HexModelGrid.angle_of_link
    ~landlab.grid.hex.HexModelGrid.downwind_links_at_node
    ~landlab.grid.hex.HexModelGrid.face_at_link
    ~landlab.grid.hex.HexModelGrid.fixed_links
    ~landlab.grid.hex.HexModelGrid.length_of_link
    ~landlab.grid.hex.HexModelGrid.link_at_face
    ~landlab.grid.hex.HexModelGrid.link_at_node_is_downwind
    ~landlab.grid.hex.HexModelGrid.link_at_node_is_upwind
    ~landlab.grid.hex.HexModelGrid.link_dirs_at_node
    ~landlab.grid.hex.HexModelGrid.links_at_node
    ~landlab.grid.hex.HexModelGrid.node_at_link_head
    ~landlab.grid.hex.HexModelGrid.node_at_link_tail
    ~landlab.grid.hex.HexModelGrid.number_of_active_links
    ~landlab.grid.hex.HexModelGrid.number_of_fixed_links
    ~landlab.grid.hex.HexModelGrid.number_of_links
    ~landlab.grid.hex.HexModelGrid.number_of_links_at_node
    ~landlab.grid.hex.HexModelGrid.resolve_values_on_active_links
    ~landlab.grid.hex.HexModelGrid.resolve_values_on_links
    ~landlab.grid.hex.HexModelGrid.status_at_link
    ~landlab.grid.hex.HexModelGrid.unit_vector_xcomponent_at_link
    ~landlab.grid.hex.HexModelGrid.unit_vector_ycomponent_at_link
    ~landlab.grid.hex.HexModelGrid.upwind_links_at_node

Information about cells
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.hex.HexModelGrid.area_of_cell
    ~landlab.grid.hex.HexModelGrid.cell_area_at_node
    ~landlab.grid.hex.HexModelGrid.cell_at_node
    ~landlab.grid.hex.HexModelGrid.core_cells
    ~landlab.grid.hex.HexModelGrid.faces_at_cell
    ~landlab.grid.hex.HexModelGrid.node_at_cell
    ~landlab.grid.hex.HexModelGrid.node_at_core_cell
    ~landlab.grid.hex.HexModelGrid.number_of_cells
    ~landlab.grid.hex.HexModelGrid.number_of_core_cells
    ~landlab.grid.hex.HexModelGrid.number_of_faces_at_cell

Information about faces
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.hex.HexModelGrid.active_faces
    ~landlab.grid.hex.HexModelGrid.face_at_link
    ~landlab.grid.hex.HexModelGrid.faces_at_cell
    ~landlab.grid.hex.HexModelGrid.link_at_face
    ~landlab.grid.hex.HexModelGrid.number_of_active_faces
    ~landlab.grid.hex.HexModelGrid.number_of_faces
    ~landlab.grid.hex.HexModelGrid.number_of_faces_at_cell
    ~landlab.grid.hex.HexModelGrid.width_of_face

Information about patches
+++++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.hex.HexModelGrid.nodes_at_patch
    ~landlab.grid.hex.HexModelGrid.number_of_patches
    ~landlab.grid.hex.HexModelGrid.patches_at_node

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

    ~landlab.grid.hex.HexModelGrid.at_node
    ~landlab.grid.hex.HexModelGrid.at_cell
    ~landlab.grid.hex.HexModelGrid.at_link
    ~landlab.grid.hex.HexModelGrid.at_face

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

    ~landlab.grid.hex.HexModelGrid.add_empty
    ~landlab.grid.hex.HexModelGrid.add_field
    ~landlab.grid.hex.HexModelGrid.add_ones
    ~landlab.grid.hex.HexModelGrid.add_zeros
    ~landlab.grid.hex.HexModelGrid.delete_field
    ~landlab.grid.hex.HexModelGrid.set_units

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
    ~landlab.grid.hex.HexModelGrid.field_units
    ~landlab.grid.hex.HexModelGrid.field_values
    ~landlab.field.grouped.ModelDataFields.groups

    i.e., call, e.g. mg.has_field('node', 'my_field_name')

    # START HERE check that all functions listed below are included above,
    # ignore ones that start with underscores(_)

Gradients, fluxes, and divergences on the grid
----------------------------------------------

Landlab is designed to easily calculate gradients in quantities across the
grid, and to construct fluxes and flux divergences from them. Because these
calculations tend to be a little more involved than property lookups, the
methods tend to start with `calc_`.

.. autosummary::
    :toctree: generated/

    ~landlab.grid.hex.HexModelGrid.calc_diff_at_link
    ~landlab.grid.hex.HexModelGrid.calc_flux_div_at_node
    ~landlab.grid.hex.HexModelGrid.calc_grad_at_link
    ~landlab.grid.hex.HexModelGrid.calc_grad_at_patch
    ~landlab.grid.hex.HexModelGrid.calc_net_flux_at_node
    ~landlab.grid.hex.HexModelGrid.calc_slope_at_node
    ~landlab.grid.hex.HexModelGrid.calc_slope_at_patch
    ~landlab.grid.hex.HexModelGrid.calc_unit_normal_at_patch

Mappers
-------

These methods allow mapping of values defined on one grid element type onto a
second, e.g., mapping upwind node values onto links, or mean link values onto
nodes.

    ~landlab.grid.hex.HexModelGrid.map_downwind_node_link_max_to_node
    ~landlab.grid.hex.HexModelGrid.map_downwind_node_link_mean_to_node
    ~landlab.grid.hex.HexModelGrid.map_link_head_node_to_link
    ~landlab.grid.hex.HexModelGrid.map_link_tail_node_to_link
    ~landlab.grid.hex.HexModelGrid.map_link_vector_to_nodes
    ~landlab.grid.hex.HexModelGrid.map_max_of_link_nodes_to_link
    ~landlab.grid.hex.HexModelGrid.map_max_of_node_links_to_node
    ~landlab.grid.hex.HexModelGrid.map_mean_of_link_nodes_to_link
    ~landlab.grid.hex.HexModelGrid.map_min_of_link_nodes_to_link
    ~landlab.grid.hex.HexModelGrid.map_min_of_node_links_to_node
    ~landlab.grid.hex.HexModelGrid.map_node_to_cell
    ~landlab.grid.hex.HexModelGrid.map_upwind_node_link_max_to_node
    ~landlab.grid.hex.HexModelGrid.map_upwind_node_link_mean_to_node
    ~landlab.grid.hex.HexModelGrid.map_value_at_downwind_node_link_max_to_node
    ~landlab.grid.hex.HexModelGrid.map_value_at_max_node_to_link
    ~landlab.grid.hex.HexModelGrid.map_value_at_min_node_to_link
    ~landlab.grid.hex.HexModelGrid.map_value_at_upwind_node_link_max_to_node

Boundary condition control
--------------------------

These are the primary properties for getting and setting the grid boundary
conditions. Changes made to :meth:`~.ModelGrid.status_at_node` and
:meth:`~.ModelGrid.status_at_node` will automatically update the conditions
defined at other grid elements automatically.

.. autosummary::
    :toctree: generated/

    ~landlab.grid.hex.HexModelGrid.active_faces
    ~landlab.grid.hex.HexModelGrid.active_links
    ~landlab.grid.hex.HexModelGrid.boundary_nodes
    ~landlab.grid.hex.HexModelGrid.closed_boundary_nodes
    ~landlab.grid.hex.HexModelGrid.core_cells
    ~landlab.grid.hex.HexModelGrid.core_nodes
    ~landlab.grid.hex.HexModelGrid.fixed_gradient_boundary_nodes
    ~landlab.grid.hex.HexModelGrid.fixed_links
    ~landlab.grid.hex.HexModelGrid.fixed_value_boundary_nodes
    ~landlab.grid.hex.HexModelGrid.node_at_core_cell
    ~landlab.grid.hex.HexModelGrid.node_is_boundary
    ~landlab.grid.hex.HexModelGrid.number_of_active_faces
    ~landlab.grid.hex.HexModelGrid.number_of_active_links
    ~landlab.grid.hex.HexModelGrid.number_of_core_cells
    ~landlab.grid.hex.HexModelGrid.number_of_core_nodes
    ~landlab.grid.hex.HexModelGrid.number_of_fixed_links
    ~landlab.grid.hex.HexModelGrid.open_boundary_nodes
    ~landlab.grid.hex.HexModelGrid.set_nodata_nodes_to_closed
    ~landlab.grid.hex.HexModelGrid.set_nodata_nodes_to_fixed_gradient
    ~landlab.grid.hex.HexModelGrid.status_at_link
    ~landlab.grid.hex.HexModelGrid.status_at_node

Identifying node subsets
------------------------

These methods are useful in identifying subsets of nodes, e.g., closest node
to a point; nodes at edges.

(None are available for this grid type)

Surface analysis
----------------

These methods permit the kinds of surface analysis that you might expect to
find in GIS software.

.. autosummary::
    :toctree: generated/

    ~landlab.grid.hex.HexModelGrid.calc_aspect_at_node
    ~landlab.grid.hex.HexModelGrid.calc_hillshade_at_node
    ~landlab.grid.hex.HexModelGrid.calc_slope_at_node

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
"""

import numpy
import six

from landlab.grid.voronoi import VoronoiDelaunayGrid


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

    def __init__(self, base_num_rows=0, base_num_cols=0, dx=1.0,
                 orientation='horizontal', shape='hex', reorient_links=True,
                 **kwds):
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
        # Set number of nodes, and initialize if caller has given dimensions
        if base_num_rows * base_num_cols > 0:
            self._initialize(base_num_rows, base_num_cols, dx, orientation,
                             shape, reorient_links)
        super(HexModelGrid, self).__init__(**kwds)

    @classmethod
    def from_dict(cls, params):
        shape = params['shape']
        spacing = params.get('spacing', 1.)

        return cls(shape[0], shape[1], spacing)

    def _initialize(self, base_num_rows, base_num_cols, dx, orientation,
                    shape, reorient_links=True):
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
        orientation : string
            Either 'horizontal' (default in __init__) or 'vertical'
        shape : string
            Either 'hex' (default in __init__) or 'rect'
        reorient_links : bool
            Whether or not to re-orient all links to point between -45 deg
            and +135 deg clockwise from "north" (i.e., along y axis)

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
            six.print_('HexModelGrid._initialize(' + str(base_num_rows) +
                       ', ' + str(base_num_cols) + ', ' + str(dx) + ')')

        # Make sure the parameter *orientation* is correct
        assert (orientation[0].lower() == 'h' or
                orientation[0].lower() == 'v'), \
            'orientation must be either "horizontal" (default) or "vertical"'

        # Make sure the parameter *shape* is correct
        assert (shape[0].lower() == 'h' or shape[0].lower() == 'r'), \
            'shape must be either "hex" (default) or "rect"'

        # Create a set of hexagonally arranged points. These will be our nodes.
        if orientation == 'horizontal' and shape == 'hex':
            pts = HexModelGrid._hex_points_with_horizontal_hex(
                base_num_rows, base_num_cols, dx)
            self.orientation = 'horizontal'
            self._nrows = base_num_rows
        elif orientation == 'horizontal' and shape == 'rect':
            pts = HexModelGrid._hex_points_with_horizontal_rect(
                base_num_rows, base_num_cols, dx)
            self.orientation = 'horizontal'
            self._nrows = base_num_rows
            self._ncols = base_num_cols
        elif orientation == 'vertical' and shape == 'hex':
            pts = HexModelGrid._hex_points_with_vertical_hex(
                base_num_rows, base_num_cols, dx)
            self.orientation = 'vertical'
            self._ncols = base_num_cols
        else:
            pts = HexModelGrid._hex_points_with_vertical_rect(
                base_num_rows, base_num_cols, dx)
            self.orientation = 'vertical'
            self._nrows = base_num_rows
            self._ncols = base_num_cols

        # Call the VoronoiDelaunayGrid constructor to triangulate/Voronoi
        # the nodes into a grid.
        super(HexModelGrid, self)._initialize(
            pts[:, 0], pts[:, 1], reorient_links)

        # Remember grid spacing
        self._dx = dx

    def _create_cell_areas_array(self):
        r"""Create an array of surface areas of hexagonal cells.

        Creates and returns an array containing the surface areas of the
        hexagonal (Voronoi) cells.

        These cells are perfect hexagons in which the apothem is dx/2. The
        formula for area is:

        .. math::
            A = 3 dx^2 / 2 \sqrt{3} \approx 0.866 dx^2
        """
        self._area_of_cell = (0.8660254 * self._dx ** 2 +
                              numpy.zeros(self.number_of_cells))
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
            npts = num_rows * base_num_cols + \
                ((num_rows - 1) // 2) * ((num_rows - 1) // 2)
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
            xshift = - half_dxh * extra_cols

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
            npts = base_num_rows * num_cols + \
                ((num_cols - 1) // 2) * ((num_cols - 1) // 2)
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
            yshift = - half_dxv * extra_rows

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
        >>> grid = HexModelGrid(5, 5, shape='rect')
        >>> grid.number_of_node_rows
        5
        """
        return self._nrows

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
        if self.orientation == 'horizontal':   # horizontal
            offsets[:, 0] = array(
                [0., apothem, apothem, 0., -apothem, -apothem])
            offsets[:, 1] = array(
                [radius, radius / 2.0, -radius / 2.0, -radius, -radius / 2.0,
                 radius / 2.0])
        else:   # vertical
            offsets[:, 0] = array(
                [radius / 2.0, radius, radius / 2.0, -radius / 2.0, -radius,
                 -radius / 2.0])
            offsets[:, 1] = array(
                [apothem, 0., -apothem, -apothem, 0., apothem])

        patches = []
        for i in range(self.number_of_nodes):
            poly_verts[:, 0] = self.node_x[i] + offsets[:, 0]
            poly_verts[:, 1] = self.node_y[i] + offsets[:, 1]
            p = Polygon(poly_verts, True)
            patches.append(p)

        self._hexplot_pc = PatchCollection(
            patches, cmap=color_map, edgecolor='none', linewidth=0.0)

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
        """
        from numpy import array, amin, amax
        import matplotlib.pyplot as plt
        import copy

        try:
            self._hexplot_configured is True
        except:
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


def from_dict(param_dict):
    """
    Create a HexModelGrid from the dictionary-like object, *param_dict*.
    Required keys of the dictionary are NUM_ROWS, NUM_COLS. Raises a KeyError
    if either of these are missing.  If GRID_SPACING is given, use it as the
    HexModelGrid *dx* parameter, otherwise default to unit spacing.
    """
    # Read and create a basic HexModelGrid
    try:
        n_rows = int(param_dict['NUM_ROWS'])
        n_cols = int(param_dict['NUM_COLS'])
        dx = float(param_dict.get('GRID_SPACING', 1.))
    except KeyError as e:
        raise
    except ValueError as e:
        raise
    else:
        hg = HexModelGrid(n_rows, n_cols, dx)

    return hg
