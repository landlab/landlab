General class methods and attributes of the `landlab.grid.base` module
----------------------------------------------------------------------

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

    ~landlab.grid.base.ModelGrid.axis_name
    ~landlab.grid.base.ModelGrid.axis_units
    ~landlab.grid.base.ModelGrid.move_origin
    ~landlab.grid.base.ModelGrid.ndim
    ~landlab.grid.base.ModelGrid.node_axis_coordinates
    ~landlab.grid.base.ModelGrid.number_of_elements
    ~landlab.grid.base.ModelGrid.size

Information about nodes
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.base.ModelGrid.active_link_dirs_at_node
    ~landlab.grid.base.ModelGrid.active_neighbors_at_node
    ~landlab.grid.base.ModelGrid.all_node_azimuths_map
    ~landlab.grid.base.ModelGrid.all_node_distances_map
    ~landlab.grid.base.ModelGrid.boundary_nodes
    ~landlab.grid.base.ModelGrid.calc_distances_of_nodes_to_point
    ~landlab.grid.base.ModelGrid.cell_area_at_node
    ~landlab.grid.base.ModelGrid.cell_at_node
    ~landlab.grid.base.ModelGrid.closed_boundary_nodes
    ~landlab.grid.base.ModelGrid.core_nodes
    ~landlab.grid.base.ModelGrid.downwind_links_at_node
    ~landlab.grid.base.ModelGrid.fixed_gradient_boundary_nodes
    ~landlab.grid.base.ModelGrid.fixed_value_boundary_nodes
    ~landlab.grid.base.ModelGrid.link_at_node_is_downwind
    ~landlab.grid.base.ModelGrid.link_at_node_is_upwind
    ~landlab.grid.base.ModelGrid.link_dirs_at_node
    ~landlab.grid.base.ModelGrid.links_at_node
    ~landlab.grid.base.ModelGrid.neighbors_at_node
    ~landlab.grid.base.ModelGrid.node_at_cell
    ~landlab.grid.base.ModelGrid.node_at_core_cell
    ~landlab.grid.base.ModelGrid.node_at_link_head
    ~landlab.grid.base.ModelGrid.node_at_link_tail
    ~landlab.grid.base.ModelGrid.node_axis_coordinates
    ~landlab.grid.base.ModelGrid.node_is_boundary
    ~landlab.grid.base.ModelGrid.node_x
    ~landlab.grid.base.ModelGrid.node_y
    ~landlab.grid.base.ModelGrid.nodes
    ~landlab.grid.base.ModelGrid.number_of_core_nodes
    ~landlab.grid.base.ModelGrid.number_of_links_at_node
    ~landlab.grid.base.ModelGrid.number_of_nodes
    ~landlab.grid.base.ModelGrid.number_of_patches_present_at_node
    ~landlab.grid.base.ModelGrid.open_boundary_nodes
    ~landlab.grid.base.ModelGrid.patches_present_at_node
    ~landlab.grid.base.ModelGrid.set_nodata_nodes_to_closed
    ~landlab.grid.base.ModelGrid.set_nodata_nodes_to_fixed_gradient
    ~landlab.grid.base.ModelGrid.status_at_node
    ~landlab.grid.base.ModelGrid.unit_vector_sum_xcomponent_at_node
    ~landlab.grid.base.ModelGrid.unit_vector_sum_ycomponent_at_node
    ~landlab.grid.base.ModelGrid.upwind_links_at_node
    ~landlab.grid.base.ModelGrid.x_of_node
    ~landlab.grid.base.ModelGrid.y_of_node

Information about links
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.base.ModelGrid.active_link_dirs_at_node
    ~landlab.grid.base.ModelGrid.active_links
    ~landlab.grid.base.ModelGrid.angle_of_link
    ~landlab.grid.base.ModelGrid.angle_of_link_about_head
    ~landlab.grid.base.ModelGrid.downwind_links_at_node
    ~landlab.grid.base.ModelGrid.face_at_link
    ~landlab.grid.base.ModelGrid.fixed_links
    ~landlab.grid.base.ModelGrid.length_of_link
    ~landlab.grid.base.ModelGrid.link_at_face
    ~landlab.grid.base.ModelGrid.link_at_node_is_downwind
    ~landlab.grid.base.ModelGrid.link_at_node_is_upwind
    ~landlab.grid.base.ModelGrid.link_dirs_at_node
    ~landlab.grid.base.ModelGrid.links_at_node
    ~landlab.grid.base.ModelGrid.node_at_link_head
    ~landlab.grid.base.ModelGrid.node_at_link_tail
    ~landlab.grid.base.ModelGrid.number_of_active_links
    ~landlab.grid.base.ModelGrid.number_of_fixed_links
    ~landlab.grid.base.ModelGrid.number_of_links
    ~landlab.grid.base.ModelGrid.number_of_links_at_node
    ~landlab.grid.base.ModelGrid.number_of_patches_present_at_link
    ~landlab.grid.base.ModelGrid.patches_present_at_link
    ~landlab.grid.base.ModelGrid.resolve_values_on_active_links
    ~landlab.grid.base.ModelGrid.resolve_values_on_links
    ~landlab.grid.base.ModelGrid.status_at_link
    ~landlab.grid.base.ModelGrid.unit_vector_xcomponent_at_link
    ~landlab.grid.base.ModelGrid.unit_vector_ycomponent_at_link
    ~landlab.grid.base.ModelGrid.upwind_links_at_node
    ~landlab.grid.base.ModelGrid.x_of_link
    ~landlab.grid.base.ModelGrid.y_of_link

Information about cells
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.base.ModelGrid.area_of_cell
    ~landlab.grid.base.ModelGrid.cell_area_at_node
    ~landlab.grid.base.ModelGrid.cell_at_node
    ~landlab.grid.base.ModelGrid.core_cells
    ~landlab.grid.base.ModelGrid.faces_at_cell
    ~landlab.grid.base.ModelGrid.node_at_cell
    ~landlab.grid.base.ModelGrid.node_at_core_cell
    ~landlab.grid.base.ModelGrid.number_of_cells
    ~landlab.grid.base.ModelGrid.number_of_core_cells
    ~landlab.grid.base.ModelGrid.number_of_faces_at_cell
    ~landlab.grid.base.ModelGrid.x_of_cell
    ~landlab.grid.base.ModelGrid.y_of_cell

Information about faces
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.base.ModelGrid.active_faces
    ~landlab.grid.base.ModelGrid.face_at_link
    ~landlab.grid.base.ModelGrid.faces_at_cell
    ~landlab.grid.base.ModelGrid.link_at_face
    ~landlab.grid.base.ModelGrid.number_of_active_faces
    ~landlab.grid.base.ModelGrid.number_of_faces
    ~landlab.grid.base.ModelGrid.number_of_faces_at_cell
    ~landlab.grid.base.ModelGrid.width_of_face
    ~landlab.grid.base.ModelGrid.x_of_face
    ~landlab.grid.base.ModelGrid.y_of_face

Information about patches
+++++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.base.ModelGrid.number_of_patches_present_at_link
    ~landlab.grid.base.ModelGrid.number_of_patches_present_at_node
    ~landlab.grid.base.ModelGrid.patches_present_at_link
    ~landlab.grid.base.ModelGrid.patches_present_at_node

Information about corners
+++++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.base.ModelGrid.number_of_corners


Data Fields in ModelGrid
------------------------
:class:`~.ModelGrid` inherits from the :class:`~.ModelDataFields` class. This
provides `~.ModelGrid`, and its subclasses, with the ability to, optionally,
store data values that are associated with the different types grid elements
(nodes, cells, etc.). In particular, as part of ``ModelGrid.__init__()``,
data field *groups* are added to the `ModelGrid` that provide containers to
put data fields into. There is one group for each of the eight grid elements
(node, cell, link, face, core_node, core_cell, active_link, and active_face).
There is an additional group at_grid that can store arrays of length one
intended as a place to store varibles global to the grid.

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
    ~landlab.grid.base.ModelGrid.at_patch
    ~landlab.grid.base.ModelGrid.at_corner
    ~landlab.grid.base.ModelGrid.at_grid

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

    ~landlab.grid.base.ModelGrid.add_empty
    ~landlab.grid.base.ModelGrid.add_field
    ~landlab.grid.base.ModelGrid.add_ones
    ~landlab.grid.base.ModelGrid.add_zeros
    ~landlab.grid.base.ModelGrid.delete_field
    ~landlab.grid.base.ModelGrid.set_units

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
    ~landlab.grid.base.ModelGrid.field_units
    ~landlab.grid.base.ModelGrid.field_values
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

    ~landlab.grid.base.ModelGrid.calc_diff_at_link
    ~landlab.grid.base.ModelGrid.calc_flux_div_at_node
    ~landlab.grid.base.ModelGrid.calc_grad_at_link
    ~landlab.grid.base.ModelGrid.calc_grad_at_patch
    ~landlab.grid.base.ModelGrid.calc_net_flux_at_node
    ~landlab.grid.base.ModelGrid.calc_slope_at_node
    ~landlab.grid.base.ModelGrid.calc_slope_at_patch
    ~landlab.grid.base.ModelGrid.calc_unit_normal_at_patch

Mappers
-------

These methods allow mapping of values defined on one grid element type onto a
second, e.g., mapping upwind node values onto links, or mean link values onto
nodes.

.. autosummary::
    :toctree: generated/

    ~landlab.grid.base.ModelGrid.map_downwind_node_link_max_to_node
    ~landlab.grid.base.ModelGrid.map_downwind_node_link_mean_to_node
    ~landlab.grid.base.ModelGrid.map_link_head_node_to_link
    ~landlab.grid.base.ModelGrid.map_link_tail_node_to_link
    ~landlab.grid.base.ModelGrid.map_link_vector_sum_to_patch
    ~landlab.grid.base.ModelGrid.map_link_vector_to_nodes
    ~landlab.grid.base.ModelGrid.map_max_of_link_nodes_to_link
    ~landlab.grid.base.ModelGrid.map_max_of_node_links_to_node
    ~landlab.grid.base.ModelGrid.map_max_of_patch_nodes_to_patch
    ~landlab.grid.base.ModelGrid.map_mean_of_link_nodes_to_link
    ~landlab.grid.base.ModelGrid.map_mean_of_patch_nodes_to_patch
    ~landlab.grid.base.ModelGrid.map_min_of_link_nodes_to_link
    ~landlab.grid.base.ModelGrid.map_min_of_node_links_to_node
    ~landlab.grid.base.ModelGrid.map_min_of_patch_nodes_to_patch
    ~landlab.grid.base.ModelGrid.map_node_to_cell
    ~landlab.grid.base.ModelGrid.map_upwind_node_link_max_to_node
    ~landlab.grid.base.ModelGrid.map_upwind_node_link_mean_to_node
    ~landlab.grid.base.ModelGrid.map_value_at_downwind_node_link_max_to_node
    ~landlab.grid.base.ModelGrid.map_value_at_max_node_to_link
    ~landlab.grid.base.ModelGrid.map_value_at_min_node_to_link
    ~landlab.grid.base.ModelGrid.map_value_at_upwind_node_link_max_to_node

Boundary condition control
--------------------------

These are the primary properties for getting and setting the grid boundary
conditions. Changes made to :meth:`~.ModelGrid.status_at_node` and
:meth:`~.ModelGrid.status_at_node` will automatically update the conditions
defined at other grid elements automatically.

.. autosummary::
    :toctree: generated/

    ~landlab.grid.base.ModelGrid.active_faces
    ~landlab.grid.base.ModelGrid.active_links
    ~landlab.grid.base.ModelGrid.active_neighbors_at_node
    ~landlab.grid.base.ModelGrid.boundary_nodes
    ~landlab.grid.base.ModelGrid.closed_boundary_nodes
    ~landlab.grid.base.ModelGrid.core_cells
    ~landlab.grid.base.ModelGrid.core_nodes
    ~landlab.grid.base.ModelGrid.fixed_gradient_boundary_nodes
    ~landlab.grid.base.ModelGrid.fixed_links
    ~landlab.grid.base.ModelGrid.fixed_value_boundary_nodes
    ~landlab.grid.base.ModelGrid.node_at_core_cell
    ~landlab.grid.base.ModelGrid.node_is_boundary
    ~landlab.grid.base.ModelGrid.number_of_active_faces
    ~landlab.grid.base.ModelGrid.number_of_active_links
    ~landlab.grid.base.ModelGrid.number_of_core_cells
    ~landlab.grid.base.ModelGrid.number_of_core_nodes
    ~landlab.grid.base.ModelGrid.number_of_fixed_links
    ~landlab.grid.base.ModelGrid.number_of_patches_present_at_link
    ~landlab.grid.base.ModelGrid.number_of_patches_present_at_node
    ~landlab.grid.base.ModelGrid.open_boundary_nodes
    ~landlab.grid.base.ModelGrid.set_nodata_nodes_to_closed
    ~landlab.grid.base.ModelGrid.set_nodata_nodes_to_fixed_gradient
    ~landlab.grid.base.ModelGrid.status_at_link
    ~landlab.grid.base.ModelGrid.status_at_node

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

    ~landlab.grid.base.ModelGrid.calc_aspect_at_node
    ~landlab.grid.base.ModelGrid.calc_slope_at_node
    ~landlab.grid.base.ModelGrid.calc_hillshade_at_node

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
>>> grid = RasterModelGrid((3, 3))
>>> groups = list(grid.groups)
>>> groups.sort()
>>> groups
['cell', 'corner', 'face', 'grid', 'link', 'node', 'patch']

Create Field Arrays
+++++++++++++++++++
If you just want to create an array but not add it to the grid, you can use
the :meth:`~.ModelGrid.ones` method.

>>> grid.ones(at='node')
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
