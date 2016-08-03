#! /usr/env/python
"""
This script auto-constructs an elementary web documentation structure & API for
Landlab grid methods, based around Sphinx and the new LLCATS type declaration
system.

THIS DOES NOT YET WORK, but illustrates the principle of how the utils work.
"""
from landlab.core.utils import get_categories_from_grid_methods
from copy import copy
import re

grid_types = ('ModelGrid', 'RasterModelGrid', 'VoronoiDelaunayGrid',
              'HexModelGrid', 'RadialModelGrid')
str_sequence = ('Raster', 'Irregular Voronoi-cell', 'Hexagonal', 'Radial')
paths = ('raster', 'voronoi', 'hex', 'radial')

autosummary = '\n\n.. autosummary::\n    :toctree: generated/\n\n'


all_methods_for_cat_allgrid = {}
all_cats_for_method_allgrid = {}
for grid_type in grid_types:
    catdict, griddict, _ = get_categories_from_grid_methods(grid_type)
    all_methods_for_cat_allgrid[grid_type] = copy(catdict)
    all_cats_for_method_allgrid[grid_type] = copy(griddict)

auto_index = (
'''

.. _auto_index:

=================================
Landlab grid method documentation
=================================

This is structured documentation built to describe methods that users will
find useful when working with Landlab's grids. It is structured first by
the kind of action you are trying to perform using the grid, then further
subdivided by grid type.

Where a grid method is relevant in multiple places, it will be listed
multiple times.

This documentation is built and maintained automatically using Landlab's
LLCATS system.


Getting Information about a Grid
================================

:ref:`These methods <auto_whole_grid_info>` give information about the grid
and its individual elements: numbers of elements present, and quantifying
their properties (size, shape, orientation).

:ref:`Getting information about a grid <auto_whole_grid_info>`

Information about nodes

Information about links

Information about patches

Information about cells

Information about faces

Information about corners


Element Connectivity
====================

These methods describe the topology of the grid: how the elements fit
together, which joins to which, and what the neighbors of each element are.

Connections to nodes

Connections to links

Connections to patches

Connections to cells

Connections to faces

Connections to corners


Data Fields in ModelGrid
========================

Landlab fields store the values of data defined on the elements of a
Landlab grid. Grid fields are the primary mechanism by which Landlab
components share data between themselves. These methods describe how to
create, interact with, and delete fields.

Creating and deleting fields
Access and modify data in an existing field


Gradients, Fluxes, and Divergences on the Grid
==============================================

These methods allow calculation of differences, gradients, fluxes,
divergences, and similar operations on data defined across a Landlab grid.


Surface analysis
================

These methods provide GIS-like methods to describe a surface defined over
the grid, for example, mean local slopes, aspects, and hillshades.


Mapping between Elements
========================

These methods allow the mapping of values defined on one Landlab element
(e.g., links) onto another (e.g., nodes). They incorporate typically useful
schemes like averaging neighboring values or upwind (maximum value) and
downwind (minimum value) mappings.


Boundary condition control
==========================

These methods allow the setting and handling of boundary conditions defined
on the various elements of Landlab grids.


Indentifying element subsets
============================

Sometimes, you need to extract subsets of elements from the grid based on
spatial criteria like:

"What is the nearest node to this point?"
"Which cell is this point inside?"
"What nodes define the left edge?"
"Is this coordinate within the grid?"

These methods will help.

''')

f = open('./auto_index.rst', "wb")
f.write(auto_index)
f.close()

grid_inf_subheads = (
    'Information about the grid as a whole', 'Information about nodes',
    'Information about links', 'Information about patches',
    'Information about cells', 'Information about faces',
    'Information about corners')

grid_info__main = (
'''

.. _auto_whole_grid_info:

Getting Information about a Grid
================================

The following attributes, properties, and methods provide data about the
grid, its geometry, and the connectivity among the various elements. Each
grid element has an ID number, which is also its position in an array that
contains information about that type of element. For example, the *x*
coordinate of node 5 would be found at `grid.node_x[5]`.

The naming of grid-element arrays is `[attribute]_at_[element]`, where
`attribute` is the name of the data in question, and `element` is the
element to which the attribute applies. For example, the property
`node_at_cell` contains the ID of the node associated with each cell. For
example, `node_at_cell[3]` contains the node ID of the node associated
with cell 3.
The attribute is singular if there is only one value per element; for
example, there is only one node associated with each cell. It is plural
when there are multiple values per element; for example, the
`faces_at_cell` array contains multiple faces for each cell. Exceptions to
these general rules are functions that return indices of a subset of
all elements of a particular type.
For example, you can obtain an array with IDs of only the core nodes using
`core_nodes`, while `active_links` provides an array of IDs of active links
(only). Finally, attributes that represent a measurement of something, such
as the length of a link or the surface area of a cell, are described using
`_of_`, as in the example `area_of_cell`.


''')

for subhead, grid in zip(grid_inf_subheads, grid_types[1:]):
    grid_info__main += subhead + '\n' + '-'*len(subhead) + '\n\n'
    for subsub in str_sequence:
        grid_info__main += '    :ref:`' + subsub + ' <GINF_' + grid + '>` \n'
    grid_info__main += '\n\n'

f = open('./auto_grid_info__main.rst', "wb")
f.write(grid_info__main)
f.close()

grid_info__whole_grid_0 = (
'''
Information about the grid as a whole
=====================================

''')
for grid, print_name, path in zip(grid_types[1:], str_sequence, paths):
    grid_info__whole_grid_0 += (
        '\n' + '.. _GINF_' + grid + ':\n\n' +
        print_name + '\n' + '-'*len(print_name) + autosummary)
    allmethsGINF = all_methods_for_cat_allgrid[grid]['GINF']
    allmethsGINF.sort()
    for meth in allmethsGINF:
        if 'DEPR' not in all_cats_for_method_allgrid[grid][meth]:
            grid_info__whole_grid_0 += (
                '    ~landlab.grid.' + path + '.' + grid + '.' + meth + '\n'
            )
    grid_info__whole_grid_0 += '\n\n'

f = open('./auto_grid_info__whole_grid_0.rst', "wb")
f.write(grid_info__whole_grid_0)
f.close()












'''
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
    ~landlab.grid.base.ModelGrid.display_grid
    ~landlab.grid.base.ModelGrid.ndim
    ~landlab.grid.base.ModelGrid.number_of_elements
    ~landlab.grid.base.ModelGrid.set_units
    ~landlab.grid.base.ModelGrid.size
    ~landlab.grid.raster.RasterModelGrid.shape

Information about nodes
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.base.ModelGrid.number_of_nodes
    ~landlab.grid.base.ModelGrid.number_of_core_nodes
    ~landlab.grid.base.ModelGrid.nodes
    ~landlab.grid.base.ModelGrid.core_nodes
    ~landlab.grid.base.ModelGrid.boundary_nodes
    ~landlab.grid.base.ModelGrid.open_boundary_nodes
    ~landlab.grid.base.ModelGrid.closed_boundary_nodes
    ~landlab.grid.base.ModelGrid.fixed_value_boundary_nodes
    ~landlab.grid.base.ModelGrid.fixed_gradient_boundary_nodes
    ~landlab.grid.base.ModelGrid.node_x
    ~landlab.grid.base.ModelGrid.node_y
    ~landlab.grid.base.ModelGrid.status_at_node
    ~landlab.grid.base.ModelGrid.cell_at_node
    ~landlab.grid.base.ModelGrid.patches_at_node
    ~landlab.grid.base.ModelGrid.links_at_node
    ~landlab.grid.base.ModelGrid.link_dirs_at_node
    ~landlab.grid.base.ModelGrid.active_link_dirs_at_node
    ~landlab.grid.base.ModelGrid.neighbors_at_node

Information about links
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.base.ModelGrid.number_of_links
    ~landlab.grid.base.ModelGrid.number_of_active_links
    ~landlab.grid.base.ModelGrid.number_of_fixed_links
    ~landlab.grid.base.ModelGrid.active_links
    ~landlab.grid.base.ModelGrid.fixed_links
    ~landlab.grid.base.ModelGrid.node_at_link_tail
    ~landlab.grid.base.ModelGrid.node_at_link_head
    ~landlab.grid.base.ModelGrid.face_at_link
    ~landlab.grid.base.ModelGrid.status_at_link

Information about cells
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.base.ModelGrid.number_of_cells
    ~landlab.grid.base.ModelGrid.number_of_core_cells
    ~landlab.grid.base.ModelGrid.core_cells
    ~landlab.grid.base.ModelGrid.node_at_cell
    ~landlab.grid.base.ModelGrid.node_at_core_cell
    ~landlab.grid.base.ModelGrid.area_of_cell

Information about faces
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.base.ModelGrid.number_of_faces
    ~landlab.grid.base.ModelGrid.number_of_active_faces
    ~landlab.grid.base.ModelGrid.active_faces
    ~landlab.grid.base.ModelGrid.link_at_face

Information about patches
+++++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.base.ModelGrid.number_of_patches
    ~landlab.grid.base.ModelGrid.nodes_at_patch

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

    # START HERE check that all functions listed below are included above, ignore ones that start with underscores(_)

Gradients, fluxes, and divergences on the grid
----------------------------------------------

Landlab is designed to easily calculate gradients in quantities across the
grid, and to construct fluxes and flux divergences from them. Because these
calculations tend to be a little more involved than property lookups, the
methods tend to start with `calc_`.

.. autosummary::
    :toctree: generated/

    ~landlab.grid.base.ModelGrid.calc_diff_at_link
    ~landlab.grid.base.ModelGrid.calc_grad_of_link
    ~landlab.grid.base.ModelGrid.calc_net_flux_at_node
    ~landlab.grid.base.ModelGrid.calc_flux_div_at_node
    ~landlab.grid.base.ModelGrid.calc_unit_normal_of_patch
    ~landlab.grid.raster.RasterModelGrid.calc_unit_normals_of_patch_subtriangles
    ~landlab.grid.base.ModelGrid.calc_grad_of_patch
    ~landlab.grid.base.ModelGrid.calc_slope_of_patch
    ~landlab.grid.base.ModelGrid.calc_slope_of_node

Mappers
-------

These methods allow mapping of values defined on one grid element type onto a
second, e.g., mapping upwind node values onto links, or mean link values onto
nodes.

...


Boundary condition control
--------------------------

These are the primary properties for getting and setting the grid boundary
conditions. Changes made to :meth:`~.ModelGrid.status_at_node` and
:meth:`~.ModelGrid.status_at_node` will automatically update the conditions
defined at other grid elements automatically.

.. autosummary::
    :toctree: generated/

    ~landlab.grid.base.ModelGrid.status_at_node
    ~landlab.grid.base.ModelGrid.status_at_link
    ~landlab.grid.base.ModelGrid.update_links_nodes_cells_to_new_BCs

...

Identifying node subsets
------------------------

These methods are useful in identifying subsets of nodes, e.g., closest node
to a point; nodes at edges.

...

Surface analysis
----------------

These methods permit the kinds of surface analysis that you might expect to
find in GIS software.

.. autosummary::
    :toctree: generated/

    ~landlab.grid.base.ModelGrid.calc_slope_of_node
    ~landlab.grid.base.ModelGrid.hillshade
    ~landlab.grid.base.ModelGrid.aspect
    ~landlab.grid.base.ModelGrid.calc_distances_of_nodes_to_point

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

'''