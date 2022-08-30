.. _api.grid.grid_summary:

====================
List of Grid Methods
====================

--------------------------------
Getting Information about a Grid
--------------------------------

The following attributes, properties, and methods provide data about the grid,
its geometry, and the connectivity among the various elements. Each grid
element has an ID number, which is also its position in an array that
contains information about that type of element. For example, the *x*
coordinate of node 5 would be found at ``grid.x_of_node[5]``.

The naming of grid-element arrays is ``<attribute>_at_<element>``, where
*attribute* is the name of the data in question, and *element* is the element
to which the attribute applies. For example, the property ``node_at_cell``
contains the ID of the node associated with each cell. For example,
``node_at_cell[3]`` contains the *ID* of the node associated with cell 3.
The *attribute* is singular if there is only one value per element; for
example, there is only one node associated with each cell. It is plural when
there are multiple values per element; for example, the `faces_at_cell` array
contains multiple faces for each cell. Exceptions to these general rules are
functions that return indices of a subset of all elements of a particular type.
For example, you can obtain an array with IDs of only the core nodes using
``core_nodes``, while `active_links` provides an array of IDs of active links
(only). Finally, attributes that represent a measurement of something, such as
the length of a link or the surface area of a cell, are described using ``_of_``,
as in the example ``area_of_cell``.


.. toctree::
  :caption: Summary of grid methods
  :glob:
  :hidden:
  
  grid_methods/*


Fields
======

:class:`~.ModelGrid` inherits from the :class:`~.GraphFields` class. This
provides :class:`~.ModelGrid`, and its subclasses, with the ability to, optionally,
store data values associated with the different types grid elements
(nodes, cells, etc.). In particular, as part of :meth:`~.ModelGrid.__init__`,
data field *groups* are added to the :class:`~.ModelGrid` that provide containers to
put data fields into. There is one group for each of the eight grid elements
(node, cell, link, face, core_node, core_cell, active_link, and active_face).

Creating Fields
===============

:class:`~.ModelGrid` inherits several useful methods for creating new data
fields and adding new data fields to a :class:`~.ModelGrid` instance. Methods to add or
create a new data array follow the ``numpy`` syntax for creating arrays. The
following methods create and, optionally, initialize new arrays. These arrays
are of the correct size but a new field will not be added to the field:

.. autosummary::
    :nosignatures:

    ~landlab.field.graph_field.GraphFields.empty
    ~landlab.field.graph_field.GraphFields.ones
    ~landlab.field.graph_field.GraphFields.zeros

Adding Fields to a grid
=======================

Unlike the equivalent ``numpy`` functions, these do not take a size argument
as the size of the returned arrays is determined from the size of the
:class:`~.ModelGrid`. However, the keyword arguments are the same as those of the ``numpy``
equivalents.

The following methods will create a new array and add a reference to that
array to the ModelGrid:

.. autosummary::
    :nosignatures:

    ~landlab.grid.raster.RasterModelGrid.add_empty
    ~landlab.grid.raster.RasterModelGrid.add_field
    ~landlab.grid.raster.RasterModelGrid.add_ones
    ~landlab.grid.raster.RasterModelGrid.add_zeros
    ~landlab.grid.raster.RasterModelGrid.delete_field

These methods operate in the same way as the previous set except that, in
addition to creating a new array, the newly-created array is added to the
ModelGrid. The calling signature is the same but with the addition of an
argument that gives the name of the new field as a string. The additional
method, :meth:`~.GraphFields.add_field`, adds a previously allocation
array to the ModelGrid. If the array is of the incorrect size it will raise
``ValueError``.

Field info
==========

Use the following methods/attributes get information about the stored data
fields:

.. autosummary::
    :nosignatures:

    ~landlab.field.graph_field.GraphFields.size
    ~landlab.field.graph_field.GraphFields.keys
    ~landlab.field.graph_field.GraphFields.has_group
    ~landlab.field.graph_field.GraphFields.has_field
    ~landlab.grid.raster.RasterModelGrid.field_units
    ~landlab.grid.raster.RasterModelGrid.field_values
    ~landlab.field.graph_field.GraphFields.groups

Example:

.. code:: python

  >>> grid.has_field("my_field_name", at="node")

