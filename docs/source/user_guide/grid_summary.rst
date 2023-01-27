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
to which the attribute applies. For example, the property :meth:`~.RasterModelGrid.node_at_cell`
is an array of length *number_of_cells* where each element of the array
is the *ID* of the node associated with that cell.
(e.g. ``node_at_cell[3]`` is the *ID* of the node associated with cell 3).
In this case the *attribute* is singular since there is only one value per element.
Sometimes there are multiple attributes associated with each element. In this
case, the attribute is plural. For example, the :meth:`~.RasterModelGrid.faces_at_cell` array
contains multiple faces for each cell. Exceptions to these general rules are
functions that return indices of a subset of all elements of a particular type.
For example, you can obtain an array with *IDs* of only the core nodes using
:meth:`~.RasterModelGrid.core_nodes`, while :meth:`~.RasterModelGrid.active_links` provides an array
of *IDs* of active links (only). Finally, attributes that represent a measurement
of something, such as the length of a link or the surface area of a cell, are
described using ``_of_`` (e.g. :meth:`~.RasterModelGrid.area_of_cell`).

------
Fields
------

:class:`~.ModelGrid` inherits from the :class:`~.GraphFields` class. This
provides :class:`~.ModelGrid`, and its subclasses, with the ability to, optionally,
store data values associated with the different types grid elements
(nodes, cells, etc.). In particular, as part of :meth:`~.ModelGrid.__init__`,
data field *groups* are added to the :class:`~.ModelGrid` that provide containers to
put data fields into. There is one group for each of the eight grid elements
(node, cell, link, face, core_node, core_cell, active_link, and active_face).


.. toctree::
  :caption: Summary of grid methods
  :glob:
  :hidden:

  grid_methods/*
