
How do I set the boundary codes for the edges of a grid?
--------------------------------------------------------

By default, the boundary nodes around the perimeter of a grid are all
open boundaries. For a raster grid, if you want to make one or more sides
closed boundaries, use the grid method :func:`~landlab.grid.base.ModelGrid.set_inactive_boundaries`.

The following code snippet sets the southern boudary nodes to be :data:`~landlab.grid.base.INACTIVE_BOUNDARY`.

  >>> import landlab
  >>> grid = landlab.RasterModelGrid(3, 4)
  >>> grid.set_inactive_boundaries(True, False, False, False)
  >>> grid.node_status
  array([4, 4, 4, 4, 1, 0, 0, 1, 1, 1, 1, 1], dtype=int8)

.. seealso::

  :func:`~landlab.grid.base.ModelGrid.set_fixed_value_boundaries`,
  :func:`~landlab.grid.base.ModelGrid.set_nodata_nodes_to_inactive`


Can I import Landlab output into ParaView or VisIt?
---------------------------------------------------

See :ref:`how_to_export_to_netcdf` below.


.. _how_to_export_to_netcdf:

How do I get netCDF output?
---------------------------

At present, Landlab can write output to a netCDF file if you are using a raster grid
(support for unstructured grids is coming later). A tutorial example is provided in
:ref:`landlab_tools_and_tricks`.  To create netCDF output, use the function
:func:`~landlab.io.netcdf.write_netcdf`. This function will write to file
  1. the grid geometry, and
  2. any data arrays that are linked to the grid
this will automatically include any arrays that you created with functions
such as :func:`~landlab.grid.base.ModelGrid.add_zeros`, as long as you
provided a name for the array as one of the arguments.


How do I assign values from nodes to links?
-------------------------------------------

Suppose you have a set of values, such as water depths, that are defined at nodes. How do
you figure out what the corresponding values would be at the links, so you can multiply
these by some other quantity (such as water-surface slope) that is defined on links? Here 
are some options:

(1) To assign the *average* 


Why are there no other FAQs besides these few?
----------------------------------------------

Because the FAQ section hasn't been finished yet.
