Frequently Asked Questions
==========================

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


Why are there no other FAQs besides these few?
----------------------------------------------

Because the FAQ section hasn't been finished yet.
