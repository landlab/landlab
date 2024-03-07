.. _faq:

Frequently Asked Questions
==========================


What is the difference between a cell and a node?
-------------------------------------------------

A cell is a polygon surrounding a node. Nodes on the grid perimeter do not have
cells. Cells have area, nodes have coordinates.

Why is my node data a 1d array? I'm using a raster...
-----------------------------------------------------

All Landlab data structures have to be compatible with both regular and
irregular grids. A 2d structure for node data might make sense for a raster,
but it wouldn't work for an irregular grid - and moreover, there is also no
sensible way to represent link data as a 2d raster either, even for a raster.
Thus, instead, Landlab stores all data in an order set by *node ID number*.
For a raster, IDs begin at zero in the bottom left corner of the grid, then
run along each row in turn. For links, the IDs also start in the bottom left
and run across then up, but in this case all vertical links are listed, then
all horizontal links. Here's a sketch summary of this scheme for a 4x5 raster::

    NODES:                                    LINKS:
    15------16------17------18------19        *--27-->*--28-->*--29-->*--30-->*
    |       |       |       |       |         ^       ^       ^       ^       ^
    |       |       |       |       |        22      23      24      25      26
    |       |       |       |       |         |       |       |       |       |
    10------11------12------13------14        *--18-->*--19-->*--20-->*--21-->*
    |       |       |       |       |         ^       ^       ^       ^       ^
    |       |       |       |       |        13      14      15      16      17
    |       |       |       |       |         |       |       |       |       |
    5-------6-------7-------8-------9         *---9-->*--10-->*--11-->*--12-->*
    |       |       |       |       |         ^       ^       ^       ^       ^
    |       |       |       |       |         4       5       6       7       8
    |       |       |       |       |         |       |       |       |       |
    0-------1-------2-------3-------4         *---0-->*---1-->*---2-->*---3-->*


How do I set the boundary codes for the edges of a grid?
--------------------------------------------------------

By default, the boundary nodes around the perimeter of a grid are all
open boundaries. For a raster grid, if you want to make one or more sides
closed boundaries, use the grid method
:py:func:`RasterModelGrid.set_closed_boundaries_at_grid_edges <RasterModelGrid.set_closed_boundaries_at_grid_edges>`.

The following code snippet sets the southern boundary nodes to be closed:

.. code-block:: python

    import landlab

    grid = landlab.RasterModelGrid(3, 4)
    grid.set_closed_boundaries_at_grid_edges(True, False, False, False)
    grid.status_at_node
    array([4, 4, 4, 4, 1, 0, 0, 1, 1, 1, 1, 1], dtype=int8)

It's also possible to set the boundary conditions "by hand", if you know the ID of the element you're trying to change:
:

.. code-block:: python

    mynodes_to_close = z < 0.0  # if z is some elevation field
    grid.status_at_node[mynodes_to_close] = grid.BC_NODE_IS_CLOSED
    my_fixed_node = mg.find_nearest_node((1.2, 2.3))
    my_fixed_node
    9
    grid.status_at_node[my_fixed_node] = (
        grid.BC_NODE_IS_FIXED_GRADIENT
    )  # to fix the node closest to (1.2, 2.3)

See also:

  - :py:func:`RasterModelGrid.set_fixed_value_boundaries_at_grid_edges <landlab.grid.raster.RasterModelGrid.set_fixed_value_boundaries_at_grid_edges>`
  - :py:func:`ModelGrid.set_nodata_nodes_to_closed <landlab.grid.base.ModelGrid.set_nodata_nodes_to_closed>`


Can I import Landlab output into ParaView or VisIt?
---------------------------------------------------

See :ref:`How do I get netCDF output? <how_do_i_get_netcdf_output>` below.

.. _how_do_i_get_netcdf_output:

How do I get netCDF output?
---------------------------

At present, Landlab can write output to a netCDF file if you are using a raster grid
(support for unstructured grids is coming later). To create netCDF output, use the function
:py:func:`landlab.io.netcdf.write_netcdf <landlab.io.netcdf.write_netcdf>`.
This function will write to file

(1) the grid geometry, and
(2) any data arrays that are linked to the grid

this will automatically include any arrays that you created with functions
such as
:py:func:`landlab.grid.base.ModelGrid.add_zeros <landlab.grid.base.ModelGrid.add_zeros>`,
as long as you provided a name for the array as one of the arguments.


How do I assign values from nodes to links?
-------------------------------------------

Suppose you have a set of values, such as water depths, that are defined at nodes. How do
you figure out what the corresponding values would be at the links, so you can multiply
these by some other quantity (such as water-surface slope) that is defined on links? Here
are some options:

(1) assign the *average*
(2) assign the upstream value
(3) assign the downstream value

Look at this
`Tutorial <https://mybinder.org/v2/gh/landlab/landlab/master?filepath=notebooks/tutorials/file=mappers/mappers.ipynb>`_
for all the Landlab mappers

How do I test whether my grid is regular or irregular?
------------------------------------------------------

There are a number of cases when designing Landlab components where you'll want to do
something one way if the grid is a raster, or another if it's a Voronoi-derived type.
The way to do this is:

.. code-block:: python

    from landlab import RasterModelGrid, VoronoiDelaunayGrid

    # ...
    if isinstance(mg, RasterModelGrid):
        print("Doing it one way")
    elif isinstance(mg, VoronoiDelaunayGrid):
        print("Doing it the other way")
    else:
        raise TypeError("Landlab did not recognize your grid type!")


How do I modify boundary conditions for part of the grid where I know the coordinates?
--------------------------------------------------------------------------------------

See `this tutorial <https://mybinder.org/v2/gh/landlab/landlab/master?filepath=notebooks/tutorials/boundary_conds/set_BCs_from_xy.ipynb>`_.

I am having trouble installing Landlab on Ubuntu without Anaconda. What is the fix?
-----------------------------------------------------------------------------------

Andy Wickert (5/16) suggests the following:

"The version of setuptools that comes standard on Ubuntu is out-of-date with respect to Landlab's Cython code. Here is the fix:"

.. code-block:: bash

    sudo apt-get install python-setuptools # if you don't have it already
    sudo easy_install pip
    sudo apt-get remove python-setuptools
    pip install setuptools # add "--upgrade" if needed

And then you can cd to landlab and this works:

.. code-block:: bash

    pip install -e .

Support: How can I ask more questions and get help?
---------------------------------------------------

File an issue at
`https://github.com/landlab/landlab <https://github.com/landlab/landlab/issues>`__
using the ``New issue`` button in the upper right.
Tell us about your issue, and we'll be in touch.

How do I keep in touch with Landlab developments?
-------------------------------------------------

There are a few ways to follow Landlab developments. You can

- follow Landlab on `Twitter <https://mobile.twitter.com/landlabtoolkit>`_  @landlabtoolkit,
- "watch" Landlab's GitHub repository,
- file a pull request or an issue at `https://github.com/landlab/landlab <https://github.com/landlab/landlab>`__,

Why are there no other FAQs besides these few?
----------------------------------------------

Because we need your questions. Please feel free to add your own questions by making a GitHub Issue.
