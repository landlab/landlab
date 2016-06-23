Base Class (`ModelGrid`)
========================

General class methods and attributes of the `landlab.grid.base` module
----------------------------------------------------------------------

.. automodule:: landlab.grid.base
    :members:
    :undoc-members:
    :show-inheritance:

Mapping data between different grid elements
--------------------------------------------

.. automodule:: landlab.grid.mappers
    :members:
    :undoc-members:
    :show-inheritance:

Gradient calculation functions
------------------------------

.. automodule:: landlab.grid.gradients
    :members:
    :undoc-members:
    :show-inheritance:

Miscellaneous functions
-----------------------

.. automodule:: landlab.grid.grid_funcs
  :members:
  :undoc-members:
  :show-inheritance:


Grid creation from a formatted input file
-----------------------------------------

.. automodule:: landlab.grid.create
    :members:
    :undoc-members:
    :show-inheritance:

Function/method decorators
--------------------------

.. automodule:: landlab.grid.decorators
    :members:
    :undoc-members:
    :show-inheritance:

Regular Rectilinear Grids (`RasterModelGrid`)
=============================================
Inherits from `ModelGrid` and adds:

  .. automodule:: landlab.grid.raster
      :members:
      :undoc-members:
      :show-inheritance:

Mapping data between different grid elements
--------------------------------------------

.. automodule:: landlab.grid.raster_mappers
    :members:
    :undoc-members:
    :show-inheritance:

Gradient calculation functions
------------------------------

.. automodule:: landlab.grid.raster_gradients
    :members:
    :undoc-members:
    :show-inheritance:

Slope-aspect calculation functions
----------------------------------

.. automodule:: landlab.grid.raster_aspect
    :members:
    :undoc-members:
    :show-inheritance:

Steepest-descent functions
--------------------------

.. automodule:: landlab.grid.raster_steepest_descent
    :members:
    :undoc-members:
    :show-inheritance:

Boundary handling for grid edges
--------------------------------

.. automodule:: landlab.grid.raster_set_status
    :members:
    :undoc-members:
    :show-inheritance:

Miscellaneous raster-grid functions
-----------------------------------

  .. automodule:: landlab.grid.raster_funcs
      :members:
      :undoc-members:
      :show-inheritance:

Voronoi-DeLaunay Grids (`VoronoiDelaunayGrid`)
==============================================
Inherits from `ModelGrid` and adds:

.. automodule:: landlab.grid.voronoi
    :members:
    :undoc-members:
    :show-inheritance:

Hex Grids (`HexModelGrid`)
==========================
Inherits from `VoronoiDelauneyGrid` and adds:

.. automodule:: landlab.grid.hex
    :members:
    :undoc-members:
    :show-inheritance:

Radial Grids (`HexModelGrid`)
=============================
Inherits from `VoronoiDelauneyGrid` and adds:

.. automodule:: landlab.grid.radial
    :members:
    :undoc-members:
    :show-inheritance:
