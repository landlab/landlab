.. _api.grid:

=============
Landlab Grids
=============

Grid types
----------

Landlab presently supports seven types of grids. The base class is ``ModelGrid``
with subclasses ``RasterModelGrid``, ``VoronoiDelaunayGrid``, and ``IcosphereGlobalGrid``.
``VoronoiDelaunayGrid`` has three further specialized subclasses: ``FramedVoronoiGrid``,
``HexModelGrid`` and ``RadialModelGrid``. A final class is ``NetworkModelGrid``.

The following is an introduction to their properties and methods:

.. toctree::
   :maxdepth: 1

   raster
   voronoi
   framed_voronoi
   hex
   radial
   network
   icosphere


Additional Methods and Properties
---------------------------------

.. toctree::
   :maxdepth: 1

   base
   create
   decorators
   diagonals
   divergence
   gradients
   grid_funcs
   link_status
   mappers
   nodestatus
   raster_aspect
   raster_funcs
   raster_gradients
   raster_mappers
   raster_set_status
   warnings

API for each grid type
----------------------
.. toctree::
   :maxdepth: 1

   base
   raster
   voronoi
   framed_voronoi
   hex
   radial
   network
   icosphere

Additional Grid Base Classes
----------------------------

.. toctree::

  unstructured
