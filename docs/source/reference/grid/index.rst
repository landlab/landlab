.. _api.grid:

=============
Landlab Grids
=============

An extensive index to the Landlab grid and its methods is found on
:py:mod:`the following page <landlab.grid>`. Below is a short hyperlinked
summary.

Grid types
----------

Landlab presently supports five types of grids:

-  :ref:`Raster <Raster>`
-  :ref:`Voronoi-Delaunay <Voronoi>`
-  :ref:`Hex <Hex>`
-  :ref:`Radial <Radial>`
-  :ref:`Network <Network>`

The base class is ``ModelGrid`` with subclasses ``RasterModelGrid``
``VoronoiDelaunayGrid``.

``VoronoiDelaunayGrid`` has two further specialized subclasses: ``HexModelGrid``
and ``RadialModelGrid``.

A final class is ``NetworkModelGrid``.

.. toctree::
   :maxdepth: 1

   raster
   voronoi
   hex
   radial
   network

Systematic Information about Grid Elements
------------------------------------------

.. toctree::
   :maxdepth: 2

   auto/index


Additional Methods and Properties
---------------------------------

.. toctree::
  :maxdepth: 1

   base
   cfuncs
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

Additional Grid Base Classes
----------------------------

.. toctree::

  structured_quad
  unstructured

  Module contents
  ---------------

.. automodule:: landlab.grid
   :members:
   :undoc-members:
   :show-inheritance:
