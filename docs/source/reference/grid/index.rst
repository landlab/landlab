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
   :maxdepth: 2

   landlab.grid.raster
   landlab.grid.voronoi
   landlab.grid.hex
   landlab.grid.radial
   landlab.grid.network


Methods and properties common to all grids
------------------------------------------

.. toctree::
   :maxdepth: 2

   landlab.grid.base
   landlab.grid.mappers
   landlab.grid.gradients
   landlab.grid.divergence
   landlab.grid.create
   landlab.grid.decorators


Additional methods
------------------

.. toctree::
   :maxdepth: 2

    landlab.grid.raster_aspect
    landlab.grid.raster_funcs
    landlab.grid.raster_mappers
    landlab.grid.raster_set_status


Additional Grid Base Classes
----------------------------

.. toctree::
   :maxdepth: 2

    landlab.grid.structured_quad
    landlab.grid.unstructured
