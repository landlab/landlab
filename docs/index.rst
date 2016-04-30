.. landlab documentation master file, created by
   sphinx-quickstart on Tue Apr 23 23:40:10 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Find Landlab's `User Guide <https://github.com/landlab/landlab/wiki/User-Guide>`_ on the `Landlab Wiki <https://github.com/landlab/landlab/wiki/User-Guide>`_

==============================
Landlab Developer API
==============================

The *Landlab Developer API* is a general reference manual for Landlab.

Grids
=======================

As of Landlab version 0.2, there are four types of Landlab grid:
 - Raster
 - Voronoi-DeLaunay
 - Hex
 - Radial

The base class is `ModelGrid` with subclasses `RasterModelGrid` and `VoronoiDelaunayGrid`.

`VoronoiDelaunayGrid` has two further specialized subclasses: `HexModelGrid` and `RadialModelGrid`.

Methods and properties common to all grids
--------------------------
.. toctree::
   :maxdepth: 4

   landlab.grid.mappers
   landlab.grid.gradients
   landlab.grid.divergence
   landlab.grid.grid_funcs
   landlab.grid.create
   landlab.grid.base
   landlab.grid.decorators

Specialized methods and properties for Rectilinear Grids 'raster grids'
--------------------------
Inherits from 'ModelGrid' and adds:

.. toctree::
   :maxdepth: 4

   landlab.grid.raster
   landlab.grid.raster_mappers
   landlab.grid.raster_gradients
   landlab.grid.raster_aspect
   landlab.grid.raster_steepest_descent
   landlab.grid.raster_set_status
   landlab.grid.raster_funcs

Specialized methods and properties for Voronoi-Delaunay grids
--------------------------

Inherits from 'ModelGrid' and adds:

.. toctree::
   :maxdepth: 4

   landlab.grid.voronoi

Specialized methods and properties for hex grids
--------------------------

Inherits from 'VoronoiDelauneyGrid' and adds:

.. toctree::
   :maxdepth: 4

   landlab.grid.hex

Specialized methods and properties for radial grids
--------------------------

Inherits from 'VoronoiDelauneyGrid' and adds:

.. toctree::
   :maxdepth: 4

   landlab.grid.radial


Components
=======================

.. toctree::
   :maxdepth: 4

   landlab.components


Input/Output (IO)
=======================

.. toctree::
   :maxdepth: 4

   landlab.io


Plotting and Visualization
=======================

.. toctree::
   :maxdepth: 4

   landlab.plot


Cellular Automata (CA)
=======================

.. toctree::
   :maxdepth: 4

   landlab.ca


Contributing to Landlab
=======================

If you're intending to make changes to the Landlab code base,
or want to develop your own components, we recommend you follow
these specialized developer install instructions.

.. toctree::
 :maxdepth: 3

 dev_guide_install
 dev_guide_components


References
==========

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
