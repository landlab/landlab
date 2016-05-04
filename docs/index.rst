.. landlab documentation master file, created by
   sphinx-quickstart on Tue Apr 23 23:40:10 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Find Landlab's `User Guide <https://github.com/landlab/landlab/wiki/User-Guide>`_ on the `Landlab Wiki <https://github.com/landlab/landlab/wiki/User-Guide>`_

==============================================
Landlab Reference Manual and API Documentation
==============================================

The *Landlab Developer API* is a general reference manual for Landlab.

Grids
=====

Grid types
-------------

As of Landlab version 0.2, there are four types of Landlab grid:
 - Raster
 - Voronoi-Delaunay
 - Hex
 - Radial

The base class is `ModelGrid` with subclasses `RasterModelGrid` and `VoronoiDelaunayGrid`.

`VoronoiDelaunayGrid` has two further specialized subclasses: `HexModelGrid` and `RadialModelGrid`.

Methods and properties common to all grids
------------------------------------------
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
-----------------------------------------------------------------------

Landlab's rectilinear grids are implemented by the class `RasterModelGrid`,
which inherits from `ModelGrid` and adds the following:

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
-------------------------------------------------------------

Landlab's Voronoi-Delaunay grids are implemented by the class `VoronoiDelaunayGrid`, which inherits from `ModelGrid` and adds the following:

.. toctree::
   :maxdepth: 4

   landlab.grid.voronoi

Specialized methods and properties for hex grids
------------------------------------------------

Landlab's hex/trigonal grids are implemented by the class `HexModelGrid`, which inherits from `VoronoiDelauneyGrid` and adds the following:

.. toctree::
   :maxdepth: 4

   landlab.grid.hex

Specialized methods and properties for radial grids
---------------------------------------------------

Landlab's radial grids are implemented by the class `RadialModelGrid`, which inherits from `VoronoiDelauneyGrid` and adds the following:

.. toctree::
   :maxdepth: 4

   landlab.grid.radial


Components
==========

This section contains documentation and API reference information for the following categories of components:

Hillslope geomorphology
---------------------------

.. toctree::
   :maxdepth: 4

   landlab.components.diffusion
   landlab.components.nonlinear_diffusion

Fluvial geomorphology
---------------------------

.. toctree::
   :maxdepth: 4

   landlab.components.stream_power
   landlab.components.transport_limited_fluvial

Flow routing
-------------------

.. toctree::
   :maxdepth: 4

   landlab.components.flow_routing

Shallow water hydrodynamics
-------------------

.. toctree::
   :maxdepth: 4

   landlab.components.overland_flow

Land surface hydrology
-------------------

.. toctree::
  :maxdepth: 4

  landlab.components.radiation
  landlab.components.PET
  landlab.components.soil_moisture

Vegetation
-------------------

.. toctree::
  :maxdepth: 4

  landlab.components.single_vegetation
  landlab.components.vegetation_ca

Precipitation
-------------------

.. toctree::
  :maxdepth: 4

  landlab.components.uniform_precip

Terrain Analysis
-------------------

.. toctree::
  :maxdepth: 4

  landlab.components.steepness_index
  landlab.components.chi_index
  landlab.components.dem_support

Glacial Processes
-------------------

.. toctree::
  :maxdepth: 4

  landlab.components.glacier_thin_ice_model

Tectonics
-------------------

.. toctree::
  :maxdepth: 4

  landlab.components.flexure
  landlab.components.gflex

Fire
-------------------

.. toctree::
  :maxdepth: 4

  landlab.components.fire_generator

Impact cratering
-------------------

.. toctree::
  :maxdepth: 4

  landlab.components.craters

Initial conditions: random field generators
-------------------

.. toctree::
  :maxdepth: 4

  landlab.components.fracture_grid


Input/Output (IO)
=================

.. toctree::
   :maxdepth: 4

   landlab.io


Plotting and Visualization
==========================

.. toctree::
   :maxdepth: 4

   landlab.plot


Cellular Automata (CA)
======================

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

* :ref:`modindex`
* :ref:`search`


Search the Index
==================

* :ref:`genindex`
