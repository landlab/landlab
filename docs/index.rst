.. landlab documentation master file, created by
   sphinx-quickstart on Tue Apr 23 23:40:10 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Find Landlab's
`User Guide <https://github.com/landlab/landlab/wiki/User-Guide>`_ on the
`Landlab Wiki <https://github.com/landlab/landlab/wiki/User-Guide>`_

==============================================
Landlab Reference Manual and API Documentation
==============================================

A guide to Landlab's classes and code.

Grids
=====

Grid types
----------

As of Landlab version 0.2, there are four types of Landlab grid:

-  :ref:`Raster <Raster>`
-  :ref:`Voronoi-Delaunay <Voronoi>`
-  :ref:`Hex <Hex>`
-  :ref:`Radial <Radial>`

The base class is `ModelGrid` with subclasses `RasterModelGrid` and
`VoronoiDelaunayGrid`.

`VoronoiDelaunayGrid` has two further specialized subclasses: `HexModelGrid`
and `RadialModelGrid`.

Methods and properties common to all grids
------------------------------------------

.. toctree::
   :maxdepth: 4

   landlab.grid.base
   landlab.grid.mappers
   landlab.grid.gradients
   landlab.grid.divergence
   landlab.grid.create
   landlab.grid.decorators

Specialized methods and properties for Rectilinear Grids 'raster grids'
-----------------------------------------------------------------------

Landlab's rectilinear grids are implemented by the class `RasterModelGrid`,
which inherits from `ModelGrid` and adds the following:

.. toctree::
   :maxdepth: 4

   landlab.grid.raster

Specialized methods and properties for Voronoi-Delaunay grids
-------------------------------------------------------------

Landlab's Voronoi-Delaunay grids are implemented by the class
`VoronoiDelaunayGrid`, which inherits from `ModelGrid` and adds the following:

.. toctree::
   :maxdepth: 4

   landlab.grid.voronoi

Specialized methods and properties for hex grids
------------------------------------------------

Landlab's hex/trigonal grids are implemented by the class `HexModelGrid`,
which inherits from `VoronoiDelauneyGrid` and adds the following:

.. toctree::
   :maxdepth: 4

   landlab.grid.hex

Specialized methods and properties for radial grids
---------------------------------------------------

Landlab's radial grids are implemented by the class `RadialModelGrid`, which
inherits from `VoronoiDelauneyGrid` and adds the following:

.. toctree::
   :maxdepth: 4

   landlab.grid.radial


Components
==========

This section contains documentation and API reference information for the
following categories of components:

Hillslope geomorphology
-----------------------

.. toctree::
   :maxdepth: 4

   landlab.components.diffusion
   landlab.components.nonlinear_diffusion

Fluvial geomorphology
---------------------

.. toctree::
   :maxdepth: 4

   landlab.components.stream_power
   landlab.components.detachment_ltd_erosion

Flow routing
------------

.. toctree::
   :maxdepth: 4

   landlab.components.flow_routing
   landlab.components.sink_fill

Shallow water hydrodynamics
---------------------------

.. toctree::
   :maxdepth: 4

   landlab.components.overland_flow

Land surface hydrology
----------------------

.. toctree::
  :maxdepth: 4

  landlab.components.radiation
  landlab.components.pet
..  landlab.components.soil_moisture

Landslides
----------

.. toctree::
  :maxdepth: 4

  landlab.components.landslides

.. Vegetation
.. ----------

.. toctree::
..  :maxdepth: 4

..  landlab.components.single_vegetation
..  landlab.components.vegetation_ca

Precipitation
-------------

.. toctree::
  :maxdepth: 4

  landlab.components.uniform_precip

Terrain Analysis
----------------

.. toctree::
  :maxdepth: 4

  landlab.components.steepness_index
  landlab.components.chi_index

Tectonics
---------

.. toctree::
  :maxdepth: 4

  landlab.components.flexure
  landlab.components.gflex

.. Fire
.. ----

.. toctree::
..  :maxdepth: 4

..  landlab.components.fire_generator

.. Initial conditions: random field generators
.. -------------------------------------------

.. toctree::
..  :maxdepth: 4

..  landlab.components.fracture_grid

The Component base class
------------------------

.. toctree::
  :maxdepth: 4

  landlab.core.model_component

Input/Output (IO)
=================

This section documents various methods you can use to bring in data and write
output to a file.

.. toctree::
  :maxdepth: 4

  landlab.io.esri_ascii
  landlab.io.netcdf
  landlab.io.vtk
  landlab.io.native_landlab
  landlab.core.model_parameter


Plotting and Visualization
==========================

.. toctree::
  :maxdepth: 4

  landlab.plot.imshow
  landlab.plot.channel_profile
  landlab.plot.colors


Utilities and Decorators
========================

.. toctree::
  :maxdepth: 4

  landlab.core.utils
  landlab.utils.decorators
..  landlab.utils.jaggedarray
..  landlab.utils.fault_facet_finder
..  landlab.utils.structured_grid


Cellular Automata (CA)
======================

.. toctree::
  :maxdepth: 4

  landlab.ca.celllab_cts
  landlab.ca.raster_cts
  landlab.ca.oriented_raster_cts
  landlab.ca.hex_cts
  landlab.ca.oriented_hex_cts
  landlab.ca.boundaries



Contributing to Landlab
=======================

If you're intending to make changes to the Landlab code base,
or want to develop your own components, we recommend you follow
these `specialized developer install instructions
<https://github.com/landlab/landlab/wiki/Installing-Landlab-from-source-code-(%22developer-install%22)>`_.

For more information about contributing to the Landlab code base (i.e. developing components), see `Developing with git and github <https://github.com/landlab/landlab/wiki/Developing-with-github-and-git>`_.



References
==========

* :ref:`modindex`
* :ref:`search`


.. Search the Index
.. ==================

.. * :ref:`genindex`
