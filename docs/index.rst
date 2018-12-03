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

Layers
======

Landlab has the ability to add layers to the grid. Two types of layers are
currently supported. First is EventLayers in which each event is preserved as
an entry into the datastructure, even if no deposition occurs. If you are
interested in chronostratigraphy, this is probably what you are interested in.
Second is MaterialLayers, in which each layer must contain some material.
If an entire layer is eroded in MaterialLayers, the layer is removed.
MaterialLayers will likely use less memory than EventLayers.

  .. toctree::
     :maxdepth: 4

     landlab.layers.eventlayers
     landlab.layers.materiallayers

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
   landlab.components.depth_dependent_diffusion
   landlab.components.transport_length_diffusion
   landlab.components.taylor_nonlinear_hillslope_flux
   landlab.components.depth_dependent_taylor_soil_creep

Fluvial geomorphology
---------------------

.. toctree::
   :maxdepth: 4

   landlab.components.stream_power
   landlab.components.detachment_ltd_erosion
   landlab.components.erosion_deposition
   landlab.components.space

Flow routing
------------

.. toctree::
   :maxdepth: 4

   landlab.components.flow_director
   landlab.components.flow_accum
   landlab.components.flow_routing
   landlab.components.lake_fill
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
  landlab.components.soil_moisture
  landlab.components.greenampt

Landslides
----------

.. toctree::
  :maxdepth: 4

  landlab.components.landslides

Vegetation
----------

.. toctree::
  :maxdepth: 4

  landlab.components.vegetation_dynamics
  landlab.components.plant_competition_ca

Precipitation
-------------

.. toctree::
  :maxdepth: 4

  landlab.components.uniform_precip
  landlab.components.spatial_precip

Weathering
----------

.. toctree::
  :maxdepth: 4

  landlab.components.weathering

Terrain Analysis
----------------

.. toctree::
  :maxdepth: 4

  landlab.components.steepness_index
  landlab.components.chi_index
  landlab.components.drainage_density

Tectonics
---------

.. toctree::
  :maxdepth: 4

  landlab.components.flexure
  landlab.components.gflex
  landlab.components.normal_fault

Fire
----

.. toctree::
  :maxdepth: 4

  landlab.components.fire_generator

Initial conditions: random field generators
-------------------------------------------

.. toctree::
  :maxdepth: 4

  landlab.components.fracture_grid


Lithology
---------
Two objects based on the EventLayers object exist to make it easier to deal
with spatially variable lithology and associated properties. The Lithology
components contain information about spatially variable lithology and connect
with the Landlab model grid so that when rock is eroded or advected upward by
rock uplift the values of rock propeties at the topographic surface are updated.

First is the Lithology component which is a generic object for variable
lithology.

  .. toctree::
     :maxdepth: 4

     landlab.components.lithology

Second is LithoLayers which makes it easy to make layered rock.

   .. toctree::
      :maxdepth: 4

      landlab.components.litholayers


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

  landlab.io


Plotting and Visualization
==========================

.. toctree::
  :maxdepth: 4

  landlab.plot


Utilities and Decorators
========================

.. toctree::
  :maxdepth: 4

  landlab.core.utils
  landlab.utils.decorators
  landlab.grid.decorators
  landlab.utils.source_tracking_algorithm


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
 dev_guide_releases
 dev_guide_components


References
==========

* :ref:`modindex`
* :ref:`search`


.. Search the Index
.. ==================

.. * :ref:`genindex`
