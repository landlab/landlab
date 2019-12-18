.. _api.components:

==========
Components
==========

This section contains documentation and API reference information for the
following categories of components:


Hillslope geomorphology
-----------------------

.. toctree::
   :maxdepth: 2

   diffusion
   nonlinear_diffusion
   depth_dependent_diffusion
   transport_length_diffusion
   taylor_nonlinear_hillslope_flux
   depth_dependent_taylor_soil_creep

Fluvial geomorphology
---------------------

.. toctree::
   :maxdepth: 2

   stream_power
   detachment_ltd_erosion
   erosion_deposition
   space

Flow routing
------------

.. toctree::
   :maxdepth: 2

   flow_director
   flow_accum
   flow_routing
   lake_fill
   sink_fill

Shallow water hydrodynamics
---------------------------

.. toctree::
   :maxdepth: 2

   overland_flow

Land surface hydrology
----------------------

.. toctree::
  :maxdepth: 2

  radiation
  pet
  soil_moisture
  greenampt

Groundwater hydrology
---------------------

.. toctree::
  :maxdepth: 2

  groundwater

Landslides
----------

.. toctree::
  :maxdepth: 2

  landslides

Vegetation
----------

.. toctree::
  :maxdepth: 2

  vegetation_dynamics
  plant_competition_ca

Precipitation
-------------

.. toctree::
  :maxdepth: 2

  uniform_precip
  spatial_precip

Weathering
----------

.. toctree::
  :maxdepth: 2

  weathering

Terrain Analysis
----------------

.. toctree::
  :maxdepth: 2

  steepness_index
  chi_index
  drainage_density
  profile
  channel_profiler
  hack_calculator

Tectonics
---------

.. toctree::
  :maxdepth: 2

  flexure
  flexure.ext
  gflex
  normal_fault

Fire
----

.. toctree::
  :maxdepth: 2

  fire_generator


Fracture Generation
-------------------

.. toctree::
  :maxdepth: 2

  fracture_grid


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
     :maxdepth: 2

     lithology

Second is LithoLayers which makes it easy to make layered rock.

   .. toctree::
      :maxdepth: 2

      litholayers


Alphabetical Listing of Modules
-------------------------------
.. toctree::

   channel_profiler
   chi_index
   depth_dependent_diffusion
   depth_dependent_taylor_soil_creep
   detachment_ltd_erosion
   diffusion
   drainage_density
   erosion_deposition
   fire_generator
   flexure
   flow_accum
   flow_director
   flow_routing
   fracture_grid
   gflex
   groundwater
   hack_calculator
   lake_fill
   landslides
   lateral_erosion
   lithology
   nonlinear_diffusion
   normal_fault
   overland_flow
   pet
   plant_competition_ca
   potentiality_flowrouting
   radiation
   sink_fill
   soil_moisture
   space
   spatial_precip
   steepness_index
   stream_power
   taylor_nonlinear_hillslope_flux
   transport_length_diffusion
   uniform_precip
   vegetation_dynamics
   weathering

Module contents
---------------

.. automodule:: landlab.components
   :members:
   :undoc-members:
   :show-inheritance:
