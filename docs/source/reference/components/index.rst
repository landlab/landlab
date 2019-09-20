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

   landlab.components.diffusion
   landlab.components.nonlinear_diffusion
   landlab.components.depth_dependent_diffusion
   landlab.components.transport_length_diffusion
   landlab.components.taylor_nonlinear_hillslope_flux
   landlab.components.depth_dependent_taylor_soil_creep

Fluvial geomorphology
---------------------

.. toctree::
   :maxdepth: 2

   landlab.components.stream_power
   landlab.components.detachment_ltd_erosion
   landlab.components.erosion_deposition
   landlab.components.space

Flow routing
------------

.. toctree::
   :maxdepth: 2

   landlab.components.flow_director
   landlab.components.flow_accum
   landlab.components.flow_routing
   landlab.components.lake_fill
   landlab.components.sink_fill

Shallow water hydrodynamics
---------------------------

.. toctree::
   :maxdepth: 2

   landlab.components.overland_flow

Land surface and groundwater hydrology
--------------------------------------

.. toctree::
  :maxdepth: 2

  landlab.components.radiation
  landlab.components.pet
  landlab.components.soil_moisture
  landlab.components.greenampt
  landlab.components.groundwater

Landslides
----------

.. toctree::
  :maxdepth: 2

  landlab.components.landslides

Vegetation
----------

.. toctree::
  :maxdepth: 2

  landlab.components.vegetation_dynamics
  landlab.components.plant_competition_ca

Precipitation
-------------

.. toctree::
  :maxdepth: 2

  landlab.components.uniform_precip
  landlab.components.spatial_precip

Weathering
----------

.. toctree::
  :maxdepth: 2

  landlab.components.weathering

Terrain Analysis
----------------

.. toctree::
  :maxdepth: 2

  landlab.components.steepness_index
  landlab.components.chi_index
  landlab.components.drainage_density
  landlab.components.channel_profile
  landlab.components.hack_calculator

Tectonics
---------

.. toctree::
  :maxdepth: 2

  landlab.components.flexure
  landlab.components.flexure.ext
  landlab.components.gflex
  landlab.components.normal_fault

Fire
----

.. toctree::
  :maxdepth: 2

  landlab.components.fire_generator


Fracture Generation
-------------------

.. toctree::
  :maxdepth: 2

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
     :maxdepth: 2

     landlab.components.lithology

Second is LithoLayers which makes it easy to make layered rock.

   .. toctree::
      :maxdepth: 2

      landlab.components.litholayers
