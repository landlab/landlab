
List of Landlab Models
======================

(pending)


List of Landlab Components
==========================

This is a (probably incomplete) list of the existing components for 
landlab. Most of these remain quite experimental, and may change in the 
future! Documentation may also be wanting in some cases.

Where a specific method is listed beneath the module class, this is the primary method of the module. This is the method you call in a loop to actually "crank the handle" of the process.


Functional modules
------------------

.. autoclass:: landlab.components.diffusion.diffusion.DiffusionComponent
	:members: diffuse


.. autoclass:: landlab.components.fire_generator.generate_fire.FireGenerator


.. autoclass:: landlab.components.flexure.flexure.Flexure


.. autoclass:: landlab.components.flow_routing.route_flow_dn.FlowRouter
	:members: route_flow


.. autoclass:: landlab.components.nonlinear_diffusion.Perron_nl_diffuse.PerronNLDiffuse
	:members: diffuse


.. autoclass:: landlab.components.overland_flow.generate_overland_flow_deAlmeida.OverlandFlow
	:members: run_one_step
	*This might not work...*


.. autoclass:: landlab.components.pet.potential_evapotranspiration_field.PotentialEvapotranspiration
	:members: update
	*This component implements a model for potential evapotranspiration.*


.. autoclass:: landlab.components.radiation.radiation_field.Radiation
	:members: update
	*This component implements a model to compute spatial radiation distribution.

.. autoclass:: landlab.components.stream_power.fastscape_stream_power.FastscapeEroder
	:members: erode


.. autoclass:: landlab.components.single_vegetation.single_vegetation_field.Vegetation
	:members: update
	*This component implements spatial vegetation dynamics.*


.. autoclass:: landlab.components.soil_moisture.soil_moisture_field.SoilMoisture
	:members: update
	*This component implements spatial soil moisture dynamics.*

.. autoclass:: landlab.components.uniform_precip.generate_uniform_precip.PrecipitationDistribution
	:members: update
	*This class generates rainfall events and interstorm intervals based
	on the Poisson-like statistical model of Eagleson (1978).*


Here be dragons!
----------------

These components are all in principle working, but either deviate from up 
to date recommended Landlab style guidelines, or are otherwise 
"confusing". 
All remain definitely experimental.
Procede at your own risk!

.. automodule:: landlab.components.fracture_grid.fracture_grid
	Note that this module runs in a script-style, and uses a syntax quite different to the "recommended" up to date landlab guidelines.

.. autoclass:: landlab.components.linkca.link_ca.LinkCellularAutomaton


Under development
-----------------

These components are currently under development for Landlab:

	* Impact cratering
	* Glacial flow, in a variety of forms
	* ...
