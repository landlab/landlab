Landlab Components
==================

For a tutorial introduction to using the component library, see `here <https://github.com/landlab/drivers/blob/master/notebooks/component_tutorial.ipynb>`_.

The Landlab Component Library
-----------------------------

Landlab offers an ever-growing library of components that aim to describe individual or closely associated suites of surface processes. Components are designed to be “plug-and-play”, and to interact with each other with the minimum of technical difficulties. Each makes use of Landlab grid fields to enable the sharing of data between the components, and we aim to have a relatively standardized way of interacting with and using each different one.

Landlab components exist as classes, and can be imported from *landlab.components*.

.. note::

    We emphasize that at the moment most Landlab components are under active development, and it is possible that changes may occur to e.g. input argument format with little or no warning. Please let the Landlab development team know if you’re making heavy use of a component, so we can avoid unnecessarily breaking your code in future Landlab releases! Components are presented as-is, and in this essentially beta release we don't make any guarantees about the stability or otherwise of the implementations.


Available Landlab components
----------------------------

The current library includes the following components in essentially full working order:

* linear diffusion: *diffusion.diffusion.LinearDiffuser*
* nonlinear hillslope diffusion: *nonlinear_diffusion.Perron_nl_diffuse.PerronNLDiffuse*
* a fire generator: *fire_generator.generate_fire.FireGenerator*
* a simple lithospheric flexure model: *flexure.flexure.FlexureComponent*
* a more complex flexural model, incorporating A. Wickert’s gFlex code: *gflex.flexure.gFlex*
* a single-direction (“D8 generalized” or “Dn”) flow router: *flow_routing.route_flow_dn.FlowRouter*
* a thin ice glacial approximation: *glacier_thin_ice_model.glacier.Glacier*
* a shallow overland flow approximation (following Bates et al., 2010): *overland_flow.generate_overland_flow_DEM.OverlandFlow*
* a potential evapotranspiration module: *pet.potential_evapotranspiration_field.PotentialEvapotranspiration*
* a solar total incident shortwave radiation calculator: *radiation.radiation_field.Radiation*
* a soil moisture module: *soil_moisture.soil_moisture_field.SoilMoisture*
* a simple stream power (pure detachment limited) fluvial component: *stream_power.stream_power.StreamPowerEroder*, or *stream_power.fastscape_stream_power.SPEroder*
* a transport-limited fluvial component: *transport_limited_fluvial.tl_fluvial_monodirectional_v3.TransportLimitedEroder*
* a sediment-flux dependent shear stress based fluvial incision component: *stream_power.sed_flux_dep_incision.SedDepEroder*
* a generator for storms over a landscape: *uniform_precip.generate_uniform_precip.PrecipitationDistribution*
* A powerful, generalized, continuous-time stochastic cellular automaton, that can be specialized for many other tasks: *cellular_automata.landlab_ca.LandlabCellularAutomaton*

Under active development are:

* a vegetation cellular automaton
* an impact cratering simulator
* divergent and mixed convergent-divergent flow routers
* a deltaic simulator


Implementing a Component
------------------------

Although the vagaries of their precise implementation vary due to their development histories, Landlab components tend to share a basic grammar for their usage:

A component class is imported from the library as 
>>> from landlab.components.[noun].[what_component_does] import ComponentClass

A component is instantiated like:
>>> compobj = ComponentClass(ModelGrid, ‘parameter_file.txt’, …)

A component has a primary method, which “does the thing”. The documentation for the component will typically tell you what is it called and how to work it (NB: you can get documentation for any object you have created in an interactive python environment by typing a “?” after it, e.g., “compobj?”; quit by pressing “q”). However, most components will run for a single timestep with a syntax something like:

>>> compobj.do_the_process(timestep, …[any other arguments])

Running one of these methods will update the fields held in common by the single grid object which you linked to all your components during component instantiation. If you look inside the grid fields having run one of these methods, you’ll see the new fields it has created and populated. The docstrings for the component should make it clear which fields the component needs to have in the grid as inputs, and which it modifies and/or creates as outputs. **ALWAYS check the documentation for a component you are about to use!**


Component Standard Properties
-----------------------------

As part of our rolling efforts to standardize and improve Landlab, we are also trying to implement a standardized set of properties that all components will have. These give automated information on fields, units, etc. For a fully compliant component, you will find you can call:

	component._name 		(a string)
	component._input_var_names 	(a set giving input field names)
	component._output_var_names	(a set giving output field names)
	component._var_units 		(a dict, with var_name keys)
	component._var_mapping		(a dict with var_name keys, giving ‘node’, ‘link’, etc)
	component._var_defs		(a dict with var_name keys, giving short descriptions)

See `the tutorials <https://github.com/landlab/drivers/blob/master/notebooks/component_tutorial.ipynb>`_ for examples of use cases with one, two and more coupled components.


Landlab Standard Naming Conventions
-----------------------------------

.. note::

    We are currently in the process of improving the Landlab-wide standardization of our naming conventions. Currently in use standard field names may change suddenly in the near future! Note we will however be making efforts to maintain backward compatibility with names currently in use.

The Landlab component library attempts to make use of a relatively standardized set of names across the various components, in order to maximize ease of component coupling. If you’re familiar with the concept of the `CSDMS standard naming conventions <http://csdms.colorado.edu/wiki/CSDMS_Standard_Names>`_, note that we have tried to strike a balance between the rigor and uniqueness of those names and a more user-friendly, succinct approach. Nonetheless, you may recognise the basic style of the names:

	**thing_descibed__what_is_described**

e.g., *topographic__elevation*, *water_surface__gradient*, *water__volume_flux*

 You can see a list of the names currently in use here: :ref:`Landlab Standard Names <standard_names>`.


Dealing with nonstandard names
++++++++++++++++++++++++++++++

The large number of developers on Landlab and historical accident have meant that despite our best efforts you’ll inevitably find instances where different components use different names for the same thing. In these cases, you need to make equivalent two fields in the grid which have different names so that two components can talk to each other. This is actually easy; you can just do:

>>> mg.add_field(‘node’, ‘second_name’, mg.at_node[‘first_name’])

Note that we are making slow progress towards truly standardizing the component library, but these kind of idiosyncrasies might yet persist for a while! 
