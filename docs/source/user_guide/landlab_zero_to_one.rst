.. _zero_to_one:

Transition from Landlab 0.x to 1.0
==================================

Landlab 1.0beta was introduced May 12, 2016. This new version—the first
official stable release—introduces some powerful new capabilities. It
also brings some standardization, which means changes to certain
function names and ways of doing things. In order to minimize the impact
of these changes on your work, we've provided here a simple guide to
some of the major changes.

Quick summary
-------------

how to make your Landlab programs work with version 1.0beta

Grid element and attribute names
````````````````````````````````

Old: ``link_length``
New: ``length_of_link``

Old: ``node_x``, ``node_y``
New: ``x_of_node``, ``y_of_node`` OR ``node_x``, ``node_y``

(NOTE: ``node_x`` and ``node_y`` are still supported, but best practice
is now to use ``x_of_node`` and ``y_of_node``)

Old:``cell_area``
New:``area_of_cell``

Old: ``face_width``
New: ``width_of_face``

Common grid functions
`````````````````````

Old: ``calculate_gradients_at_active_links``,
``calculate_gradients_at_links``

New: ``calc_grad_at_link`` (input is node-based array or field name;
output is a link-based array [**not** just active links])

Old: ``calculate_flux_divergence_at_nodes``

New: ``calc_flux_div_at_node`` (input is a *link-based array* [**not**
just active links])

New standardized function usage and syntax
------------------------------------------

Syntax of grid properties and functions
---------------------------------------

Data associated with Landlab grids include:

-  arrays that contain information about connections among grid elements

-  arrays that contain numerical attributes

-  grid functions

Arrays that contain information about how various grid elements now use
a naming scheme that looks like this:

``<element(s)>_at_<element>``

Examples: ``link_at_face``, ``node_at_cell``, ``node_at_link_tail``,
``links_at_node``, ``faces_at_cell``

The word before ``_at_`` indicates what is being returned: the IDs of
links, or nodes, or faces, etc. A plural word before ``_at_``, such as
``links_at_node``, indicates that there are multiple values for each
element. For example, there are multiple links connected to every grid
node (in a raster grid, most nodes have four links).

The word(s) after ``_at_`` indicate the dimension of the array returned.
For example, ``grid.link_at_face`` contains an array of the IDs of a
link associated with a face. The array length is equal to the number of
faces. Element 0 contains the ID of the link that crosses face 0, and so
on.

In general, arrays containing a numerical attribute associated with
elements, such as the length of links, are named with the pattern:

``<attribute>_of_<element>``

Examples: ``length_of_link``, ``width_of_face``, ``area_of_cell``

Grid functions that return arrays of information have names that hint at
the complexity of the underlying computation. In general, names that
simply describe a property or attribute, such as ``length_of_link``,
have minimum computational cost. Often they either simply provide an
existing array, or compute the values once but retain them thereafter so
that no further cost is incurred on future calls (though creating the
array may have an extra memory cost). More computationally expensive
functions, such that require one or more multiple floating-point
operations at each grid element, have the prefix ``calc_``. For example,
``calc_grad_at_link`` performs a subtraction and division operation at
each grid link.

Gradient and flux divergence functions
--------------------------------------

New functionality in Landlab's built-in gradient and flux-divergence
methods include:

1 - The primary gradient-calculation function (for scalar quantities
defined on nodes) is now ``calc_grad_at_link``. This replaces the
previous functions ``calculate_gradients_at_active_links`` and
``calculates_gradients_at_links``.

2 - ``calc_grad_at_link`` returns gradient values at **all** links, not
just active links. In general, Landlab's grid functions now only operate
with "complete" sets of grid elements. They will neither accept nor
return partial sets of a given element type (such as active or inactive
links, core or boundary nodes, etc.).

3 - The primary divergence operation function is now
``calc_flux_div_at_node``. This replaces
``calculate_flux_divergence_at_nodes``. Note that it accepts an array
(or grid field name) with flux values defined on **all links** (not just
active links). Alternatively, if you provide an array (or grid field
name) with a length equal to the number of cell faces, it assumes that
these represent flux vectors aligned with the links that cross the
faces.

4 - In general, previous function names that began with ``calculate_``
have been shortened to ``calc_``. Similarly, names involving
``_gradient_`` are shortened to ``_grad_``, and those involving
``_divergence_`` are shortened to ``_div_``.

3 - Names are singular. For example, ``calc_grad_at_link`` replaces
``calculate_gradients_at_links``, even though the function returns a
value for *every* link.

4 - Landlab now provides a *net flux* function:
``calc_net_flux_at_node`` adds up fluxes around the perimeter of each
cell, essentially performing a numerical line integral. Values at
perimeter nodes (which do not have corresponding cells) are either set
to zero, or (if the user passes in an array to hold the result) are left
unchanged. The function ``calc_net_flux_at_node`` is identical to
\`calc_flux_div_at_node except that the net flux is not divided by cell
area.

5 - Stay tuned for equivalent functions that calculates the gradient
along **faces**, and net flux and flux divergence within **patches**
(with the resulting values tied to **corners**). These will appear in a
future version, and will be named: ``calc_grad_at_face``,
``calc_net_flux_at_corner``, and ``calc_flux_div_at_corner``.

For more information and examples about the new **grad** and **div**
functions, see the entries in the Reference manual for gradient
functions and divergence functions.

Changes in Component Usage
--------------------------

Components are now standardised across Landlab, eliminating the numerous
inconsistencies that existed before. This standardisation is as follows:

1 - All component classes can now be imported directly from
``landlab.components``, e.g.:

``from landlab.components import FlowRouter``

2 - All components are now instantiated using the same signature. The
first argument is the ModelGrid. Subsequent arguments are the keywords
to provide the input parameters needed by the component, and any other
flags.

Components no longer take an input file, per se. Rather, in keeping with
the spirit of numpy, the necessary input parameters are passed as
keywords. This has several advantages, including allowing explicit
default values to be present, and clear to the user, and also allowing
dynamic Python objects (e.g., an existing array of values) to be passed
in as arguments. Note however that it is still possible (and indeed
encouraged) to use an input file, but now you will need to turn it into
a Python dictionary before passing it to the component (see below). The
recommended way to do this is with the ``load_params`` function, which
performs typing of arguments automatically and can read a variety of
file types. However, the older ways of using the Landlab
``ModelParameterDictionary`` will also still work (though are
deprecated).

This construction format will be listed explicitly in the documentation.
Try ``help(MyComponent)`` in an interactive session to see it, or look
it up online.

For the moment, many components are back compatible with the old ways of
doing things, but this is deprecated functionality and no longer
documented. It may disappear entirely in future releases.

All this means that all of the following are possible ways to
instantiate a component:

::

   >>> from landlab.components import FastscapeEroder
   >>> from landlab import RasterModelGrid, load_params, ModelParameterDictionary
   >>> mg = RasterModelGrid((4, 5), 1.0)
   >>> sp1 = FastscapeEroder(
   ...     mg, K_sp=1.0e-6
   ... )  # the minimum information needed, passed by hand, OR
   >>> sp2 = FastscapeEroder(
   ...     mg,
   ...     K_sp=np.random.rand(20.0),
   ...     m_sp=0.5,
   ...     n_sp=1.0,
   ...     threshold_sp=0.0,
   ...     rainfall_intensity=1.0,
   ... )  # note the array, OR
   >>> myparamdict1 = load_params("my_input_file.txt")
   >>> sp3 = FastscapeEroder(
   ...     mg, **myparamdict1
   ... )  # note the "**". Necessary args come from the dict, OR
   >>> myparamdict2 = ModelParameterDictionary("my_input_file.txt", auto_type=True)
   >>> sp4 = FastscapeEroder(
   ...     mg, **myparamdict2
   ... )  # ...but it's best practice to use load_params instead
   >>> sp5 = FastscapeEroder(
   ...     mg, "my_input_file.txt"
   ... )  # still works in many cases, but DEPRECATED

3 - All components now have a "run method" with the standardised name
``run_one_step``. The first argument is always the timestep, dt, if
needed by the component. Subsequent arguments may be present as flags to
control run behaviour. As an example:

::

   >>> sp = FastscapeEroder(mg, K_sp=1.0e-6)
   >>> dt = 1000.0
   >>> for i in range(100):  # 100 ka of erosion
   ...     sp.run_one_step(dt)
   ...

The old run methods still exist inside many components, but we encourage
migration to this new standardised format.

4 - ``run_one_step()`` never returns anything. There is no need; the
grid object will already have been updated as necessary.

5 - All components should now have comprehensive and up-to-date
documentation. View it on the website, or in an interactive Python
session use either ``help(MyComponent)`` or ``MyComponent?``.

Standardisation of Component Standard Field Names
-------------------------------------------------

In the interests of internal self consistency and repeatability, the
currently in-use component standard field names have been overhauled.
This is likely to break quite a bit of code, but a search-and-replace
will fix things very fast.

The following represents a (hopefully) almost complete list of the name
substitutions:

-  'channel_bed_shear_stress' → 'channel__bed_shear_stress'
-  'channel_depth' → 'channel__depth'
-  'channel_discharge' → 'channel__discharge'
-  'channel_width' → 'channel__width'
-  'drainage_area' –> We're keeping it BOOM
-  'effective_fluvial_diffusivity' → field removed
-  'elevation' –> 'topographic__elevation'
-  'flow_receiver' –> 'flow__receiver_node'
-  'flow_sinks' –> 'flow__sink_flag'
-  'fluvial_sediment_flux_into_node' →
   'channel_sediment__volumetric_flux'
-  'Fluvial_sediment_transport_capacity' →
   'channel_sediment__volumetric_transport_capacity'
-  'Links_to_flow_receiver' →flow__link_to_receiver_node'
-  'lithosphere__elevation' –> 'lithosphere_surface__elevation'
-  'lithosphere__elevation_increment' →
   'lithosphere_surface__elevation_increment'
-  'planet_surface_sediment__deposition_increment' –>
   'sediment__deposition_increment'
-  'potentiality_field' –> 'flow__potential'
-  'relative_sediment_flux' → 'channel_sediment__relative_flux'
-  'shear_stress' –> 'channel__bed_shear_stress'
-  'slope_at_nodes' –> 'topographic__steepest_slope' (slope === downhill
   gradient)
-  'stream_power_erosion' → field removed
-  'surface_gradient' –> 'topographic__slope'
-  'upstream_ID_order' –> 'flow__upstream_node_order'
-  'Upstream_node_order' –> 'flow__upstream_node_order'
-  'water__volume_flux' → 'water__discharge'
-  'water__volume_flux_in' → 'water__unit_flux_in' (special case in flow
   router)
-  'water__volume_flux_magnitude', → 'water__discharge'
-  'water__volume_flux_xcomponent', → 'water__discharge_x_component'
-  'water__volume_flux_ycomponent', → 'water__discharge_y_component'
-  'water_depth' –> 'water__depth'
-  'water_discharge' –> 'water__discharge'
-  'water_discharge_at_nodes' –> 'water__discharge'
-  'water_surface_slope_at_nodes' –> 'water_surface__gradient'

These changes are likely to occur in components probably not released as
part of LL1.0, but will have likely occurred once the components return
in a future release:

-  'ActualEvapotranspiration' –> surface__evapotranspiration_rate
-  'CumulativeWaterStress' –> vegetation__cumulative_water_stress
-  'DeadBiomass' –> vegetation__dead_biomass
-  'DeadLeafAreaIndex –> vegetation__dead_leaf_area_index
-  'Drainage' –> duplicate of 'drainage_area' ?
-  'Elevation' –> duplicate of topographic__elevation, or needs to be
   more specific
-  'LiveBiomass' –> vegetation__live_biomass
-  'LiveLeafAreaIndex' –> vegetation__live_leaf_area_index
-  'NetLongWaveRadiation' –> radiation__net_longwave
-  'NetRadiation' –> radiation__net
-  'NetShortWaveRadiation' –> radiation__net_shortwave
-  'PlantAge' –> plant__age
-  'PlantLiveIndex' –> plant__live_index
-  'PotentialEvapotranspiration' –>
   surface__potential_evapotranspiration_rate
-  'RadiationFactor' –> radiation__ratio_to_flat_surface
-  'Runoff' –> I think I'm OK with runoff__rate
-  'SaturationFraction' –> soil_moisture__saturation_fraction
-  'TotalShortWaveRadiation' –> radiation__incoming_shortwave
-  'VegetationCover', –> vegetation__cover_fraction
-  'VegetationType' –> vegetation__type
-  'WaterStress' –> soil_moisture__water_stress
