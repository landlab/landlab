.. _standard_name_changes:

Changes to standard names from Landlab 0.x to 1.x
=================================================

As part of our push to version 1 of Landlab, the standard names have been overhauled to enhance
internal consistency. Most of this work happened before our beta launch at the CSDMS meeting, so
should not cause too many problems. However, if in doubt interrogate the most current input and
output names for the component you're currently using with `[component].input_var_names` and
`[component].output_var_names`.

However, a few standard names have had to change since the version 1 beta. To our best knowledge
most of these were not widely used or public-facing. The list is as follows::

    'water__discharge' is now 'surface_water__discharge'
    'water__depth' is now 'surface_water__depth'
    'unit_flux' is now 'hillslope_sediment__unit_volume_flux'
    'lithosphere__vertical_displacement' is now 'lithosphere_surface__elevation_increment'
    'rainfall__daily' is now 'rainfall__daily_depth'

Of these, `'water__depth'` is most likely to impact people, as it formed an input to the
`StreamPowerEroder`. However, for back compatibility, you should still find that that component
is still able to handle both the old and new names.
