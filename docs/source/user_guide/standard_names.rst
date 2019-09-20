`Landlab <http://landlab.github.io>`_ |
[[ About | About ]] |
[[ Examples | Examples ]] |
[[ User Guide | User-Guide ]] |
`Reference Manual <http://landlab.readthedocs.org/en/latest/#developer-documentation>`_ |
[[ Tutorials| Tutorials ]] |
[[ FAQs | FAQs ]]

[[ ← Previous topic: Components | Components ]]



Standard naming in Landlab
------------------------------------
See `this page <https://github.com/landlab/landlab/wiki/Components#landlab-standard-naming-conventions>`_
for information on the naming conventions.

.. _standard_name_list:
List of Standard names in Landlab
------------------------------------

+---------------------------------------------------+------------------------------+------------------------------------+ 
| Name                                              | Used by                      | Provided by                        | 
+===================================================+==============================+====================================+ 
| area_coefficient                                  | 'DrainageDensity'            |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| area_exponent                                     | 'DrainageDensity'            |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| bedrock__elevation                                |                              | 'DepthDependentDiffuser'           |
|                                                   |                              | 'DepthDependentCubicDiffuser'      | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| channel__bed_shear_stress                         |                              | 'SedDepEroder'                     | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| channel__chi_index                                |                              | 'ChiFinder'                        | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| channel__depth                                    |                              | 'SedDepEroder'                     | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| channel__discharge                                |                              | 'SedDepEroder'                     | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| channel__mask                                     | 'DrainageDensity'            |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| channel__steepness_index                          |                              |'SteepnessFinder'                   | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| channel__width                                    |                              | 'SedDepEroder'                     | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
|channel_sediment__relative_flux                    |                              | 'SedDepEroder'                     | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| channel_sediment__volumetric_flux                 |                              | 'SedDepEroder'                     | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| channel_sediment__volumetric_transport_capacity   |                              | 'SedDepEroder'                     | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| channelization_threshold                          | 'DrainageDensity'            |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| depression__depth                                 |                              | 'DepressionFinderAndRouter'        | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| depression__outlet_node                           |                              | 'DepressionFinderAndRouter'        | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| drainage_area                                     | 'ChiFinder'                  | 'FlowRouter'                       | 
|                                                   | 'StreamPowerEroder'          | 'FlowAccumulator'                  |
|                                                   | 'StreamPowerSmooth...        |                                    |
|                                                   | ...ThresholdEroder'          |                                    |
|                                                   | 'FastscapeEroder'            |                                    |
|                                                   | 'SedDepEroder'               |                                    |
|                                                   | 'SteepnessFinder'            |                                    |
|                                                   | 'Space'                      |                                    |
|                                                   | 'ErosionDeposition'          |                                    |
+---------------------------------------------------+------------------------------+------------------------------------+ 
| flow__data_structure_D                            |                              | 'FlowAccumulator'                  | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| flow__data_structure_delta                        |                              | 'FlowAccumulator'                  | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| flow__link_to_receiver_node                       | 'ChiFinder'                  | 'FlowRouter'                       | 
|                                                   | 'StreamPowerEroder'          | 'FlowDirectorD8'                   |
|                                                   | 'StreamPowerSmooth...        | 'FlowDirectorSteepest'             |
|                                                   | ...ThresholdEroder'          | 'FlowDirectorMFD'                  |
|                                                   | 'FastscapeEroder'            | 'FlowDirectorDINF'                 |
|                                                   | 'SedDepEroder'               |                                    |
|                                                   | 'SteepnessFinder'            |                                    |
|                                                   | 'DrainageDensity'.           |                                    |
|                                                   | 'Space'                      |                                    |
|                                                   | 'ErosionDeposition'          |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| flow__nodes_not_in_stack                          |                              | 'FlowAccumulator'                  | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| flow__potential                                   |                              | 'PotentialityFlowRouter'           | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| flow__receiver_node                               | 'ChiFinder'                  | 'FlowRouter'                       | 
|                                                   | 'StreamPowerEroder'          | 'FlowDirectorD8'                   |
|                                                   | 'StreamPowerSmooth...        | 'FlowDirectorSteepest'             |
|                                                   | ...ThresholdEroder'          | 'FlowDirectorMFD'                  |
|                                                   | 'FastscapeEroder'            | 'FlowDirectorDINF'                 |
|                                                   | 'SedDepEroder'               |                                    |
|                                                   | 'SteepnessFinder'            |                                    |
|                                                   | 'DrainageDensity'            |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| flow__receiver_nodes                              |                              | 'FlowDirectorMFD'                  |
|                                                   |                              | 'FlowDirectorDINF'                 | 
+---------------------------------------------------+------------------------------+------------------------------------+
| flow__receiver_proportionsflow...                 |                              | 'FlowDirectorMFD'                  |
| ...__link_to_receiver_nodes                       |                              | 'FlowDirectorDINF'                 | 
+---------------------------------------------------+------------------------------+------------------------------------+
| flow__sink_flag.                                  |                              | 'FlowRouter'                       | 
|                                                   |                              | 'FlowDirectorD8'                   |
|                                                   |                              | 'FlowDirectorSteepest'             |
|                                                   |                              | 'FlowDirectorMFD'                  |
|                                                   |                              | 'FlowDirectorDINF'                 |
+---------------------------------------------------+------------------------------+------------------------------------+
| flow__upstream_node_order                         | 'ChiFinder'                  | 'FlowRouter'                       | 
|                                                   | 'StreamPowerEroder'          | 'FlowAccumulator'                  |
|                                                   | 'StreamPowerSmooth...        |                                    |
|                                                   | ...ThresholdEroder'          |                                    |
|                                                   | 'FastscapeEroder'            |                                    |
|                                                   | 'SedDepEroder'               |                                    |
|                                                   | 'SteepnessFinder'            |                                    |
|                                                   | 'Space'                      |                                    |
|                                                   | 'ErosionDeposition'          |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+   
| hillslope_sediment__unit_volume_flux              |                              | 'LinearDiffuser'                   | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| landslide__probability_of_failure                 |                              | 'LandslideProbability'             | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| lithosphere__overlying_pressure_increment         | 'Flexure'                    |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| lithosphere_surface__elevation_increment          |                              | 'Flexure'                          | 
|                                                   |                              | 'gFlex'                            | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| plant__age                                        |                              | 'VegCA'                            | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| plant__live_index                                 |                              | 'VegCA'                            | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| radiation__incoming_shortwave_flux                |                              | 'PotentialEvapotranspiration'      |
|                                                   |                              | 'Radiation'                        | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| radiation__net_flux                               |                              | 'PotentialEvapotranspiration'      | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| radiation__net_longwave_flux                      |                              | 'PotentialEvapotranspiration'      | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| radiation__net_shortwave_flux                     |                              | 'PotentialEvapotranspiration'      |
|                                                   |                              | 'Radiation'                        | 
+---------------------------------------------------+------------------------------+------------------------------------+
| radiation__ratio_to_flat_surface                  | 'PotentialEvapotranspiration'| 'Radiation'                        |
+---------------------------------------------------+------------------------------+------------------------------------+  
| rainfall__daily_depth                             | 'SoilMoisture'               |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| sediment_fill__depth                              |                              | 'SinkFiller'                       | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| slope_coefficient                                 | 'DrainageDensity'            |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| slope_exponent                                    | 'DrainageDensity'            |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| soil__density                                     | 'LandslideProbability'       |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+
| soil__depth                                       | 'ExponentialWeatherer'       | 'DepthDependentDiffuser'           |
|                                                   | 'DepthDependentDiffuser'     | 'DepthDependentCubicDiffuser'      |
|                                                   | 'Space'                      |                                    |
|                                                   | 'DepthDependentCubicDiffuser'|                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+
| soil__flux                                        |                              | 'DepthDependentDiffuser'           |
|                                                   |                              | 'CubicNonLinearDiffuser'           |
|                                                   |                              | 'DepthDependentCubicDiffuser'      |
+---------------------------------------------------+------------------------------+------------------------------------+
| soil__internal_friction_angle                     | 'LandslideProbability'       |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+
| soil__maximum_total_cohesion                      | 'LandslideProbability'       |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+
| soil__mean_relative_wetness                       |                              | 'LandslideProbability'             | 
+---------------------------------------------------+------------------------------+------------------------------------+
| soil__minimum_total_cohesion                      | 'LandslideProbability'       |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+
| soil__mode_total_cohesion                         | 'LandslideProbability'       |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+
| soil__probability_of_saturation                   |                              | 'LandslideProbability'             | 
+---------------------------------------------------+------------------------------+------------------------------------+
| soil__saturated_hydraulic_conductivity            | 'LandslideProbability'       |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+
| soil__thickness                                   | 'LandslideProbability'       |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+
| soil__transmissivity                              | 'LandslideProbability'       |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+
| soil_moisture__initial_saturation_fraction        | 'SoilMoisture'               |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+
| soil_moisture__root_zone_leakage                  | 'SoilMoisture'               |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+
| soil_moisture__saturation_fraction                |                              | 'SoilMoisture'                     | 
+---------------------------------------------------+------------------------------+------------------------------------+
| soil_production__rate                             | 'DepthDependentDiffuser'     | 'ExponentialWeatherer'             |
|                                                   | 'DepthDependentCubicDiffuser'|                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+
| soil_water_infiltration__depth                    | 'SoilInfiltrationGreenAmpt'  | 'SoilInfiltrationGreenAmpt'        | 
+---------------------------------------------------+------------------------------+------------------------------------+
| surface__evapotranspiration                       | 'Vegetation'                 | 'SoilMoisture'                     | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| surface__potential_evapotranspiration_30day_mean  | 'Vegetation'                 |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+
| surface__potential_evapotranspiration_rate        | 'SoilMoisture'               | 'PotentialEvapotranspiration'      |
|                                                   | 'Vegetation'                 |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| surface__runoff                                   |                              | 'SoilMoisture'                     | 
+---------------------------------------------------+------------------------------+------------------------------------+
| surface_load__stress                              | 'gFlex'                      |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| surface_to_channel__minimum_distance              |                              | 'DrainageDensity'                  | 
+---------------------------------------------------+------------------------------+------------------------------------+
| surface_water__depth                              | 'OverlandFlowBates'          | 'OverlandFlowBates'                |
|                                                   | 'OverlandFlow'               | 'OverlandFlow'                     |
|                                                   | 'SoilInfiltrationGreenAmpt'  | 'KinwaveImplicitOverlandFlow'      |
|                                                   | 'DepthSlopeProductErosion'   | 'PotentialityFlowRouter'           | 
|                                                   |                              | 'SoilInfiltrationGreenAmpt'        | 
+---------------------------------------------------+------------------------------+------------------------------------+ 
| surface_water__discharge                          | 'DetachmentLtdErosion'       | 'FlowRouter'                       |
|                                                   |                              | 'OverlandFlowBates'                |
|                                                   |                              | 'OverlandFlow'                     |
|                                                   |                              | 'PotentialityFlowRouter'           | 
|                                                   |                              | 'FlowAccumulator'                  | 
+---------------------------------------------------+------------------------------+------------------------------------+
| surface_water_inflow__discharge                   |                              | 'KinwaveImplicitOverlandFlow'      | 
+---------------------------------------------------+------------------------------+------------------------------------+
| topographic__elevation                            | 'ChiFinder'                  | 'LinearDiffuser'                   |
|                                                   | 'LinearDiffuser'             | 'PerronNLDiffuse'                  |
|                                                   | 'FlowRouter'                 | 'SinkFiller'                       |
|                                                   | 'DepressionFinderAndRouter'  | 'StreamPowerEroder'                |
|                                                   | 'PerronNLDiffuse'            | 'StreamPowerSmooth...              |
|                                                   | 'OverlandFlowBates'          | ...ThresholdEroder'                |
|                                                   | 'OverlandFlow'               | 'FastscapeEroder'                  |
|                                                   | 'KinwaveImplicitOverlandFlow'| 'SedDepEroder'                     |
|                                                   | 'PotentialityFlowRouter'     | 'DetachmentLtdErosion'             |
|                                                   | 'Radiation'                  | 'gFlex'                            |
|                                                   | 'SinkFiller'                 | 'DepthDependentDiffuser'           |
|                                                   | 'StreamPowerEroder'          | 'CubicNonLinearDiffuser'           |
|                                                   | 'StreamPowerSmooth...        | 'DepthSlopeProductErosion'         |
|                                                   | ...ThresholdEroder'          | 'DepthDependentCubicDiffuser'      |
|                                                   | 'FastscapeEroder'            |                                    |
|                                                   | 'SedDepEroder'               |                                    |
|                                                   | 'SteepnessFinder'            |                                    |
|                                                   | 'DetachmentLtdErosion'       |                                    |
|                                                   | 'DepthDependentDiffuser'     |                                    |
|                                                   | 'CubicNonLinearDiffuser'     |                                    |
|                                                   | 'DepthSlopeProductErosion'   |                                    |
|                                                   | 'FlowDirectorD8'             |                                    |
|                                                   | 'FlowDirectorSteepest'       |                                    |
|                                                   | 'FlowDirectorMFD'            |                                    |
|                                                   | 'FlowDirectorDINF'           |                                    |
|                                                   | 'FlowAccumulator'            |                                    |
|                                                   | 'DepthDependentCubicDiffuser'|                                    |
+---------------------------------------------------+------------------------------+------------------------------------+
| topographic__gradient                             |                              | LinearDiffuser'                    |
|                                                   |                              | 'KinwaveImplicitOverlandFlow'      | 
+---------------------------------------------------+------------------------------+------------------------------------+
| topographic__slope                                | 'DetachmentLtdErosion'       | 'DepthDependentDiffuser'           |
|                                                   | 'DepthSlopeProductErosion'   | 'CubicNonLinearDiffuser'           |
|                                                   | 'LandslideProbability'       | 'DepthDependentCubicDiffuser'      | 
+---------------------------------------------------+------------------------------+------------------------------------+
| topographic__specific_contributing_area           | 'LandslideProbability'       |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+
| topographic__steepest_slope                       | 'ChiFinder'                  | 'FlowRouter'                       |
|                                                   | 'StreamPowerEroder'          | 'FlowDirectorD8'                   |
|                                                   | 'SedDepEroder'               | 'FlowDirectorSteepest'             |
|                                                   | 'SteepnessFinder'            | 'FlowDirectorMFD'                  |
|                                                   | 'DrainageDensity'            | 'FlowDirectorDINF'                 |
|                                                   | 'Space'                      |                                    |
|                                                   | 'ErosionDeposition'          |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+
| vegetation__cover_fraction                        | 'SoilMoisture'               | 'Vegetation'                       | 
+---------------------------------------------------+------------------------------+------------------------------------+
| vegetation__cumulative_water_stress               | 'VegCA'                      |                                    | 
+---------------------------------------------------+------------------------------+------------------------------------+
| vegetation__dead_biomass                          |                              | 'Vegetation'                       | 
+---------------------------------------------------+------------------------------+------------------------------------+
| vegetation__dead_leaf_area_index                  |                              | 'Vegetation'                       |
+---------------------------------------------------+------------------------------+------------------------------------+
| vegetation__live_biomass                          |                              | 'Vegetation'                       |
+---------------------------------------------------+------------------------------+------------------------------------+  
| vegetation__live_leaf_area_index                  | 'SoilMoisture'               | 'Vegetation'                       |
+---------------------------------------------------+------------------------------+------------------------------------+  
| vegetation__plant_functional_type                 | 'SoilMoisture'               |                                    |
|                                                   | 'Vegetation'                 |                                    |
|                                                   | 'VegCA'                      |                                    |
+---------------------------------------------------+------------------------------+------------------------------------+  
| vegetation__water_stress                          | 'Vegetation'                 | 'SoilMoisture'                     |
+---------------------------------------------------+------------------------------+------------------------------------+
| water__unit_flux_in                               | 'FlowRouter'                 |                                    | 
|                                                   | 'PotentialityFlowRouter'     |                                    |
|                                                   | 'FlowAccumulator'            |                                    |
+---------------------------------------------------+------------------------------+------------------------------------+ 
| water_surface__gradient                           |                              | 'OverlandFlowBates'                |
|                                                   |                              | 'OverlandFlow'                     |
+---------------------------------------------------+------------------------------+------------------------------------+    




.. _standard_name_changes:

Changes to standard names in Landlab
------------------------------------

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

[[ ← Previous topic: Components | Components ]]