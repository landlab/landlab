.. _standard_name_list:

Landlab Standard Names
======================

The table below gives a list of all the standard names in use in Landlab as of version
1.0. It may be out of date, so always check your component documentation! The table is
also available in Hobley et al., ESurf, in prep.

Several names have changed since the beta release; these are:

* 'water__discharge' is now 'surface_water__discharge'
* 'water__depth' is now 'surface_water__depth'
* 'unit_flux' is now 'hillslope_sediment__unit_volume_flux'
* 'lithosphere__vertical_displacement' is now 'lithosphere_surface__elevation_increment'
* 'rainfall__daily' is now 'rainfall__daily_depth'

Hopefully these will not impact anyone.

You can use the landlab command line interface to get information about field use in
Landlab as a whole; try::

    > (landlab provided_by && landlab used_by) | sort | uniq

...for an alphabetised list of all names currently in use. For more information on the
available command line options, just call `landlab` on its own from the bash prompt.


.. htmlonly::

+--------------------------------------------------+-----------------------------+-----------------------------+
| Field name                                       | Used by                     | Provided by                 |
+==================================================+=============================+=============================+
| channel__bed_shear_stress                        |                             | SedDepEroder                |
+--------------------------------------------------+-----------------------------+-----------------------------+
| channel__chi_index                               |                             | ChiFinder                   |
+--------------------------------------------------+-----------------------------+-----------------------------+
| channel__depth                                   |                             | SedDepEroder                |
+--------------------------------------------------+-----------------------------+-----------------------------+
| channel__discharge                               |                             | SedDepEroder                |
+--------------------------------------------------+-----------------------------+-----------------------------+
| channel__steepness_index                         |                             | SteepnessFinder             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| channel__width                                   |                             | SedDepEroder                |
+--------------------------------------------------+-----------------------------+-----------------------------+
| channel_sediment__relative_flux                  |                             | SedDepEroder                |
+--------------------------------------------------+-----------------------------+-----------------------------+
| channel_sediment__volumetric_flux                |                             | SedDepEroder                |
+--------------------------------------------------+-----------------------------+-----------------------------+
| channel_sediment__volumetric_transport_capacity  |                             | SedDepEroder                |
+--------------------------------------------------+-----------------------------+-----------------------------+
| depression__depth                                | DepressionFinderAndRouter   |                             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| depression__outlet_node                          | DepressionFinderAndRouter   |                             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| drainage_area                                    | ChiFinder                   | FlowRouter                  |
|                                                  | FastscapeEroder             |                             |
|                                                  | SedDepEroder                |                             |
|                                                  | SteepnessFinder             |                             |
|                                                  | StreamPowerEroder           |                             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| flow__link_to_receiver_node                      | ChiFinder                   | FlowRouter                  |
|                                                  | FastscapeEroder             |                             |
|                                                  | SedDepEroder                |                             |
|                                                  | SteepnessFinder             |                             |
|                                                  | StreamPowerEroder           |                             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| flow__receiver_node                              | ChiFinder                   | FlowRouter                  |
|                                                  | FastscapeEroder             |                             |
|                                                  | SedDepEroder                |                             |
|                                                  | SteepnessFinder             |                             |
|                                                  | StreamPowerEroder           |                             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| flow__sink_flag                                  |                             | FlowRouter                  |
+--------------------------------------------------+-----------------------------+-----------------------------+
| flow__upstream_node_order                        | ChiFinder                   | FlowRouter                  |
|                                                  | FastscapeEroder             |                             |
|                                                  | SedDepEroder                |                             |
|                                                  | SteepnessFinder             |                             |
|                                                  | StreamPowerEroder           |                             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| lithosphere__overlying_pressure_increment        | Flexure                     |                             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| lithosphere_surface__elevation_increment         |                             | Flexure                     |
|                                                  |                             | gFlex                       |
+--------------------------------------------------+-----------------------------+-----------------------------+
| plant__age                                       |                             | VegCA                       |
+--------------------------------------------------+-----------------------------+-----------------------------+
| plant__live_index                                |                             | VegCA                       |
+--------------------------------------------------+-----------------------------+-----------------------------+
| radiation__incoming_shortwave_flux               |                             | PotentialEvapotranspiration |
|                                                  |                             | Radiation                   |
+--------------------------------------------------+-----------------------------+-----------------------------+
| radiation__net_flux                              |                             | PotentialEvapotranspiration |
+--------------------------------------------------+-----------------------------+-----------------------------+
| radiation__net_longwave_flux                     |                             | PotentialEvapotranspiration |
+--------------------------------------------------+-----------------------------+-----------------------------+
| radiation__net_shortwave_flux                    |                             | PotentialEvapotranspiration |
|                                                  |                             | Radiation                   |
+--------------------------------------------------+-----------------------------+-----------------------------+
| radiation__ratio_to_flat_surface                 | PotentialEvapotranspiration | Radiation                   |
+--------------------------------------------------+-----------------------------+-----------------------------+
| rainfall__daily_depth                            | SoilMoisture                |                             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| sediment_fill__depth                             |                             | SinkFiller                  |
+--------------------------------------------------+-----------------------------+-----------------------------+
| soil_moisture__initial_saturation_fraction       | SoilMoisture                |                             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| soil_moisture__root_zone_leakage                 |                             | SoilMoisture                |
+--------------------------------------------------+-----------------------------+-----------------------------+
| soil_moisture__saturation_fraction               |                             | SoilMoisture                |
+--------------------------------------------------+-----------------------------+-----------------------------+
| soil_water_infiltration__depth                   | SoilInfiltrationGreenAmpt   | SoilInfiltrationGreenAmpt   |
+--------------------------------------------------+-----------------------------+-----------------------------+
| surface__evapotranspiration                      | Vegetation                  | SoilMoisture                |
+--------------------------------------------------+-----------------------------+-----------------------------+
| surface__potential_evapotranspiration_30day_mean | Vegetation                  |                             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| surface__potential_evapotranspiration_rate       | SoilMoisture                | PotentialEvapotranspiration |
|                                                  | Vegetation                  |                             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| surface__runoff                                  |                             | SoilMoisture                |
+--------------------------------------------------+-----------------------------+-----------------------------+
| surface_load__stress                             | gFlex                       |                             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| surface_water__depth                             | KinematicWaveRengers        | KinematicWaveRengers        |
|                                                  | OverlandFlow                | OverlandFlow                |
|                                                  | OverlandFlowBates           | OverlandFlowBates           |
|                                                  | SoilInfiltrationGreenAmpt   | SoilInfiltrationGreenAmpt   |
+--------------------------------------------------+-----------------------------+-----------------------------+
| surface_water__discharge                         | DetachmentLtdErosion        | FlowRouter                  |
|                                                  |                             | KinematicWaveRengers        |
|                                                  |                             | OverlandFlow                |
|                                                  |                             | OverlandFlowBates           |
+--------------------------------------------------+-----------------------------+-----------------------------+
| surface_water__velocity                          |                             | KinematicWaveRengers        |
+--------------------------------------------------+-----------------------------+-----------------------------+
| topographic__elevation                           | ChiFinder                   | DetachmentLtdErosion        |
|                                                  | DepressionFinderAndRouter   | FastscapeEroder             |
|                                                  | DetachmentLtdErosion        | gFlex                       |
|                                                  | FastscapeEroder             | LinearDiffuser              |
|                                                  | FlowRouter                  | PerronNLDiffuse             |
|                                                  | KinematicWaveRengers        | SedDepEroder                |
|                                                  | LinearDiffuser              | SinkFiller                  |
|                                                  | OverlandFlow                | StreamPowerEroder           |
|                                                  | OverlandFlowBates           |                             |
|                                                  | PerronNLDiffuse             |                             |
|                                                  | Radiation                   |                             |
|                                                  | SedDepEroder                |                             |
|                                                  | SinkFiller                  |                             |
|                                                  | SteepnessFinder             |                             |
|                                                  | StreamPowerEroder           |                             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| topographic__gradient                            |                             | LinearDiffuser              |
+--------------------------------------------------+-----------------------------+-----------------------------+
| topographic__slope                               | DetachmentLtdErosion        |                             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| topographic__steepest_slope                      | ChiFinder                   | FlowRouter                  |
|                                                  | SedDepEroder                |                             |
|                                                  | SteepnessFinder             |                             |
|                                                  | StreamPowerEroder           |                             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| hillslope_sediment__unit_volume_flux             |                             | LinearDiffuser              |
+--------------------------------------------------+-----------------------------+-----------------------------+
| vegetation__cover_fraction                       | SoilMoisture                |                             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| vegetation__cumulative_water_stress              | VegCA                       | Vegetation                  |
+--------------------------------------------------+-----------------------------+-----------------------------+
| vegetation__dead_biomass                         |                             | Vegetation                  |
+--------------------------------------------------+-----------------------------+-----------------------------+
| vegetation__dead_leaf_area_index                 |                             | Vegetation                  |
+--------------------------------------------------+-----------------------------+-----------------------------+
| vegetation__live_biomass                         |                             | Vegetation                  |
+--------------------------------------------------+-----------------------------+-----------------------------+
| vegetation__live_leaf_area_index                 | SoilMoisture                | Vegetation                  |
+--------------------------------------------------+-----------------------------+-----------------------------+
| vegetation__plant_functional_type                | SoilMoisture                |                             |
|                                                  | VegCA                       |                             |
|                                                  | Vegetation                  |                             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| vegetation__water_stress                         | Vegetation                  | SoilMoisture                |
+--------------------------------------------------+-----------------------------+-----------------------------+
| water__unit_flux_in                              | FlowRouter                  |                             |
+--------------------------------------------------+-----------------------------+-----------------------------+
| water_surface__gradient                          |                             | OverlandFlow                |
|                                                  |                             | OverlandFlowBates           |
+--------------------------------------------------+-----------------------------+-----------------------------+