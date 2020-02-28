.. _standard_name_mapping:

Landlab Standard Name Field-Component Mapping
=============================================

+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| Field Name                                       | Provided By                                                             | Used By                                                                 |
+==================================================+=========================================================================+=========================================================================+
| aquifer__thickness                               |                                                                         | :py:class:`~landlab.components.GroundwaterDupuitPercolator` (node)      |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| aquifer_base__elevation                          | :py:class:`~landlab.components.GroundwaterDupuitPercolator` (node)      |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| aquifer_base__gradient                           |                                                                         | :py:class:`~landlab.components.GroundwaterDupuitPercolator` (link)      |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| area_coefficient                                 | :py:class:`~landlab.components.DrainageDensity` (node)                  |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| area_exponent                                    | :py:class:`~landlab.components.DrainageDensity` (node)                  |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| bedrock__elevation                               |                                                                         | :py:class:`~landlab.components.DepthDependentDiffuser` (node)           |
|                                                  |                                                                         | :py:class:`~landlab.components.DepthDependentTaylorDiffuser` (node)     |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| channel__bed_shear_stress                        |                                                                         | :py:class:`~landlab.components.SedDepEroder` (node)                     |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| channel__chi_index                               |                                                                         | :py:class:`~landlab.components.ChiFinder` (node)                        |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| channel__depth                                   |                                                                         | :py:class:`~landlab.components.SedDepEroder` (node)                     |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| channel__discharge                               |                                                                         | :py:class:`~landlab.components.SedDepEroder` (node)                     |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| channel__mask                                    | :py:class:`~landlab.components.DrainageDensity` (node)                  |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| channel__steepness_index                         |                                                                         | :py:class:`~landlab.components.SteepnessFinder` (node)                  |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| channel__width                                   |                                                                         | :py:class:`~landlab.components.SedDepEroder` (node)                     |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| channel_sediment__relative_flux                  |                                                                         | :py:class:`~landlab.components.SedDepEroder` (node)                     |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| channel_sediment__volumetric_flux                |                                                                         | :py:class:`~landlab.components.SedDepEroder` (node)                     |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| channel_sediment__volumetric_transport_capacity  |                                                                         | :py:class:`~landlab.components.SedDepEroder` (node)                     |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| channelization_threshold                         | :py:class:`~landlab.components.DrainageDensity` (node)                  |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| depression__depth                                |                                                                         | :py:class:`~landlab.components.DepressionFinderAndRouter` (node)        |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| depression__outlet_node                          |                                                                         | :py:class:`~landlab.components.DepressionFinderAndRouter` (node)        |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| distance_to_divide                               |                                                                         | :py:class:`~landlab.components.HackCalculator` (node)                   |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| drainage_area                                    | :py:class:`~landlab.components.ChannelProfiler` (node)                  | :py:class:`~landlab.components.FlowAccumulator` (node)                  |
|                                                  | :py:class:`~landlab.components.ChiFinder` (node)                        | :py:class:`~landlab.components.LakeMapperBarnes` (node)                 |
|                                                  | :py:class:`~landlab.components.FastscapeEroder` (node)                  | :py:class:`~landlab.components.LossyFlowAccumulator` (node)             |
|                                                  | :py:class:`~landlab.components.HackCalculator` (node)                   |                                                                         |
|                                                  | :py:class:`~landlab.components.LakeMapperBarnes` (node)                 |                                                                         |
|                                                  | :py:class:`~landlab.components.LateralEroder` (node)                    |                                                                         |
|                                                  | :py:class:`~landlab.components.SedDepEroder` (node)                     |                                                                         |
|                                                  | :py:class:`~landlab.components.SteepnessFinder` (node)                  |                                                                         |
|                                                  | :py:class:`~landlab.components.StreamPowerEroder` (node)                |                                                                         |
|                                                  | :py:class:`~landlab.components.StreamPowerSmoothThresholdEroder` (node) |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| flood_status_code                                |                                                                         | :py:class:`~landlab.components.DepressionFinderAndRouter` (node)        |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| flow__data_structure_delta                       | :py:class:`~landlab.components.LakeMapperBarnes` (node)                 | :py:class:`~landlab.components.FlowAccumulator` (node)                  |
|                                                  |                                                                         | :py:class:`~landlab.components.LakeMapperBarnes` (node)                 |
|                                                  |                                                                         | :py:class:`~landlab.components.LossyFlowAccumulator` (node)             |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| flow__link_direction                             |                                                                         | :py:class:`~landlab.components.FlowDirectorSteepest` (link)             |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| flow__link_to_receiver_node                      | :py:class:`~landlab.components.ChannelProfiler` (node)                  | :py:class:`~landlab.components.FlowDirectorD8` (node)                   |
|                                                  | :py:class:`~landlab.components.ChiFinder` (node)                        | :py:class:`~landlab.components.FlowDirectorDINF` (node)                 |
|                                                  | :py:class:`~landlab.components.DrainageDensity` (node)                  | :py:class:`~landlab.components.FlowDirectorMFD` (node)                  |
|                                                  | :py:class:`~landlab.components.ErosionDeposition` (node)                | :py:class:`~landlab.components.FlowDirectorSteepest` (node)             |
|                                                  | :py:class:`~landlab.components.FastscapeEroder` (node)                  | :py:class:`~landlab.components.LakeMapperBarnes` (node)                 |
|                                                  | :py:class:`~landlab.components.HackCalculator` (node)                   |                                                                         |
|                                                  | :py:class:`~landlab.components.LakeMapperBarnes` (node)                 |                                                                         |
|                                                  | :py:class:`~landlab.components.SedDepEroder` (node)                     |                                                                         |
|                                                  | :py:class:`~landlab.components.Space` (node)                            |                                                                         |
|                                                  | :py:class:`~landlab.components.SteepnessFinder` (node)                  |                                                                         |
|                                                  | :py:class:`~landlab.components.StreamPowerEroder` (node)                |                                                                         |
|                                                  | :py:class:`~landlab.components.StreamPowerSmoothThresholdEroder` (node) |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| flow__potential                                  |                                                                         | :py:class:`~landlab.components.DischargeDiffuser` (node)                |
|                                                  |                                                                         | :py:class:`~landlab.components.PotentialityFlowRouter` (node)           |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| flow__receiver_node                              | :py:class:`~landlab.components.ChannelProfiler` (node)                  | :py:class:`~landlab.components.FlowDirectorD8` (node)                   |
|                                                  | :py:class:`~landlab.components.ChiFinder` (node)                        | :py:class:`~landlab.components.FlowDirectorDINF` (node)                 |
|                                                  | :py:class:`~landlab.components.DrainageDensity` (node)                  | :py:class:`~landlab.components.FlowDirectorMFD` (node)                  |
|                                                  | :py:class:`~landlab.components.ErosionDeposition` (node)                | :py:class:`~landlab.components.FlowDirectorSteepest` (node)             |
|                                                  | :py:class:`~landlab.components.FastscapeEroder` (node)                  | :py:class:`~landlab.components.LakeMapperBarnes` (node)                 |
|                                                  | :py:class:`~landlab.components.HackCalculator` (node)                   |                                                                         |
|                                                  | :py:class:`~landlab.components.LakeMapperBarnes` (node)                 |                                                                         |
|                                                  | :py:class:`~landlab.components.LateralEroder` (node)                    |                                                                         |
|                                                  | :py:class:`~landlab.components.SedDepEroder` (node)                     |                                                                         |
|                                                  | :py:class:`~landlab.components.Space` (node)                            |                                                                         |
|                                                  | :py:class:`~landlab.components.SteepnessFinder` (node)                  |                                                                         |
|                                                  | :py:class:`~landlab.components.StreamPowerEroder` (node)                |                                                                         |
|                                                  | :py:class:`~landlab.components.StreamPowerSmoothThresholdEroder` (node) |                                                                         |
|                                                  | :py:class:`~landlab.components.TransportLengthHillslopeDiffuser` (node) |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| flow__receiver_proportions                       |                                                                         | :py:class:`~landlab.components.FlowDirectorDINF` (node)                 |
|                                                  |                                                                         | :py:class:`~landlab.components.FlowDirectorMFD` (node)                  |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| flow__sink_flag                                  | :py:class:`~landlab.components.LakeMapperBarnes` (node)                 | :py:class:`~landlab.components.FlowDirectorD8` (node)                   |
|                                                  |                                                                         | :py:class:`~landlab.components.FlowDirectorDINF` (node)                 |
|                                                  |                                                                         | :py:class:`~landlab.components.FlowDirectorMFD` (node)                  |
|                                                  |                                                                         | :py:class:`~landlab.components.FlowDirectorSteepest` (node)             |
|                                                  |                                                                         | :py:class:`~landlab.components.LakeMapperBarnes` (node)                 |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| flow__upstream_node_order                        | :py:class:`~landlab.components.ChiFinder` (node)                        | :py:class:`~landlab.components.FlowAccumulator` (node)                  |
|                                                  | :py:class:`~landlab.components.DrainageDensity` (node)                  | :py:class:`~landlab.components.LakeMapperBarnes` (node)                 |
|                                                  | :py:class:`~landlab.components.ErosionDeposition` (node)                | :py:class:`~landlab.components.LossyFlowAccumulator` (node)             |
|                                                  | :py:class:`~landlab.components.FastscapeEroder` (node)                  |                                                                         |
|                                                  | :py:class:`~landlab.components.HackCalculator` (node)                   |                                                                         |
|                                                  | :py:class:`~landlab.components.LakeMapperBarnes` (node)                 |                                                                         |
|                                                  | :py:class:`~landlab.components.LateralEroder` (node)                    |                                                                         |
|                                                  | :py:class:`~landlab.components.SedDepEroder` (node)                     |                                                                         |
|                                                  | :py:class:`~landlab.components.Space` (node)                            |                                                                         |
|                                                  | :py:class:`~landlab.components.SteepnessFinder` (node)                  |                                                                         |
|                                                  | :py:class:`~landlab.components.StreamPowerEroder` (node)                |                                                                         |
|                                                  | :py:class:`~landlab.components.StreamPowerSmoothThresholdEroder` (node) |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| fracture_at_node                                 |                                                                         | :py:class:`~landlab.components.FractureGridGenerator` (node)            |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| groundwater__specific_discharge                  |                                                                         | :py:class:`~landlab.components.GroundwaterDupuitPercolator` (link)      |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| groundwater__velocity                            |                                                                         | :py:class:`~landlab.components.GroundwaterDupuitPercolator` (link)      |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| hillslope_sediment__unit_volume_flux             |                                                                         | :py:class:`~landlab.components.LinearDiffuser` (link)                   |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| hydraulic__gradient                              |                                                                         | :py:class:`~landlab.components.GroundwaterDupuitPercolator` (link)      |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| is_pit                                           |                                                                         | :py:class:`~landlab.components.DepressionFinderAndRouter` (node)        |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| landslide__probability_of_failure                |                                                                         | :py:class:`~landlab.components.LandslideProbability` (node)             |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| lateral_erosion__depth_increment                 |                                                                         | :py:class:`~landlab.components.LateralEroder` (node)                    |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| lithosphere__increment_of_overlying_pressure     | :py:class:`~landlab.components.Flexure1D` (node)                        |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| lithosphere__overlying_pressure_increment        | :py:class:`~landlab.components.Flexure` (node)                          |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| lithosphere_surface__elevation_increment         |                                                                         | :py:class:`~landlab.components.Flexure` (node)                          |
|                                                  |                                                                         | :py:class:`~landlab.components.gFlex` (node)                            |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| lithosphere_surface__increment_of_elevation      |                                                                         | :py:class:`~landlab.components.Flexure1D` (node)                        |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| plant__age                                       |                                                                         | :py:class:`~landlab.components.VegCA` (cell)                            |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| plant__live_index                                |                                                                         | :py:class:`~landlab.components.VegCA` (cell)                            |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| radiation__incoming_shortwave_flux               |                                                                         | :py:class:`~landlab.components.PotentialEvapotranspiration` (cell)      |
|                                                  |                                                                         | :py:class:`~landlab.components.Radiation` (cell)                        |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| radiation__net_flux                              |                                                                         | :py:class:`~landlab.components.PotentialEvapotranspiration` (cell)      |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| radiation__net_longwave_flux                     |                                                                         | :py:class:`~landlab.components.PotentialEvapotranspiration` (cell)      |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| radiation__net_shortwave_flux                    |                                                                         | :py:class:`~landlab.components.PotentialEvapotranspiration` (cell)      |
|                                                  |                                                                         | :py:class:`~landlab.components.Radiation` (cell)                        |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| radiation__ratio_to_flat_surface                 | :py:class:`~landlab.components.PotentialEvapotranspiration` (cell)      | :py:class:`~landlab.components.Radiation` (cell)                        |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| rainfall__daily_depth                            | :py:class:`~landlab.components.SoilMoisture` (cell)                     |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| rainfall__flux                                   |                                                                         | :py:class:`~landlab.components.PrecipitationDistribution` (grid)        |
|                                                  |                                                                         | :py:class:`~landlab.components.SpatialPrecipitationDistribution` (node) |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| rainfall__total_depth_per_year                   |                                                                         | :py:class:`~landlab.components.SpatialPrecipitationDistribution` (node) |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| sediment__deposition_coeff                       |                                                                         | :py:class:`~landlab.components.TransportLengthHillslopeDiffuser` (node) |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| sediment__deposition_rate                        |                                                                         | :py:class:`~landlab.components.TransportLengthHillslopeDiffuser` (node) |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| sediment__discharge_in                           | :py:class:`~landlab.components.DischargeDiffuser` (node)                |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| sediment__erosion_rate                           |                                                                         | :py:class:`~landlab.components.TransportLengthHillslopeDiffuser` (node) |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| sediment__flux                                   |                                                                         | :py:class:`~landlab.components.ErosionDeposition` (node)                |
|                                                  |                                                                         | :py:class:`~landlab.components.LateralEroder` (node)                    |
|                                                  |                                                                         | :py:class:`~landlab.components.Space` (node)                            |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| sediment__flux_in                                |                                                                         | :py:class:`~landlab.components.TransportLengthHillslopeDiffuser` (node) |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| sediment__flux_out                               |                                                                         | :py:class:`~landlab.components.TransportLengthHillslopeDiffuser` (node) |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| sediment__transfer_rate                          |                                                                         | :py:class:`~landlab.components.TransportLengthHillslopeDiffuser` (node) |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| sediment_fill__depth                             |                                                                         | :py:class:`~landlab.components.SinkFiller` (node)                       |
|                                                  |                                                                         | :py:class:`~landlab.components.SinkFillerBarnes` (node)                 |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| slope_coefficient                                | :py:class:`~landlab.components.DrainageDensity` (node)                  |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| slope_exponent                                   | :py:class:`~landlab.components.DrainageDensity` (node)                  |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| soil__density                                    | :py:class:`~landlab.components.LandslideProbability` (node)             |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| soil__depth                                      | :py:class:`~landlab.components.DepthDependentDiffuser` (node)           | :py:class:`~landlab.components.DepthDependentDiffuser` (node)           |
|                                                  | :py:class:`~landlab.components.DepthDependentTaylorDiffuser` (node)     | :py:class:`~landlab.components.DepthDependentTaylorDiffuser` (node)     |
|                                                  | :py:class:`~landlab.components.ExponentialWeatherer` (node)             | :py:class:`~landlab.components.Space` (node)                            |
|                                                  | :py:class:`~landlab.components.Space` (node)                            |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| soil__flux                                       |                                                                         | :py:class:`~landlab.components.DepthDependentDiffuser` (link)           |
|                                                  |                                                                         | :py:class:`~landlab.components.DepthDependentTaylorDiffuser` (link)     |
|                                                  |                                                                         | :py:class:`~landlab.components.TaylorNonLinearDiffuser` (link)          |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| soil__internal_friction_angle                    | :py:class:`~landlab.components.LandslideProbability` (node)             |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| soil__maximum_total_cohesion                     | :py:class:`~landlab.components.LandslideProbability` (node)             |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| soil__mean_relative_wetness                      |                                                                         | :py:class:`~landlab.components.LandslideProbability` (node)             |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| soil__minimum_total_cohesion                     | :py:class:`~landlab.components.LandslideProbability` (node)             |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| soil__mode_total_cohesion                        | :py:class:`~landlab.components.LandslideProbability` (node)             |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| soil__probability_of_saturation                  |                                                                         | :py:class:`~landlab.components.LandslideProbability` (node)             |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| soil__saturated_hydraulic_conductivity           | :py:class:`~landlab.components.LandslideProbability` (node)             |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| soil__thickness                                  | :py:class:`~landlab.components.LandslideProbability` (node)             |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| soil__transmissivity                             | :py:class:`~landlab.components.LandslideProbability` (node)             |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| soil_moisture__initial_saturation_fraction       | :py:class:`~landlab.components.SoilMoisture` (cell)                     |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| soil_moisture__root_zone_leakage                 |                                                                         | :py:class:`~landlab.components.SoilMoisture` (cell)                     |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| soil_moisture__saturation_fraction               |                                                                         | :py:class:`~landlab.components.SoilMoisture` (cell)                     |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| soil_production__rate                            | :py:class:`~landlab.components.DepthDependentDiffuser` (node)           | :py:class:`~landlab.components.ExponentialWeatherer` (node)             |
|                                                  | :py:class:`~landlab.components.DepthDependentTaylorDiffuser` (node)     |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| soil_water_infiltration__depth                   | :py:class:`~landlab.components.SoilInfiltrationGreenAmpt` (node)        | :py:class:`~landlab.components.SoilInfiltrationGreenAmpt` (node)        |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| surface__evapotranspiration                      | :py:class:`~landlab.components.Vegetation` (cell)                       | :py:class:`~landlab.components.SoilMoisture` (cell)                     |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| surface__potential_evapotranspiration_30day_mean | :py:class:`~landlab.components.Vegetation` (cell)                       |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| surface__potential_evapotranspiration_rate       | :py:class:`~landlab.components.SoilMoisture` (cell)                     | :py:class:`~landlab.components.PotentialEvapotranspiration` (cell)      |
|                                                  | :py:class:`~landlab.components.Vegetation` (cell)                       |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| surface__runoff                                  |                                                                         | :py:class:`~landlab.components.SoilMoisture` (cell)                     |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| surface_load__stress                             | :py:class:`~landlab.components.gFlex` (node)                            |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| surface_to_channel__minimum_distance             |                                                                         | :py:class:`~landlab.components.DrainageDensity` (node)                  |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| surface_water__depth                             | :py:class:`~landlab.components.DepthSlopeProductErosion` (node)         | :py:class:`~landlab.components.KinwaveImplicitOverlandFlow` (node)      |
|                                                  | :py:class:`~landlab.components.OverlandFlow` (node)                     | :py:class:`~landlab.components.KinwaveOverlandFlowModel` (node)         |
|                                                  | :py:class:`~landlab.components.OverlandFlowBates` (node)                | :py:class:`~landlab.components.OverlandFlow` (node)                     |
|                                                  | :py:class:`~landlab.components.SoilInfiltrationGreenAmpt` (node)        | :py:class:`~landlab.components.OverlandFlowBates` (node)                |
|                                                  |                                                                         | :py:class:`~landlab.components.PotentialityFlowRouter` (node)           |
|                                                  |                                                                         | :py:class:`~landlab.components.SoilInfiltrationGreenAmpt` (node)        |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| surface_water__discharge                         | :py:class:`~landlab.components.DetachmentLtdErosion` (node)             | :py:class:`~landlab.components.DischargeDiffuser` (node)                |
|                                                  | :py:class:`~landlab.components.ErosionDeposition` (node)                | :py:class:`~landlab.components.FlowAccumulator` (node)                  |
|                                                  | :py:class:`~landlab.components.LakeMapperBarnes` (node)                 | :py:class:`~landlab.components.LakeMapperBarnes` (node)                 |
|                                                  | :py:class:`~landlab.components.Space` (node)                            | :py:class:`~landlab.components.LossyFlowAccumulator` (node)             |
|                                                  |                                                                         | :py:class:`~landlab.components.OverlandFlow` (link)                     |
|                                                  |                                                                         | :py:class:`~landlab.components.OverlandFlowBates` (link)                |
|                                                  |                                                                         | :py:class:`~landlab.components.PotentialityFlowRouter` (node)           |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| surface_water__discharge_loss                    |                                                                         | :py:class:`~landlab.components.LossyFlowAccumulator` (node)             |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| surface_water__specific_discharge                |                                                                         | :py:class:`~landlab.components.GroundwaterDupuitPercolator` (node)      |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| surface_water_inflow__discharge                  |                                                                         | :py:class:`~landlab.components.KinwaveImplicitOverlandFlow` (node)      |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| topographic__elevation                           | :py:class:`~landlab.components.ChiFinder` (node)                        | :py:class:`~landlab.components.DepthDependentDiffuser` (node)           |
|                                                  | :py:class:`~landlab.components.DepressionFinderAndRouter` (node)        | :py:class:`~landlab.components.DepthDependentTaylorDiffuser` (node)     |
|                                                  | :py:class:`~landlab.components.DepthDependentDiffuser` (node)           | :py:class:`~landlab.components.DepthSlopeProductErosion` (node)         |
|                                                  | :py:class:`~landlab.components.DepthDependentTaylorDiffuser` (node)     | :py:class:`~landlab.components.DetachmentLtdErosion` (node)             |
|                                                  | :py:class:`~landlab.components.DepthSlopeProductErosion` (node)         | :py:class:`~landlab.components.DischargeDiffuser` (node)                |
|                                                  | :py:class:`~landlab.components.DetachmentLtdErosion` (node)             | :py:class:`~landlab.components.ErosionDeposition` (node)                |
|                                                  | :py:class:`~landlab.components.DischargeDiffuser` (node)                | :py:class:`~landlab.components.FastscapeEroder` (node)                  |
|                                                  | :py:class:`~landlab.components.ErosionDeposition` (node)                | :py:class:`~landlab.components.gFlex` (node)                            |
|                                                  | :py:class:`~landlab.components.FastscapeEroder` (node)                  | :py:class:`~landlab.components.LakeMapperBarnes` (node)                 |
|                                                  | :py:class:`~landlab.components.FlowAccumulator` (node)                  | :py:class:`~landlab.components.LateralEroder` (node)                    |
|                                                  | :py:class:`~landlab.components.FlowDirectorD8` (node)                   | :py:class:`~landlab.components.LinearDiffuser` (node)                   |
|                                                  | :py:class:`~landlab.components.FlowDirectorDINF` (node)                 | :py:class:`~landlab.components.NormalFault` (node)                      |
|                                                  | :py:class:`~landlab.components.FlowDirectorMFD` (node)                  | :py:class:`~landlab.components.PerronNLDiffuse` (node)                  |
|                                                  | :py:class:`~landlab.components.FlowDirectorSteepest` (node)             | :py:class:`~landlab.components.SedDepEroder` (node)                     |
|                                                  | :py:class:`~landlab.components.GroundwaterDupuitPercolator` (node)      | :py:class:`~landlab.components.SinkFiller` (node)                       |
|                                                  | :py:class:`~landlab.components.HackCalculator` (node)                   | :py:class:`~landlab.components.SinkFillerBarnes` (node)                 |
|                                                  | :py:class:`~landlab.components.KinwaveImplicitOverlandFlow` (node)      | :py:class:`~landlab.components.Space` (node)                            |
|                                                  | :py:class:`~landlab.components.KinwaveOverlandFlowModel` (node)         | :py:class:`~landlab.components.StreamPowerEroder` (node)                |
|                                                  | :py:class:`~landlab.components.LakeMapperBarnes` (node)                 | :py:class:`~landlab.components.StreamPowerSmoothThresholdEroder` (node) |
|                                                  | :py:class:`~landlab.components.LateralEroder` (node)                    | :py:class:`~landlab.components.TaylorNonLinearDiffuser` (node)          |
|                                                  | :py:class:`~landlab.components.LinearDiffuser` (node)                   | :py:class:`~landlab.components.TransportLengthHillslopeDiffuser` (node) |
|                                                  | :py:class:`~landlab.components.LossyFlowAccumulator` (node)             |                                                                         |
|                                                  | :py:class:`~landlab.components.NormalFault` (node)                      |                                                                         |
|                                                  | :py:class:`~landlab.components.OverlandFlow` (node)                     |                                                                         |
|                                                  | :py:class:`~landlab.components.OverlandFlowBates` (node)                |                                                                         |
|                                                  | :py:class:`~landlab.components.PerronNLDiffuse` (node)                  |                                                                         |
|                                                  | :py:class:`~landlab.components.PotentialityFlowRouter` (node)           |                                                                         |
|                                                  | :py:class:`~landlab.components.Radiation` (node)                        |                                                                         |
|                                                  | :py:class:`~landlab.components.SedDepEroder` (node)                     |                                                                         |
|                                                  | :py:class:`~landlab.components.SinkFiller` (node)                       |                                                                         |
|                                                  | :py:class:`~landlab.components.SinkFillerBarnes` (node)                 |                                                                         |
|                                                  | :py:class:`~landlab.components.Space` (node)                            |                                                                         |
|                                                  | :py:class:`~landlab.components.SpatialPrecipitationDistribution` (node) |                                                                         |
|                                                  | :py:class:`~landlab.components.SteepnessFinder` (node)                  |                                                                         |
|                                                  | :py:class:`~landlab.components.StreamPowerEroder` (node)                |                                                                         |
|                                                  | :py:class:`~landlab.components.StreamPowerSmoothThresholdEroder` (node) |                                                                         |
|                                                  | :py:class:`~landlab.components.TaylorNonLinearDiffuser` (node)          |                                                                         |
|                                                  | :py:class:`~landlab.components.TransportLengthHillslopeDiffuser` (node) |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| topographic__gradient                            | :py:class:`~landlab.components.KinwaveOverlandFlowModel` (link)         | :py:class:`~landlab.components.KinwaveImplicitOverlandFlow` (link)      |
|                                                  |                                                                         | :py:class:`~landlab.components.LinearDiffuser` (link)                   |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| topographic__slope                               | :py:class:`~landlab.components.DepthSlopeProductErosion` (node)         | :py:class:`~landlab.components.DepthDependentDiffuser` (link)           |
|                                                  | :py:class:`~landlab.components.DetachmentLtdErosion` (node)             | :py:class:`~landlab.components.DepthDependentTaylorDiffuser` (link)     |
|                                                  | :py:class:`~landlab.components.LandslideProbability` (node)             | :py:class:`~landlab.components.TaylorNonLinearDiffuser` (link)          |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| topographic__specific_contributing_area          | :py:class:`~landlab.components.LandslideProbability` (node)             |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| topographic__steepest_slope                      | :py:class:`~landlab.components.ChiFinder` (node)                        | :py:class:`~landlab.components.FlowDirectorD8` (node)                   |
|                                                  | :py:class:`~landlab.components.DrainageDensity` (node)                  | :py:class:`~landlab.components.FlowDirectorDINF` (node)                 |
|                                                  | :py:class:`~landlab.components.ErosionDeposition` (node)                | :py:class:`~landlab.components.FlowDirectorMFD` (node)                  |
|                                                  | :py:class:`~landlab.components.LateralEroder` (node)                    | :py:class:`~landlab.components.FlowDirectorSteepest` (node)             |
|                                                  | :py:class:`~landlab.components.SedDepEroder` (node)                     |                                                                         |
|                                                  | :py:class:`~landlab.components.Space` (node)                            |                                                                         |
|                                                  | :py:class:`~landlab.components.SteepnessFinder` (node)                  |                                                                         |
|                                                  | :py:class:`~landlab.components.TransportLengthHillslopeDiffuser` (node) |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| vegetation__cover_fraction                       | :py:class:`~landlab.components.SoilMoisture` (cell)                     | :py:class:`~landlab.components.Vegetation` (cell)                       |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| vegetation__cumulative_water_stress              | :py:class:`~landlab.components.VegCA` (cell)                            |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| vegetation__dead_biomass                         |                                                                         | :py:class:`~landlab.components.Vegetation` (cell)                       |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| vegetation__dead_leaf_area_index                 |                                                                         | :py:class:`~landlab.components.Vegetation` (cell)                       |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| vegetation__live_biomass                         |                                                                         | :py:class:`~landlab.components.Vegetation` (cell)                       |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| vegetation__live_leaf_area_index                 | :py:class:`~landlab.components.SoilMoisture` (cell)                     | :py:class:`~landlab.components.Vegetation` (cell)                       |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| vegetation__plant_functional_type                | :py:class:`~landlab.components.SoilMoisture` (cell)                     |                                                                         |
|                                                  | :py:class:`~landlab.components.VegCA` (cell)                            |                                                                         |
|                                                  | :py:class:`~landlab.components.Vegetation` (cell)                       |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| vegetation__water_stress                         | :py:class:`~landlab.components.Vegetation` (cell)                       | :py:class:`~landlab.components.SoilMoisture` (cell)                     |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| volume__lateral_erosion                          |                                                                         | :py:class:`~landlab.components.LateralEroder` (node)                    |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| water__discharge_in                              | :py:class:`~landlab.components.DischargeDiffuser` (node)                |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| water__specific_discharge                        |                                                                         | :py:class:`~landlab.components.KinwaveOverlandFlowModel` (link)         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| water__unit_flux_in                              | :py:class:`~landlab.components.FlowAccumulator` (node)                  |                                                                         |
|                                                  | :py:class:`~landlab.components.LossyFlowAccumulator` (node)             |                                                                         |
|                                                  | :py:class:`~landlab.components.PotentialityFlowRouter` (node)           |                                                                         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| water__velocity                                  |                                                                         | :py:class:`~landlab.components.KinwaveOverlandFlowModel` (link)         |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| water_surface__gradient                          |                                                                         | :py:class:`~landlab.components.OverlandFlow` (link)                     |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| water_table__elevation                           |                                                                         | :py:class:`~landlab.components.GroundwaterDupuitPercolator` (node)      |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
| water_table__velocity                            |                                                                         | :py:class:`~landlab.components.GroundwaterDupuitPercolator` (node)      |
+--------------------------------------------------+-------------------------------------------------------------------------+-------------------------------------------------------------------------+
