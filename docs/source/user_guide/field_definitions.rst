.. _standard_name_definitions:


Landlab Standard Name Definitions
=================================

+--------------------------------------------------+------------------------------------------+
| Field Name                                       | Definition                               |
+==================================================+==========================================+
| aquifer__thickness                               | thickness of saturated zone              |
+--------------------------------------------------+------------------------------------------+
| aquifer_base__elevation                          | elevation of impervious layer            |
+--------------------------------------------------+------------------------------------------+
| aquifer_base__gradient                           | gradient of the aquifer base in the link |
|                                                  | directio                                 |
+--------------------------------------------------+------------------------------------------+
| area_coefficient                                 | Area coefficient to define channels.     |
+--------------------------------------------------+------------------------------------------+
| area_exponent                                    | Area exponent to define channels.        |
+--------------------------------------------------+------------------------------------------+
| bedrock__elevation                               | elevation of the bedrock surface         |
+--------------------------------------------------+------------------------------------------+
| channel__bed_shear_stress                        | Shear exerted on the bed of the channel, |
|                                                  | assuming all discharge travels along a   |
|                                                  | single, self-formed channel              |
+--------------------------------------------------+------------------------------------------+
| channel__chi_index                               | the local steepness index                |
+--------------------------------------------------+------------------------------------------+
| channel__depth                                   | Depth of the a single channel carrying   |
|                                                  | all runoff through the node              |
+--------------------------------------------------+------------------------------------------+
| channel__discharge                               | Volumetric water flux of the a single    |
|                                                  | channel carrying all runoff through the  |
|                                                  | node                                     |
+--------------------------------------------------+------------------------------------------+
| channel__mask                                    | Logical map of at which grid nodes       |
|                                                  | channels are present                     |
+--------------------------------------------------+------------------------------------------+
| channel__steepness_index                         | the local steepness index                |
+--------------------------------------------------+------------------------------------------+
| channel__width                                   | Width of the a single channel carrying   |
|                                                  | all runoff through the node              |
+--------------------------------------------------+------------------------------------------+
| channel_sediment__relative_flux                  | The fluvial_sediment_flux_into_node      |
|                                                  | divided by the                           |
|                                                  | fluvial_sediment_transport_capacity      |
+--------------------------------------------------+------------------------------------------+
| channel_sediment__volumetric_flux                | Total volumetric fluvial sediment flux   |
|                                                  | brought into the node from upstream      |
+--------------------------------------------------+------------------------------------------+
| channel_sediment__volumetric_transport_capacity  | Volumetric transport capacity of a       |
|                                                  | channel carrying all runoff through the  |
|                                                  | node, assuming the Meyer-Peter Muller    |
|                                                  | transport equation                       |
+--------------------------------------------------+------------------------------------------+
| channelization_threshold                         | Channelization threshold for use with    |
|                                                  | area and slope coefficients and          |
|                                                  | exponents.                               |
+--------------------------------------------------+------------------------------------------+
| depression__depth                                | Depth of depression below its spillway   |
|                                                  | point                                    |
+--------------------------------------------------+------------------------------------------+
| depression__outlet_node                          | If a depression, the id of the outlet    |
|                                                  | node for that depression, otherwise      |
|                                                  | BAD_INDEX_VALUE                          |
+--------------------------------------------------+------------------------------------------+
| distance_to_divide                               | Distance from drainage divide.           |
+--------------------------------------------------+------------------------------------------+
| drainage_area                                    | Upstream accumulated surface area        |
|                                                  | contributing to the node's discharge     |
+--------------------------------------------------+------------------------------------------+
| flood_status_code                                | Map of flood status (_PIT,               |
|                                                  | _CURRENT_LAKE, _UNFLOODED, or _FLOODED). |
+--------------------------------------------------+------------------------------------------+
| flow__data_structure_delta                       | Node array containing the elements       |
|                                                  | delta[1:] of the data structure 'delta'  |
|                                                  | used for construction of the downstream- |
|                                                  | to-upstream node array                   |
+--------------------------------------------------+------------------------------------------+
| flow__link_direction                             | Direction of flow on link. A value of -1 |
|                                                  | indicates that water flow goes from head |
|                                                  | node to tail node, while a value of 1    |
|                                                  | indicates that water flow goes from tail |
|                                                  | node to head node.                       |
+--------------------------------------------------+------------------------------------------+
| flow__link_to_receiver_node                      | ID of link downstream of each node,      |
|                                                  | which carries the discharge              |
+--------------------------------------------------+------------------------------------------+
| flow__potential                                  | Value of the hypothetical field 'K',     |
|                                                  | used to force water flux to flow         |
|                                                  | downhill                                 |
+--------------------------------------------------+------------------------------------------+
| flow__receiver_node                              | Node array of receivers (node that       |
|                                                  | receives flow from current node)         |
+--------------------------------------------------+------------------------------------------+
| flow__receiver_proportions                       | Node array of proportion of flow sent to |
|                                                  | each receiver.                           |
+--------------------------------------------------+------------------------------------------+
| flow__sink_flag                                  | Boolean array, True at local lows        |
+--------------------------------------------------+------------------------------------------+
| flow__upstream_node_order                        | Node array containing downstream-to-     |
|                                                  | upstream ordered list of node IDs        |
+--------------------------------------------------+------------------------------------------+
| fracture_at_node                                 | presence (1) or absence (0) of fracture  |
+--------------------------------------------------+------------------------------------------+
| groundwater__specific_discharge                  | discharge per width in link dir          |
+--------------------------------------------------+------------------------------------------+
| groundwater__velocity                            | velocity of groundwater in link          |
|                                                  | direction                                |
+--------------------------------------------------+------------------------------------------+
| hillslope_sediment__unit_volume_flux             | Volume flux per unit width along links   |
+--------------------------------------------------+------------------------------------------+
| hydraulic__gradient                              | gradient of water table in link          |
|                                                  | direction                                |
+--------------------------------------------------+------------------------------------------+
| is_pit                                           | Boolean flag indicating whether a node   |
|                                                  | is a pit.                                |
+--------------------------------------------------+------------------------------------------+
| landslide__probability_of_failure                | number of times FS is <=1 out of number  |
|                                                  | of iterations user selected              |
+--------------------------------------------------+------------------------------------------+
| lateral_erosion__depth_increment                 | Change in elevation at each node from    |
|                                                  | lateral erosion during time step         |
+--------------------------------------------------+------------------------------------------+
| lithosphere__increment_of_overlying_pressure     | Applied pressure to the lithosphere over |
|                                                  | a time step                              |
+--------------------------------------------------+------------------------------------------+
| lithosphere__overlying_pressure_increment        | Applied pressure to the lithosphere over |
|                                                  | a time step                              |
+--------------------------------------------------+------------------------------------------+
| lithosphere_surface__elevation_increment         | The change in elevation of the top of    |
|                                                  | the lithosphere (the land surface) in    |
|                                                  | one timestep                             |
+--------------------------------------------------+------------------------------------------+
| lithosphere_surface__increment_of_elevation      | The change in elevation of the top of    |
|                                                  | the lithosphere (the land surface) in    |
|                                                  | one timestep                             |
+--------------------------------------------------+------------------------------------------+
| plant__age                                       | Age of plant                             |
+--------------------------------------------------+------------------------------------------+
| plant__live_index                                | 1 - vegetation__cumulative_water_stress  |
+--------------------------------------------------+------------------------------------------+
| radiation__incoming_shortwave_flux               | total incident shortwave radiation over  |
|                                                  | the time step                            |
+--------------------------------------------------+------------------------------------------+
| radiation__net_flux                              | net total radiation over the time step   |
+--------------------------------------------------+------------------------------------------+
| radiation__net_longwave_flux                     | net incident longwave radiation over the |
|                                                  | time step                                |
+--------------------------------------------------+------------------------------------------+
| radiation__net_shortwave_flux                    | net incident shortwave radiation over    |
|                                                  | the time step                            |
+--------------------------------------------------+------------------------------------------+
| radiation__ratio_to_flat_surface                 | ratio of total incident shortwave        |
|                                                  | radiation on sloped surface to flat      |
|                                                  | surface                                  |
+--------------------------------------------------+------------------------------------------+
| rainfall__daily_depth                            | Rain in (mm) as a field, allowing        |
|                                                  | spatio-temporal soil moisture saturation |
|                                                  | analysis.                                |
+--------------------------------------------------+------------------------------------------+
| rainfall__flux                                   | Depth of water delivered per unit time   |
|                                                  | in each storm                            |
+--------------------------------------------------+------------------------------------------+
| rainfall__total_depth_per_year                   | Depth of water delivered in total in     |
|                                                  | each model year                          |
+--------------------------------------------------+------------------------------------------+
| sediment__deposition_coeff                       | Fraction of incoming sediment that is    |
|                                                  | deposited on the node                    |
+--------------------------------------------------+------------------------------------------+
| sediment__deposition_rate                        | Deposition rate on node                  |
+--------------------------------------------------+------------------------------------------+
| sediment__discharge_in                           | Sediment discharge into a node.          |
+--------------------------------------------------+------------------------------------------+
| sediment__erosion_rate                           | Erosion rate on node                     |
+--------------------------------------------------+------------------------------------------+
| sediment__flux                                   | Sediment flux (volume per unit time of   |
|                                                  | sediment entering each node)             |
+--------------------------------------------------+------------------------------------------+
| sediment__flux_in                                | Incoming sediment rate on node (=qs/dx)  |
+--------------------------------------------------+------------------------------------------+
| sediment__flux_out                               | Outgoing sediment rate on node =         |
|                                                  | sediment eroded on node + sediment       |
|                                                  | transported across node from upstream    |
+--------------------------------------------------+------------------------------------------+
| sediment__transfer_rate                          | Rate of transferred sediment across a    |
|                                                  | node (incoming sediment - deposited      |
|                                                  | sediment on node)                        |
+--------------------------------------------------+------------------------------------------+
| sediment_fill__depth                             | Depth of sediment added at eachnode      |
+--------------------------------------------------+------------------------------------------+
| slope_coefficient                                | Slope coefficient to define channels.    |
+--------------------------------------------------+------------------------------------------+
| slope_exponent                                   | Slope exponent to define channels.       |
+--------------------------------------------------+------------------------------------------+
| soil__density                                    | wet bulk density of soil                 |
+--------------------------------------------------+------------------------------------------+
| soil__depth                                      | Depth of soil or weathered bedrock       |
+--------------------------------------------------+------------------------------------------+
| soil__flux                                       | flux of soil in direction of link        |
+--------------------------------------------------+------------------------------------------+
| soil__internal_friction_angle                    | critical angle just before failure due   |
|                                                  | to friction between particles            |
+--------------------------------------------------+------------------------------------------+
| soil__maximum_total_cohesion                     | maximum of combined root and soil        |
|                                                  | cohesion at node                         |
+--------------------------------------------------+------------------------------------------+
| soil__mean_relative_wetness                      | Indicator of soil wetness; relative      |
|                                                  | depth perched water table within the     |
|                                                  | soil layer                               |
+--------------------------------------------------+------------------------------------------+
| soil__minimum_total_cohesion                     | minimum of combined root and soil        |
|                                                  | cohesion at node                         |
+--------------------------------------------------+------------------------------------------+
| soil__mode_total_cohesion                        | mode of combined root and soil cohesion  |
|                                                  | at node                                  |
+--------------------------------------------------+------------------------------------------+
| soil__probability_of_saturation                  | number of times relative wetness is >=1  |
|                                                  | out of number of iterations user         |
|                                                  | selected                                 |
+--------------------------------------------------+------------------------------------------+
| soil__saturated_hydraulic_conductivity           | mode rate of water transmitted through   |
|                                                  | soil - provided if transmissivity is NOT |
|                                                  | provided to calculate tranmissivity      |
|                                                  | with soil depth                          |
+--------------------------------------------------+------------------------------------------+
| soil__thickness                                  | soil depth to restrictive layer          |
+--------------------------------------------------+------------------------------------------+
| soil__transmissivity                             | mode rate of water transmitted through a |
|                                                  | unit width of saturated soil - either    |
|                                                  | provided or calculated with Ksat and     |
|                                                  | soil depth                               |
+--------------------------------------------------+------------------------------------------+
| soil_moisture__initial_saturation_fraction       | initial                                  |
|                                                  | soil_moisture__saturation_fraction       |
+--------------------------------------------------+------------------------------------------+
| soil_moisture__root_zone_leakage                 | leakage of water into deeper portions of |
|                                                  | the soil not accessible to the plant     |
+--------------------------------------------------+------------------------------------------+
| soil_moisture__saturation_fraction               | relative volumetric water content        |
|                                                  | (theta) - limits=[0,1]                   |
+--------------------------------------------------+------------------------------------------+
| soil_production__rate                            | rate of soil production at nodes         |
+--------------------------------------------------+------------------------------------------+
| soil_water_infiltration__depth                   | Water column height above the surface    |
|                                                  | previously absorbed into the soil. Note  |
|                                                  | that this is NOT the actual depth of the |
|                                                  | wetted front, which also depends on the  |
|                                                  | porosity.                                |
+--------------------------------------------------+------------------------------------------+
| surface__evapotranspiration                      | actual sum of evaporation and plant      |
|                                                  | transpiration                            |
+--------------------------------------------------+------------------------------------------+
| surface__potential_evapotranspiration_30day_mean | 30 day mean of                           |
|                                                  | surface__potential_evapotranspiration    |
+--------------------------------------------------+------------------------------------------+
| surface__potential_evapotranspiration_rate       | potential sum of evaporation and         |
|                                                  | potential transpiration                  |
+--------------------------------------------------+------------------------------------------+
| surface__runoff                                  | runoff from ground surface               |
+--------------------------------------------------+------------------------------------------+
| surface_load__stress                             | Magnitude of stress exerted by surface   |
|                                                  | load                                     |
+--------------------------------------------------+------------------------------------------+
| surface_to_channel__minimum_distance             | Distance from each node to the nearest   |
|                                                  | channel                                  |
+--------------------------------------------------+------------------------------------------+
| surface_water__depth                             | Depth of water on the surface            |
+--------------------------------------------------+------------------------------------------+
| surface_water__discharge                         | Volumetric discharge of surface water    |
+--------------------------------------------------+------------------------------------------+
| surface_water__discharge_loss                    | Total volume of water per second lost    |
|                                                  | during all flow out of the node          |
+--------------------------------------------------+------------------------------------------+
| surface_water__specific_discharge                | rate of seepage to surface               |
+--------------------------------------------------+------------------------------------------+
| surface_water_inflow__discharge                  | water volume inflow rate to the cell     |
|                                                  | around each node                         |
+--------------------------------------------------+------------------------------------------+
| topographic__elevation                           | Land surface topographic elevation       |
+--------------------------------------------------+------------------------------------------+
| topographic__gradient                            | Gradient of the ground surface           |
+--------------------------------------------------+------------------------------------------+
| topographic__slope                               | gradient of the ground surface           |
+--------------------------------------------------+------------------------------------------+
| topographic__specific_contributing_area          | specific contributing (upslope area/cell |
|                                                  | face ) that drains to node               |
+--------------------------------------------------+------------------------------------------+
| topographic__steepest_slope                      | The steepest *downhill* slope            |
+--------------------------------------------------+------------------------------------------+
| vegetation__cover_fraction                       | fraction of land covered by vegetation   |
+--------------------------------------------------+------------------------------------------+
| vegetation__cumulative_water_stress              | cumulative vegetation__water_stress over |
|                                                  | the growing season                       |
+--------------------------------------------------+------------------------------------------+
| vegetation__dead_biomass                         | weight of dead organic mass per unit     |
|                                                  | area - measured in terms of dry matter   |
+--------------------------------------------------+------------------------------------------+
| vegetation__dead_leaf_area_index                 | one-sided dead leaf area per unit ground |
|                                                  | surface area                             |
+--------------------------------------------------+------------------------------------------+
| vegetation__live_biomass                         | weight of green organic mass per unit    |
|                                                  | area - measured in terms of dry matter   |
+--------------------------------------------------+------------------------------------------+
| vegetation__live_leaf_area_index                 | one-sided green leaf area per unit       |
|                                                  | ground surface area                      |
+--------------------------------------------------+------------------------------------------+
| vegetation__plant_functional_type                | classification of plants (int), grass=0, |
|                                                  | shrub=1, tree=2, bare=3,                 |
|                                                  | shrub_seedling=4, tree_seedling=5        |
+--------------------------------------------------+------------------------------------------+
| vegetation__water_stress                         | parameter that represents nonlinear      |
|                                                  | effects of water deficit on plants       |
+--------------------------------------------------+------------------------------------------+
| volume__lateral_erosion                          | Array tracking volume eroded at each     |
|                                                  | node from lateral erosion                |
+--------------------------------------------------+------------------------------------------+
| water__discharge_in                              | Incoming water discharge at node.        |
+--------------------------------------------------+------------------------------------------+
| water__specific_discharge                        | flow discharge component in the          |
|                                                  | direction of the link                    |
+--------------------------------------------------+------------------------------------------+
| water__unit_flux_in                              | External volume water per area per time  |
|                                                  | input to each node (e.g., rainfall rate) |
+--------------------------------------------------+------------------------------------------+
| water__velocity                                  | flow velocity component in the direction |
|                                                  | of the link                              |
+--------------------------------------------------+------------------------------------------+
| water_surface__gradient                          | Downstream gradient of the water         |
|                                                  | surface.                                 |
+--------------------------------------------------+------------------------------------------+
| water_table__elevation                           | elevation of water table                 |
+--------------------------------------------------+------------------------------------------+
| water_table__velocity                            | rate of change of water table elevation  |
+--------------------------------------------------+------------------------------------------+
