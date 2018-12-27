#! /usr/bin/env python

STANDARD_NAME = {
    "channel__bed_shear_stress": "channel_bottom_water__shear_stress",
    "channel__chi_index": None,
    "channel__depth": "channel_x-section__mean_of_depth",
    "channel__discharge": "channel_water_x-section__volume_flow_rate",
    "channel__steepness_index": "land_surface__steepness_index",
    "channel__width": "channel_x-section_top__width",
    "channel_sediment__relative_flux": "channel_water_sediment~suspended__time_max_normalized_volume_flux_ratio",
    "channel_sediment__volumetric_flux": "channel_water_sediment~suspended__volume_flux",  # (flux is area/time) or is it flow_rate?
    "channel_sediment__volumetric_transport_capacity": "channel_water_sediment~suspended__potential_volume_flow_rate",  # ?
    "depression__depth": "land_depression__max-fill_depth",
    "depression__outlet_node": "model_grid_node_land_depression_pour-point__id",  # necessary? model-specific? does anything need this?
    "drainage_area": "basin__total_contributing_area",  # model_grid_cell__total_contributing_area
    "flow__link_to_receiver_node": "model_grid_node_link~downstream__index",
    "flow__potential": None,
    "flow__receiver_node": None,
    "flow__sink_flag": None,  # model_*__flag boolean?
    "flow__upstream_node_order": None,
    "lithosphere__overlying_pressure_increment": "lithosphere_top_surface__increment_of_static_pressure",
    "lithosphere_surface__elevation_increment": "lithosphere_top_surface__increment_of_elevation",
    "lithosphere_surface__elevation_increment": "lithosphere_top_surface__increment_of_elevation",
    "plant__age": "plant__age",
    "plant__live_index": None,
    "radiation__incoming_shortwave_flux": "earth_surface_radiation~incoming~shortwave__energy_flux",
    "radiation__net_flux": "earth_surface_radiation~net~total__energy_flux",
    "radiation__net_longwave_flux": "earth_surface_radiation~net~longwave__energy_flux",
    "radiation__net_shortwave_flux": "earth_surface_radiation~net~shortwave__energy_flux",
    "radiation__ratio_to_flat_surface": "earth_surface_radiation~incoming~shortwave-to-flat_surface__energy_flux",
    "rainfall__daily_depth": "atmosphere_water__rainfall_volume_flux",
    "sediment_fill__depth": "land_surface_sediment__deposition_depth",
    "soil_moisture__cumulative_water_stress": "land_vegetation__time_integral_of_water_stress",
    "soil_moisture__initial_saturation_fraction": "soil_water__initial_saturated_volume_fraction",
    "soil_moisture__root_zone_leakage_rate": "soil_water_root-zone__volume_flux",
    "soil_moisture__saturation_fraction": "soil_water__saturated_volume_fraction",
    "soil_moisture__water_stress": "land_vegetation__water_stress",
    "surface__evapotranspiration_rate": "land_surface_water__evaptranspiration_volume_flux",
    "surface__potential_evapotranspiration_30day_mean": "land_surface_water__time_mean_of_potential_evaptranspiration_volume_flux",
    "surface__potential_evapotranspiration_rate": "land_surface_water__potential_evaptranspiration_volume_flux",
    "surface__runoff_rate": "land_surface_water__runoff_volume_flux",
    "surface_load__stress": "lithosphere_surface__normal_component_of_stress",
    "topographic__elevation": "land_surface__elevation",
    "topographic__gradient": "land_surface__slope",
    "topographic__slope": "land_surface__slope_angle",
    "topographic__steepest_slope": "model_grid_cell__max_of_d8_slope",
    "hillslope_sediment__unit_volume_flux": None,
    "vegetation__cover_fraction": "land_vegetation__area_fraction",
    "vegetation__dead_biomass": "land_vegetation~biomass~dead__volume_flux",
    "vegetation__dead_leaf_area_index": "land_vegetation~dead__leaf-area_index",
    "vegetation__live_biomass": "land_vegetation~biomass~live__volume_flux",
    "vegetation__live_leaf_area_index": "land_vegetation~live__leaf-area_index",
    "vegetation__plant_functional_type": None,
    "surface_water__depth": "land_surface_water__depth",
    "surface_water__discharge": "land_surface_water__volume_flow_rate",
    "water__unit_flux_in": "model_grid_cell_water~incoming__volume_flow_rate",
    "water_surface__gradient": "land_surface_water_surface__slope",
}


LANDLAB_NAME = dict((value, key) for key, value in STANDARD_NAME.items() if key)
