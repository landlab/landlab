
=======================
Information about nodes
=======================


.. _NINF_ModelGrid:

Base class
----------

.. autosummary::
    :toctree: generated/

    :py:meth:`~landlab.grid.base.ModelGrid.active_adjacent_nodes_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.active_link_dirs_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.active_neighbors_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.adjacent_nodes_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.all_node_azimuths_map`
    :py:meth:`~landlab.grid.base.ModelGrid.all_node_distances_map`
    :py:meth:`~landlab.grid.base.ModelGrid.boundary_nodes`
    :py:meth:`~landlab.grid.base.ModelGrid.calc_distances_of_nodes_to_point`
    :py:meth:`~landlab.grid.base.ModelGrid.cell_area_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.cell_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.closed_boundary_nodes`
    :py:meth:`~landlab.grid.base.ModelGrid.core_nodes`
    :py:meth:`~landlab.grid.base.ModelGrid.downwind_links_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.fixed_gradient_boundary_nodes`
    :py:meth:`~landlab.grid.base.ModelGrid.fixed_value_boundary_nodes`
    :py:meth:`~landlab.grid.base.ModelGrid.link_at_node_is_downwind`
    :py:meth:`~landlab.grid.base.ModelGrid.link_at_node_is_upwind`
    :py:meth:`~landlab.grid.base.ModelGrid.link_dirs_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.links_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.neighbors_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.node_at_cell`
    :py:meth:`~landlab.grid.base.ModelGrid.node_at_core_cell`
    :py:meth:`~landlab.grid.base.ModelGrid.node_at_link_head`
    :py:meth:`~landlab.grid.base.ModelGrid.node_at_link_tail`
    :py:meth:`~landlab.grid.base.ModelGrid.node_axis_coordinates`
    :py:meth:`~landlab.grid.base.ModelGrid.node_has_boundary_neighbor`
    :py:meth:`~landlab.grid.base.ModelGrid.node_is_boundary`
    :py:meth:`~landlab.grid.base.ModelGrid.node_x`
    :py:meth:`~landlab.grid.base.ModelGrid.node_y`
    :py:meth:`~landlab.grid.base.ModelGrid.nodes`
    :py:meth:`~landlab.grid.base.ModelGrid.nodes_at_link`
    :py:meth:`~landlab.grid.base.ModelGrid.number_of_core_nodes`
    :py:meth:`~landlab.grid.base.ModelGrid.number_of_links_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.number_of_nodes`
    :py:meth:`~landlab.grid.base.ModelGrid.number_of_patches_present_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.open_boundary_nodes`
    :py:meth:`~landlab.grid.base.ModelGrid.patches_present_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.set_nodata_nodes_to_closed`
    :py:meth:`~landlab.grid.base.ModelGrid.set_nodata_nodes_to_fixed_gradient`
    :py:meth:`~landlab.grid.base.ModelGrid.status_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.unit_vector_sum_xcomponent_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.unit_vector_sum_ycomponent_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.upwind_links_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.x_of_node`
    :py:meth:`~landlab.grid.base.ModelGrid.xy_of_node`
    :py:meth:`~landlab.grid.base.ModelGrid.y_of_node`



.. _NINF_RasterModelGrid:

Raster
------

.. autosummary::
    :toctree: generated/

    :py:meth:`~landlab.grid.raster.RasterModelGrid.active_adjacent_nodes_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.active_link_dirs_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.active_neighbors_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.adjacent_nodes_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.all_node_azimuths_map`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.all_node_distances_map`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.boundary_nodes`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.calc_distances_of_nodes_to_point`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.cell_area_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.cell_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.closed_boundary_nodes`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.core_nodes`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.d8s_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.diagonal_adjacent_nodes_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.diagonals_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.downwind_links_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.find_nearest_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.fixed_gradient_boundary_nodes`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.fixed_value_boundary_nodes`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.grid_coords_to_node_id`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.link_at_node_is_downwind`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.link_at_node_is_upwind`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.link_dirs_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.links_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.neighbors_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.node_at_cell`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.node_at_core_cell`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.node_at_link_head`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.node_at_link_tail`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.node_axis_coordinates`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.node_has_boundary_neighbor`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.node_is_boundary`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.node_is_core`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.node_vector_to_raster`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.node_x`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.node_y`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.nodes`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.nodes_are_all_core`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.nodes_around_point`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.nodes_at_bottom_edge`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.nodes_at_corners_of_grid`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.nodes_at_edge`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.nodes_at_left_edge`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.nodes_at_link`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.nodes_at_patch`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.nodes_at_right_edge`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.nodes_at_top_edge`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.number_of_cell_columns`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.number_of_core_nodes`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.number_of_interior_nodes`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.number_of_links_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.number_of_node_columns`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.number_of_node_rows`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.number_of_nodes`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.number_of_patches_present_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.open_boundary_nodes`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.patches_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.patches_present_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.roll_nodes_ud`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.set_nodata_nodes_to_closed`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.set_nodata_nodes_to_fixed_gradient`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.shape`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.status_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.unit_vector_sum_xcomponent_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.unit_vector_sum_ycomponent_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.upwind_links_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.x_of_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.xy_of_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.y_of_node`



.. _NINF_VoronoiDelaunayGrid:

Irregular Voronoi-cell
----------------------

.. autosummary::
    :toctree: generated/

    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.active_adjacent_nodes_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.active_link_dirs_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.active_neighbors_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.adjacent_nodes_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.all_node_azimuths_map`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.all_node_distances_map`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.boundary_nodes`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.calc_distances_of_nodes_to_point`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.cell_area_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.cell_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.closed_boundary_nodes`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.core_nodes`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.downwind_links_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.fixed_gradient_boundary_nodes`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.fixed_value_boundary_nodes`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.link_at_node_is_downwind`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.link_at_node_is_upwind`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.link_dirs_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.links_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.neighbors_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.node_at_cell`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.node_at_core_cell`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.node_at_link_head`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.node_at_link_tail`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.node_axis_coordinates`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.node_has_boundary_neighbor`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.node_is_boundary`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.node_x`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.node_y`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.nodes`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.nodes_at_link`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.nodes_at_patch`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_core_nodes`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_links_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_nodes`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_patches_present_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.open_boundary_nodes`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.patches_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.patches_present_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.set_nodata_nodes_to_closed`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.set_nodata_nodes_to_fixed_gradient`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.status_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.unit_vector_sum_xcomponent_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.unit_vector_sum_ycomponent_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.upwind_links_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.x_of_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.xy_of_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.y_of_node`



.. _NINF_HexModelGrid:

Hexagonal
---------

.. autosummary::
    :toctree: generated/

    :py:meth:`~landlab.grid.hex.HexModelGrid.active_adjacent_nodes_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.active_link_dirs_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.active_neighbors_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.adjacent_nodes_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.all_node_azimuths_map`
    :py:meth:`~landlab.grid.hex.HexModelGrid.all_node_distances_map`
    :py:meth:`~landlab.grid.hex.HexModelGrid.boundary_nodes`
    :py:meth:`~landlab.grid.hex.HexModelGrid.calc_distances_of_nodes_to_point`
    :py:meth:`~landlab.grid.hex.HexModelGrid.cell_area_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.cell_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.closed_boundary_nodes`
    :py:meth:`~landlab.grid.hex.HexModelGrid.core_nodes`
    :py:meth:`~landlab.grid.hex.HexModelGrid.downwind_links_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.fixed_gradient_boundary_nodes`
    :py:meth:`~landlab.grid.hex.HexModelGrid.fixed_value_boundary_nodes`
    :py:meth:`~landlab.grid.hex.HexModelGrid.link_at_node_is_downwind`
    :py:meth:`~landlab.grid.hex.HexModelGrid.link_at_node_is_upwind`
    :py:meth:`~landlab.grid.hex.HexModelGrid.link_dirs_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.links_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.neighbors_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.node_at_cell`
    :py:meth:`~landlab.grid.hex.HexModelGrid.node_at_core_cell`
    :py:meth:`~landlab.grid.hex.HexModelGrid.node_at_link_head`
    :py:meth:`~landlab.grid.hex.HexModelGrid.node_at_link_tail`
    :py:meth:`~landlab.grid.hex.HexModelGrid.node_axis_coordinates`
    :py:meth:`~landlab.grid.hex.HexModelGrid.node_has_boundary_neighbor`
    :py:meth:`~landlab.grid.hex.HexModelGrid.node_is_boundary`
    :py:meth:`~landlab.grid.hex.HexModelGrid.node_x`
    :py:meth:`~landlab.grid.hex.HexModelGrid.node_y`
    :py:meth:`~landlab.grid.hex.HexModelGrid.nodes`
    :py:meth:`~landlab.grid.hex.HexModelGrid.nodes_at_bottom_edge`
    :py:meth:`~landlab.grid.hex.HexModelGrid.nodes_at_left_edge`
    :py:meth:`~landlab.grid.hex.HexModelGrid.nodes_at_link`
    :py:meth:`~landlab.grid.hex.HexModelGrid.nodes_at_patch`
    :py:meth:`~landlab.grid.hex.HexModelGrid.nodes_at_right_edge`
    :py:meth:`~landlab.grid.hex.HexModelGrid.nodes_at_top_edge`
    :py:meth:`~landlab.grid.hex.HexModelGrid.number_of_core_nodes`
    :py:meth:`~landlab.grid.hex.HexModelGrid.number_of_links_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.number_of_node_columns`
    :py:meth:`~landlab.grid.hex.HexModelGrid.number_of_node_rows`
    :py:meth:`~landlab.grid.hex.HexModelGrid.number_of_nodes`
    :py:meth:`~landlab.grid.hex.HexModelGrid.number_of_patches_present_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.open_boundary_nodes`
    :py:meth:`~landlab.grid.hex.HexModelGrid.patches_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.patches_present_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.set_nodata_nodes_to_closed`
    :py:meth:`~landlab.grid.hex.HexModelGrid.set_nodata_nodes_to_fixed_gradient`
    :py:meth:`~landlab.grid.hex.HexModelGrid.status_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.unit_vector_sum_xcomponent_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.unit_vector_sum_ycomponent_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.upwind_links_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.x_of_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.xy_of_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.y_of_node`



.. _NINF_RadialModelGrid:

Radial
------

.. autosummary::
    :toctree: generated/

    :py:meth:`~landlab.grid.radial.RadialModelGrid.active_adjacent_nodes_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.active_link_dirs_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.active_neighbors_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.adjacent_nodes_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.all_node_azimuths_map`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.all_node_distances_map`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.boundary_nodes`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.calc_distances_of_nodes_to_point`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.cell_area_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.cell_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.closed_boundary_nodes`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.core_nodes`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.downwind_links_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.fixed_gradient_boundary_nodes`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.fixed_value_boundary_nodes`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.link_at_node_is_downwind`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.link_at_node_is_upwind`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.link_dirs_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.links_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.neighbors_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.node_at_cell`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.node_at_core_cell`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.node_at_link_head`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.node_at_link_tail`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.node_axis_coordinates`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.node_has_boundary_neighbor`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.node_is_boundary`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.node_x`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.node_y`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.nodes`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.nodes_at_link`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.nodes_at_patch`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.number_of_core_nodes`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.number_of_links_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.number_of_nodes`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.number_of_nodes_in_shell`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.number_of_patches_present_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.open_boundary_nodes`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.patches_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.patches_present_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.radius_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.set_nodata_nodes_to_closed`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.set_nodata_nodes_to_fixed_gradient`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.status_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.unit_vector_sum_xcomponent_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.unit_vector_sum_ycomponent_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.upwind_links_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.x_of_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.xy_of_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.y_of_node`


