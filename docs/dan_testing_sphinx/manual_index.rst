==============================================
An index for key functions in the Landlab grid
==============================================

This document is designed to present a more ordered list of methods available in
the Landlab grids. It is broken down primarily into functional groups according to 
the type of task you may want to achieve:

* `Grid creation`_
* `Create data in the grid fields`_
* `Resolve, project, and move data between grid element types`_
* `Access information about the grid`_
* `Access the grid geometry`_
* `Link coordinates and distances to nodes`_
* `Derive offsets, gradients, flux divergences, and steepest descents`_
* `Control boundary conditions`_
* `Manipulate arrays for plotting and display`_


At the moment, this document only covers the base grid and the raster. Voronoi support
coming soon.


Grid creation
=============

.. autoclass:: landlab.grid.raster.RasterModelGrid
    :members: __init__

These should be broken out as this format:
.. automethod:: landlab.grid.base.ModelGrid.resolve_values_on_links


Create data in the grid fields
==============================

    .. autoclass:: landlab.grid.base.ModelGrid
        :members: create_node_array_zeros, create_active_link_array_zeros, zeros, empty, \
                ones, set_nodata_nodes_to_inactive


Resolve, project, and move data between grid element types
==========================================================

    .. autoclass:: landlab.grid.base.ModelGrid
        :members: resolve_values_on_links, resolve_values_on_active_links, \
                assign_upslope_vals_to_active_links, max_of_link_end_node_values


Access information about the grid
=================================

    .. autoclass:: landlab.grid.base.ModelGrid
        :members: ndim, number_of_nodes, number_of_cells, number_of_links, number_of_faces, \
                number_of_active_nodes, number_of_active_cells, number_of_active_links, \
                number_of_active_faces, number_of_elements, axis_units, axis_name, cell_areas, \
                forced_cell_areas, active_link_length, link_length, min_active_link_length, \
                max_active_link_length, calculate_link_length, calculate_numbers_of_node_neighbors

    .. autoclass:: landlab.grid.raster.RasterModelGrid
        :members: shape, dx, get_grid_xdimension, get_grid_ydimension, number_of_interior_nodes, \
                number_of_nodes, number_of_node_columns, number_of_node_rows, node_spacing, \
                min_active_link_length, max_active_link_length, link_length, d8_active_links
                           

Access the grid geometry
========================

    .. autoclass:: landlab.grid.base.ModelGrid
        :members: node_index_at_cells, node_boundary_status, open_nodes, open_boundary_nodes, \
                closed_boundary_nodes, active_links, node_index_at_active_cells, \
                active_cell_index_at_nodes, active_cell_index, node_index_at_link_head, \
                node_index_at_link_tail, face_index_at_links, get_interior_nodes, get_node_status, \
                node_x, node_y, node_axis_coordinates, get_active_cell_node_ids, \
                get_active_link_connecting_node_pair, active_link_length, link_length, \
                calculate_numbers_of_node_neighbors, is_boundary, get_boundary_nodes
                

    .. autoclass:: landlab.grid.raster.RasterModelGrid
        :members: node_links, active_node_links, cell_faces, face_links, link_faces, \
                corner_nodes, get_neighbor_list, create_neighbor_list, has_boundary_neighbor, \
                get_diagonal_list, create_diagonal_list, is_interior, are_all_interior, \
                get_face_connecting_cell_pair, get_link_connecting_node_pair, \
                get_active_link_connecting_node_pair, top_edge_node_ids, bottom_edge_node_ids, \
                left_edge_node_ids, right_edge_node_ids, grid_coords_to_node_id


Link coordinates and distances to nodes
=======================================

    .. autoclass:: landlab.grid.base.ModelGrid
        :members: node_x, node_y, node_axis_coordinates, get_distances_of_nodes_to_point, \
                build_all_node_distances_azimuths_maps

    .. autoclass:: landlab.grid.raster.RasterModelGrid
        :members: is_point_on_grid, get_nodes_around_point, find_nearest_node, \
                grid_coords_to_node_id
                

Derive offsets, gradients, flux divergences, and steepest descents
==================================================================

    .. autoclass:: landlab.grid.base.ModelGrid
        :members: calculate_diff_at_links, calculate_diff_at_active_links, \
                calculate_gradients_at_links, calculate_gradients_at_active_links, \
                calculate_flux_divergence_at_active_cells, calculate_flux_divergence_at_nodes
    
    .. autoclass:: landlab.grid.raster.RasterModelGrid
        :members: calculate_gradient_across_cell_faces, calculate_gradient_across_cell_corners, \
                calculate_steepest_descent_across_cell_faces, calculate_steepest_descent_across_cell_corners, \
                calculate_steepest_descent_across_adjacent_cells, calculate_flux_divergence_at_nodes
                
                
Control boundary conditions
===========================

    .. autoclass:: landlab.grid.base.ModelGrid
        :members: node_boundary_status, open_nodes, open_boundary_nodes, closed_boundary_nodes, \
                get_node_status, set_nodata_nodes_to_inactive, is_boundary, get_boundary_nodes, \
                set_inactive_boundaries, set_inactive_nodes
        
    .. autoclass:: landlab.grid.raster.RasterModelGrid
        :members: set_inactive_boundaries, set_looped_boundaries, set_fixed_gradient_boundaries, \
                force_boundaries_from_gradients, has_boundary_neighbor, is_interior, are_all_interior, \
                get_boundary_code, top_edge_node_ids, bottom_edge_node_ids, left_edge_node_ids, \
                right_edge_node_ids
                
                
Manipulate arrays for plotting and display
==========================================

NB: these slope and aspect methods were devised for display purposes, and not 
intended or tested for quantitative use. But if you wish to explore their uses for such,
have at it!

    .. autoclass:: landlab.grid.base.ModelGrid
        :members: display_grid

    .. autoclass:: landlab.grid.raster.RasterModelGrid
        :members: node_vector_to_raster, cell_vector_to_raster, calculate_aspect_at_nodes_bestFitPlane, \
                calculate_slope_at_nodes_bestFitPlane, calculate_slope_aspect_at_nodes_bestFitPlane, \
                hillshade