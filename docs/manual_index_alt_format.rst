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


At the moment, the raster grid is much better developed than the Voronoi grid (and its
derived unstructured grids). Most calls to ModelGrid methods below will work for such
unstructured grids, but some may result in errors. Please report such missing
functionality through the github page if you need it urgently!


Grid creation
=============

*These methods are used to create grids.*

.. _RMG_link:
.. automethod:: landlab.grid.raster.RasterModelGrid.__init__
.. _VMG_link:
.. automethod:: landlab.grid.voronoi.VoronoiDelaunayGrid.__init__
.. _Radial_link:
.. automethod:: landlab.grid.radial.RadialModelGrid.__init__
.. _Hex_link:
.. automethod:: landlab.grid.hex.HexModelGrid.__init__


Create data in the grid fields
==============================

*These methods are used to associate data - e.g., elevations, discharge values, flow
depths - with the grid object. Several can also be used to create new arrays of lengths
related to numbers of grid elements (e.g., number of active nodes), but not to link them
to the grid.*
*Data stored inside the grid object allows you to pass this information between components
more easily. The data can be accessed by, e.g., mygrid.at_node('my_data_name'), or a
number of alternative methods (see below).*

# .. automethod:: landlab.grid.base.ModelGrid.create_active_link_array_zeros
# .. automethod:: landlab.grid.base.ModelGrid.create_node_array_zeros
.. automethod:: landlab.grid.base.ModelGrid.empty
.. automethod:: landlab.grid.base.ModelGrid.ones
.. automethod:: landlab.grid.base.ModelGrid.set_nodata_nodes_to_inactive
.. automethod:: landlab.grid.base.ModelGrid.zeros


DEM and NetCDF input/output
===========================

*The landlab/io folder contains the various methods that allow Landlab to ingest DEMs,
and to import and export NetCDF files. i/o with vtk files is also possible, but not
detailed here.*

.. autofunction:: landlab.io.esri_ascii.read_esri_ascii
.. autofunction:: landlab.io.netcdf.read.read_netcdf
.. autofunction:: landlab.io.netcdf.write.write_netcdf

Access data in the grid fields
==============================

*Once you've created the fields, these methods can be used to access and modify the data
stored in them. Many other methods are available, see the docstring of base.py.*

    **grid.at_node['my_data_name']**

    **grid['node']['my_data_name']**

    **grid.at_node.keys()**
        Get the names of the data fields already stored on nodes in the grid.

(see also entry for ModelGrid.create_node_array_zeros, and the docstrings of
the base.py module.)


Resolve, project, and move data between grid element types
==========================================================

*These methods allow manipulation of data between links and nodes. e.g., if I have some
data defined on the nodes, but I need to use it on the links, these methods may help.*

.. automethod:: landlab.grid.base.ModelGrid.assign_upslope_vals_to_active_links
.. automethod:: landlab.grid.base.ModelGrid.max_of_link_end_node_values
.. automethod:: landlab.grid.base.ModelGrid.resolve_values_on_links
.. automethod:: landlab.grid.base.ModelGrid.resolve_values_on_active_links


Access information about the grid
=================================

*These methods allow you to access descriptive data about the grid itself. Note that many
of them are properties, so are accessed like grid.my_property, not grid.my_method().*

.. autoattribute:: landlab.grid.base.ModelGrid.active_link_length
.. autoattribute:: landlab.grid.base.ModelGrid.axis_name
.. autoattribute:: landlab.grid.base.ModelGrid.axis_units
.. automethod:: landlab.grid.base.ModelGrid.calculate_numbers_of_node_neighbors
.. autoattribute:: landlab.grid.base.ModelGrid.cell_areas
.. autoattribute:: landlab.grid.base.ModelGrid.forced_cell_areas
.. autoattribute:: landlab.grid.base.ModelGrid.link_length
.. autoattribute:: landlab.grid.base.ModelGrid.ndim
.. autoattribute:: landlab.grid.base.ModelGrid.number_of_active_cells
.. autoattribute:: landlab.grid.base.ModelGrid.number_of_core_cells
.. autoattribute:: landlab.grid.base.ModelGrid.number_of_active_faces
.. autoattribute:: landlab.grid.base.ModelGrid.number_of_active_links
.. autoattribute:: landlab.grid.base.ModelGrid.number_of_active_nodes
.. autoattribute:: landlab.grid.base.ModelGrid.number_of_core_nodes
.. autoattribute:: landlab.grid.base.ModelGrid.number_of_cells
.. automethod:: landlab.grid.base.ModelGrid.number_of_elements
.. autoattribute:: landlab.grid.base.ModelGrid.number_of_faces
.. autoattribute:: landlab.grid.base.ModelGrid.number_of_links
.. autoattribute:: landlab.grid.base.ModelGrid.number_of_nodes

.. automethod:: landlab.grid.raster.RasterModelGrid.d8_active_links
.. autoattribute:: landlab.grid.raster.RasterModelGrid.dx
.. autoattribute:: landlab.grid.raster.RasterModelGrid.extent
.. autoattribute:: landlab.grid.raster.RasterModelGrid.link_length
.. autoattribute:: landlab.grid.raster.RasterModelGrid.node_spacing
.. autoattribute:: landlab.grid.raster.RasterModelGrid.number_of_interior_nodes
.. autoattribute:: landlab.grid.raster.RasterModelGrid.number_of_nodes
.. autoattribute:: landlab.grid.raster.RasterModelGrid.number_of_node_columns
.. autoattribute:: landlab.grid.raster.RasterModelGrid.number_of_node_rows
.. autoattribute:: landlab.grid.raster.RasterModelGrid.shape


Access the grid geometry
========================

*These methods give you data about the structural elements of the grid, i.e., nodes,
cells, links, faces. This might include their statuses (active/inactive; core/boundary;
open/closed; interior/perimeter), their lengths and sizes, and their positions.*

.. autoattribute:: landlab.grid.base.ModelGrid.node_inlink_matrix
.. autoattribute:: landlab.grid.base.ModelGrid.node_outlink_matrix
.. autoattribute:: landlab.grid.base.ModelGrid.active_links
.. autoattribute:: landlab.grid.base.ModelGrid.active_link_length
.. automethod:: landlab.grid.base.ModelGrid.calculate_numbers_of_node_neighbors
.. autoattribute:: landlab.grid.base.ModelGrid.closed_boundary_nodes
.. autoattribute:: landlab.grid.base.ModelGrid.cell_index
.. autoattribute:: landlab.grid.base.ModelGrid.cell_index_at_nodes
.. autoattribute:: landlab.grid.base.ModelGrid.core_cell_index
.. autoattribute:: landlab.grid.base.ModelGrid.core_cell_index_at_nodes
.. autoattribute:: landlab.grid.base.ModelGrid.core_nodes
.. autoattribute:: landlab.grid.base.ModelGrid.face_index_at_links
.. automethod:: landlab.grid.base.ModelGrid.active_link_connecting_node_pair
.. autoattribute:: landlab.grid.base.ModelGrid.link_length
.. automethod:: landlab.grid.base.ModelGrid.is_boundary
.. automethod:: landlab.grid.base.ModelGrid.node_axis_coordinates
.. autoattribute:: landlab.grid.base.ModelGrid.node_boundary_status
.. autoattribute:: landlab.grid.base.ModelGrid.node_index_at_core_cells
.. autoattribute:: landlab.grid.base.ModelGrid.node_index_at_cells
.. autoattribute:: landlab.grid.base.ModelGrid.node_index_at_link_head
.. autoattribute:: landlab.grid.base.ModelGrid.node_index_at_link_tail
.. autoattribute:: landlab.grid.base.ModelGrid.node_x
.. autoattribute:: landlab.grid.base.ModelGrid.node_y
.. autoattribute:: landlab.grid.base.ModelGrid.open_boundary_nodes

.. automethod:: landlab.grid.raster.RasterModelGrid.are_all_core
.. autoattribute:: landlab.grid.raster.RasterModelGrid.corner_nodes
.. automethod:: landlab.grid.raster.RasterModelGrid.create_diagonal_list
.. automethod:: landlab.grid.raster.RasterModelGrid.create_neighbor_list
.. automethod:: landlab.grid.raster.RasterModelGrid.active_link_connecting_node_pair
.. automethod:: landlab.grid.raster.RasterModelGrid.get_diagonal_list
.. automethod:: landlab.grid.raster.RasterModelGrid.get_face_connecting_cell_pair
.. automethod:: landlab.grid.raster.RasterModelGrid.get_link_connecting_node_pair
.. automethod:: landlab.grid.raster.RasterModelGrid.grid_coords_to_node_id
.. automethod:: landlab.grid.raster.RasterModelGrid.has_boundary_neighbor
.. automethod:: landlab.grid.raster.RasterModelGrid.is_core


Link coordinates and distances to nodes
=======================================

*These methods are focused on the specifically spatial relationships between grid
elements. e.g., Where in x,y space is my element? How far is it from one node to another?*

.. automethod:: landlab.grid.base.ModelGrid.build_all_node_distances_azimuths_maps
.. automethod:: landlab.grid.base.ModelGrid.get_distances_of_nodes_to_point
.. automethod:: landlab.grid.base.ModelGrid.node_axis_coordinates
.. autoattribute:: landlab.grid.base.ModelGrid.node_x
.. autoattribute:: landlab.grid.base.ModelGrid.node_y

.. automethod:: landlab.grid.raster.RasterModelGrid.find_nearest_node
.. automethod:: landlab.grid.raster.RasterModelGrid.get_nodes_around_point
.. automethod:: landlab.grid.raster.RasterModelGrid.grid_coords_to_node_id
.. automethod:: landlab.grid.raster.RasterModelGrid.is_point_on_grid


Derive offsets, gradients, flux divergences, and steepest descents
==================================================================

*These methods allow calculation of derivatives, divergences, and steepest paths through
data defined on the grid.*

.. automethod:: landlab.grid.base.ModelGrid.calculate_diff_at_links
.. automethod:: landlab.grid.base.ModelGrid.calculate_diff_at_active_links
.. automethod:: landlab.grid.base.ModelGrid.calculate_flux_divergence_at_core_nodes
.. automethod:: landlab.grid.base.ModelGrid.calculate_flux_divergence_at_nodes
.. automethod:: landlab.grid.base.ModelGrid.calculate_gradients_at_links
.. automethod:: landlab.grid.base.ModelGrid.calculate_gradients_at_active_links

.. automethod:: landlab.grid.raster.RasterModelGrid.calculate_flux_divergence_at_nodes
.. automethod:: landlab.grid.raster.RasterModelGrid.calculate_gradient_across_cell_faces
.. automethod:: landlab.grid.raster.RasterModelGrid.calculate_gradient_across_cell_corners
.. automethod:: landlab.grid.raster.RasterModelGrid.calculate_steepest_descent_across_adjacent_cells
.. automethod:: landlab.grid.raster.RasterModelGrid.calculate_steepest_descent_across_cell_corners
.. automethod:: landlab.grid.raster.RasterModelGrid.calculate_steepest_descent_across_cell_faces


Control boundary conditions
===========================

*These methods allow explicit control of the boundary nodes in the grid, and their
properties.*
*Note that boundary condition handling may change somewhat in future development, in
particular improving functionality for Voronoi grids and rasters with non-perimeter
boundary nodes.*

.. autoattribute:: landlab.grid.base.ModelGrid.closed_boundary_nodes
.. autoattribute:: landlab.grid.base.ModelGrid.core_nodes
.. automethod:: landlab.grid.base.ModelGrid.is_boundary
.. autoattribute:: landlab.grid.base.ModelGrid.node_boundary_status
.. autoattribute:: landlab.grid.base.ModelGrid.open_boundary_nodes
.. automethod:: landlab.grid.base.ModelGrid.set_closed_nodes
.. automethod:: landlab.grid.base.ModelGrid.set_nodata_nodes_to_closed
.. automethod:: landlab.grid.base.ModelGrid.update_links_nodes_cells_to_new_BCs

.. automethod:: landlab.grid.raster.RasterModelGrid.are_all_core
.. automethod:: landlab.grid.raster.RasterModelGrid.has_boundary_neighbor
.. automethod:: landlab.grid.raster.RasterModelGrid.is_core
.. automethod:: landlab.grid.raster.RasterModelGrid.set_fixed_value_boundaries_at_grid_edges
.. automethod:: landlab.grid.raster.RasterModelGrid.set_closed_boundaries_at_grid_edges
.. automethod:: landlab.grid.raster.RasterModelGrid.set_looped_boundaries


Manipulate arrays for plotting and display
==========================================

*These methods are intended to make plotting up visualizations of the grid easier.
NB: these slope and aspect methods were devised for display purposes, and not
intended or tested for quantitative use. But if you wish to explore their uses for such,
have at it!*

.. automethod:: landlab.grid.base.ModelGrid.display_grid

.. automethod:: landlab.grid.raster.RasterModelGrid.calculate_aspect_at_nodes_bestFitPlane
.. automethod:: landlab.grid.raster.RasterModelGrid.calculate_slope_at_nodes_bestFitPlane
.. automethod:: landlab.grid.raster.RasterModelGrid.cell_vector_to_raster
.. automethod:: landlab.grid.raster.RasterModelGrid.hillshade
.. automethod:: landlab.grid.raster.RasterModelGrid.node_vector_to_raster
