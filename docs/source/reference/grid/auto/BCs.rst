
==========================
Boundary condition control
==========================


.. _BC_ModelGrid:

Base class
----------

.. autosummary::
    :toctree: generated/

    :py:meth:`~landlab.grid.base.ModelGrid.active_adjacent_nodes_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.active_faces`
    :py:meth:`~landlab.grid.base.ModelGrid.active_links`
    :py:meth:`~landlab.grid.base.ModelGrid.active_neighbors_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.boundary_nodes`
    :py:meth:`~landlab.grid.base.ModelGrid.closed_boundary_nodes`
    :py:meth:`~landlab.grid.base.ModelGrid.core_cells`
    :py:meth:`~landlab.grid.base.ModelGrid.core_nodes`
    :py:meth:`~landlab.grid.base.ModelGrid.fixed_gradient_boundary_nodes`
    :py:meth:`~landlab.grid.base.ModelGrid.fixed_links`
    :py:meth:`~landlab.grid.base.ModelGrid.fixed_value_boundary_nodes`
    :py:meth:`~landlab.grid.base.ModelGrid.node_at_core_cell`
    :py:meth:`~landlab.grid.base.ModelGrid.node_has_boundary_neighbor`
    :py:meth:`~landlab.grid.base.ModelGrid.node_is_boundary`
    :py:meth:`~landlab.grid.base.ModelGrid.number_of_active_faces`
    :py:meth:`~landlab.grid.base.ModelGrid.number_of_active_links`
    :py:meth:`~landlab.grid.base.ModelGrid.number_of_core_cells`
    :py:meth:`~landlab.grid.base.ModelGrid.number_of_core_nodes`
    :py:meth:`~landlab.grid.base.ModelGrid.number_of_fixed_links`
    :py:meth:`~landlab.grid.base.ModelGrid.number_of_patches_present_at_link`
    :py:meth:`~landlab.grid.base.ModelGrid.number_of_patches_present_at_node`
    :py:meth:`~landlab.grid.base.ModelGrid.open_boundary_nodes`
    :py:meth:`~landlab.grid.base.ModelGrid.set_nodata_nodes_to_closed`
    :py:meth:`~landlab.grid.base.ModelGrid.set_nodata_nodes_to_fixed_gradient`
    :py:meth:`~landlab.grid.base.ModelGrid.status_at_link`
    :py:meth:`~landlab.grid.base.ModelGrid.status_at_node`



.. _BC_RasterModelGrid:

Raster
------

.. autosummary::
    :toctree: generated/

    :py:meth:`~landlab.grid.raster.RasterModelGrid.active_adjacent_nodes_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.active_faces`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.active_links`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.active_neighbors_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.boundary_nodes`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.closed_boundary_nodes`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.core_cells`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.core_nodes`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.fixed_gradient_boundary_nodes`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.fixed_links`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.fixed_value_boundary_nodes`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.node_at_core_cell`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.node_has_boundary_neighbor`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.node_is_boundary`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.node_is_core`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.nodes_are_all_core`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.nodes_at_bottom_edge`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.nodes_at_edge`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.nodes_at_left_edge`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.nodes_at_right_edge`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.nodes_at_top_edge`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.number_of_active_faces`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.number_of_active_links`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.number_of_core_cells`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.number_of_core_nodes`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.number_of_fixed_links`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.number_of_patches_present_at_link`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.number_of_patches_present_at_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.open_boundary_nodes`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.second_ring_looped_neighbors_at_cell`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.set_closed_boundaries_at_grid_edges`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.set_fixed_link_boundaries_at_grid_edges`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.set_fixed_value_boundaries_at_grid_edges`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.set_looped_boundaries`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.set_nodata_nodes_to_closed`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.set_nodata_nodes_to_fixed_gradient`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.set_open_nodes_disconnected_from_watershed_to_closed`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.set_status_at_node_on_edges`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.set_watershed_boundary_condition`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.set_watershed_boundary_condition_outlet_coords`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.set_watershed_boundary_condition_outlet_id`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.status_at_link`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.status_at_node`



.. _BC_VoronoiDelaunayGrid:

Irregular Voronoi-cell
----------------------

.. autosummary::
    :toctree: generated/

    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.active_adjacent_nodes_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.active_faces`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.active_links`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.active_neighbors_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.boundary_nodes`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.closed_boundary_nodes`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.core_cells`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.core_nodes`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.fixed_gradient_boundary_nodes`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.fixed_links`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.fixed_value_boundary_nodes`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.node_at_core_cell`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.node_has_boundary_neighbor`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.node_is_boundary`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_active_faces`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_active_links`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_core_cells`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_core_nodes`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_fixed_links`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_patches_present_at_link`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_patches_present_at_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.open_boundary_nodes`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.set_nodata_nodes_to_closed`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.set_nodata_nodes_to_fixed_gradient`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.status_at_link`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.status_at_node`



.. _BC_HexModelGrid:

Hexagonal
---------

.. autosummary::
    :toctree: generated/

    :py:meth:`~landlab.grid.hex.HexModelGrid.active_adjacent_nodes_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.active_faces`
    :py:meth:`~landlab.grid.hex.HexModelGrid.active_links`
    :py:meth:`~landlab.grid.hex.HexModelGrid.active_neighbors_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.boundary_nodes`
    :py:meth:`~landlab.grid.hex.HexModelGrid.closed_boundary_nodes`
    :py:meth:`~landlab.grid.hex.HexModelGrid.core_cells`
    :py:meth:`~landlab.grid.hex.HexModelGrid.core_nodes`
    :py:meth:`~landlab.grid.hex.HexModelGrid.fixed_gradient_boundary_nodes`
    :py:meth:`~landlab.grid.hex.HexModelGrid.fixed_links`
    :py:meth:`~landlab.grid.hex.HexModelGrid.fixed_value_boundary_nodes`
    :py:meth:`~landlab.grid.hex.HexModelGrid.node_at_core_cell`
    :py:meth:`~landlab.grid.hex.HexModelGrid.node_has_boundary_neighbor`
    :py:meth:`~landlab.grid.hex.HexModelGrid.node_is_boundary`
    :py:meth:`~landlab.grid.hex.HexModelGrid.nodes_at_bottom_edge`
    :py:meth:`~landlab.grid.hex.HexModelGrid.nodes_at_left_edge`
    :py:meth:`~landlab.grid.hex.HexModelGrid.nodes_at_right_edge`
    :py:meth:`~landlab.grid.hex.HexModelGrid.nodes_at_top_edge`
    :py:meth:`~landlab.grid.hex.HexModelGrid.number_of_active_faces`
    :py:meth:`~landlab.grid.hex.HexModelGrid.number_of_active_links`
    :py:meth:`~landlab.grid.hex.HexModelGrid.number_of_core_cells`
    :py:meth:`~landlab.grid.hex.HexModelGrid.number_of_core_nodes`
    :py:meth:`~landlab.grid.hex.HexModelGrid.number_of_fixed_links`
    :py:meth:`~landlab.grid.hex.HexModelGrid.number_of_patches_present_at_link`
    :py:meth:`~landlab.grid.hex.HexModelGrid.number_of_patches_present_at_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.open_boundary_nodes`
    :py:meth:`~landlab.grid.hex.HexModelGrid.set_nodata_nodes_to_closed`
    :py:meth:`~landlab.grid.hex.HexModelGrid.set_nodata_nodes_to_fixed_gradient`
    :py:meth:`~landlab.grid.hex.HexModelGrid.set_watershed_boundary_condition`
    :py:meth:`~landlab.grid.hex.HexModelGrid.set_watershed_boundary_condition_outlet_id`
    :py:meth:`~landlab.grid.hex.HexModelGrid.status_at_link`
    :py:meth:`~landlab.grid.hex.HexModelGrid.status_at_node`



.. _BC_RadialModelGrid:

Radial
------

.. autosummary::
    :toctree: generated/

    :py:meth:`~landlab.grid.radial.RadialModelGrid.active_adjacent_nodes_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.active_faces`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.active_links`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.active_neighbors_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.boundary_nodes`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.closed_boundary_nodes`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.core_cells`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.core_nodes`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.fixed_gradient_boundary_nodes`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.fixed_links`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.fixed_value_boundary_nodes`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.node_at_core_cell`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.node_has_boundary_neighbor`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.node_is_boundary`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.number_of_active_faces`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.number_of_active_links`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.number_of_core_cells`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.number_of_core_nodes`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.number_of_fixed_links`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.number_of_patches_present_at_link`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.number_of_patches_present_at_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.open_boundary_nodes`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.set_nodata_nodes_to_closed`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.set_nodata_nodes_to_fixed_gradient`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.status_at_link`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.status_at_node`


