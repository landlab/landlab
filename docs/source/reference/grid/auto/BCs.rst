
==========================
Boundary condition control
==========================


.. _BC_ModelGrid:

Base class
----------

.. currentmodule:: landlab 

.. autosummary::

    ~landlab.grid.base.ModelGrid.active_adjacent_nodes_at_node
    ~landlab.grid.base.ModelGrid.active_faces
    ~landlab.grid.base.ModelGrid.active_links
    ~landlab.grid.base.ModelGrid.boundary_nodes
    ~landlab.grid.base.ModelGrid.closed_boundary_nodes
    ~landlab.grid.base.ModelGrid.core_cells
    ~landlab.grid.base.ModelGrid.core_nodes
    ~landlab.grid.base.ModelGrid.fixed_gradient_boundary_nodes
    ~landlab.grid.base.ModelGrid.fixed_links
    ~landlab.grid.base.ModelGrid.fixed_value_boundary_nodes
    ~landlab.grid.base.ModelGrid.node_at_core_cell
    ~landlab.grid.base.ModelGrid.node_has_boundary_neighbor
    ~landlab.grid.base.ModelGrid.node_is_boundary
    ~landlab.grid.base.ModelGrid.number_of_active_faces
    ~landlab.grid.base.ModelGrid.number_of_active_links
    ~landlab.grid.base.ModelGrid.number_of_core_cells
    ~landlab.grid.base.ModelGrid.number_of_core_nodes
    ~landlab.grid.base.ModelGrid.number_of_fixed_links
    ~landlab.grid.base.ModelGrid.number_of_patches_present_at_link
    ~landlab.grid.base.ModelGrid.number_of_patches_present_at_node
    ~landlab.grid.base.ModelGrid.open_boundary_nodes
    ~landlab.grid.base.ModelGrid.set_nodata_nodes_to_closed
    ~landlab.grid.base.ModelGrid.set_nodata_nodes_to_fixed_gradient
    ~landlab.grid.base.ModelGrid.status_at_link
    ~landlab.grid.base.ModelGrid.status_at_node



.. _BC_RasterModelGrid:

Raster
------

.. currentmodule:: landlab 

.. autosummary::

    ~landlab.grid.raster.RasterModelGrid.active_adjacent_nodes_at_node
    ~landlab.grid.raster.RasterModelGrid.active_faces
    ~landlab.grid.raster.RasterModelGrid.active_links
    ~landlab.grid.raster.RasterModelGrid.boundary_nodes
    ~landlab.grid.raster.RasterModelGrid.closed_boundary_nodes
    ~landlab.grid.raster.RasterModelGrid.core_cells
    ~landlab.grid.raster.RasterModelGrid.core_nodes
    ~landlab.grid.raster.RasterModelGrid.fixed_gradient_boundary_nodes
    ~landlab.grid.raster.RasterModelGrid.fixed_links
    ~landlab.grid.raster.RasterModelGrid.fixed_value_boundary_nodes
    ~landlab.grid.raster.RasterModelGrid.node_at_core_cell
    ~landlab.grid.raster.RasterModelGrid.node_has_boundary_neighbor
    ~landlab.grid.raster.RasterModelGrid.node_is_boundary
    ~landlab.grid.raster.RasterModelGrid.number_of_active_faces
    ~landlab.grid.raster.RasterModelGrid.number_of_active_links
    ~landlab.grid.raster.RasterModelGrid.number_of_core_cells
    ~landlab.grid.raster.RasterModelGrid.number_of_core_nodes
    ~landlab.grid.raster.RasterModelGrid.number_of_fixed_links
    ~landlab.grid.raster.RasterModelGrid.number_of_patches_present_at_link
    ~landlab.grid.raster.RasterModelGrid.number_of_patches_present_at_node
    ~landlab.grid.raster.RasterModelGrid.open_boundary_nodes
    ~landlab.grid.raster.RasterModelGrid.second_ring_looped_neighbors_at_cell
    ~landlab.grid.raster.RasterModelGrid.set_closed_boundaries_at_grid_edges
    ~landlab.grid.raster.RasterModelGrid.set_fixed_value_boundaries_at_grid_edges
    ~landlab.grid.raster.RasterModelGrid.set_looped_boundaries
    ~landlab.grid.raster.RasterModelGrid.set_nodata_nodes_to_closed
    ~landlab.grid.raster.RasterModelGrid.set_nodata_nodes_to_fixed_gradient
    ~landlab.grid.raster.RasterModelGrid.set_open_nodes_disconnected_from_watershed_to_closed
    ~landlab.grid.raster.RasterModelGrid.set_status_at_node_on_edges
    ~landlab.grid.raster.RasterModelGrid.set_watershed_boundary_condition
    ~landlab.grid.raster.RasterModelGrid.set_watershed_boundary_condition_outlet_coords
    ~landlab.grid.raster.RasterModelGrid.set_watershed_boundary_condition_outlet_id
    ~landlab.grid.raster.RasterModelGrid.status_at_link
    ~landlab.grid.raster.RasterModelGrid.status_at_node



.. _BC_VoronoiDelaunayGrid:

Irregular Voronoi-cell
----------------------

.. currentmodule:: landlab 

.. autosummary::

    ~landlab.grid.voronoi.VoronoiDelaunayGrid.active_adjacent_nodes_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.active_faces
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.active_links
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.boundary_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.closed_boundary_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.core_cells
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.core_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.fixed_gradient_boundary_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.fixed_links
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.fixed_value_boundary_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.node_at_core_cell
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.node_has_boundary_neighbor
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.node_is_boundary
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_active_faces
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_active_links
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_core_cells
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_core_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_fixed_links
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_patches_present_at_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.number_of_patches_present_at_node
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.open_boundary_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.perimeter_nodes
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.set_nodata_nodes_to_closed
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.set_nodata_nodes_to_fixed_gradient
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.status_at_link
    ~landlab.grid.voronoi.VoronoiDelaunayGrid.status_at_node



.. _BC_HexModelGrid:

Hexagonal
---------

.. currentmodule:: landlab 

.. autosummary::

    ~landlab.grid.hex.HexModelGrid.active_adjacent_nodes_at_node
    ~landlab.grid.hex.HexModelGrid.active_faces
    ~landlab.grid.hex.HexModelGrid.active_links
    ~landlab.grid.hex.HexModelGrid.boundary_nodes
    ~landlab.grid.hex.HexModelGrid.closed_boundary_nodes
    ~landlab.grid.hex.HexModelGrid.core_cells
    ~landlab.grid.hex.HexModelGrid.core_nodes
    ~landlab.grid.hex.HexModelGrid.fixed_gradient_boundary_nodes
    ~landlab.grid.hex.HexModelGrid.fixed_links
    ~landlab.grid.hex.HexModelGrid.fixed_value_boundary_nodes
    ~landlab.grid.hex.HexModelGrid.node_at_core_cell
    ~landlab.grid.hex.HexModelGrid.node_has_boundary_neighbor
    ~landlab.grid.hex.HexModelGrid.node_is_boundary
    ~landlab.grid.hex.HexModelGrid.number_of_active_faces
    ~landlab.grid.hex.HexModelGrid.number_of_active_links
    ~landlab.grid.hex.HexModelGrid.number_of_core_cells
    ~landlab.grid.hex.HexModelGrid.number_of_core_nodes
    ~landlab.grid.hex.HexModelGrid.number_of_fixed_links
    ~landlab.grid.hex.HexModelGrid.number_of_patches_present_at_link
    ~landlab.grid.hex.HexModelGrid.number_of_patches_present_at_node
    ~landlab.grid.hex.HexModelGrid.open_boundary_nodes
    ~landlab.grid.hex.HexModelGrid.set_nodata_nodes_to_closed
    ~landlab.grid.hex.HexModelGrid.set_nodata_nodes_to_fixed_gradient
    ~landlab.grid.hex.HexModelGrid.set_watershed_boundary_condition
    ~landlab.grid.hex.HexModelGrid.set_watershed_boundary_condition_outlet_id
    ~landlab.grid.hex.HexModelGrid.status_at_link
    ~landlab.grid.hex.HexModelGrid.status_at_node



.. _BC_RadialModelGrid:

Radial
------

.. currentmodule:: landlab 

.. autosummary::

    ~landlab.grid.radial.RadialModelGrid.active_adjacent_nodes_at_node
    ~landlab.grid.radial.RadialModelGrid.active_faces
    ~landlab.grid.radial.RadialModelGrid.active_links
    ~landlab.grid.radial.RadialModelGrid.boundary_nodes
    ~landlab.grid.radial.RadialModelGrid.closed_boundary_nodes
    ~landlab.grid.radial.RadialModelGrid.core_cells
    ~landlab.grid.radial.RadialModelGrid.core_nodes
    ~landlab.grid.radial.RadialModelGrid.fixed_gradient_boundary_nodes
    ~landlab.grid.radial.RadialModelGrid.fixed_links
    ~landlab.grid.radial.RadialModelGrid.fixed_value_boundary_nodes
    ~landlab.grid.radial.RadialModelGrid.node_at_core_cell
    ~landlab.grid.radial.RadialModelGrid.node_has_boundary_neighbor
    ~landlab.grid.radial.RadialModelGrid.node_is_boundary
    ~landlab.grid.radial.RadialModelGrid.number_of_active_faces
    ~landlab.grid.radial.RadialModelGrid.number_of_active_links
    ~landlab.grid.radial.RadialModelGrid.number_of_core_cells
    ~landlab.grid.radial.RadialModelGrid.number_of_core_nodes
    ~landlab.grid.radial.RadialModelGrid.number_of_fixed_links
    ~landlab.grid.radial.RadialModelGrid.number_of_patches_present_at_link
    ~landlab.grid.radial.RadialModelGrid.number_of_patches_present_at_node
    ~landlab.grid.radial.RadialModelGrid.open_boundary_nodes
    ~landlab.grid.radial.RadialModelGrid.perimeter_nodes
    ~landlab.grid.radial.RadialModelGrid.set_nodata_nodes_to_closed
    ~landlab.grid.radial.RadialModelGrid.set_nodata_nodes_to_fixed_gradient
    ~landlab.grid.radial.RadialModelGrid.status_at_link
    ~landlab.grid.radial.RadialModelGrid.status_at_node


