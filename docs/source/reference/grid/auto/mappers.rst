
=======
Mappers
=======


.. _MAP_ModelGrid:

Base class
----------

.. autosummary::
    :toctree: generated/

    :py:meth:`~landlab.grid.base.ModelGrid.map_downwind_node_link_max_to_node`
    :py:meth:`~landlab.grid.base.ModelGrid.map_downwind_node_link_mean_to_node`
    :py:meth:`~landlab.grid.base.ModelGrid.map_link_head_node_to_link`
    :py:meth:`~landlab.grid.base.ModelGrid.map_link_tail_node_to_link`
    :py:meth:`~landlab.grid.base.ModelGrid.map_link_vector_sum_to_patch`
    :py:meth:`~landlab.grid.base.ModelGrid.map_link_vector_to_nodes`
    :py:meth:`~landlab.grid.base.ModelGrid.map_max_of_link_nodes_to_link`
    :py:meth:`~landlab.grid.base.ModelGrid.map_max_of_node_links_to_node`
    :py:meth:`~landlab.grid.base.ModelGrid.map_max_of_patch_nodes_to_patch`
    :py:meth:`~landlab.grid.base.ModelGrid.map_mean_of_link_nodes_to_link`
    :py:meth:`~landlab.grid.base.ModelGrid.map_mean_of_patch_nodes_to_patch`
    :py:meth:`~landlab.grid.base.ModelGrid.map_min_of_link_nodes_to_link`
    :py:meth:`~landlab.grid.base.ModelGrid.map_min_of_node_links_to_node`
    :py:meth:`~landlab.grid.base.ModelGrid.map_min_of_patch_nodes_to_patch`
    :py:meth:`~landlab.grid.base.ModelGrid.map_node_to_cell`
    :py:meth:`~landlab.grid.base.ModelGrid.map_upwind_node_link_max_to_node`
    :py:meth:`~landlab.grid.base.ModelGrid.map_upwind_node_link_mean_to_node`
    :py:meth:`~landlab.grid.base.ModelGrid.map_value_at_downwind_node_link_max_to_node`
    :py:meth:`~landlab.grid.base.ModelGrid.map_value_at_max_node_to_link`
    :py:meth:`~landlab.grid.base.ModelGrid.map_value_at_min_node_to_link`
    :py:meth:`~landlab.grid.base.ModelGrid.map_value_at_upwind_node_link_max_to_node`



.. _MAP_RasterModelGrid:

Raster
------

.. autosummary::
    :toctree: generated/

    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_downwind_node_link_max_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_downwind_node_link_mean_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_link_head_node_to_link`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_link_tail_node_to_link`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_link_vector_sum_to_patch`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_link_vector_to_nodes`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_max_of_inlinks_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_max_of_link_nodes_to_link`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_max_of_node_links_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_max_of_outlinks_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_max_of_patch_nodes_to_patch`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_mean_of_horizontal_active_links_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_mean_of_horizontal_links_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_mean_of_inlinks_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_mean_of_link_nodes_to_link`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_mean_of_links_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_mean_of_outlinks_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_mean_of_patch_nodes_to_patch`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_mean_of_vertical_active_links_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_mean_of_vertical_links_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_min_of_inlinks_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_min_of_link_nodes_to_link`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_min_of_node_links_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_min_of_outlinks_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_min_of_patch_nodes_to_patch`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_node_to_cell`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_sum_of_inlinks_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_sum_of_outlinks_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_upwind_node_link_max_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_upwind_node_link_mean_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_value_at_downwind_node_link_max_to_node`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_value_at_max_node_to_link`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_value_at_min_node_to_link`
    :py:meth:`~landlab.grid.raster.RasterModelGrid.map_value_at_upwind_node_link_max_to_node`



.. _MAP_VoronoiDelaunayGrid:

Irregular Voronoi-cell
----------------------

.. autosummary::
    :toctree: generated/

    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_downwind_node_link_max_to_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_downwind_node_link_mean_to_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_link_head_node_to_link`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_link_tail_node_to_link`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_link_vector_sum_to_patch`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_link_vector_to_nodes`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_max_of_link_nodes_to_link`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_max_of_node_links_to_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_max_of_patch_nodes_to_patch`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_mean_of_link_nodes_to_link`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_mean_of_patch_nodes_to_patch`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_min_of_link_nodes_to_link`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_min_of_node_links_to_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_min_of_patch_nodes_to_patch`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_node_to_cell`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_upwind_node_link_max_to_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_upwind_node_link_mean_to_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_value_at_downwind_node_link_max_to_node`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_value_at_max_node_to_link`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_value_at_min_node_to_link`
    :py:meth:`~landlab.grid.voronoi.VoronoiDelaunayGrid.map_value_at_upwind_node_link_max_to_node`



.. _MAP_HexModelGrid:

Hexagonal
---------

.. autosummary::
    :toctree: generated/

    :py:meth:`~landlab.grid.hex.HexModelGrid.map_downwind_node_link_max_to_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_downwind_node_link_mean_to_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_link_head_node_to_link`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_link_tail_node_to_link`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_link_vector_sum_to_patch`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_link_vector_to_nodes`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_max_of_link_nodes_to_link`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_max_of_node_links_to_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_max_of_patch_nodes_to_patch`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_mean_of_link_nodes_to_link`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_mean_of_patch_nodes_to_patch`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_min_of_link_nodes_to_link`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_min_of_node_links_to_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_min_of_patch_nodes_to_patch`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_node_to_cell`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_upwind_node_link_max_to_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_upwind_node_link_mean_to_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_value_at_downwind_node_link_max_to_node`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_value_at_max_node_to_link`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_value_at_min_node_to_link`
    :py:meth:`~landlab.grid.hex.HexModelGrid.map_value_at_upwind_node_link_max_to_node`



.. _MAP_RadialModelGrid:

Radial
------

.. autosummary::
    :toctree: generated/

    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_downwind_node_link_max_to_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_downwind_node_link_mean_to_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_link_head_node_to_link`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_link_tail_node_to_link`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_link_vector_sum_to_patch`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_link_vector_to_nodes`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_max_of_link_nodes_to_link`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_max_of_node_links_to_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_max_of_patch_nodes_to_patch`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_mean_of_link_nodes_to_link`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_mean_of_patch_nodes_to_patch`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_min_of_link_nodes_to_link`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_min_of_node_links_to_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_min_of_patch_nodes_to_patch`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_node_to_cell`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_upwind_node_link_max_to_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_upwind_node_link_mean_to_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_value_at_downwind_node_link_max_to_node`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_value_at_max_node_to_link`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_value_at_min_node_to_link`
    :py:meth:`~landlab.grid.radial.RadialModelGrid.map_value_at_upwind_node_link_max_to_node`


