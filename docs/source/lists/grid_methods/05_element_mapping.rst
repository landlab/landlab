.. _api.grid.grid_summary.mappers:

Mapping between elements
========================

These methods allow mapping of values defined on one grid element onto a
second, e.g., mapping upwind node values onto links, or mean link values onto
nodes.

.. jinja:: llcats

  .. currentmodule:: landlab

  {% for grid, label in [('RasterModelGrid', 'Raster'), ('HexModelGrid', 'Hex'), ('RadialModelGrid', 'Radial'), ('VoronoiDelaunayGrid', 'Voronoi'), ('FramedVoronoiGrid', 'FramedVoronoi')] %}

  .. tab:: {{ label }}

    .. autosummary::
      :nosignatures:

      {% for func in grids[grid]['map'] %}
        ~{{func}}
      {% endfor %}
  {% endfor %}
