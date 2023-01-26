Surface analysis
================

These methods permit the kinds of surface analysis that you might expect to
find in GIS software.

.. jinja:: llcats

  .. currentmodule:: landlab

  {% for grid, label in [('RasterModelGrid', 'Raster'), ('HexModelGrid', 'Hex'), ('RadialModelGrid', 'Radial'), ('VoronoiDelaunayGrid', 'Voronoi'), ('FramedVoronoiGrid', 'FramedVoronoi')] %}

  .. tab:: {{ label }}

    .. autosummary::
      :nosignatures:

      {% for func in grids[grid]['surface'] %}
        ~{{func}}
      {% endfor %}
  {% endfor %}
