Subsets of elements
===================

These methods are useful in identifying subsets of grid elements, e.g., closest node
to a point; nodes at edges.

.. jinja:: llcats

  .. currentmodule:: landlab

  {% for grid, label in [('RasterModelGrid', 'Raster'), ('HexModelGrid', 'Hex'), ('RadialModelGrid', 'Radial'), ('VoronoiDelaunayGrid', 'Voronoi'), ('FramedVoronoiGrid', 'FramedVoronoi')] %}

  .. tab:: {{ label }}

    .. autosummary::
      :nosignatures:

      {% for func in grids[grid]['subset'] %}
        ~{{func}}
      {% endfor %}
  {% endfor %}
