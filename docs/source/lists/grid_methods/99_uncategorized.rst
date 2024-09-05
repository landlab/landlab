Uncategorized or Deprecated
===========================

The following functions are either *uncategorized* or *deprecated*. Uncategorized functions are simply
those to which we have not yet assigned a category (but we will as we continue to improve the
documentation).
Although functions marked as deprecated are currently still available, **they will be removed** in
a future *Landlab* release and so their use is discouraged.

.. jinja:: llcats

  .. currentmodule:: landlab

  {% for grid, label in [('RasterModelGrid', 'Raster'), ('HexModelGrid', 'Hex'), ('RadialModelGrid', 'Radial'), ('VoronoiDelaunayGrid', 'Voronoi'), ('FramedVoronoiGrid', 'FramedVoronoi')] %}

  .. tab:: {{ label }}

    {% for cat, label in [('uncategorized', 'Uncategorized'), ('deprecated', 'Deprecated')] %}

      .. tab:: {{label}}

        .. autosummary::
          :nosignatures:

          {% for func in grids[grid][cat] %}
            ~{{func}}
          {% endfor %}
    {% endfor %}
  {% endfor %}
