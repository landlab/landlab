Nodes, Links, and Patches
=========================

.. jinja:: llcats

  .. currentmodule:: landlab

  {% for grid, label in [('RasterModelGrid', 'Raster'), ('HexModelGrid', 'Hex'), ('RadialModelGrid', 'Radial'), ('VoronoiDelaunayGrid', 'Voronoi'), ('FramedVoronoiGrid', 'FramedVoronoi')] %}

  .. tab:: {{ label }}

    {% for cat, label in [('info-node', 'Nodes'), ('info-link', 'Links'), ('info-patch', 'Patches')] %}

      .. tab:: {{label}}

        .. autosummary::
          :nosignatures:

          {% for func in grids[grid][cat] %}
            ~{{func}}
          {% endfor %}
    {% endfor %}
  {% endfor %}
