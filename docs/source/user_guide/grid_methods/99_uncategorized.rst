Uncategorized and Deprecated
============================

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
  
