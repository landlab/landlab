.. _api.grid.grid_summary.bc:

Boundary conditions
===================

These are the primary properties for getting and setting the grid boundary
conditions. Changes made to :meth:`~.ModelGrid.status_at_node` automatically
update the conditions defined at other grid elements.

.. jinja:: llcats
  
  .. currentmodule:: landlab
    
  {% for grid, label in [('RasterModelGrid', 'Raster'), ('HexModelGrid', 'Hex'), ('RadialModelGrid', 'Radial'), ('VoronoiDelaunayGrid', 'Voronoi')] %}
  
  .. tab:: {{ label }}
    
    {% for cat, label in [('boundary-condition', 'Boundary Conditions')] %}
    
      .. tab:: {{label}}
      
        .. autosummary::
          :nosignatures:
        
          {% for func in grids[grid][cat] %}
            ~{{func}}      
          {% endfor %}
    {% endfor %}
  {% endfor %}