(api-grid-grid-summary-bc)=

# Boundary conditions

These are the primary properties for getting and setting the grid boundary
conditions. Changes made to {meth}`~.ModelGrid.status_at_node` automatically
update the conditions defined at other grid elements.

```{eval-rst}
.. jinja:: llcats

  .. currentmodule:: landlab

  {% for grid in ['RasterModelGrid'] + grids
    | reject('equalto', 'RasterModelGrid')
    | sort
  %}
  {% set label = grid.replace('ModelGrid', '').replace('Grid', '') %}

  .. tab:: {{ label }}

    .. autosummary::
      :nosignatures:

      {% for func in grids[grid]['boundary-condition'] %}
        ~{{func}}
      {% endfor %}
  {% endfor %}
```
