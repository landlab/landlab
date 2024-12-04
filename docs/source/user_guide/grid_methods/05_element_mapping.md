(api-grid-grid-summary-mappers)=

# Mapping between elements

These methods allow mapping of values defined on one grid element onto a
second, e.g., mapping upwind node values onto links, or mean link values onto
nodes.

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

      {% for func in grids[grid]['map'] %}
        ~{{func | replace("landlab.", "")}}
      {% endfor %}
  {% endfor %}
```
