# Surface analysis

These methods permit the kinds of surface analysis that you might expect to
find in GIS software.

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

      {% for func in grids[grid]['surface'] %}
        ~{{func}}
      {% endfor %}
  {% endfor %}
```
