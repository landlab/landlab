# Subsets of elements

These methods are useful in identifying subsets of grid elements, e.g., closest node
to a point; nodes at edges.

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

      {% for func in grids[grid]['subset'] %}
        ~{{func}}
      {% endfor %}
  {% endfor %}
```
