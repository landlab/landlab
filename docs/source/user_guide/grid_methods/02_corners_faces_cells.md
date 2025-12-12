# Corners, Faces, and Cells

```{eval-rst}
.. jinja:: llcats

  .. currentmodule:: landlab

  {% for grid in ['RasterModelGrid'] + grids
    | reject('equalto', 'RasterModelGrid')
    | sort
  %}
  {% set label = grid.replace('ModelGrid', '').replace('Grid', '') %}

  .. tab:: {{ label }}

    {% for cat, label in [('info-corner', 'Corners'), ('info-face', 'Faces'), ('info-cell', 'Cells')] %}

      .. tab:: {{label}}

        .. autosummary::
          :nosignatures:

          {% for func in grids[grid][cat] %}
            ~{{func | replace("landlab.", "")}}
          {% endfor %}
    {% endfor %}
  {% endfor %}
```
