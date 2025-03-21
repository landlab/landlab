# Gradients, fluxes, and divergences

Landlab is designed to easily calculate gradients in quantities across the
grid, and to construct fluxes and flux divergences from them. Because these
calculations tend to be a little more involved than property lookups, the
methods tend to start with `calc_`.

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

      {% for func in grids[grid]['gradient'] %}
        ~{{func | replace("landlab.", "")}}
      {% endfor %}
  {% endfor %}
```
