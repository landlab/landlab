# Uncategorized or Deprecated

The following functions are either *uncategorized* or *deprecated*. Uncategorized functions are simply
those to which we have not yet assigned a category (but we will as we continue to improve the
documentation).
Although functions marked as deprecated are currently still available, **they will be removed** in
a future *Landlab* release and so their use is discouraged.

```{eval-rst}
.. jinja:: llcats

  .. currentmodule:: landlab

  {% for grid in ['RasterModelGrid'] + grids
    | reject('equalto', 'RasterModelGrid')
    | sort
  %}
  {% set label = grid.replace('ModelGrid', '').replace('Grid', '') %}

  .. tab:: {{ label }}

    {% for cat, label in [('uncategorized', 'Uncategorized'), ('deprecated', 'Deprecated')] %}

      .. tab:: {{label}}

        .. autosummary::
          :nosignatures:

          {% for func in grids[grid][cat] %}
            ~{{func | replace("landlab.", "")}}
          {% endfor %}
    {% endfor %}
  {% endfor %}
```
