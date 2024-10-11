(api-grid)=

# Landlab Grids

## Grid types

Landlab presently supports multiple grid types. The base class for all of these
grid types is {class}`~.ModelGrid`.

The following is an introduction to their properties and methods:

```{jinja} llcats
{% for name in grids|sort %}
* {{"{"}}class{{"}"}}`~.{{name}}`
{% endfor %}
```

## Additional Methods and Properties

```{toctree}
:maxdepth: 1

/generated/api/landlab.grid.base
/generated/api/landlab.grid.create
/generated/api/landlab.grid.decorators
/generated/api/landlab.grid.diagonals
/generated/api/landlab.grid.divergence
/generated/api/landlab.grid.gradients
/generated/api/landlab.grid.grid_funcs
/generated/api/landlab.grid.linkstatus
/generated/api/landlab.grid.mappers
/generated/api/landlab.grid.nodestatus
/generated/api/landlab.grid.raster_aspect
/generated/api/landlab.grid.raster_funcs
/generated/api/landlab.grid.raster_gradients
/generated/api/landlab.grid.raster_mappers
/generated/api/landlab.grid.raster_set_status
/generated/api/landlab.grid.warnings
```

## API for each grid type

```{toctree}
:maxdepth: 1

/generated/api/landlab.grid.base
/generated/api/landlab.grid.raster
/generated/api/landlab.grid.voronoi
/generated/api/landlab.grid.framed_voronoi
/generated/api/landlab.grid.hex
/generated/api/landlab.grid.radial
/generated/api/landlab.grid.network
/generated/api/landlab.grid.icosphere
```

## Additional Grid Base Classes

```{toctree}
/generated/api/landlab.grid.unstructured
```
