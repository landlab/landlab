(api-grid)=

# Landlab Grids

## Grid types

Landlab presently supports multiple grid types. The base class for all of these
grid types is {class}`~.ModelGrid`.

The following is an introduction to their properties and methods:

```{jinja} llcats
{% for name in grids|sort %}
* {{"{"}}class{{"}"}}`landlab.{{name}}`
{% endfor %}
```

## Additional Methods and Properties

* {mod}`landlab.grid.create`
* {mod}`landlab.grid.decorators`
* {mod}`landlab.grid.diagonals`
* {mod}`landlab.grid.divergence`
* {mod}`landlab.grid.gradients`
* {mod}`landlab.grid.grid_funcs`
* {mod}`landlab.grid.linkstatus`
* {mod}`landlab.grid.mappers`
* {mod}`landlab.grid.nodestatus`
* {mod}`landlab.grid.raster_aspect`
* {mod}`landlab.grid.raster_funcs`
* {mod}`landlab.grid.raster_gradients`
* {mod}`landlab.grid.raster_mappers`
* {mod}`landlab.grid.raster_set_status`
* {mod}`landlab.grid.warnings`

## API for each grid type

* {mod}`landlab.grid.base`
* {mod}`landlab.grid.raster`
* {mod}`landlab.grid.voronoi`
* {mod}`landlab.grid.framed_voronoi`
* {mod}`landlab.grid.hex`
* {mod}`landlab.grid.radial`
* {mod}`landlab.grid.network`
* {mod}`landlab.grid.icosphere`

## Additional Grid Base Classes

* {mod}`landlab.grid.unstructured`
