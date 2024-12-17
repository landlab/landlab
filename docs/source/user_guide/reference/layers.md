(api-layers)=

# Layers

Landlab has the ability to add layers to the grid. Two types of layers are
currently supported. First is {class}`~landlab.layers.eventlayers.EventLayers`
in which each event is preserved as an entry into the datastructure, even if no
deposition occurs. If you are interested in chronostratigraphy, this is probably what
you are interested in. Second is {class}`~landlab.layers.materiallayers.MaterialLayers`,
in which each layer must contain some material. If an entire layer is eroded in
{class}`~landlab.layers.materiallayers.MaterialLayers`, the layer is removed.
{class}`~landlab.layers.materiallayers.MaterialLayers` will likely use less memory
than {class}`~landlab.layers.eventlayers.EventLayers`.

* {class}`landlab.layers.eventlayers.EventLayers`
* {class}`landlab.layers.materiallayers.MaterialLayers`

## Lithology

Two objects based on the {class}`~landlab.layers.eventlayers.EventLayers` object exist
to make it easier to deal with spatially variable lithology and associated properties.
The {mod}`~landlab.components.lithology` components contain information about spatially
variable lithology and connect with the Landlab model grid so that when rock is eroded
or advected upward by rock uplift the values of rock propeties at the topographic
surface are updated.

First is the {class}`~landlab.components.lithology.lithology.Lithology` component,
which is a generic object for variable lithology.

* {class}`landlab.components.lithology.lithology.Lithology`

Second is {class}`~landlab.components.lithology.litholayers.LithoLayers` which makes
it easy to make layered rock.

* {class}`landlab.components.lithology.litholayers.LithoLayers`
