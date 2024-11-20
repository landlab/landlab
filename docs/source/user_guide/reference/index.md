(api)=

# API reference

This page gives an overview of all public Landlab objects, functions and
methods.

## Grids

```{toctree}
:maxdepth: 2

grid
```

## Layers

```{toctree}
:maxdepth: 2

layers
```

## Lithology

Two objects based on the EventLayers object exist to make it easier to deal
with spatially variable lithology and associated properties. The Lithology
components contain information about spatially variable lithology and connect
with the Landlab model grid so that when rock is eroded or advected upward by
rock uplift the values of rock propeties at the topographic surface are updated.

First is the Lithology component which is a generic object for variable
lithology.

```{toctree}
:maxdepth: 1

/generated/api/landlab.components.lithology.lithology
```

Second is LithoLayers which makes it easy to make layered rock.

```{toctree}
:maxdepth: 1

/generated/api/landlab.components.lithology.litholayers
```

## Values

```{toctree}
:maxdepth: 2

values
```

## Components

```{toctree}
:maxdepth: 2

components
```

## References

- {ref}`modindex`
- {ref}`search`
