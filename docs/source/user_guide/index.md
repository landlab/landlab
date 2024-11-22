(user-guide)=

# User Guide

The User Guide describes Landlab by topic area.

Users brand-new to Landlab should start with {ref}`10min`.

Further information on any specific method can be obtained in the
{ref}`api`.

The Landlab project creates an environment in which scientists can build a
numerical surface process model without having to code all of the individual
components. Surface process models compute flows of mass, such as water, sediment,
glacial ice, volcanic material, or landslide debris, across a gridded terrain
surface. Surface process models have a number of commonalities, such as operating
on a grid of points and routing material across the grid. Scientists who want
to use a surface process model often build their own unique model from the ground
up, re-coding the basic building blocks of their surface process model rather than
taking advantage of codes that have already been written.

A list of papers and presentations using Landlab can be found {ref}`here <papers>`.

```{toctree}
:caption: Grids and Components
:hidden: true
:maxdepth: 2

Grid & Component reference <reference/index>
component_list
field_definitions
grid_summary
```

```{toctree}
:caption: The Landlab Grid
:hidden: true
:maxdepth: 2

grid
```

```{toctree}
:caption: Model with Landlab and Components
:hidden: true
:maxdepth: 2

components
build_a_model
```

```{toctree}
:caption: Landlab and Units
:hidden: true
:maxdepth: 2

units
```

```{toctree}
:caption: Additional resources
:hidden: true
:maxdepth: 2

time_steps
faq
```

```{toctree}
:caption: Overland flow User Guide
:hidden: true
:maxdepth: 2

overland_flow_user_guide
```

```{toctree}
:caption: CellLab-CTS User Guide
:hidden: true
:maxdepth: 2

cell_lab_user_guide
```
