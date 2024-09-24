(tutorials)=

# The Landlab tutorials

Two papers have been written describing Landlab, both are open access:
* [Creative computing with Landlab: an open-source toolkit for building,
  coupling, and exploring two-dimensional numerical models of Earth-surface
  dynamics][esurf-2017]
* [Short communication: Landlab v2.0: a software package for Earth
  surface dynamics][esurf-2020]

[esurf-2017]: https://esurf.copernicus.org/articles/5/21/2017/
[esurf-2020]: https://esurf.copernicus.org/articles/8/379/2020/
[reference-guide]: /user_guide/reference/index
[tutorial-bc-perimeter]: boundary_conditions/set_BCs_on_raster_perimeter.ipynb
[tutorial-bc-voronoi]: boundary_conditions/set_BCs_on_voronoi.ipynb
[tutorial-bc-watersheds]: boundary_conditions/set_watershed_BCs_raster.ipynb
[tutorial-bc-xy]: boundary_conditions/set_BCs_from_xy.ipynb
[tutorial-components]: component_tutorial/component_tutorial.ipynb
[tutorial-div-grad]: gradient_and_divergence/gradient_and_divergence.ipynb
[tutorial-fault-scarp]: fault_scarp/landlab-fault-scarp.ipynb
[tutorial-fields]: fields/working_with_fields.ipynb
[tutorial-grid-objects]: grids/grid_object_demo.ipynb
[tutorial-mappers]: mappers/mappers.ipynb
[tutorial-plotting]: plotting/landlab-plotting.ipynb
[tutorial-reading-dem]: reading_dem_into_landlab/reading_dem_into_landlab.ipynb
[user-guide]: /user_guide/index
[user-guide-components]: /user_guide/components
[user-guide-grids]: /user_guide/reference/grid


```{note}
We highly recommend reading both before starting on the steps below.
```

##  1. Format and Outline

You will alternate between reading documentation on the [User Guide][user-guide]
finding information in the [Reference Manual][reference-guide], and working
through the tutorials.

The tutorials are presented as Jupyter notebooks, which contain a mixture of text,
images, and code blocks. When you look at the tutorials, don't just read them. Start
by clearing the results by selecting ``Kernel`` → ``Restart & Clear Output``, then go
ahead and try running each code block as you come to it.

### 1.1 A motivating example

- [Notebook: Introduction to Landlab: example model of fault-scarp
  degradation][tutorial-fault-scarp]

### 1.2 Using the Documentation

The Landlab Reference Manual contains documentation for most functions in the
Landlab package. It is the comprehensive counterpart to the anecdotal tutorials.

- Spend some time clicking around in the [User Guide][user-guide] and
[Reference Manual][reference-guide] to get a sense for what is there. Tip: to
find a particular command, click on Index and use your browser's search function
to search for a command by name or keyword. For example, look at the
{class}`~.LinearDiffuser`, which you just used in the prior tutorial.

### 1.3 Introduction to the Landlab Grid and Fields

First, lets look at the [User Guide page on Landlab grids][user-guide-grids].

- [Notebook: Introduction to the model grid object][tutorial-grid-objects]
  Grid topology; how landlab represents data; connectivity of grid elements.
- [Notebook: Introduction to Landlab data fields][tutorial-fields]
How Landlab stores spatial data on the grid; a little on naming conventions.

Extra credit: Go back to the [Hobley et al. 2017 publication][esurf-2017]
and identify the ordering conventions of nodes, links, and other grid elements.

### 1.4 Working with Digital Elevtion Models (DEMs)

- [Notebook: Reading DEMs into Landlab][tutorial-reading-dem] Getting an ARC
  ESRI ASCII into Landlab; getting the boundary conditions set right.

### 1.5 Plotting

- [Notebook: Introduction to plotting output with Landlab][tutorial-plotting]
  The basics of plotting with Landlab; combining matplotlib and out plots; the
  {func}`~.imshow_grid` function.

### 1.6 Boundary conditions

- Setting boundary conditions on Landlab grids (several tutorials): How Landlab
  conceptualises boundary conditions; various ways to interact and work with them.
  - [Notebook: Raster perimeter][tutorial-bc-perimeter]
  - [Notebook: Based on X-Y values][tutorial-bc-xy]
  - [Notebook: Watersheds][tutorial-bc-watersheds]
  - [Notebook: Voronoi][tutorial-bc-voronoi]

### 1.7 Introduction to Components

- Read the [Task: Component page in the User Guide][user-guide-components]
- [Notebook: Introduction to using the Landlab component library][tutorial-components]
  The basics of working with and coupling components, using {class}`~.LinearDiffuser`,
  {class}`~.FastscapeEroder`, and {class}`~.PrecipitationDistribution`.

### 1.8 Advanced Grid and Fields: Gradients, Flux-Divergence, Mapping

In addition to having lots of important information about adjacency of nodes, links,
and other grid elements, the Landlab Grid object has a number of built-in functions
for calculating quantities like gradients and flux-divergence, and for mapping
quantities from nodes to links and so forth. Work through these tutorials to get a
sense of this functionality:

- [Notebook: Using the gradient and flux-divergence functions][tutorial-div-grad]
  Landlab as solving environment for staggered grid finite difference differential
  approximations; functions available to help you do this.
- [Notebook: Mapping values from nodes to links][tutorial-mappers] Options
  for getting data on links to nodes, nodes to links, etc.; min, max, and mean;
  upwinding and downwinding schemes; one-to-one, one-to-many, and many-to-one mappings.


```{toctree}
:caption: Gallery
:hidden:
:glob:

Gallery </generated/tutorials/index>
```
