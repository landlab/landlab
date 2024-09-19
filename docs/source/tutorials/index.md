(tutorials)=

# Recommended order for Landlab tutorials

Note, a paper has been written describing Landlab. It is Open Access, and a link to the PDF is [here](https://www.earth-surf-dynam.net/5/21/2017/esurf-5-21-2017.pdf)

**We highly recommend reading it before starting on the steps below.**

##  1. Format and Outline
You will alternate between reading documentation on the [User Guide](https://landlab.readthedocs.io/en/latest/user_guide/index.html), finding information in the [Reference Manual](https://landlab.readthedocs.io/en/latest/reference/index.html), and working through the tutorials.

The tutorials are presented as Jupyter notebooks, which contain a mixture of text, images, and code blocks. When you look at the tutorials, don't just read them. Start by clearing the results by selecting "Kernel ==> Restart & Clear Output," then go ahead and try running each code block as you come to it.

### 1.1 A motivating example

- [Notebook: Introduction to Landlab: example model of fault-scarp degradation](fault_scarp/landlab-fault-scarp.ipynb)

### 1.2 Using the Documentation

The Landlab Reference Manual contains documentation for most functions in the Landlab package. It is the comprehensive counterpart to the anecdotal tutorials.

- Spend some time clicking around in the [User Guide](https://landlab.readthedocs.io/en/latest/user_guide/index.html) and [Reference Manual](https://landlab.readthedocs.io/en/latest/reference/index.html) to get a sense for what is there. Tip: to find a particular command, click on Index and use your browser's search function to search for a command by name or keyword. For example, look at the [documentation for the LinearDiffuser](https://landlab.readthedocs.io/en/latest/reference/components/diffusion.html) which you just used in the prior tutorial.

### 1.3 Introduction to the Landlab Grid and Fields

First, lets look at the [User Guide page on Landlab grids](https://landlab.readthedocs.io/en/latest/user_guide/grid.html).

- [Notebook: Introduction to the model grid object](grids/grid_object_demo.ipynb) Grid topology; how landlab represents data; connectivity of grid elements.
- [Notebook: Introduction to Landlab data fields](fields/working_with_fields.ipynb) How Landlab stores spatial data on the grid; a little on naming conventions.

Extra credit: Go back to the [Hobley et al. 2017 publication](https://www.earth-surf-dynam.net/5/21/2017/esurf-5-21-2017.html) and identify the ordering conventions of nodes, links, and other grid elements.

### 1.4 Working with Digital Elevtion Models (DEMs)

- [Notebook: Reading DEMs into Landlab](reading_dem_into_landlab/reading_dem_into_landlab.ipynb) Getting an ARC ESRI ASCII into Landlab; getting the boundary conditions set right.

### 1.5 Plotting

- [Notebook: Introduction to plotting output with Landlab](plotting/landlab-plotting.ipynb) The basics of plotting with Landlab; combining matplotlib and out plots; the all-powerful [``imshow_grid()``](https://landlab.readthedocs.io/en/latest/reference/plot/index.html#landlab.plot.imshow_grid) function.

### 1.6 Boundary conditions

- Setting boundary conditions on Landlab grids (several tutorials): How Landlab conceptualises boundary conditions; various ways to interact and work with them.
  - [Notebook: Raster perimeter](boundary_conditions/set_BCs_on_raster_perimeter.ipynb)
  - [Notebook: Based on X-Y values](boundary_conditions/set_BCs_from_xy.ipynb)
  - [Notebook: Watersheds](boundary_conditions/set_watershed_BCs_raster.ipynb)
  - [Notebook: Voronoi](boundary_conditions/set_BCs_on_voronoi.ipynb)

### 1.7 Introduction to Components

- Read the [Task: Component page in the User Guide](https://landlab.readthedocs.io/en/latest/user_guide/components.html)
- [Notebook: Introduction to using the Landlab component library](component_tutorial/component_tutorial.ipynb) The basics of working with and coupling components, using *diffusion*, *stream power*, and a *storm generator* as examples.

### 1.8 Advanced Grid and Fields: Gradients, Flux-Divergence, Mapping

In addition to having lots of important information about adjacency of nodes, links, and other grid elements, the Landlab Grid object has a number of built-in functions for calculating quantities like gradients and flux-divergence, and for mapping quantities from nodes to links and so forth. Work through these tutorials to get a sense of this functionality:

- [Notebook: Using the gradient and flux-divergence functions](gradient_and_divergence/gradient_and_divergence.ipynb) Landlab as solving environment for staggered grid finite difference differential approximations; functions available to help you do this.
- [Notebook: Mapping values from nodes to links](mappers/mappers.ipynb) Options for getting data on links to nodes, nodes to links, etc.; min, max, and mean; upwinding and downwinding schemes; one-to-one, one-to-many, and many-to-one mappings.


```{toctree}
:caption: Tutorials
:hidden:
:glob:

/generated/tutorials/**/_index
```
