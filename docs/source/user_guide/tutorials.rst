.. _tutorials:

Tutorials
=========

The `Landlab tutorials <https://github.com/landlab/tutorials>`__
repository contains `IPython
notebooks <https://ipython.org/notebook.html>`__ that are both
unexpanded (so you can run them in an IPython notebook viewer) and
expanded (so you can read them like a regular text tutorial) along with
code examples. Landlab’s tutorials repo can also be accessed from
`nbviewer <https://nbviewer.jupyter.org/github/landlab/tutorials>`__.

IPython notebook tutorials
~~~~~~~~~~~~~~~~~~~~~~~~~~

Instructions on how to run an IPython notebook can be found
`here <https://github.com/landlab/tutorials/blob/release/README.md>`__.

A short IPython notebook tutorial along with a screencast can be found
`here <http://www.randalolson.com/2012/05/12/a-short-demo-on-how-to-use-ipython-notebook-as-a-research-notebook/>`__
(the tutorial uses an example with statistics, but you can substitute
Landlab!)

`Click here to download all the
tutorials <https://github.com/landlab/tutorials/archive/release.zip>`__

A suggested introduction to Landlab follows roughly this order:

-  `Introduction to Python and
   NumPy <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/python_intro/Python_intro.ipynb>`__.
   *Learn about:* The very basics of Python.
-  `Introduction to Landlab: example model of fault-scarp
   degradation <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/fault_scarp/landlab-fault-scarp.ipynb>`__.
   A short overview of some of the things Landlab can do.
-  `Introduction to the model grid
   object <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/grid_object_demo/grid_object_demo.ipynb>`__.
   Grid topology; how landlab represents data; connectivity of grid
   elements.
-  `Introduction to Landlab data
   fields <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/fields/working_with_fields.ipynb>`__.
   How Landlab stores spatial data on the grid; a little on naming
   conventions.
-  `Introduction to plotting output with
   Landlab <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/plotting/landlab-plotting.ipynb>`__.
   The basics of plotting with Landlab; combining matplotlib and out
   plots; the all-powerful ``imshow_grid()`` function.
-  `Introduction to using the Landlab component
   library <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/component_tutorial/component_tutorial.ipynb>`__.
   The basics of working with and coupling components, using
   *diffusion*, *stream power*, and a *storm generator* as examples.
-  `Using the gradient and flux-divergence
   functions <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/gradient_and_divergence/gradient_and_divergence.ipynb>`__.
   Landlab as solving environment for staggered grid finite difference
   differential approximations; functions available to help you do this.
-  `Mapping values from nodes to
   links <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/mappers/mappers.ipynb>`__.
   Options for getting data on links to nodes, nodes to links, etc.;
   min, max, and mean; upwinding and downwinding schemes; one-to-one,
   one-to-many, and many-to-one mappings.
-  `Setting boundary conditions on Landlab grids (several
   tutorials) <https://nbviewer.jupyter.org/github/landlab/tutorials/tree/release/boundary_conds/>`__
   How Landlab conceptualises boundary conditions; various ways to
   interact and work with them.
-  `Reading DEMs into
   Landlab <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/reading_dem_into_landlab/reading_dem_into_landlab.ipynb>`__
   Getting an ARC ESRI ASCII into Landlab; getting the boundary
   conditions set right.
-  `How to write a Landlab
   component <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/making_components/making_components.ipynb>`__
   What makes up a Landlab Component Standard Interface; how to make one
   for your process model.

Notebook tutorials on Landlab’s components include: \* Flow Direction
and Accumulation \* `Introduction to the FlowDirector
Components <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/flow_direction_and_accumulation/the_FlowDirectors.ipynb>`__
\* `Introduction to the FlowAccumulator
Component <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/flow_direction_and_accumulation/the_FlowAccumulator.ipynb>`__
\* `Comparison of FlowDirector
Components <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/flow_direction_and_accumulation/compare_FlowDirectors.ipynb>`__
\*
`Flexure <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/flexure/lots_of_loads.ipynb>`__
\* `Overland
flow <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/overland_flow/overland_flow_driver.ipynb>`__
\* `Diffusion, stream power, and the storm
generator <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/component_tutorial/component_tutorial.ipynb>`__
\* `Ecohydrology Model on Flat
Domain <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/ecohydrology/cellular_automaton_vegetation_flat_surface/cellular_automaton_vegetation_flat_domain.ipynb>`__
\* `Ecohydrology Model on Actual
Landscape <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/ecohydrology/cellular_automaton_vegetation_DEM/cellular_automaton_vegetation_DEM.ipynb>`__

Tutorial-like examples
~~~~~~~~~~~~~~~~~~~~~~

-  `How to import topography data from a grid-based digital elevation
   model
   (DEM) <https://github.com/landlab/landlab/wiki/Grid#importing-a-dem>`__

-  `OverlandFlow component User
   Guide <https://github.com/landlab/landlab/wiki/OverlandFlow-Component-Users-Manual>`__

Landlab clinics and workshops
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For more examples and tutorials, see also our `Clinics & workshops
page <https://github.com/landlab/landlab/wiki/Landlab-Clinics-and-Workshops>`__

If you write a Landlab tutorial or Gist, please add a link to it here:
