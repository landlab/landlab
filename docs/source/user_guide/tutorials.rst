.. _tutorials:

Tutorials
=========

The `Landlab tutorials <https://github.com/landlab/tutorials>`_
repository contains `IPython
notebooks <https://ipython.org/notebook.html>`_ that are both
unexpanded (so you can run them in an IPython notebook viewer) and
expanded (so you can read them like a regular text tutorial) along with
code examples. Landlab's tutorials repo can also be accessed from
`nbviewer <https://mybinder.org/v2/gh/landlab/tutorials/release>`_.

IPython notebook tutorials
~~~~~~~~~~~~~~~~~~~~~~~~~~

Instructions on how to run an IPython notebook can be found
`here <https://github.com/landlab/tutorials/blob/release/README.md>`__.

A short IPython notebook tutorial along with a screencast can be found
`here <http://www.randalolson.com/2012/05/12/a-short-demo-on-how-to-use-ipython-notebook-as-a-research-notebook/>`__
(the tutorial uses an example with statistics, but you can substitute
Landlab!)

`Click here to download all the
tutorials <https://codeload.github.com/landlab/tutorials/zip/release>`_

A suggested introduction to Landlab follows roughly this order:

-  `Introduction to Python and
   NumPy <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=python_intro/Python_intro.ipynb>`_.
   *Learn about:* The very basics of Python.
-  `Introduction to Landlab: example model of fault-scarp
   degradation <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=fault_scarp/landlab-fault-scarp.ipynb>`_.
   A short overview of some of the things Landlab can do.
-  `Introduction to the model grid
   object <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=grid_object_demo/grid_object_demo.ipynb>`_.
   Grid topology; how landlab represents data; connectivity of grid
   elements.
-  `Introduction to Landlab data
   fields <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=fields/working_with_fields.ipynb>`_.
   How Landlab stores spatial data on the grid; a little on naming
   conventions.
-  `Introduction to plotting output with
   Landlab <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=plotting/landlab-plotting.ipynb>`_.
   The basics of plotting with Landlab; combining matplotlib and out
   plots; the all-powerful ``imshow_grid()`` function.
-  `Introduction to using the Landlab component
   library <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=component_tutorial/component_tutorial.ipynb>`_.
   The basics of working with and coupling components, using
   *diffusion*, *stream power*, and a *storm generator* as examples.
-  `Using the gradient and flux-divergence
   functions <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=gradient_and_divergence/gradient_and_divergence.ipynb>`_.
   Landlab as solving environment for staggered grid finite difference
   differential approximations; functions available to help you do this.
-  `Mapping values from nodes to
   links <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=mappers/mappers.ipynb>`_.
   Options for getting data on links to nodes, nodes to links, etc.;
   min, max, and mean; upwinding and downwinding schemes; one-to-one,
   one-to-many, and many-to-one mappings.
-  `Setting boundary conditions on Landlab grids (several
   tutorials) <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=boundary_conds/>`_
   How Landlab conceptualises boundary conditions; various ways to
   interact and work with them.
-  `Reading DEMs into
   Landlab <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=reading_dem_into_landlab/reading_dem_into_landlab.ipynb>`_
   Getting an ARC ESRI ASCII into Landlab; getting the boundary
   conditions set right.
-  `How to write a Landlab
   component <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=making_components/making_components.ipynb>`_
   What makes up a Landlab Component Standard Interface; how to make one
   for your process model.

Notebook tutorials on Landlab's components include: \* Flow Direction
and Accumulation \* `Introduction to the FlowDirector
Components <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=flow_direction_and_accumulation/the_FlowDirectors.ipynb>`_
\* `Introduction to the FlowAccumulator
Component <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=flow_direction_and_accumulation/the_FlowAccumulator.ipynb>`_
\* `Comparison of FlowDirector
Components <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=flow_direction_and_accumulation/compare_FlowDirectors.ipynb>`_
\*
`Flexure <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=flexure/lots_of_loads.ipynb>`_
\* `Overland
flow <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=overland_flow/overland_flow_driver.ipynb>`_
\* `Diffusion, stream power, and the storm
generator <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=component_tutorial/component_tutorial.ipynb>`_
\* `Ecohydrology Model on Flat
Domain <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=ecohydrology/cellular_automaton_vegetation_flat_surface/cellular_automaton_vegetation_flat_domain.ipynb>`_
\* `Ecohydrology Model on Actual
Landscape <https://mybinder.org/v2/gh/landlab/tutorials/release?filepath=ecohydrology/cellular_automaton_vegetation_DEM/cellular_automaton_vegetation_DEM.ipynb>`_

Tutorial-like examples
~~~~~~~~~~~~~~~~~~~~~~

-  :ref:`How to import topography data from a grid-based digital elevation
   model
   (DEM) <importing_a_dem>`

-  :ref:`OverlandFlow component User
   Guide <overland_flow_manual>`

Landlab clinics and workshops
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For more examples and tutorials, see also our :ref:`Clinics & workshops
page <clinics_workshops>`.

If you write a Landlab tutorial or Gist, please add a link to it here:
