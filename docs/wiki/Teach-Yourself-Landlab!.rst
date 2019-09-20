Self-study course for learning to use Landlab in your research
==============================================================

WARNING! WARNING! THIS SELF STUDY COURSE IS UNDER DEVELOPMENT AND
INITIAL TESTING. THIS WILL BE REMOVED WHEN TESTING IS COMPLETE. UNTIL
THEN… LINKS MAY BE BROKEN, ORDER OF ACTIVITIES CHANGED, INSTRUCTIONS
WRONG, COMMENTS INCOMPLETE AND MISSPELLED.

Introduction and scope
----------------------

The following set of steps is intended as a self-guided journey in which
you will learn how to use Landlab for your research.

You will : - Get a comprehensive overview of Landlab’s model grid,
process components, and utilities. - Make different models. - Practice
using Landlab to build a model on your own.

You will not : - Learn the details of using Git to collaborate. - Learn
how to modify the Landlab source code. - Get a comprehensive
introduction to landscape evolution theory.

After completion, you may consider the following continuation: Learning
to develop with Landlab [Link once this exists].

This is a self-study course, so you should work through it at your own
pace. As with most self-studies, you will get out of it what you put in.

Finally, a paper has been written describing Landlab. It is Open Access,
and a link to the PDF is
`here <https://www.earth-surf-dynam.net/5/21/2017/esurf-5-21-2017.pdf>`__.
**We highly recommend reading it before starting on the steps below.**

Basic Outline
-------------

1. Installation and Updating Landlab
2. Basics of Python and Numpy
3. Introduction to Landlab
4. Capstone tutorial

1. Installing and Updating Landlab
----------------------------------

Follow [[these instructions for a standard installation \|
Installing-Landlab ]].

It’s a good idea to keep your version up to date. Once you have
installed Landlab, you can update it to the latest release version by
entering the following on the command line:

[ADD CONDA UPDATE COMMAND HERE]

2. Basics of Python, Numpy, and Matplotlib
------------------------------------------

Landlab is build with Python, so it is essential that you feel
comfortable with programming in Python. Unless you use the python
packages NumPy and SciPy regularly, we highly recommend that you review
the resources listed below. Within these pages are links to additional
resources for learning Python. We recommend you look at them as part of
this self-study course.

First, read this page on the Landlab User Guide: `Python NumPy, SciPy,
Cython <https://github.com/landlab/landlab/wiki/Python%2C-NumPy%2C-SciPy%2C-Cython>`__

Then do this tutorial: \* `Introduction to Python and
NumPy <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/master/python_intro/Python_intro.ipynb>`__.
*Learn about:* The very basics of Python.

3. Introduction to Landlab
--------------------------

This is the core of this self-study course. You will alternate between
reading documentation on the User Guide, finding information in the
Reference Manual, and working through the tutorials.

The tutorials are presented as Jupyter notebooks, which contain a
mixture of text, images, and code blocks. When you look at the
tutorials, don’t just read them. Start by clearing the results by
selecting “Kernel ==> Restart & Clear Output,” then go ahead and try
running each code block as you come to it.

A motivating example
~~~~~~~~~~~~~~~~~~~~

-  `Introduction to Landlab: example model of fault-scarp
   degradation <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/master/fault_scarp/landlab-fault-scarp.ipynb>`__.
   A short overview of some of the things Landlab can do.

Introduction to the Landlab Grid and Fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, read the `User Guide page on Landlab
grids <https://github.com/landlab/landlab/wiki/Grid>`__. Then work
through the following tutorials:

-  `Introduction to the model grid
   object <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/master/grid_object_demo/grid_object_demo.ipynb>`__.
   Grid topology; how landlab represents data; connectivity of grid
   elements.
-  `Introduction to Landlab data
   fields <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/master/fields/working_with_fields.ipynb>`__.
   How Landlab stores spatial data on the grid; a little on naming
   conventions.

Now, go back to the [Hobley et al. 2017 publication] and identify the
ordering conventions of nodes, links, and other grid elements.

Plotting
~~~~~~~~

-  `Introduction to plotting output with
   Landlab <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/master/plotting/landlab-plotting.ipynb>`__.
   The basics of plotting with Landlab; combining matplotlib and out
   plots; the all-powerful ``imshow_grid()`` function.

How to Model with Landlab
~~~~~~~~~~~~~~~~~~~~~~~~~

Read the `Build-A-Model page in the User
Guide <https://github.com/landlab/landlab/wiki/Build-a-Model>`__.

Introduction to Components
~~~~~~~~~~~~~~~~~~~~~~~~~~

Read the `Component page in the User
Guide <https://github.com/landlab/landlab/wiki/Components>`__.

-  `Introduction to using the Landlab component
   library <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/master/component_tutorial/component_tutorial.ipynb>`__.
   The basics of working with and coupling components, using
   *diffusion*, *stream power*, and a *storm generator* as examples.

How to Use the Reference Manual
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Landlab Reference Manual contains documentation for most functions
in the Landlab package. It is the comprehensive counterpart to the
anecdotal tutorials.

Look at the `documentation for the
LinearDiffuser <http://landlab.readthedocs.io/en/latest/landlab.components.diffusion.html>`__,
which you just used in the prior tutorial.

Then spend some time (we recommend at least 30 minutes) clicking around
in the rest of the `Reference
Manual <http://landlab.readthedocs.io/en/latest/#developer-documentation>`__
getting a sense for what is there. Tip: to find a particular command,
click on Index and use your browser’s search function to search for a
command by name or keyword.

An important thing to appreciate about Components is that they often
make new fields that are used by other components. A very common example
of this is the FlowAccumulator and FlowDirector components.

Task: Before moving on to the next section identify two grid fields each
that the FlowAccumulator and FlowDirectors make.

Tutorials about specific Components
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some components are more sophisticated than others. Tutorials are
provided for many of these more elaborate components. You can find them
`near the bottom of the tutorials
page <https://github.com/landlab/landlab/wiki/Tutorials>`__ on the User
Guide.

Look at all tutorials on Flow Direction and Accumulation, and at least
one other component-specific tutorial based on your interests.

Interacting with the Landlab Developers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| You may find yourself with have a question to which you can’t find the
  answer in the User Guide or Reference Manual pages. Maybe you’ve
  stumbled on a bug (heaven forbid!). Or you might have a feature
  request. For such matters, the best way to communicate with the
  Landlab Developer Team
| is through `GitHub
  Issues <https://github.com/landlab/landlab/issues>`__: post an issue,
  and we’ll do our best to get back to you within 48 hours.

Task: Make an issue describing an improvement that you think should be
made to the Reference Manual Documentation based on your experience
reading it.

Advanced Grid and Fields: Gradients, Flux-Divergence, Mapping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to having lots of important information about adjacency of
nodes, links, and other grid elements, the Landlab Grid object has a
number of built-in functions for calculating quantities like gradients
and flux-divergence, and for mapping quantities from nodes to links and
so forth. Work through these tutorials to get a sense of this
functionality:

-  `Using the gradient and flux-divergence
   functions <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/master/gradient_and_divergence/gradient_and_divergence.ipynb>`__.
   Landlab as solving environment for staggered grid finite difference
   differential approximations; functions available to help you do this.
-  `Mapping values from nodes to
   links <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/master/mappers/mappers.ipynb>`__.
   Options for getting data on links to nodes, nodes to links, etc.;
   min, max, and mean; upwinding and downwinding schemes; one-to-one,
   one-to-many, and many-to-one mappings.

Boundary conditions
~~~~~~~~~~~~~~~~~~~

-  `Setting boundary conditions on Landlab grids (several
   tutorials) <https://nbviewer.jupyter.org/github/landlab/tutorials/tree/master/boundary_conds/>`__
   How Landlab conceptualizes boundary conditions; various ways to
   interact and work with them.

Working with Digital Elevtion Models (DEMs)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `Reading DEMs into
   Landlab <https://nbviewer.jupyter.org/github/landlab/tutorials/blob/master/reading_dem_into_landlab/reading_dem_into_landlab.ipynb>`__
   Getting a DEM in ESRI ASCII format into Landlab; getting the boundary
   conditions set right.

4. Capstone tutorial
--------------------

[There will eventually be an assignment here: an empty notebook with
instructions for the student to work through. Probably something like
Weathering, Depth Dependent HS transport, Stream Power, Accumulation.]
