In the previous section, we showed you most of the core functionality of the Landlab grid. In this section, we introduce you to how to actually use it.

Using Landlab requires that you build a Python script to import, instantiate, and then run your landscape model which you are implementing with Landlab. We describe such a script as a driver. It’s also possible to do the same set of processes on the fly in an interactive Python environment like iPython.)

Typically, a driver file will consist of six distinct sections:

Import the Python and Landlab libraries you’ll need to run your model
Instantiate the Landlab elements (grid and, if using them, components)
Load any necessary data into the grid fields
Set the boundary conditions
Run the model, typically by creating a for loop or using a Landlab generator (see below)
Finalize and handle the data (e.g., plot, export)

Beyond the driver, if you’re using Landlab components, you’ll probably also need a parameter file. This file supplies the components with the additional parameter and setup information they need. Landlab parameter files are text files (.txt), have fixed format, and for convenience (so you only have to specify the minimum of path information in the file name) should be placed in the same folder as the driver file. Find out more about driver files here (XXX LINK to under component section). However, if you’re not using components, there’s little need to create a parameter file; you can just directly other parameters to the grid in the driver. 


What is a Landlab component?

A key strength of Landlab is that not only is it designed to make implementing your own process simulations as simple as possible, it also offers an off-the-shelf library of pre-designed process descriptions that you can use in your drivers. We call these process simulators Landlab components. The intention is that each component be:

Plug-and-play
Interoperable with all other components
Implementable in your driver in only one or two lines of code

Note that by no means is using the component library necessary or even always desirable when working with Landlab. However, we hope that their availability and ease of use will dramatically reduce the time investment needed to implement a wide variety of modelling scenarios. In particular, components should make production of models coupling more than one process significantly easier, as existing, off-the-shelf components can be slotted in alongside novel process descriptions.

At the time of writing, the existing library contains the following components:

Detachment limited/stream power channel incision
Transport limited channel incision
Linear (hillslope) diffusion
Nonlinear hillslope diffusion
D8 /Dn (D8 generalized to Voronoi) flow routing
Simple and not-so-simple crustal flexure
A simple fire simulator
A thin-ice glacial approximation
Precipitation/evapotranspiration
A relatively simple router for overland flow
A radiative intensity calculator
A generalized framework for cellular automata
A vegetation cellular automaton
A hillslope particle cellular automaton

Note that not all components will run under all conditions, but that any limitations should be made clear in the inline documentation associated with that component (access help either through the indices you can find on this site (XXX LINK) or by typing “[component or method]?” in an interactive Python session). In particular, some components may demand you are running on a regular grid. It should probably also be emphasised that most of these components are still under active development within this beta release of Landlab, and may behave in idiosyncratic ways or be subject to sudden changes with little or no warning. In all cases, we’d recommend contacting the original coder of the component to let them know they have external users to think about before setting out on any major research challenges using it!


Implementing a Landlab driver

As noted above, the process of creating a driver is essentially equivalent whether you want to implement Landlab components, purely your own code, or some mixture of the two. Here we take a closer look at the various steps.

1. Import the libraries and functions you need

Landlab handles a lot like numpy, and like numpy you’ll need to import the various libraries and functions that you’ll want to use. At the very least, we suspect you’ll need from outside Landlab:
numpy itself
rudimentary pylab plotting routines: plot, show, figure
Also useful can be:
the Python module time, to time various parts of your code
elements from SciPy, the scientific computing library. Lots of useful methods (e.g., matrix solutions, curve fitting) can be found in here, to avoid reinventing the wheel.

From inside Landlab, you’ll also need
A grid class. Choose from RasterModelGrid, VoronoiDelaunayGrid, or some of the more specialized classes
Any components you want to run
Any Landlab utilities you need, like the plotters (imshow_node_grid) or io functions

A specific example might be:

>>> import numpy as np
>>> from pylab import show, figure, plot
>>> import time
>>> from landlab import RasterModelGrid
>>> from landlab.components.flow_routing import route_flow_dn
>>> from landlab.plot.imshow import imshow_node_grid


2. Instantiate objects

As noted above, Landlab is coded in an object-oriented (XXX LINK) style. This means that we need to “instantiate” the various Landlab objects that we will use to store data and run the model. The grid and all the components are objects, so we need to instantiate them next.

Note that most components require the grid object be passed to them as one of their arguments during instantiation, so the first thing you’ll want to instantiate will be the grid.

Check the docstrings for each class (grid, component) you want to instantiate for a detailed description of what you need to supply as arguments. For a RasterModelGrid, this will be (number_of_node_rows, number_of_node_columns, node_spacing(optional)). For a VoronoiDelaunayGrid, it will be (array_of_node_x_coords, array_of_node_y_coords). For a generic component, it will typically be (ModelGrid, ‘path_to_parameter_file.txt’), though there may be some variation, and optional inputs may also be available.

Give each object you instantiate a variable name. We like “mg” for ModelGrid objects, and some appropriate abbreviation for a component.

An example might be:
>>> mg = RasterModelGrid(10,10,1.) #100 nodes, spacing of 1.
>>> fr = route_flow_dn(mg, ‘./params.txt’) #this assumes params.txt is in the current directory


3. Load/create data in fields

(see this section (XXX LINK) if you don’t know what a Landlab field is)

Now we need some data to work with.
