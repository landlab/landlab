==============================================
Using Landlab's coupled component capabilities
==============================================

Landlab is designed as a modular modelling tool. Existing component modules describing
individual geological, geomorphological, and climatic processes can be loaded into a
landlab driver script and run with the minimum of fuss. This document briefly illustrates
how this works in practice, using examples of linear and nonlinear diffusion, dN flow 
routing, and a stream power based fluvial erosion component.

Note that the advantages of many of these modules is that they implement implicit methods,
or are otherwise accelerated above what could be quickly coded up by hand using the
grid methods. In particular, 
:class:`landlab.components.flow_routing.route_flow_dn.FlowRouter` 
is considerably faster than the equivalent grid method 
:method:`landlab.RasterModelGrid.calculate_steepest_descent_across_adjacent_cells`,
exploiting the Braun-Willett (2013) "Fastscape" flow routing algorithms.

Typically, components require you to pass them the model grid object as an input, and
often also a variable representing time elapsed. They normally also demand that you are
storing data on the grid as fields, and might require you to use specific naming 
conventions for those fields. This should all be documented in the docstrings for any
given component.


Example 1: Using the diffusion components
=========================================

A typical Landlab driver has the structure:

    - Import Python packages and Landlab components
    - Import parameters needed to run the model from a .txt file using the
        ModelParameterDictionary
    - Create an instance of the grid
    - Set up that grid, e.g., boundary conditions, data on the grid
    - Instantiate the component objects
    - Create a loop to actually run the processes, and update the BCs (e.g., uplift)
    - Output, e.g., save data, prep data for plotting, and plot

Open up **landlab/examples/diffusion_driver_2.py** and have a look at it to see how it
follows this structure. It's heavily commented for clarity. It's set up to allow you
to alternate between linear and nonlinear diffusive models, hence the fact we import
both components in lines 1 & 2.

Run the model with linear diffusion, then run it with nonlinear diffusion by switching
which of lines 55 and 56 is commented out. Note that both of these modules are set up to
include uplift internally, so it doesn't need to be manually added in the loop. Feel free
to play around with the parameters in the input file (which we identify in line 10).

It's also worth taking a look at the plotting routines seen from lines 59 onwards:

.. literalinclude:: landlab.examples.diffusion_driver_2
    :linenos:
    :language: python
    :lines: 48, 55-56, 58-82
    
Note the ways we can plot both cross-sections and raster views of the data with minimal
code, the former exploiting the power of Python's `slicing syntax 
<http://wiki.scipy.org/Tentative_NumPy_Tutorial#head-864862d3f2bb4c32f04260fac61eb4ef34788c4c>`_.

