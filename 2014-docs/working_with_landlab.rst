.. _drive_a_model:

What goes into a Landlab Model?
===============================

In the previous section, we showed you most of the core functionality of the Landlab grid. In this section, we introduce you to how to actually use it, to build models and work with the Landlab component library.

Using Landlab requires that you build a Python script to import, instantiate, and then run your landscape model. We describe such a script as a **driver**.  It’s also possible to do the same set of processes on the fly in an interactive Python environment like iPython.

Typically, a driver file will consist of six distinct sections:

* **Import** the Python and Landlab libraries you’ll need to run your model
* **Instantiate** the Landlab elements (grid and, if using them, components)
* **Load** any necessary data into the grid fields
* Set the **boundary conditions**
* **Run** the model, typically by creating a for loop or using a Landlab generator (see below)
* **Finalize** and handle the data (e.g., plot, export)

Beyond the driver, if you’re using Landlab components, you’ll probably also need a **parameter file**. This file supplies the components with the additional parameter and setup information they need. Landlab parameter files are text files (.txt), have fixed format, and for convenience (so you only have to specify the minimum of path information in the file name) should be placed in the same folder as the driver file. Find out more about parameter files :ref:`here <input_files>`. However, if you’re not using components, there’s little need to create a parameter file; you can just directly other parameters to the grid in the driver.


A Brief Introduction to Components
----------------------------------

A key strength of Landlab is that not only is it designed to make implementing your own process simulations as simple as possible, it also offers an off-the-shelf library of pre-designed process descriptions that you can use in your drivers. We call these process simulators Landlab **components**. The intention is that each component be:

* Plug-and-play
* Interoperable with all other components
* Implementable in your driver in only one or two lines of code

Note that by no means is using the component library necessary or even always desirable when working with Landlab. However, we hope that their availability and ease of use will dramatically reduce the time investment needed to implement a wide variety of modelling scenarios. In particular, components should make production of models coupling more than one process significantly easier, as existing, off-the-shelf components can be slotted in alongside novel process descriptions.

At the time of writing, the existing library contains the following components:

* Detachment limited/stream power channel incision
* Transport limited channel incision
* Linear (hillslope) diffusion
* Nonlinear hillslope diffusion
* D8 /Dn (D8 generalized to Voronoi) flow routing
* Simple and not-so-simple crustal flexure
* A simple fire simulator
* A thin-ice glacial approximation
* Precipitation/evapotranspiration
* A relatively simple router for overland flow
* A radiative intensity calculator
* A generalized framework for cellular automata
* A hillslope particle cellular automaton

And under active development are:

* A vegetation cellular automaton
* An impact cratering simulator
* Divergent and mixed convergent-divergent flow routers
* A deltaic simulator

Note that not all components will run under all conditions, but that any limitations should be made clear in the inline documentation associated with that component (access help either through the indices you can find `on this site <http://landlab.readthedocs.org/en/latest/users_guide.html#list-of-landlab-components>`_ or by typing “[component or method]?” in an interactive Python session). In particular, some components may demand you are running on a regular grid. It should probably also be emphasised that most of these components are still under active development within this beta release of Landlab, and may behave in idiosyncratic ways or be subject to sudden changes with little or no warning. In all cases, we’d recommend contacting the original coder of the component to let them know they have external users to think about before setting out on any major research challenges using it!


Implementing a Landlab driver
-----------------------------

As noted above, the process of creating a driver is essentially equivalent whether you want to implement Landlab components, purely your own code, or some mixture of the two. Here we take a closer look at the various steps.

1. Import the libraries and functions you need
++++++++++++++++++++++++++++++++++++++++++++++

Landlab handles a lot like numpy, and like numpy you’ll need to import the various modules and functions that you’ll want to use. At the very least, we suspect you’ll need from outside Landlab:

* *numpy* itself
* rudimentary pylab plotting routines: *plot*, *show*, *figure*

Also useful can be:

* the Python module *time*, to time various parts of your code
* elements from *SciPy*, the scientific computing library. Lots of useful methods (e.g., matrix solutions, curve fitting) can be found in here, to avoid reinventing the wheel.

From inside Landlab, you’ll also need:

* A **grid** class. Choose from :class:`~landlab.grid.raster.RasterModelGrid`, :class:`~landlab.grid.voronoi.VoronoiDelaunayGrid`, or some of the more specialized classes
* Any **components** you want to run
* Any Landlab **utilities** you need, like the plotters (:func:`~landlab.plot.imshow.imshow_node_grid`) or **io** functions

A specific example might be:

.. code-block:: python

    import numpy as np
    from pylab import show, figure, plot
    import time
    from landlab import RasterModelGrid
    from landlab.components.flow_routing import FlowRouter
    from landlab.plot.imshow import imshow_node_grid


2. Instantiate objects
++++++++++++++++++++++

As noted in previous sections, Landlab is coded in an `object-oriented style <http://code.tutsplus.com/articles/python-from-scratch-object-oriented-programming--net-21476>`_. This means that we need to “instantiate” the various Landlab objects that we will use to store data and run the model. The grid and all the components are objects, so we need to instantiate them next.

Note that most components require the grid object be passed to them as one of their arguments during instantiation, so the first thing you’ll want to instantiate will be the grid.

Check the docstrings for each class (grid, component) you want to instantiate for a detailed description of what you need to supply as arguments. For a RasterModelGrid, this will be a tuple, (i, j), where i is the number of columns and j the number of rows. A second float (or tuple) is optional, specifying the node spacing - using a tuple to specify (dy, dx) if you want them to be different.
[Landlab also recognises an older style of RasterModelGrid signature, which looks like (number_of_node_rows, number_of_node_columns, node_spacing(optional)), and is clever enough to work out this is what you're doing if your arguments are of this form.]
For a VoronoiDelaunayGrid, it will be (array_of_node_x_coords, array_of_node_y_coords). For a generic component, it will typically be (ModelGrid, ‘path_to_parameter_file.txt’), though there may be some variation, and optional inputs may also be available.

Give each object you instantiate a variable name. We like “mg” for ModelGrid objects, and some appropriate abbreviation for a component.

An example might be:

.. code-block:: python

    mg = RasterModelGrid((10, 10), (1., 2.))  # 100 nodes, dy=1., dx=2.
    fr = FlowRouter(mg, './params.txt')
    # ...this assumes params.txt is in the current directory


3. Load/create data in fields
+++++++++++++++++++++++++++++

(:ref:`See this section <fields>` if you don’t know what a Landlab field is)

Now we need some data to work with. Here we’ll assume that you’re going to be working with a DEM-style elevation map across the nodes of the grid, but similar considerations would apply for any other type of data.

You will likely be in one of two situations regarding the initial data you want to put on the grid - either you will have some external data source that you want to load in and use as your initial conditions (e.g., a DEM of some basin, or some other real topography), or you want to set up some simple analytical initial condition like a flat surface with noise or an inclined surface.

In both cases, we advocate a two step process: creating a numpy array of the data, then loading it into the grid as a field. We can illustrate for both of the above cases:

.. code-block:: python

    mg = RasterModelGrid((10, 10), 1.)  # make a grid
    z = np.zeros(100, dtype=float)  # make a flat surface, elev 0
    # or…
    z = mg.node_y*0.01  # a flat surface dipping shallowly south
    # add a little noise to the surface:
    z += np.random.rand(100.)/10000.
    # create the field:
    mg.add_field('node', 'topographic__elevation', z, units='m')

Alternatively, we can use the specialized Landlab function :func:`~landlab.io.esri_ascii.read_esri_ascii`
to import an ascii raster that can be output from ARC. Note this function both creates the grid for you and loads the data as a field if you provide ‘name’. If not, you’ll have to load the data output (*z*, below) manually::

    >>> from landlab.io import read_esri_ascii
    >>> mg, z = read_esri_ascii('my_ARC_output.asc', name='topographic__elevation')
    >>> np.all(mg.at_node['topographic__elevation'] == z)
    True

Note that if you don’t want to use any Landlab components, you can continue to work with data as “free floating” numpy arrays, and can ignore the fields (e.g., see the :ref:`simple tutorial <getting_started>` at the start of this guide).


4. Set the boundary conditions
++++++++++++++++++++++++++++++

Once you have a grid and the initial condition data you’ll need, it’s time to set the boundary conditions. If you’re working with a raster, or some pre-existing imported data, this is very straightforward using the built in RasterModelGrid functions. For a raster where only the edges are to be boundary nodes::

    >>> mg.set_fixed_value_boundaries_at_grid_edges(False, True, False, True)
    >>> mg.set_closed_boundaries_at_grid_edges(True, False, True, False)

This will give a grid with fixed value boundaries at the left and right edges, and closed boundaries at the top and bottom.

If you’re working with, say, an ARC imported array with a null value on the closed nodes
(e.g., -9999), you can do this::

    >>> mg.set_nodata_nodes_to_closed(mg.at_node['topographic__elevation'], -9999)

(Note though that you’re still likely to have to reopen an outlet node manually! In which case you’ll also need to follow the instructions below.)

If you’re working with individual nodes’ boundary statuses, you’ll need to set the BCs by hand. This means individually modifying the boundary condition status of each node or link that you want to be of the new type. Fortunately, Landlab uses some Python magic to make sure that when you update, for example, the status of a node, the statuses of attached links and cells change concomitantly. For example::

    >>> # import the BC values we'll need:
    >>> from landlab import FIXED_LINK, FIXED_GRADIENT_BOUNDARY
    >>> # find the ID of the lowest elevation core node.
    >>> # we'll make this a fixed gradient outlet:
    >>> outlet_id = mg.core_nodes[np.argmin(
                    mg.at_node['topographic__elevation'][mg.core_nodes])]
    >>> # show there are no FIXED_LINK boundary conditions in the grid yet:
    >>> np.any(mg.status_at_link==FIXED_LINK)
    False
    >>> # update the outlet node:
    >>> mg.status_at_node[outlet_id] = FIXED_VALUE_BOUNDARY
    >>> np.any(mg.status_at_link==FIXED_LINK)
    True
    >>> # the corresponding link has been automatically updated.


5. Run the model
++++++++++++++++

We’re now ready to actually implement a run of our model! Most things you might want to do with Landlab are probably time-sensitive, so in almost all cases, you’ll probably be placing the guts of your simulation inside a loop of some sort. In simple cases, you can just use some variation on a simple for loop or while statement, either:

.. code-block:: python

    dt = 10.
    for tstep in xrange(100):
        #...do the thing for one timestep dt

or:

.. code-block:: python

    dt = 10.
    accumulated_time = 0.
    while accumulated_time<1000.:
        #...do the thing for one timestep dt
        accumulated_time += dt

Both produce 1000 time units of run, with an explicit timestep of 10. Notice that the latter technique is particularly amenable to situations where your explicit timestep is varying (e.g., a storm sequence).

Landlab also however has a built in storm generator component,
:class:`~landlab.components.uniform_precip.generate_uniform_precip.PrecipitationDistribution`,
which acts as a true `Python generator <http://www.python-course.eu/generators.php>`_ (see its documentation!). The main method is
:func:`~landlab.components.uniform_precip.generate_uniform_precip.PrecipitationDistribution.yield_storm_interstorm_duration_intensity`.
This means producing a storm series in Landlab is also very easy:

.. code-block:: python

    from landlab.components.uniform_precip import PrecipitationDistribution
    time_to_run = 500000.
    precip_perturb = PrecipitationDistribution(input_file=input_file_string,
                                               total_t=time_to_run)
    for (interval_duration, rainfall_rate) in \
            precip_perturb.yield_storm_interstorm_duration_intensity():
        if rainfall_rate != 0.:
            # ...do the thing, making sure to pass it the current
            # interval_duration and rainfall_rate

Notice that the advantage of the generator is that it just stops when the desired number of events/time duration has expired! See the end of `this tutorial <http://nbviewer.ipython.org/github/landlab/drivers/blob/master/notebooks/component_tutorial.ipynb>`_ for an example of this generator in action.


What exactly “...do the thing” consists of is up to you. You can either design your own operations to do in the loop for yourself, or you can implement processes from Landlab’s component library. See :ref:`here <landlab_components_page>` for more information on using the components.


6. Finalize and handle the data
+++++++++++++++++++++++++++++++

Once the looping is complete, the model is effectively finished. However, you will still need to output the data somehow! Some obvious options are:

Save or export the data
^^^^^^^^^^^^^^^^^^^^^^^

If you’re using a raster grid, you can easily save your grid output to either ESRI ascii (i.e. ARCmap) or open source netCDF formats. netCDF in particular is a powerful format, and allows easy subsequent re-loading of a Landlab modelgrid and all its fields. Save your raster like this:

.. code-block:: python

    rmg.save('my_savename.asc', names=['field1','field2'])
    # for esri ascii, only saving the fields 1 and 2

or:

.. code-block:: python

    rmg.save('my_savename.nc')
    # save as netCDF3, saving all fields by default

The former way will give two save files, my_savename_field1.asc and my_savename_field2.asc. The latter will just give ‘my_savename.nc’.

To reload a netCDF file, use the Landlab io function
:func:`~landlab.io.netcdf.read_netcdf`::

    >>> from landlab.io.netcdf import read_netcdf
    >>> mg = read_netcdf('my_savename.nc')

Note all the original fields you had will automatically be repopulated.

If you’re using an irregular grid, the simple grid save function is not yet operational (though is under development). Instead, we recommend using Pickle, a native Python way of saving (“pickling”) any Python object. It works like this::

    >>> import cPickle as pickle  # cPickle is a lot faster than normal pickle
    >>> pickle.dump( mg, open('my_savename.pickle', 'wb') )
    # ...save the grid, and all its fields
    >>> mg = pickle.load( open('my_savename.pickle', 'rb') )
    # ...load the grid and fields back into a grid object

Unfortunately, the power of pickle comes somewhat at the expense of both disk space and speed. Saves this way can be slow and, if the grid is big, memory expensive (e.g., ~1 Gb for millions of nodes).

You can also use lower level, numpy save routines to preserve just your data (rather than the whole grid object). The numpy methods **save** and **savetxt** and **load** and **loadtxt** can be called on any numpy array, including those saved as fields. Save and load use the numpy specific .npy file format; savetxt and loadtxt use textfiles. Use them like this::

    >>> np.save('savename.npy', mg.at_node['my_field'])
    >>> mg.at_node['my_field'] = np.load('savename.npy')
    >>> np.savetxt('savename.txt', mg.at_node['my_field'])
    >>> mg.at_node['my_field'] = np.loadtxt('savename.txt')


Plot the data
^^^^^^^^^^^^^

Landlab has a fairly comprehensive suite of built in plotting functions; read more about them :ref:`here <plotting>`.

You also of course have the option of using the `matplotlib plotting library <http://matplotlib.org>`_ of Python for things like cross-sections.

If you’re careful, you can also build plotting functions into the body of a run loop for your model, so you can see how your output evolves through time. Note however that all Python save and plot functions are considerably time expensive, so it would probably be a bad idea to do this kind of thing every timestep. Instead, you can try something like:

.. code-block:: python

    import plot
    dt = 10.
    accumulated_time = 0.
    last_accumulated_time_remainder = 0.
    while accumulated_time<1000.:
        #...do the thing for one timestep dt
        accumulated_time += dt
        if last_accumulated_time_remainder < accumulated_time%100.:  # output every 100.
            plot(mg.node_vector_to_raster(z)[mg.number_of_node_rows//2,:])  # a cross section
            last_accumulated_time_remainder = accumulated_time%100.
    show()

Note that if you’re running inside an interactive Python session like iPython, all the variables and objects (both grid and component) that you’ve used in your model will still be available in the environment. Thus, you can play with your data for as long as you want!
