Tutorial: Modeling Diffusion on a Raster Grid
=============================================

The following is a simple example in which we use **ModelGrid** to build an explicit, 
finite-volume, staggered-grid model of diffusion. The mathematics of diffusion describe 
several different phenomena, including heat conduction in solids, chemical diffusion 
of solutes, transport of momentum in a viscous shear flow, and transport of 
soil on hillslopes. To make this example concrete, we will use the hillslope evolution as 
our working case study, though in fact the solution could apply to any of these systems.

To work through this example, you can type in and run the code below, or you can download
:download:`the example script <../model_grid_guide/diffusion_with_model_grid.py>`_.
The complete source code for the diffusion model is listed 
below. Line numbers are included to make it easier to refer to particular lines of code 
(of course, these numbers are not part of the source code). After the listing, we will take 
a closer look at each piece of the code in turn. Output from the the diffusion model is shown in 
:ref:`Figure 3 <diff1>`.

.. code-block:: python

	#! /usr/env/python
	"""

	2D numerical model of diffusion, implemented using Landlab's ModelGrid module.
	Provides a simple tutorial example of ModelGrid functionality.

	Last updated GT May 2014

	"""

	from landlab import RasterModelGrid
	import pylab, time

	def main():
		"""
		In this simple tutorial example, the main function does all the work: 
		it sets the parameter values, creates and initializes a grid, sets up 
		the state variables, runs the main loop, and cleans up.
		"""
	
		# INITIALIZE
	
		# User-defined parameter values
		numrows = 20          # number of rows in the grid
		numcols = 30          # number of columns in the grid
		dx = 10.0             # grid cell spacing
		kd = 0.01             # diffusivity coefficient, in m2/yr
		uplift_rate = 0.0001  # baselevel/uplift rate, in m/yr
		num_time_steps = 10000 # number of time steps in run
	
		# Derived parameters
		dt = 0.1*dx**2 / kd    # time-step size set by CFL condition
	
		# Create and initialize a raster model grid
		mg = RasterModelGrid(numrows, numcols, dx)
	
		# Set the boundary conditions
		mg.set_closed_boundaries_at_grid_edges(False, False, True, True)

		# Set up scalar values
		z = mg.add_zeros('node', 'Elevation')            # node elevations
	
		# Get a list of the core cells
		core_cells = mg.get_core_cell_node_ids()

		# Display a message, and record the current clock time
		print( 'Running diffusion_with_model_grid.py' )
		print( 'Time-step size has been set to ' + str( dt ) + ' years.' )
		start_time = time.time()

		# RUN
	
		# Main loop
		for i in range(0, num_time_steps):
		
			# Calculate the gradients and sediment fluxes
			g = mg.calculate_gradients_at_active_links(z)
			qs = -kd*g
		
			# Calculate the net deposition/erosion rate at each node
			dqsds = mg.calculate_flux_divergence_at_nodes(qs)
		
			# Calculate the total rate of elevation change
			dzdt = uplift_rate - dqsds
			
			# Update the elevations
			z[core_cells] = z[core_cells] + dzdt[core_cells] * dt


		# FINALIZE

		# Get a 2D array version of the elevations
		zr = mg.node_vector_to_raster(z)
	
		# Create a shaded image
		pylab.close()  # clear any pre-existing plot
		im = pylab.imshow(zr, cmap=pylab.cm.RdBu, extent=[0,numcols*dx,0,numrows*dx],
						  origin='lower')
		# add contour lines with labels
		cset = pylab.contour(zr, extent=[0,numcols*dx,numrows*dx,0], hold='on',
							 origin='image')
		pylab.clabel(cset, inline=True, fmt='%1.1f', fontsize=10)
	
		# add a color bar on the side
		cb = pylab.colorbar(im)
		cb.set_label('Elevation in meters')
	
		# add a title and axis labels
		pylab.title('Simulated topography with uplift and diffusion')
		pylab.xlabel('Distance (m)')
		pylab.ylabel('Distance (m)')

		# Display the plot
		pylab.show()
		print('Run time = '+str(time.time()-start_time)+' seconds')

	if __name__ == "__main__":
		main()
	   

.. _diff1:

.. figure:: images/basic_diffusion_example.png
    :figwidth: 80 %
    :align: center

    Figure 4: Output from the hillslope diffusion model.
    
Below we explore how the code works line-by-line.

Importing Packages
>>>>>>>>>>>>>>>>>>

.. code-block:: python

	from landlab import RasterModelGrid
	import pylab, time

We start by importing the grid class ``RasterModelGrid`` from the ``landlab`` package (note that the ``landlab`` package must first be installed; see instructions under :ref:`Installing Landlab <install>`). We'll also import ``pylab`` so we can plot the results, and ``time`` so we can report the time it takes to run the model.

Setting the User-Defined Parameters
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

.. code-block:: python

		# User-defined parameter values
		numrows = 20          # number of rows in the grid
		numcols = 30          # number of columns in the grid
		dx = 10.0             # grid cell spacing
		kd = 0.01             # diffusivity coefficient, in m2/yr
		uplift_rate = 0.0001  # baselevel/uplift rate, in m/yr
		num_time_steps = 10000 # number of time steps in run

The first thing we'll do in the ``main()`` function is set a group of user-defined parameters. The size of the grid is set by ``numrows`` and ``numcols``, with cell spacing ``dx``. In this example, we have a 20 by 30 grid with 10 m grid spacing, so our domain represents a 200 by 300 m rectangular patch of land. The diffusivity coefficient ``kd`` describes the efficiency of soil creep, while the ``uplift_rate`` indicates how fast the land is rising relative to base level along its boundaries. Finally, we set how many time steps we want to compute.

Note that the code for our simple program lives inside a ``main()`` function. This isn't strictly necessary---we could have put the code in the file without a ``main()`` function and it would work just fine when we run it---but it is good Python practice, and will be helpful later on.

Calculating Derived Parameters
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

.. code-block:: python

		# Derived parameters
		dt = 0.1*dx**2 / kd    # time-step size set by CFL condition
		
Next, we calculate the values of parameters that are derived from the user-defined parameters. In this case, we have just one: the time-step size, which is set by the Courant-Friedrichs-Lewy condition for an explicit, finite-difference solution to the diffusion equation (to be on the safe side, we multiply the ratio :math:`\Delta x^2 / k_d` by 0.1 instead of the theoretical limit of 1/2). With the parameter values above, :math:`\Delta t = 1000` years, so our total run duration will be one million years. Remember, though, that the same code could be used for any diffusion application with a source term. For instance, we could model conductive heat flow, with :math:`k_d` representing thermal diffusivity and ``uplift_rate`` representing heat input by, for example, radioactive decay in the earth's crust.


Creating and Configuring the Grid
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

.. code-block:: python

    # Create and initialize a raster model grid
    mg = RasterModelGrid(numrows, numcols, dx)
    
    # Set the boundary conditions
    mg.set_closed_boundaries_at_grid_edges(False, False, True, True)

Our model grid is created with a call to ``RasterModelGrid()``. (Object-oriented programmers will recognize this as the syntax for creating a new object---in this case a
raster model grid.) The variable ``mg`` now contains a ``RasterModelGrid`` object that has
been configured with 20 rows and 30 columns.

For our boundary conditions, we would like to keep the nodes along the bottom and right edges of the grid fixed at zero elevation. We also want to have the top and left boundaries represent ridge-lines with a fixed horizontal position and no flow of sediment in or out. To accomplish this, we call the ``set_closed_boundaries_at_grid_edges()`` method. (Note: the term *method* is object-oriented parlance for a function that belongs to a particular class of object, and is always called in reference to a particular object). The method takes four boolean arguments, which indicate whether there should be closed boundary condition on the bottom, right, top, and left sides of the grid. Here we have set the flag to ``True`` for the top and left sides. This means that the links connecting the interior nodes to the perimeter nodes along these two sides will be flagged as inactive, just as illustrated (with a smaller grid) in :ref:`Figure 3 <raster4x5openclosed>`. As we'll see in a moment, we will simply not bother to calculate any mass flux across these closed boundaries.


Creating Data
>>>>>>>>>>>>>

.. code-block:: python

    # Set up scalar values
    z = mg.add_zeros('node', 'Elevation')            # node elevations

Our state variable, :math:`z(x,y,t)`, represents the land surface elevation. One of the unique aspects of ModelGrid is that grid-based variables like :math:`z` are represented as 1D rather than 2D Numpy arrays. Why do it this way, if we have a regular grid that naturally lends itself to 2D arrays? The answer is that we might want to have an irregular, unstructured grid, which is much easier to handle with 1D arrays of values. By using 1D arrays for all types of ModelGrid, we allow the user to switch seamlessly between structured and unstructured grids.

We create our data structure for :math:`z` values with  ``add_zeros()``, a ModelGrid method that creates and returns a 1D Numpy array filled with zeros (behind the scenes, it also "attaches" the array to the grid; we'll see later why this is useful). The length of the array is equal to the number of nodes in the grid (:math:`20\times 30=600`), which makes sense because we want to have an elevation value associated with every node in the grid.

When we update elevation values, we will want to operate only on the core nodes. To help with this, we use the ``core_nodes`` property (a *property* in Python is essentially a variable that belongs to an object, which you can access but not modify directly). This property contains a 1D numpy array of integers that represent the node ID numbers associated with all of the core nodes (of which there are :math:`18\times 28 = 504`). Finally, we display a message to tell the user that we're about to run and with what time step size.

Main Loop
>>>>>>>>>

.. code-block:: python

		# Main loop
		for i in range(0, num_time_steps):

Our model implements a finite-volume solution to the diffusion equation. The idea here is that we calculate sediment fluxes around the perimeter of each cell. We then integrate these fluxes forward in time to calculate the net change in volume, which is divided by the cell's surface area to obtain an equivalent change in height. The numerical solution is given by:

.. math::

	\begin{equation}
	\frac{d z_i}{dt} \approx \frac{z^{T+1}_i-z^T_i}{\Delta t}
	= - \frac{1}{\Lambda_i} \sum_{j=1}^{N_i} \mathbf{q}_{Sij}^T \lambda_{ij}.
	\label{eq:dzdt}
	\end{equation}

Here, :math:`z_i^T` is the elevation at node :math:`i` at time step :math:`T`, :math:`t` is time, :math:`\Lambda_i` is the surface area of cell :math:`i`, :math:`N_i` is the number of nodes adjacent to :math:`i` (called the node's *neighbors*), :math:`\mathbf{q}_{Sij}^T` is the sediment flux per unit face width from cell :math:`i` to cell :math:`j`, and :math:`\lambda_{ij}` is the width of the face between cells :math:`i` and :math:`j`. The flux between a pair of adjacent cells is the product of the slope (positive upward) between their associated nodes, :math:`\mathbf{S}_{ij}`, and a transport coefficient, :math:`k_d`,

.. math::

	\begin{equation}
	\mathbf{q}_{Sij} = - k_d \mathbf{S}_{ij} = - k_d \frac{z_j-z_i}{L_{ij}}
	\end{equation}

where :math:`L_{ij}` is the length of the link connecting nodes :math:`i` and :math:`j`. Notice that elevation values (which are scalars) are associated with nodes, while slopes and sediment fluxes (which are vectors) are associated with links and faces. If we want to think of the slopes and fluxes as being calculated at a particular point, that point is the junction between a link and its corresponding face :ref:`Figure 1 <grid>`.

Because we are using a regular (raster) grid with node spacing :math:`\Delta x`, the face width and link length are both equal to :math:`\Delta x` everywhere, and the cell area :math:`\Lambda=\Delta x^2`. This would not be true, however, for an unstructured grid.

Calculating gradients and sediment fluxes
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

.. code-block:: python
		
			# Calculate the gradients and sediment fluxes
			g = mg.calculate_gradients_at_active_links(z)
			qs = -kd*g

In order to calculate new elevation values, the first quantity we need to know is the gradient (slope) values between all the node pairs. We can calculate this in a single line of code using ModelGrid's ``calculate_gradients_at_active_links()`` method. This method takes a single argument: a 1D numpy array of scalar values associated with nodes. The length of this array must be the same as the number of nodes in the grid. The method calculates the gradients in ``z`` between each pair of nodes. It returns a 1D numpy array, ``g`` (for gradient), the size of which is the same as the number of active links in the grid (the difference between active and inactive links is illustrated in :ref:`Figure 2 <raster4x5>` and :ref:`3 <raster4x5openclosed>`). The sign of each value of ``g`` is positive when the slope runs uphill from a link's *from-node* to its *to-node*, and negative otherwise.

To calculate the sediment fluxes, we multiply each gradient value by the transport coefficient ``kd``. The minus sign simply means that the sediment goes downhill: where the gradient is negative, the flux should be positive, and vice versa. Here, we are taking advantage of numpy's ability to perform mathematical operations on entire arrays in a single line of code, rather than having to write out a ``for`` loop. The line ``qs = -kd*g`` in our code multiplies ``ks`` by every value of ``g``, and returns the result as a numpy array the same size as ``g``.

Calculating net fluxes in and out of cells
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

.. code-block:: python
		
			# Calculate the net deposition/erosion rate at each node
			dqsds = mg.calculate_flux_divergence_at_nodes(qs)
		
Now that we know the unit flux associated with each link and its corresponding cell face, the next thing we need to do is add up the total flux around the perimeter of each cell. In other words, we need to calculate the summation in equation above. Landlab allows us to do this in one line of code, by calling the ``calculate_flux_divergence_at_nodes()`` method. This method takes a single argument: a 1D Numpy array containing the flux per unit width at each face in the grid. The method multiplies each unit flux by its corresponding face width, adds up the fluxes across each face for each cell, and divides the result by the surface area of the cell. It returns a 1D Numpy array that contains the net rate of change of volume per unit cell area. The length of this array is the same as the number of nodes in the grid. We will store the result in ``dqsds``.

If the boundary nodes around the grid perimeter do not have associated cells, why do we bother calculating net fluxes for them? In fact, we do not need to; we could have called the method ``calculate_flux_divergence_at_core_cells()`` instead. This would have given us an array the length of the number of core cells, not nodes (there is one every core node has a corresponding core cell). There are two reasons to do the net flux calculation at all nodes. The first is simply that the node-based method is slightly faster than the cell-based version. The second is that by using nodes, we retain some information about the flow of mass into the open boundary nodes. This could be useful in testing whether our model correctly balances mass (though we do not actually use that capability in this example).

Updating elevations
>>>>>>>>>>>>>>>>>>>

.. code-block:: python
		
			# Calculate the total rate of elevation change
			dzdt = uplift_rate - dqsds
			
			# Update the elevations
			z[core_cells] = z[core_cells] + dzdt[core_cells] * dt

When we calculated flux divergence, we got back an array of numbers, ``dqsds``, that represents the rate of gain or loss of sediment volume per unit area at each node. Now we need to combine this information with the source term---representing vertical motion of the soil relative to the base level at the model's fixed boundaries---in order to calculate the total rate of elevation change at the nodes. Once we've calculated rates of change, we update all node elevations by simply multiplying ``dzdt`` by our time step size. We do not want to change the elevations of the boundary nodes, however, and so we perform the update only on the interior cells. Because we are using numpy arrays, we can isolate the core nodes simply by putting our array of node IDs for core nodes inside square brackets. 


Plotting the Results
>>>>>>>>>>>>>>>>>>>>

.. code-block:: python

		# Get a 2D array version of the elevations
		zr = mg.node_vector_to_raster(z)
	
		# Create a shaded image
		pylab.close()  # clear any pre-existing plot
		im = pylab.imshow(zr, cmap=pylab.cm.RdBu, extent=[0,numcols*dx,0,numrows*dx],
						  origin='lower')
		# add contour lines with labels
		cset = pylab.contour(zr, extent=[0,numcols*dx,numrows*dx,0], hold='on',
							 origin='image')
		pylab.clabel(cset, inline=True, fmt='%1.1f', fontsize=10)
	
		# add a color bar on the side
		cb = pylab.colorbar(im)
		cb.set_label('Elevation in meters')
	
		# add a title and axis labels
		pylab.title('Simulated topography with uplift and diffusion')
		pylab.xlabel('Distance (m)')
		pylab.ylabel('Distance (m)')

		# Display the plot
		pylab.show()
		print('Run time = '+str(time.time()-start_time)+' seconds')

The last section of the ``main`` function plots the result of our calculation. We do this using pylab's ``imshow`` and ``contour`` functions to create a colored image of topography overlain by contours. To use these functions, we need our elevations to be ordered in a 2D array. We obtain a 2D array version of our ``z`` values through a call to RasterModelGrid's ``node_vector_to_raster()`` method.

Running the ``main()`` function
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

.. code-block:: python

	if __name__ == "__main__":
		main()

The last two lines of code are standard Python syntax. They will execute the ``main`` function when the code is run, but not when the code is simply imported as a module.

That's it. The 2D diffusion code is less than 100 lines long. In fact, only about 20 of these actually implement the diffusion calculation; the rest handle plotting, comments, blank lines, etc.

Checking against the analytical solution
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

To test the diffusion model against an analytical solution, we can change the setup to have closed boundaries on two opposite sides, by modifying line 39 to read:

.. code:: python

	mg.set_closed_boundaries_at_grid_edges(True, False, True, False)

This change makes the solution identical in the `y` direction, so that we can compare it with a 1D analytical solution. For a 1D steady state configuration with a constant source term (baselevel lowering) and two fixed boundaries, the elevation field is a parabola:

.. math::

	z(x') = \frac{U}{2K_d} \left( L^2 - x'^2 \right),

where :math:`L` is the half-length of the domain and :math:`x'` is a transformed :math:`x` coordinate such that :math:`x'=0` at the ridge crest. The numerical solution fits the analytical solution quite well (:ref:`Figure 5 <diffan>`).

.. _diffan:

.. figure:: images/diffusion_raster_with_analytical.png
    :scale: 50 %
    :align: center

    Figure 5: Output from the hillslope diffusion model, compared with the analytical solution (right, red curve).
