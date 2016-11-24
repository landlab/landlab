====================================
Tutorial: Overland Flow across a DEM
====================================

In the next example, we create a version of the storage-cell overland-flow model that uses a digital elevation model (DEM) for the topography, and has the flow fed by rain rather than by a boundary input. In walking through the code, we'll focus only on those aspects that are new. The code is set up to run for 40 minutes (2400 seconds) of flow, which takes about 78 seconds to execute on a 2.7 Ghz Intel Core i7 processor. The complete code listing is below and available for download :download:` here <../model_grid_guide/overland_flow_with_model_grid_dem.py>`. Note that to run this code you will also need to download the DEM which is available :download:` here <../model_grid_guide/ExampleDEM/west_bijou_gully.asc>`. Output is shown in :ref:`Figure 7 <olflowdem>`.

.. _olflowdem:

.. figure:: images/overland_flow_dem.png
    :scale: 40%
    :align: center

    Figure 7: Output from a model of overland flow run on a DEM. Left: images showing
    topography, and water depth at end of run. Right: hydrograph at catchment outlet.

.. code-block:: python

	#! /usr/env/python
	"""
	2D numerical model of shallow-water flow over topography read from a DEM, using
	the Bates et al. (2010) algorithm for storage-cell inundation modeling.

	Last updated GT May 2014
	"""

	from landlab.io import read_esri_ascii
	import time
	import os
	import pylab
	import numpy as np


	def main():
		"""
		In this simple tutorial example, the main function does all the work:
		it sets the parameter values, creates and initializes a grid, sets up
		the state variables, runs the main loop, and cleans up.
		"""

		# INITIALIZE

		# User-defined parameter values
		dem_name = 'ExampleDEM/west_bijou_gully.asc'
		outlet_row = 82
		outlet_column = 38
		next_to_outlet_row = 81
		next_to_outlet_column = 38
		n = 0.06              # roughness coefficient (Manning's n)
		h_init = 0.001        # initial thin layer of water (m)
		g = 9.8               # gravitational acceleration (m/s2)
		alpha = 0.2           # time-step factor (ND; from Bates et al., 2010)
		run_time = 2400       # duration of run, seconds
		rainfall_mmhr = 100   # rainfall rate, in mm/hr
		rain_duration = 15*60 # rainfall duration, in seconds

		# Derived parameters
		rainfall_rate = (rainfall_mmhr/1000.)/3600.  # rainfall in m/s
		ten_thirds = 10./3.   # pre-calculate 10/3 for speed
		elapsed_time = 0.0    # total time in simulation
		report_interval = 5.  # interval to report progress (seconds)
		next_report = time.time()+report_interval   # next time to report progress
		DATA_FILE = os.path.join(os.path.dirname(__file__), dem_name)

		# Create and initialize a raster model grid by reading a DEM
		print('Reading data from "'+str(DATA_FILE)+'"')
		(mg, z) = read_esri_ascii(DATA_FILE)
		print('DEM has ' + str(mg.number_of_node_rows) + ' rows, ' +
			  str(mg.number_of_node_columns) + ' columns, and cell size ' + str(mg.dx)) + ' m'

		# Modify the grid DEM to set all nodata nodes to inactive boundaries
		mg.set_nodata_nodes_to_closed(z, 0) # set nodata nodes to inactive bounds

		# Set the open boundary (outlet) cell. We want to remember the ID of the
		# outlet node and the ID of the interior node adjacent to it. We'll make
		# the outlet node an open boundary.
		outlet_node = mg.grid_coords_to_node_id(outlet_row, outlet_column)
		node_next_to_outlet = mg.grid_coords_to_node_id(next_to_outlet_row, next_to_outlet_column)
		mg.set_fixed_value_boundaries(outlet_node)

		# Set up state variables
		h = mg.add_zeros('node', 'Water_depth') + h_init     # water depth (m)
		q = mg.create_active_link_array_zeros()       # unit discharge (m2/s)

		# Get a list of the core nodes
		core_nodes = mg.core_nodes

		# To track discharge at the outlet through time, we create initially empty
		# lists for time and outlet discharge.
		q_outlet = []
		t = []
		q_outlet.append(0.)
		t.append(0.)
		outlet_link = mg.active_link_connecting_node_pair(outlet_node, node_next_to_outlet)

		# Display a message
		print( 'Running ...' )
		start_time = time.time()

		# RUN

		# Main loop
		while elapsed_time < run_time:

			# Report progress
			if time.time()>=next_report:
				print('Time = '+str(elapsed_time)+' ('
					  +str(100.*elapsed_time/run_time)+'%)')
				next_report += report_interval

			# Calculate time-step size for this iteration (Bates et al., eq 14)
			dtmax = alpha*mg.dx/np.sqrt(g*np.amax(h))

			# Calculate the effective flow depth at active links. Bates et al. 2010
			# recommend using the difference between the highest water-surface
			# and the highest bed elevation between each pair of cells.
			zmax = mg.max_of_link_end_node_values(z)
			w = h+z   # water-surface height
			wmax = mg.max_of_link_end_node_values(w)
			hflow = wmax - zmax

			# Calculate water-surface slopes
			water_surface_slope = mg.calculate_gradients_at_active_links(w)

			# Calculate the unit discharges (Bates et al., eq 11)
			q = (q-g*hflow*dtmax*water_surface_slope)/ \
				(1.+g*hflow*dtmax*n*n*abs(q)/(hflow**ten_thirds))

			# Calculate water-flux divergence at nodes
			dqds = mg.calculate_flux_divergence_at_nodes(q)

			# Update rainfall rate
			if elapsed_time > rain_duration:
				rainfall_rate = 0.

			# Calculate rate of change of water depth
			dhdt = rainfall_rate-dqds

			# Second time-step limiter (experimental): make sure you don't allow
			# water-depth to go negative
			if np.amin(dhdt) < 0.:
				shallowing_locations = np.where(dhdt<0.)
				time_to_drain = -h[shallowing_locations]/dhdt[shallowing_locations]
				dtmax2 = alpha*np.amin(time_to_drain)
				dt = np.min([dtmax, dtmax2])
			else:
				dt = dtmax

			# Update the water-depth field
			h[core_nodes] = h[core_nodes] + dhdt[core_nodes]*dt
			h[outlet_node] = h[node_next_to_outlet]

			# Update current time
			elapsed_time += dt

			# Remember discharge and time
			t.append(elapsed_time)
			q_outlet.append(q[outlet_link])


		# FINALIZE

		# Set the elevations of the nodata cells to the minimum active cell
		# elevation (convenient for plotting)
		z[np.where(z<=0.)] = 9999            # temporarily change their elevs ...
		zmin = np.amin(z)                    # ... so we can find the minimum ...
		z[np.where(z==9999)] = zmin          # ... and assign them this value.

		# Get a 2D array version of the water depths and elevations
		hr = mg.node_vector_to_raster(h)
		zr = mg.node_vector_to_raster(z)

		# Clear previous plots
		pylab.figure(1)
		pylab.close()
		pylab.figure(2)
		pylab.close()

		# Plot discharge vs. time
		pylab.figure(1)
		pylab.plot(np.array(t), np.array(q_outlet)*mg.dx)
		pylab.xlabel('Time (s)')
		pylab.ylabel('Q (m3/s)')
		pylab.title('Outlet discharge')

		# Plot topography
		pylab.figure(2)
		pylab.subplot(121)
		im = pylab.imshow(zr, cmap=pylab.cm.RdBu,
						  extent=[0, mg.number_of_node_columns * mg.dx,
								  0, mg.number_of_node_rows * mg.dx])
		cb = pylab.colorbar(im)
		cb.set_label('Elevation (m)', fontsize=12)
		pylab.title('Topography')

		# Plot water depth
		pylab.subplot(122)
		im2 = pylab.imshow(hr, cmap=pylab.cm.RdBu,
						   extent=[0, mg.number_of_node_columns * mg.dx,
								   0, mg.number_of_node_rows * mg.dx])
		pylab.clim(0, 0.25)
		cb = pylab.colorbar(im2)
		cb.set_label('Water depth (m)', fontsize=12)
		pylab.title('Water depth')

		# Display the plots
		pylab.show()
		print('Done.')
		print('Total run time = '+str(time.time()-start_time)+' seconds.')


	if __name__ == "__main__":
		main()

Loading modules
>>>>>>>>>>>>>>>

.. code-block:: python

	from landlab.io import read_esri_ascii
	import time
	import os
	import pylab
	import numpy as np

In order to import the DEM, we will use Landlab's ``read_esri_ascii`` function, so we need to import this. We also need the ``time`` module for timekeeping, ``os`` for manipulating path names, ``pylab`` for plotting, and ``numpy`` for numerical operations.

User-defined variables
>>>>>>>>>>>>>>>>>>>>>>

.. code-block:: python

    # User-defined parameter values
    dem_name = 'ExampleDEM/west_bijou_gully.asc'
    outlet_row = 82
    outlet_column = 38
    next_to_outlet_row = 81
    next_to_outlet_column = 38
    n = 0.06              # roughness coefficient (Manning's n)
    h_init = 0.001        # initial thin layer of water (m)
    g = 9.8               # gravitational acceleration (m/s2)
    alpha = 0.2           # time-step factor (ND; from Bates et al., 2010)
    run_time = 2400       # duration of run, seconds
    rainfall_mmhr = 100   # rainfall rate, in mm/hr
    rain_duration = 15*60 # rainfall duration, in seconds

We will obtain topography from a 3-m resolution digital elevation model (DEM) of a small gully watershed in the West Bijou Creek drainage basin, east-central Colorado, USA. The drainage area of this catchment is about one hectare. The topography derives from airborne lidar data. The DEM is contained in an ArcInfo-format ascii file called *west_bijou_gully.asc*, located in the *ExampleDEM* folder.

In this example, we will allow flow through a single outlet cell, which we need to flag as a fixed-value boundary. We will also monitor discharge at the outlet. To accomplish these tasks, we need the row and column of the cell that will be used as the outlet and the cell next to it.

Our run will apply water as rainfall, with a rate given by ``rainfall_mmhr`` and a duration given by ``rain_duration``. In fact, in this simple model, we won't allow any infiltration, so the rainfall rate is actually a runoff rate.

Derived parameters
>>>>>>>>>>>>>>>>>>

.. code-block:: python

    # Derived parameters
    rainfall_rate = (rainfall_mmhr/1000.)/3600.  # rainfall in m/s
    ten_thirds = 10./3.   # pre-calculate 10/3 for speed
    elapsed_time = 0.0    # total time in simulation
    report_interval = 5.  # interval to report progress (seconds)
    next_report = time.time()+report_interval   # next time to report progress
    DATA_FILE = os.path.join(os.path.dirname(__file__), dem_name)

In this block of code, we convert the rainfall rate from millimeters per hour to meters per second. We also find the full path name of the input DEM by combining the pathname of the python code file (which is stored in ``__file__``) with the specified DEM file name. We take advantage of the ``dirname`` and ``join`` functions in the OS module.

Reading and initializing the DEM
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

.. code-block:: python

    # Create and initialize a raster model grid by reading a DEM
    print('Reading data from "'+str(DATA_FILE)+'"')
    (mg, z) = read_esri_ascii(DATA_FILE)
    print('DEM has ' + str(mg.number_of_node_rows) + ' rows, ' +
          str(mg.number_of_node_columns) + ' columns, and cell size ' + str(mg.dx)) + ' m'

    # Modify the grid DEM to set all nodata nodes to inactive boundaries
    mg.set_nodata_nodes_to_closed(z, 0) # set nodata nodes to inactive bounds

Landlab's IO module allows us to read an ArcInfo ascii-format DEM with a call to the ``read_esri_ascii()`` method. The function creates and returns a ``RasterModelGrid`` of the correct size and resolution, as well as a Numpy array of node elevation values. In this example, we know that the DEM contains elevations for a small watershed; nodes outside the watershed have a no-data value of zero. We don't want any flow to cross the watershed perimeter except at a single outlet cell. The call to the ModelGrid method ``set_nodata_nodes_to_closed()`` accomplishes this by identifying all nodes for which the corresponding value in ``z`` equals the specified no-data value of zero.

Setting up the watershed outlet
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

.. code-block:: python

    # Set the open boundary (outlet) cell. We want to remember the ID of the
    # outlet node and the ID of the interior node adjacent to it. We'll make
    # the outlet node an open boundary.
    outlet_node = mg.grid_coords_to_node_id(outlet_row, outlet_column)
    node_next_to_outlet = mg.grid_coords_to_node_id(next_to_outlet_row,
                                                    next_to_outlet_column)
    mg.set_fixed_value_boundaries(outlet_node)

We will handle the outlet by keeping the water-surface slope the same as the bed-surface slope along the link that leads to the outlet boundary node. To accomplish this, the first thing we need to do is find the ID of the outlet node and the core node adjacent to it. We already know what the row and column numbers of these nodes are; to obtain the corresponding node ID, we use ModelGrid's ``grid_coords_to_node_id`` method. We then convert the outlet node to a fixed-value (i.e., open) boundary with the ``set_fixed_value_boundaries`` method. (Note that in doing this, we've converted what was a core node into a fixed boundary; had we converted a no-data node, we would end up with a waterfall at the outlet because the no-data nodes all have zero elevation, while the interior nodes all have elevations above 1600 m).

Preparing to track discharge at the outlet
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

.. code-block:: python

    # To track discharge at the outlet through time, we create initially empty
    # lists for time and outlet discharge.
    q_outlet = []
    t = []
    q_outlet.append(0.)
    t.append(0.)
    outlet_link = mg.active_link_connecting_node_pair(outlet_node,
                                                          node_next_to_outlet)

For this model, it would be nice to track discharge through time at the watershed outlet. To do this, we create two new lists: one for the time corresponding to each iteration, and one for the outlet discharge. Using lists will be slightly slower than using pre-defined Numpy arrays, but avoids forcing us to guess how many iterations there will be (recall that time-step size depends on the flow conditions in any given iteration). We append zeros to each list to represent the starting condition. To find out which active link represents the watershed outlet, we use ModelGrid's ``active_link_connecting_node_pair()`` method. This method takes a pair of node IDs as arguments. If the nodes are connected by an active link, it returns the ID of that active link; otherwise, it returns ``ModelGrid.BAD_INDEX_VALUE``.

Main loop
>>>>>>>>>

.. code-block:: python

        # Update rainfall rate
        if elapsed_time > rain_duration:
            rainfall_rate = 0.

        # Calculate rate of change of water depth
        dhdt = rainfall_rate-dqds

Most of the main loop is identical to what we saw in Example 2, and here we will only list the parts that are new or different. One difference is that we now have a source term that represents rainfall and runoff. The code listed above sets the rainfall rate to zero when the elapsed time is greater than the rainfall duration. It also adds ``rainfall_rate`` as a source term when computing :math:`dh/dt`.

.. code-block:: python

        # Update the water-depth field
        h[core_nodes] = h[core_nodes] + dhdt[core_nodes]*dt
        h[outlet_node] = h[node_next_to_outlet]

After updating water depth values for the core nodes, we also need to update the water depth at the outlet boundary so that it matches the depth at the adjacent node.

.. code-block:: python

        # Remember discharge and time
        t.append(elapsed_time)
        q_outlet.append(q[outlet_link])

The last few lines in the main loop keep track of discharge at the outlet by appending the current time and discharge to their respective lists.

Plotting the result
>>>>>>>>>>>>>>>>>>>

The plotting section is similar to what we saw in the previous two examples. One difference is that we now use two figures: one for the topography and water depth, and one for outlet discharge over time. We also use Pylab's sub-plot capability to place images of topography and water depth side by side.
