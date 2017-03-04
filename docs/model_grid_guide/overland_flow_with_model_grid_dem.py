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
from landlab.plot import imshow_grid


def main():
    """
    In this simple tutorial example, the main function does all the work:
    it sets the parameter values, creates and initializes a grid, sets up
    the state variables, runs the main loop, and cleans up.
    """

    # INITIALIZE

    # User-defined parameter values
    dem_name = 'ExampleDEM/west_bijou_gully.asc'
    outlet_row = 6
    outlet_column = 38
    next_to_outlet_row = 7
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
    node_next_to_outlet = mg.grid_coords_to_node_id(next_to_outlet_row,
                                                    next_to_outlet_column)
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
        q_outlet.append(abs(q[outlet_link]))


    # FINALIZE

    # Set the elevations of the nodata cells to the minimum active cell
    # elevation (convenient for plotting)
    z[np.where(z<=0.)] = 9999            # temporarily change their elevs ...
    zmin = np.amin(z)                    # ... so we can find the minimum ...
    z[np.where(z==9999)] = zmin          # ... and assign them this value.

    # Get a 2D array version of the water depths and elevations

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
    imshow_grid(mg, z, allow_colorbar=False)
    pylab.xlabel(None)
    pylab.ylabel(None)
    im = pylab.set_cmap('RdBu')
    cb = pylab.colorbar(im)
    cb.set_label('Elevation (m)', fontsize=12)
    pylab.title('Topography')

    # Plot water depth
    pylab.subplot(122)
    imshow_grid(mg, h, allow_colorbar=False)
    im2 = pylab.set_cmap('RdBu')
    pylab.clim(0, 0.25)
    cb = pylab.colorbar(im2)
    cb.set_label('Water depth (m)', fontsize=12)
    pylab.title('Water depth')
    #
    ## Display the plots
    pylab.show()
    print('Done.')
    print('Total run time = '+str(time.time()-start_time)+' seconds.')


if __name__ == "__main__":
    main()
