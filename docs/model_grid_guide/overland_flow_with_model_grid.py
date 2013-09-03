#! /usr/env/python
"""

2D numerical model of shallow-water flow over topography, using the
Bates et al. (2010) algorithm for storage-cell inundation modeling.

Last updated GT July 2013

"""

import landlab
import pylab, time
import numpy as np

def main():
    """
    In this simple tutorial example, the main function does all the work: 
    it sets the parameter values, creates and initializes a grid, sets up 
    the state variables, runs the main loop, and cleans up.
    """
    
    # INITIALIZE
    
    # User-defined parameter values
    numrows = 20
    numcols = 100
    dx = 50.
    n = 0.03              # roughness coefficient
    run_time = 1800       # duration of run, seconds
    h_init = 0.001        # initial thin layer of water (m)
    h_boundary = 2.5      # water depth at left side (m) 
    g = 9.8
    alpha = 0.2           # time-step factor (ND; from Bates et al., 2010)
    
    # Derived parameters
    ten_thirds = 10./3.   # pre-calculate 10/3 for speed
    elapsed_time = 0.0    # total time in simulation
    report_interval = 2.  # interval to report progress (seconds)
    next_report = time.time()+report_interval   # next time to report progress
    
    # Create and initialize a raster model grid
    mg = landlab.RasterModelGrid(numrows, numcols, dx)
    
    # Set up boundaries. We'll have the right and left sides open, the top and
    # bottom closed. The water depth on the left will be 5 m, and on the right 
    # just 1 mm.
    mg.set_inactive_boundaries(True, False, True, False)
    
    # Set up scalar values
    z = mg.create_node_dvector()   # land elevation
    h = mg.create_node_dvector() + h_init     # water depth (m)
    q = mg.create_active_link_dvector()  # unit discharge (m2/s)
    dhdt = mg.create_node_dvector()  # rate of water-depth change
    
    # Left side has deep water
    leftside = mg.left_edge_node_ids()
    h[leftside] = h_boundary
    
    # Get a list of the interior cells
    interior_cells = mg.get_active_cell_node_ids()

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
        zmax = mg.active_link_max(z)
        w = h+z   # water-surface height
        wmax = mg.active_link_max(w)
        hflow = wmax - zmax
        
        # Calculate water-surface slopes
        water_surface_slope = mg.calculate_gradients_at_active_links(w)
       
        # Calculate the unit discharges (Bates et al., eq 11)
        q = (q-g*hflow*dtmax*water_surface_slope)/ \
            (1.+g*hflow*dtmax*n*n*abs(q)/(hflow**ten_thirds))
        
        # Calculate water-flux divergence and time rate of change of water depth
        # at nodes
        dhdt = -mg.calculate_flux_divergence_at_nodes(q)
        
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
        h[interior_cells] = h[interior_cells] + dhdt[interior_cells]*dt
        
        # Update current time
        elapsed_time += dt

      
    # FINALIZE
    
    # Get a 2D array version of the elevations
    hr = mg.node_vector_to_raster(h, flip_vertically=True)
    
    # Create a shaded image
    pylab.close()  # clear any pre-existing plot
    image_extent = [0, 0.001*dx*numcols, 0, 0.001*dx*numrows] # in km
    im = pylab.imshow(hr, cmap=pylab.cm.RdBu, extent=image_extent)
    pylab.xlabel('Distance (km)', fontsize=12)
    pylab.ylabel('Distance (km)', fontsize=12)
    
    # add contour lines with labels
    cset = pylab.contour(hr, extent=image_extent)
    pylab.clabel(cset, inline=True, fmt='%1.1f', fontsize=10)
    
    # add a color bar on the side
    cb = pylab.colorbar(im)
    cb.set_label('Water depth (m)', fontsize=12)
    
    # add a title
    pylab.title('Simulated inundation')

    # Display the plot
    pylab.show()
    print('Done.')
    print('Total run time = '+str(time.time()-start_time)+' seconds.')

if __name__ == "__main__":
    main()
