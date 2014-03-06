#! /usr/env/python
"""

2D numerical model of diffusion, implemented using ModelGrid.
Provides example of an unstructured grid.

Last updated GT August 2013

"""

from landlab import HexModelGrid
import pylab
import time

def main():
    """
    In this simple tutorial example, the main function does all the work: 
    it sets the parameter values, creates and initializes a grid, sets up 
    the state variables, runs the main loop, and cleans up.
    """

    # INITIALIZE
    
    # User-defined parameter values
    numrows = 7          # number of rows in the grid
    basenumcols = 6          # number of columns in the grid
    dx = 10.0             # grid cell spacing
    kd = 0.01             # diffusivity coefficient, in m2/yr
    uplift_rate = 0.001   # baselevel/uplift rate, in m/yr
    num_time_steps = 1000 # number of time steps in run
    output_name = 'hex_grid_example'  # base name for output files
    output_interval = 100 # interval between outputs, in time steps
    
    # Derived parameters
    dt = 0.1*dx**2 / kd    # time-step size set by CFL condition
    next_output = output_interval  # time of next output to file
    
    # Create and initialize a raster model grid
    mg = HexModelGrid(numrows, basenumcols, dx)
    
    # Set up scalar values
    z = mg.create_node_array_zeros()            # node elevations
    dzdt = mg.create_node_array_zeros()  # node rate of elevation change
    
    # Get a list of the interior cells
    interior_cells = mg.get_active_cell_node_ids()

    # Display a message
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
        z[interior_cells] = z[interior_cells] + dzdt[interior_cells] * dt
        
   
    # FINALIZE
    
    import numpy

    # Test solution for a 1-cell hex grid
    if mg.number_of_nodes==7:  # single cell with 6 boundaries
        perimeter = dx*6*(numpy.sqrt(3.)/3.)
        flux = kd*(numpy.amax(z)/dx)
        total_outflux = perimeter*flux
        total_influx = mg.active_cell_areas*uplift_rate
        print 'total influx=',total_influx,'total outflux=',total_outflux
    print('Run time = '+str(time.time()-start_time))
    
    # Plot the points, colored by elevation
    pylab.figure()
    maxelev = numpy.amax(z)
    for i in range(mg.number_of_nodes):
        mycolor = str(z[i]/maxelev)
        pylab.plot(mg.node_x[i], mg.node_y[i], 'o', color=mycolor, ms=50)
    pylab.show()
    
    mg.display_grid()
    

if __name__ == "__main__":
    main()