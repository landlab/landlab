#! /usr/env/python
"""

2D numerical model of diffusion, implemented using ModelGrid.
Provides example of a radial grid.

Last updated GT August 2013

"""

from landlab import RadialModelGrid
import pylab

def main():
    """
    In this simple tutorial example, the main function does all the work: 
    it sets the parameter values, creates and initializes a grid, sets up 
    the state variables, runs the main loop, and cleans up.
    """
    
    # INITIALIZE
    
    # User-defined parameter values
    num_shells=10         # number of radial "shells" in the grid
    #numcols = 30         # not needed for a radial model grid
    dr = 10.0             # grid cell spacing
    kd = 0.01             # diffusivity coefficient, in m2/yr
    uplift_rate = 0.0001  # baselevel/uplift rate, in m/yr
    num_time_steps = 1000 # number of time steps in run
    
    # Derived parameters
    dt = 0.1*dr**2 / kd    # time-step size set by CFL condition
    
    # Create and initialize a radial model grid
    mg = RadialModelGrid()
    mg.initialize(num_shells, dr)
    
    # Set up scalar values
    z = mg.create_node_dvector()            # node elevations
    dzdt = mg.create_node_dvector()  # node rate of elevation change
    
    # Get a list of the interior cells
    interior_cells = mg.get_active_cell_node_ids()

    # Display a message
    print( 'Running diffusion_with_radial_model_grid.py' )
    print( 'Time-step size has been set to ' + str( dt ) + ' years.' )

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
    
    # Plot the points, colored by elevation
    import numpy
    maxelev = numpy.amax(z)
    for i in range(mg.num_nodes):
        mycolor = str(z[i]/maxelev)
        pylab.plot(mg.node_x[i], mg.node_y[i], 'o', color=mycolor, ms=10)
    
    mg.display_grid()
    
    # Plot the points from the side, with analytical solution
    pylab.figure(3)
    L = num_shells*dr
    xa = numpy.arange(-L, L+dr, dr)
    z_analytical = (uplift_rate/(4*kd))*(L*L-xa*xa)
    pylab.plot(mg.node_x, z, 'o')
    pylab.plot(xa, z_analytical, 'r-')
    pylab.xlabel('Distance from center (m)')
    pylab.ylabel('Height (m)')
    
    pylab.show()
    

if __name__ == "__main__":
    main()