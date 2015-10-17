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
