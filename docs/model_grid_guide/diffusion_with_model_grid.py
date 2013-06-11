#! /usr/env/python
"""

2D numerical model of diffusion, implemented using ModelGrid.
Provides a simple tutorial example of ModelGrid functionality.

Last updated GT June 2013

"""

from landlab import model_grid
import pylab

def main():
    """
    In this simple tutorial example, the main function does all the work: 
    it sets the parameter values, creates and initializes a grid, 
    sets up the state variables, runs the main loop, and cleans up.
    """
    
    # INITIALIZE
    
    # User-defined parameter values
    numrows = 10          # number of rows in the grid
    numcols = 15          # number of columns in the grid
    dx = 10.0             # grid cell spacing
    kd = 0.01             # diffusivity coefficient, in m2/yr
    uplift_rate = 0.001   # baselevel/uplift rate, in m/yr
    num_time_steps = 1000 # number of time steps in run
    
    # Derived parameters
    dt = 0.1*dx**2 / kd    # time-step size set by CFL condition
    
    # Create and initialize a raster model grid
    mg = model_grid.RasterModelGrid()
    mg.initialize(numrows, numcols, dx)
    
    # Set the boundary conditions
    mg.set_noflux_boundaries(False, True, False, True)

    # Set up scalar values
    z = mg.create_cell_dvector()     # node/cell elevations
    dzdt = mg.create_cell_dvector()  # node/cell rate of elevation change
    
    # Get a list of the interior cells
    interior_cells = mg.get_interior_cells()

    # Display a message
    print( 'Running diffusion_with_model_grid.py' )
    print( 'Time-step size has been set to' + str( dt ) + 'years.' )
    

    # RUN
    
    # Main loop
    for i in range(0, num_time_steps):
        
        # Calculate the gradients and sediment fluxes
        g = mg.calculate_face_gradients(z)
        qs = -kd*g
        
        # Calculate the net deposition/erosion rate in each cell
        dqsds = mg.calculate_flux_divergences(qs)
        
        # Calculate the total rate of elevation change in the interior cells
        for c in interior_cells:
            dzdt[c] = uplift_rate - dqsds[c]
            
        # Update the elevations
        z = z + dzdt * dt
        
        # Update the boundaries to maintain no-flux condition where needed
        z = mg.update_noflux_boundaries(z)
   
   
    # FINALIZE
    
    # Get a 2D array version of the elevations
    zr = mg.cell_vector_to_raster(z)
    
    # Create a shaded image
    im = pylab.imshow(zr, cmap=pylab.cm.RdBu)  # display a colored image
    
    # add contour lines with labels
    cset = pylab.contour(zr)
    pylab.clabel(cset, inline=True, fmt='%1.1f', fontsize=10)
    
    # add a color bar on the side
    pylab.colorbar(im)
    
    # add a title
    pylab.title('Simulated topography with uplift and diffusion')

    # Display the plot
    pylab.show()
   

if __name__ == "__main__":
    main()