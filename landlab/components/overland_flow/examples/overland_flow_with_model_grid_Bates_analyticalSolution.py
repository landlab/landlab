#! /usr/env/python
"""

2D numerical model of shallow-water flow over topography, using the
de Almeida et al. (2012) algorithm for storage-cell inundation modeling.

Last updated Jordan Adams July 2015

"""

from landlab import RasterModelGrid
import time
import matplotlib.pyplot as plt
import numpy as np
from landlab.plot import imshow_grid
from landlab.grid.structured_quad import links

"""
In this simple tutorial example, the main function does all the work: 
it sets the parameter values, creates and initializes a grid, sets up 
the state variables, runs the main loop, and cleans up.
"""

# INITIALIZE

# User-defined parameter values
numrows = 32              # number of rows in the raster grid
numcols = 240             # number of columns in the raster grid
dx = 25.                  # grid spacing, (m)
n = 0.01                  # roughness coefficient, (s/m^(1/3))
run_time = 9005        # duration of run, (s)
h_init = 0.001            # initial thin layer of water (m)
g = 9.8                   # gravity (m/s^2) 
alpha = 0.7               # time-step factor (nondimensional; from Bates et al., 2010)
u = 0.4                   # constant velocity (m/s, de Almeida et al., 2012)
theta = 0.8               # weighting factor (nondimensional; de Almeida 2012)
seven_over_three = 7./3.  # precalculated fractions for speed
three_over_seven = 3./7.
ten_thirds = 10./3.

elapsed_time = 1.0   # start time for simulation

# Create and initialize a raster model grid
mg = RasterModelGrid(numrows, numcols, dx)
    
# Set up boundaries. We'll have the right and left sides open, the top and
# bottom closed. The water depth on the left will be 5 m, and on the right 
# just 1 mm.
mg.set_closed_boundaries_at_grid_edges(True, True, True, True)

# Set up scalar values
z = mg.add_zeros('node', 'topographic_elevation')    # land elevation
h = mg.add_zeros('node', 'water_depth') + h_init     # water depth (m)
q = mg.add_zeros('link', 'water_discharge')          # unit discharge (m2/s)
slope = mg.add_zeros('link', 'topographic_slope')    # dimensionless slope
h_links = mg.add_zeros('link', 'water_depth')+h_init # water depth (m) on links 
dhdt = mg.add_zeros('node', 'water_depth_time_derivative')  # rate of water-depth change

# Left side has deep water
leftside = mg.left_edge_node_ids()
leftside = leftside+1                     # One column in to prevent issues with BC

# Get a list of the core nodes
core_nodes = mg.core_nodes
active_links = mg.active_links

# Display a message
print( 'Running ...' )
start_time = time.time()


# Main loop
while elapsed_time <= run_time:

    # Calculate time-step size for this iteration (Bates et al., eq 14)
    dtmax = alpha*mg.dx/np.sqrt(g*np.amax(h))
    
    # First we calculate our updated boundary water depth
    h_boundary = ((seven_over_three)*n*n*u*u*u*elapsed_time)**(three_over_seven)      # water depth at left side (m) 
    
    # And now we add it to the second column, in all rows that are not boundary rows.
    h[(leftside)[1:len(leftside)-1]] = h_boundary
    
    # Calculate the effective flow depth at active links. Bates et al. 2010
    # and de Almeida et al. 2012 both recommend using the the difference 
    # between the highest water-surface and the highest bed elevation
    # between each pair of nodes.
    zmax = mg.max_of_link_end_node_values(z)
    w = h+z   # water-surface height
    wmax = mg.max_of_link_end_node_values(w)
    hflow = wmax - zmax

    # Calculate water-surface slopes
    water_surface_slope = mg.calculate_gradients_at_active_links(w)

    # Calculate discharge using Eq. 11 from Bates et al., 2010
    q[active_links] = (q[active_links]-g*hflow*dtmax*water_surface_slope)/(1.+g*hflow*dtmax*n*n*abs(q[active_links])/(hflow**ten_thirds))
    
    # Calculate water-flux divergence and time rate of change of water depth at nodes
    dhdt = -mg.calculate_flux_divergence_at_nodes(q[active_links])

    # Update the water-depth field
    h[core_nodes] = h[core_nodes] + dhdt[core_nodes]*dtmax
    
    # Now we update our boundary condition, adding water to the second column
    # to prevent issues with boundary conditions in the first column.
    # First we calculate our updated boundary water depth
    h_boundary = ((seven_over_three)*n*n*u*u*u*elapsed_time)**(three_over_seven)      # water depth at left side (m) 
    
    # And now we add it to the second column, in all rows that are not boundary rows.
    h[(leftside)[1:len(leftside)-1]] = h_boundary

    # Print and update current time
    print elapsed_time
    elapsed_time += dtmax

    
# FINALIZE

# Plot the water depth values across the grid at the end of the run
plt.figure(1)
imshow_grid(mg, h, show=True)

# Compare the wave front to the analytical solution
# First, we will get an array of the distance
x = np.arange(0, ((numcols)*dx), dx)

# Now we will solve the analytical solution
h_analytical = (-seven_over_three*n*n*u*u*(x-(u*9000)))

# We can only solve the analytical solution where water depth is positive... so we weed out negative
# values to avoid NaN errors. 
h_analytical[np.where(h_analytical>0)] = h_analytical[np.where(h_analytical>0)]**three_over_seven
h_analytical[np.where(h_analytical<0)] = 0.0

# And we will reshape our depth array we solved for in the above loop to access one row for plotting.
# We will also remove the first (boundary) cell from this array, while also appending a [0] value at the
# end to keep it the same size as the 'x' array.
h_deAlmeida = h.reshape(mg.shape)
h_deAlmeida = h_deAlmeida[1][1:]
h_deAlmeida = np.append(h_deAlmeida,[0])

plt.figure(2)
plt.plot(x,h_analytical)
plt.plot(x, h_deAlmeida)
plt.show()