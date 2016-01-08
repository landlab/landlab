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
numrows = 32           # number of rows in the raster grid
numcols = 240          # number of columns in the raster grid
dx = 25.               # grid spacing, (m)
n = 0.01               # roughness coefficient, (s/m^(1/3))
run_time = 9000        # duration of run, (s)
h_init = 0.001         # initial thin layer of water (m)
g = 9.8                # gravity (m/s^2)
alpha = 0.7            # time-step factor (nondimensional, Bates et al., 2010)
u = 0.4                # constant velocity (m/s, de Almeida et al., 2012)
theta = 0.8            # weighting factor (nondimensional; de Almeida 2012)
seven_over_three = 7. / 3.  # precalculated fraction for speed
three_over_seven = 3. / 7.

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

# water depth (m) on links
h_links = mg.add_zeros('link', 'water_depth') + h_init

# rate water-depth change
dhdt = mg.add_zeros('node', 'water_depth_time_derivative')

# Left side has deep water
leftside = mg.left_edge_node_ids()
leftside = leftside + 1               # One column in to prevent issues with BC

# Get a list of the core nodes
core_nodes = mg.core_nodes
active_links = mg.active_links

# Display a message
print( 'Running ...' )
start_time = time.time()


# These functions find link neighbors for horizontal and vertical active links.
# First, we find all active links in the raster grid.
active_ids = links.active_link_ids(mg.shape, mg.status_at_node)

# Then we find all horizontal link ids...
horizontal_ids = links.horizontal_link_ids(mg.shape)

# Get ids of left-most link ids for boundary issues...
left_inactive_ids = horizontal_ids[:,0]

# Then we flatten our array so we can use it elsewhere.
horizontal_ids = horizontal_ids.flatten()

# ... and narrow them down to the active horizontal link ids.
horizontal_active_link_ids = (links.horizontal_active_link_ids(mg.shape,
                                                               active_ids))

# Here we actually identify, for each link, the id of its W and E neighbor.
# For the de Almeida solution, horizontal neighbors are west and east. Only
# active link ids are given, any inactive link id is replaced with a '-1'.
west_neighbors = links.horizontal_west_link_neighbor(mg.shape, horizontal_ids)
east_neighbors = (links.horizontal_east_link_neighbor(mg.shape,
                                                horizontal_active_link_ids))

# Now we do the same with all vertical link ids. First, we get ALL vertical ids.
vertical_ids = links.vertical_link_ids(mg.shape).flatten()

# And then we narrow them down to just the active vertical link ids.
vertical_active_link_ids = links.vertical_active_link_ids(mg.shape, active_ids)

# For the de Almeida solution, we only need N and S neighbors for vertical
# active links, so for each link, the N and S neighbor ids are given by these
# two function calls. Any inactive link id is replaced with an index of '-1'.
north_neighbors = (links.vertical_north_link_neighbor(mg.shape,
                                                    vertical_active_link_ids))
south_neighbors = (links.vertical_south_link_neighbor(mg.shape,
                                                    vertical_active_link_ids))


# First we calculate our updated boundary water depth
h_boundary = (((seven_over_three) * n * n * u * u * u * elapsed_time) **
    (three_over_seven))      # water depth at left side (m)

# And now we add it to the second column, in all rows that are not boundary
# rows.
h[(leftside)[1:len(leftside)-1]] = h_boundary

# Main loop
while elapsed_time <= run_time:

    # Calculate time-step size for this iteration (Bates et al., eq 14)
    dtmax = alpha * mg.dx / np.sqrt(g * np.amax(h))
    # First we calculate our updated boundary water depth
    h_boundary = (((seven_over_three) * n * n * u * u * u * elapsed_time) **
        (three_over_seven))      # water depth at left side (m)

    # And now we add it to the second column, in all rows that are not boundary
    # rows.
    h[(leftside)[1:len(leftside)-1]] = h_boundary
    # Calculate the effective flow depth at active links. Bates et al. 2010
    # and de Almeida et al. 2012 both recommend using the the difference
    # between the highest water-surface and the highest bed elevation
    # between each pair of nodes.
    zmax = mg.max_of_link_end_node_values(z)
    w = h + z   # water-surface height
    wmax = mg.max_of_link_end_node_values(w)
    hflow = wmax - zmax

    # Putting our active link flow depths back into an array of len(links)
    h_links[active_links] = hflow

    # Getting active link flow depths on the horizontal links...
    h_horizontal = h_links[horizontal_ids]
    # ...and the vertical links.
    h_vertical = h_links[vertical_ids]

    # Calculate water-surface slopes
    water_surface_slope = mg.calculate_gradients_at_active_links(w)
    # Put the active link slopes back into an array of len(links)
    slope[active_links] = water_surface_slope

    # Add in a neighbor boundary condition for the left-most links
    q[left_inactive_ids] = q[left_inactive_ids + 1]
    # Now a little trick to help with our '-1' indices...
    # We append a '0' to the end of our discharge array so that anywhere there
    # is an index of '-1', it will give a value of '0'.
    q = np.append(q, [0])
    # Calculate discharge on our horizontal links using Eq. 41 from
    # de Almeida et al., 2012.
    q[horizontal_ids] = ((((theta * q[horizontal_ids]) + (((1 - theta) / 2) *
        (q[west_neighbors] + q[east_neighbors]))) - (g *
        h_links[horizontal_ids] * (dtmax) * slope[horizontal_ids])) / (1 + g *
        dtmax * n * n * abs(q[horizontal_ids]) / (h_links[horizontal_ids] **
        seven_over_three)))

    # Calculate discharge on our vertical links using Eq. 41 from
    # de Almeida et al., 2012.
    q[vertical_ids] = ((((theta * q[vertical_ids]) + ((1 - theta) / 2) *
        (q[north_neighbors] + q[south_neighbors])) - (g * h_links[vertical_ids]
        * (dtmax) * slope[vertical_ids])) / (1 + g * dtmax * n * n *
        abs(q[vertical_ids]) / (h_links[vertical_ids] ** seven_over_three)))

    # Now we delete the extra value of '0' from the end of our discharge array,
    # removing our little trick that handled '-1' indices.
    q = np.delete(q, len(q) - 1)

    # To create an entire discharge array with len(links), we combine the
    # arrays with vertical and horizontal links.

    # Calculate water-flux divergence and time rate of change of water depth
    # at nodes
    dhdt = -mg.calculate_flux_divergence_at_nodes(q[active_links])

    # Update the water-depth field
    h[core_nodes] = h[core_nodes] + dhdt[core_nodes]*dtmax

    # Now we update our boundary condition, adding water to the second column
    # to prevent issues with boundary conditions in the first column.

    # First we calculate our updated boundary water depth
    h_boundary = (((seven_over_three) * n * n * u * u * u * elapsed_time) **
        (three_over_seven))      # water depth at left side (m)

    # And now we add it to the second column, in all rows that are not boundary
    # rows.
    h[(leftside)[1:len(leftside)-1]] = h_boundary


    # Print and update current time
    #print elapsed_time
    elapsed_time += dtmax


# FINALIZE

# Plot the water depth values across the grid at the end of the run
plt.figure(1)
imshow_grid(mg, h, show=True)

# Compare the wave front to the analytical solution
# First, we will get an array of the distance
x = np.arange(0, ((numcols) * dx), dx)

# Now we will solve the analytical solution
h_analytical = (-seven_over_three * n * n * u * u * (x - (u * elapsed_time)))

# We can only solve the analytical solution where water depth is positive...
# so we weed out negative values to avoid NaN errors.
h_analytical[np.where(h_analytical > 0)] = (h_analytical[
    np.where(h_analytical > 0)] ** three_over_seven)
h_analytical[np.where(h_analytical < 0)] = 0.0

# And we will reshape our depth array we solved for in the above loop to
# access one row for plotting. We will also remove the first (boundary) cell
# from this array, while also appending a [0] value at the end to keep it the
# same size as the 'x' array.
h_deAlmeida = h.reshape(mg.shape)
h_deAlmeida = h_deAlmeida[1][1:]
h_deAlmeida = np.append(h_deAlmeida,[0])

plt.figure(2)
plt.plot(x,h_analytical)
plt.plot(x, h_deAlmeida)
plt.show()
