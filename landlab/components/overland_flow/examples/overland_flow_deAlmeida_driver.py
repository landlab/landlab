# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 15:51:25 2015

Simple driver showing how to use the de Almeida's overland flow component in
Landlab. This driver sets up their test case from de Almeida et al., 2012,
where output from the Bates et al., 2010, de Almeida et al., 2012 and
analytical solutions are compared.

This is a componentized version of the driver example
overland_flow_with_model_grid_deAlmeida_analyticalSolution.py

.. codeauthor:: Jordan Adams <jadams15@tulane.edu>
"""

from __future__ import print_function

from landlab.components.overland_flow import OverlandFlow
from landlab.grid.structured_quad import links
from landlab.plot.imshow import imshow_grid
from landlab import RasterModelGrid

import numpy as np
import pylab
from time import time

numrows = 32              # number of rows in the raster grid
numcols = 240             # number of columns in the raster grid
dx = 25.                  # grid spacing, (m)
run_time = 9000           # duration of run, (s)
h_init = 0.001            # initial thin layer of water (m)
n = 0.01                  # roughness coefficient, (s/m^(1/3))
g = 9.81                  # gravity (m/s^2)
alpha = 0.7               # time-step factor (nondimensional; from Bates, 2010)
u = 0.4                   # constant velocity (m/s, de Almeida et al., 2012)

seven_over_three = 7. / 3.  # precalculated fractions for speed
three_over_seven = 3. / 7.

# Elapsed time starts at 1 second. This prevents errors when setting our
# boundary conditions
elapsed_time = 1.0

# Now we create our grid using the parameters set above.
rmg = RasterModelGrid(numrows, numcols, dx)

# Set our boundaries to closed to prevent water from flowing out of the  plane
rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)

# Create fields in the grid for topographic elevation, water depth, discharge.
rmg.add_zeros('topographic__elevation', at='node') # topographic elevation (m)
rmg.add_zeros('surface_water__depth', at='node') # water depth (m)
rmg.add_zeros('surface_water__discharge', at='link') # unit discharge (m2/s)

# Add our initial thin layer of water to the field of water depth.
#rmg['node']['surface_water__depth'] += h_init

# Now we'll identify our leftmost, but interior, column and the IDs of those
# nodes. One column in to prevent issues with BC.
inside_left_edge = rmg.nodes[1: -1, 1]

# Initializing our class...
of = OverlandFlow(rmg, mannings_n=n, theta=0.8, h_init=0.001)

# Now, we need to set a fixed value on the left edge, so we find the
# link neighbor arrays...
of.set_up_neighbor_arrays()

# ... and get a list of all horizonal ids, not just active ids (which is what
# the deAlmeida solution uses)
all_horizontal_ids = links.horizontal_link_ids(rmg.shape)

# from there, we are going to reset our west neighbor array...
of.west_neighbors = (links.horizontal_west_link_neighbor(rmg.shape,
                                                         all_horizontal_ids))

# and find the ids of the arrays along the west edge of the grid. We actually
# will set the discharge values here at every time step in the loop.
left_inactive_ids = links.left_edge_horizontal_ids(rmg.shape)

# Let's see how long this run takes...
starttime = time()

while elapsed_time < run_time:

    # Now we are going to set the left edge horizontal links to their
    # neighboring discharge value

    rmg['link']['surface_water__discharge'][left_inactive_ids] =   (rmg['link'][
        'surface_water__discharge'][left_inactive_ids + 1])

    # Now, we can generate overland flow.
    of.overland_flow()

    # Recalculate water depth at the boundary ...

    # water depth at left side (m)
    h_boundary = (seven_over_three * n * n * u * u * u *
                  elapsed_time) ** three_over_seven

    # And now we input that water depth along the left-most interior column,
    # in all rows that are not boundary rows.
    rmg.at_node['surface_water__depth'][inside_left_edge] = h_boundary

    # Increased elapsed time
    elapsed_time += of.dt

# End time...
endtime = time()

# Now we calculate how long the run actually took and print out the value.
totaltime = endtime - starttime
print("Total time: ", totaltime, " seconds")

# Plotting

# Our first figure will be the wave front on the horizontal plane
pylab.plt.figure(1)
imshow_grid(rmg, 'surface_water__depth', cmap="Blues", grid_units=("m", "m"))

# The second figure will compare the depth profiles of the analytical solution
# and our modeled solution.

# Here is the total plane distance
x = np.arange(0, ((numcols) * dx), dx)
x = np.insert(x, 0, 0)

# Now we will solve the analytical solution
h_analytical = (-seven_over_three * n * n * u * u * (x - (u * run_time)))

# We can only solve the analytical solution where water depth is positive...
# so we weed out negative values to avoid NaN errors.
h_analytical[np.where(h_analytical > 0)] = (h_analytical[np.where(
    h_analytical > 0)] ** three_over_seven)
h_analytical[np.where(h_analytical < 0)] = 0.0

# And we will reshape our depth array we solved for in the above loop to access
# one row for plotting. We will also remove the first (boundary) cell from this
# array, while also appending a [0] value at the # end to keep it the same size
# as the 'x' array.
h_deAlmeida = rmg['node']['surface_water__depth'].reshape(rmg.shape)
h_deAlmeida = h_deAlmeida[1][:]
h_deAlmeida = np.append(h_deAlmeida,[0])

# Now we actually plot the depth profiles...
pylab.plt.figure(2)
pylab.plt.title("Depth Profiles - Analytical and deAlmeida Solutions")
pylab.plt.plot(x,h_analytical, 'k-')
pylab.plt.plot(x, h_deAlmeida, 'b--')
pylab.legend(('Analytical Solution', 'Modeled Solution - deAlmeida'))
pylab.xlabel("Distance, (m)")
pylab.ylabel("Water Depth, (m)")
pylab.plt.show()
