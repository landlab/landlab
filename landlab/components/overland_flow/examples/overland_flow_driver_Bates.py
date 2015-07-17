# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 15:51:25 2015

Simple driver showing how to use the Bates' overland flow component in Landlab.
This driver sets up their test case from de Almeida et al., 2012, where output
from the Bates et al., 2010, de Almeida et al., 2012 and analytical solutions 
are compared. 

This is a componentized version of the driver example
overland_flow_with_model_grid_Bates_analyticalSolution.py

@author: Jordan Adams
"""

from __future__ import print_function

from landlab.components.overland_flow.generate_overland_flow_Bates import OverlandFlow
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
g = 9.8                   # gravity (m/s^2) 
alpha = 0.7               # time-step factor (nondimensional; from Bates et al., 2010)
u = 0.4                   # constant velocity (m/s, de Almeida et al., 2012)

seven_over_three = 7./3.  # precalculated fractions for speed
three_over_seven = 3./7.
ten_thirds = 10./3.

# Elapsed time starts at 1 second. This prevents errors when setting our
# boundary conditions
elapsed_time = 1.0

# Now we create our grid using the parameters set above.
rmg = RasterModelGrid(numrows, numcols, dx)

# Set our boundaries to closed to prevent water from flowing out of the study plane
rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)

# Create fields in the grid for topographic elevation, water depth and discharge.
rmg.create_node_array_zeros('topographic__elevation')        # topographic elevation (m)
rmg.create_node_array_zeros('water_depth')                   # water depth (m)
rmg.add_zeros('active_link', 'water_discharge')              # unit discharge (m2/s)

# Add our initial thin layer of water to the field of water depth.
rmg['node']['water_depth'] += h_init 

# Now we'll identify our leftmost, but interior, column and the IDs of those nodes.
leftside = rmg.left_edge_node_ids()
leftside = leftside+1                     # One column in to prevent issues with BC

# Initializing our class...
of = OverlandFlow(rmg)

# Let's see how long this run takes...
starttime = time()
while elapsed_time < run_time:
    # First, we calculate our time step.
    dt = of.gear_time_step(rmg)

    # Now, we generate overland flow.
    of.overland_flow(rmg, dt)
    
    # Recalculate water depth at the boundary ...
    h_boundary = ((seven_over_three)*n*n*u*u*u*elapsed_time)**(three_over_seven)      # water depth at left side (m) 
    
    # And now we input that water depth along the left-most interior column, in all rows that are not boundary rows.
    rmg['node']['water_depth'][(leftside)[1:len(leftside)-1]] = h_boundary
    
    # Print time
    #print(elapsed_time)

    # Increased elapsed time 
    elapsed_time += dt
    
# End time...
endtime = time()

# Now we calculate how long the run actually took and print out the total value.
totaltime = endtime-starttime
print("Total time: ", totaltime, " seconds")

# Plotting

# Our first figure will be the wave front on the horizontal plane
pylab.plt.figure(1)
imshow_grid(rmg, 'water_depth', cmap="Purples", grid_units=("m", "m"))

# The second figure will compare the depth profiles of the analytical solution and our modeled solution
# Here is the total plane distance
x = np.arange(0, ((numcols)*dx), dx)

# Now we will solve the analytical solution
h_analytical = (-seven_over_three*n*n*u*u*(x-(u*run_time)))

# We can only solve the analytical solution where water depth is positive... so we weed out negative
# values to avoid NaN errors. 
h_analytical[np.where(h_analytical>0)] = h_analytical[np.where(h_analytical>0)]**three_over_seven
h_analytical[np.where(h_analytical<0)] = 0.0

# And we will reshape our depth array we solved for in the above loop to access one row for plotting.
# We will also remove the first (boundary) cell from this array, while also appending a [0] value at the
# end to keep it the same size as the 'x' array.
h_Bates = rmg['node']['water_depth'].reshape(rmg.shape)
h_Bates = h_Bates[1][1:]
h_Bates = np.append(h_Bates,[0])

# Now we actually plot the depth profiles...
pylab.plt.figure(2)
pylab.plt.title("Depth Profiles - Analytical and Bates Solutions")
pylab.plt.plot(x,h_analytical, 'k--')
pylab.plt.plot(x, h_Bates, 'm-')
pylab.legend(('Analytical Solution', 'Modeled Solution - Bates'))
pylab.xlabel("Distance, (m)")
pylab.ylabel("Water Depth, (m)")
pylab.plt.show()