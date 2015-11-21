from __future__ import print_function

from landlab.components.nonlinear_diffusion.Perron_nl_diffuse import PerronNLDiffuse
# ...the two different diffusion formulations
from landlab.components.diffusion.diffusion import LinearDiffuser
from landlab import ModelParameterDictionary  # handles input from the input file

from landlab import RasterModelGrid  # the grid object
from landlab.plot.imshow import imshow_node_grid
import numpy as np
import pylab

# get the needed properties to build the grid:
input_file = './diffusion_test_params.txt'
# initialize an object that will supply the parameters:
inputs = ModelParameterDictionary(input_file)
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
leftmost_elev = inputs.read_float('leftmost_elevation')
initial_slope = inputs.read_float('initial_slope')
uplift_rate = inputs.read_float('uplift_rate')
runtime = inputs.read_float('total_time')
dt = inputs.read_float('dt')
# "//" means "divide and truncate" ("%" means "division remainder")
nt = int(runtime // dt)
uplift_per_step = uplift_rate * dt

# Instantiate the grid object
# We know which parameters are needed for input by inspecting the function where it lives, in landlab.grid.raster
# We could also look at the documentation for landlab found online
# (http://the-landlab.readthedocs.org)
mg = RasterModelGrid(nrows, ncols, dx)

# create the elevation field in the grid:
# create the field
mg.create_node_array_zeros('topographic__elevation')
# in our case, slope is zero, so the leftmost_elev is the mean elev
z = mg.create_node_array_zeros() + leftmost_elev
# put these values plus roughness into that field
mg['node']['topographic__elevation'] = z + np.random.rand(len(z)) / 100000.

# set up its boundary conditions (bottom, left, top, right)
# The mechanisms for this are all automated within the grid object
mg.set_fixed_value_boundaries_at_grid_edges(True, True, True, True)

# Display a message
print('Running ...')

# instantiate the components:
diffuse = PerronNLDiffuse(mg, input_file)
lin_diffuse = LinearDiffuser(grid=mg, input_stream=input_file)
# lin_diffuse.initialize(input_file)

# Perform the loops.
# First, we do the nonlinear diffusion:

# We're going to perform a block uplift of the interior of the grid, but
# leave the boundary nodes at their original elevations.
# access this function of the grid, and store the output with a local name
uplifted_nodes = mg.get_core_nodes()
#(Note: Node numbering runs across from the bottom left of the grid.)

# nt is the number of timesteps we calculated above, i.e., loop nt times.
# We never actually use i within the loop, but we could do.
for i in xrange(nt):
    #("xrange" is a clever memory-saving way of producing consecutive integers to govern a loop)
    # this colon-then-tab-in arrangement is what Python uses to delineate connected blocks of text, instead of brackets or parentheses
    # This line performs the actual functionality of the component:
    # mg = lin_diffuse.diffuse(dt) #linear diffusion
    mg = diffuse.diffuse(mg, i * dt)  # nonlinear diffusion
    #...swap around which line is commented out to switch between formulations of diffusion

    # now plot a N-S cross section from this stage in the run onto figure 1.
    # The sections will all be superimposed, as show() hasn't yet been called
    pylab.figure(1)
    # turn the 1-D array of elevation values into a spatially accurate 2-D
    # gridded format, for plotting
    elev_r = mg.node_vector_to_raster(mg['node']['topographic__elevation'])
    # square brackets denote a subset of nodes to use.
    im = pylab.plot(mg.dx * np.arange(nrows), elev_r[:, int(ncols // 2)])
    #...this kind of data extraction from a larger data structure ("slicing", or "fancy indexing") is extremely useful and powerful, and is one of the appeals of Python
    #...this is plot(x, y).
    # x is the distance north up the grid.
    # y is the elevation along all the rows, but only the 50th column (more
    # slicing!), i.e., halfway along the grid and N-S

    print('Completed loop ', i)

print('Completed the simulation. Plotting...')

# Finalize and plot:
# put a title on figure 1
pylab.figure(1)
pylab.title('N-S cross_section, nonlinear diffusion')
pylab.xlabel('Distance')
pylab.ylabel('Elevation')

# figure 2 is the map of the final elevations
pylab.figure(2)
# display a colored image
imshow_node_grid(mg, 'topographic__elevation')

pylab.figure(3)
# turn the 1-D array of elevation values into a spatially accurate 2-D
# gridded format, for plotting
elev_r = mg.node_vector_to_raster(mg['node']['topographic__elevation'])
im = pylab.plot(mg.dx * np.arange(nrows), elev_r[:, int(ncols // 2)])

print('Done.')

# now do the linear diffusion:

# Reset the elevation field in the grid:
mg['node']['topographic__elevation'] = z + np.random.rand(len(z)) / 100000.

# Display a message
print('Running ...')

for i in xrange(nt):
    # This line performs the actual functionality of the component:
    #***NB: the nonlinear diffuser contains an "automatic" element of uplift. If you instead use the linear diffuser, you need to add the uplift manually...
    mg['node']['topographic__elevation'][uplifted_nodes] += uplift_per_step
    mg = lin_diffuse.diffuse(dt)  # linear diffusion

    pylab.figure(4)
    elev_r = mg.node_vector_to_raster(mg['node']['topographic__elevation'])
    im = pylab.plot(mg.dx * np.arange(nrows), elev_r[:, int(ncols // 2)])

    print('Completed loop ', i)

print('Completed the simulation. Plotting...')

# Finalize and plot:
# put a title on figure 4
pylab.figure(4)
pylab.title('N-S cross_section, linear diffusion')
pylab.xlabel('Distance')
pylab.ylabel('Elevation')

# figure 5 is the map of the final elevations
pylab.figure(5)
imshow_node_grid(mg, 'topographic__elevation')

# superpose this final form onto figure 3:
pylab.figure(3)
# turn the 1-D array of elevation values into a spatially accurate 2-D
# gridded format, for plotting
elev_r = mg.node_vector_to_raster(mg['node']['topographic__elevation'])
im = pylab.plot(mg.dx * np.arange(nrows), elev_r[:, int(ncols // 2)])
pylab.xlabel('Distance')
pylab.ylabel('Elevation')

pylab.show()  # this line displays all of the figures you've issued plot commands for, since you last called show()
