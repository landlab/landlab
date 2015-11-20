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

# Perform the loops.

for i in xrange(nt):
    # This line performs the actual functionality of the component:
    #***NB: both diffusers contain an "automatic" element of uplift.
    # You can suppress this for the linear diffuser with the *internal_uplift* keyword, =False
    # See the docstrings for both classes for more details.

    # Switch these lines to switch between diffusion styles:
    # mg = diffuse.diffuse(mg, i*dt) #nonlinear diffusion
    lin_diffuse.diffuse(dt)  # linear diffusion

    # Plot a Xsection north-south through the middle of the data, once per loop
    pylab.figure(1)
    elev_r = mg.node_vector_to_raster(mg['node']['topographic__elevation'])
    im = pylab.plot(mg.dx * np.arange(nrows), elev_r[:, int(ncols // 2)])

    print('Completed loop ', i)

print('Completed the simulation. Plotting...')

# Finalize and plot:
# put a title on figure 1
pylab.figure(1)
pylab.title('N-S cross_section')
pylab.xlabel('Distance')
pylab.ylabel('Elevation')

# figure 2 is the map of the final elevations
pylab.figure(2)
imshow_node_grid(mg, 'topographic__elevation')

pylab.show()  # this line displays all of the figures you've issued plot commands for, since you last called show()
