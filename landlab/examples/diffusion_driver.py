from landlab.components.nonlinear_diffusion.Perron_nl_diffuse import PerronNLDiffuse
from landlab.components.diffusion.diffusion import DiffusionComponent #...the two different diffusion formulations
from landlab import ModelParameterDictionary #handles input from the input file

from landlab import RasterModelGrid #the grid object
from landlab.plot.imshow import imshow_node_grid
import numpy as np
import pylab

#get the needed properties to build the grid:
input_file = './diffusion_test_params.txt'
#initialize an object that will supply the parameters:
inputs = ModelParameterDictionary(input_file)
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
leftmost_elev = inputs.read_float('leftmost_elevation')
initial_slope = inputs.read_float('initial_slope')
uplift_rate = inputs.read_float('uplift_rate')
runtime = inputs.read_float('total_time')
dt = inputs.read_float('dt')
nt = int(runtime//dt) #"//" means "divide and truncate" ("%" means "division remainder")
uplift_per_step = uplift_rate * dt

#Instantiate the grid object
#We know which parameters are needed for input by inspecting the function where it lives, in landlab.grid.raster
#We could also look at the documentation for landlab found online (http://the-landlab.readthedocs.org)
mg = RasterModelGrid(nrows, ncols, dx)

##create the elevation field in the grid:
#create the field
mg.create_node_array_zeros('topographic_elevation')
z = mg.create_node_array_zeros() + leftmost_elev #in our case, slope is zero, so the leftmost_elev is the mean elev
#put these values plus roughness into that field
mg['node'][ 'topographic_elevation'] = z + np.random.rand(len(z))/100000.

#set up its boundary conditions (bottom, left, top, right)
#The mechanisms for this are all automated within the grid object
mg.set_fixed_value_boundaries_at_grid_edges(True, True, True, True)

# Display a message
print 'Running ...' 

#instantiate the components:
diffuse = PerronNLDiffuse(mg, input_file)
lin_diffuse = DiffusionComponent(grid=mg, input_stream=input_file)
#lin_diffuse.initialize(input_file)

#Perform the loops.
#First, we do the nonlinear diffusion:

#We're going to perform a block uplift of the interior of the grid, but leave the boundary nodes at their original elevations.
uplifted_nodes = mg.get_core_nodes() #access this function of the grid, and store the output with a local name
#(Note: Node numbering runs across from the bottom left of the grid.)

for i in xrange(nt): #nt is the number of timesteps we calculated above, i.e., loop nt times. We never actually use i within the loop, but we could do.
    #("xrange" is a clever memory-saving way of producing consecutive integers to govern a loop)
    #this colon-then-tab-in arrangement is what Python uses to delineate connected blocks of text, instead of brackets or parentheses
    #This line performs the actual functionality of the component:
    #mg = lin_diffuse.diffuse(mg, dt) #linear diffusion
    mg = diffuse.diffuse(mg, i*dt) #nonlinear diffusion
    #...swap around which line is commented out to switch between formulations of diffusion

    #now plot a N-S cross section from this stage in the run onto figure 1. The sections will all be superimposed, as show() hasn't yet been called
    pylab.figure(1)
    elev_r = mg.node_vector_to_raster(mg['node']['topographic_elevation']) #turn the 1-D array of elevation values into a spatially accurate 2-D gridded format, for plotting
    im = pylab.plot(mg.dx*np.arange(nrows), elev_r[:,int(ncols//2)]) #square brackets denote a subset of nodes to use.
    #...this kind of data extraction from a larger data structure ("slicing", or "fancy indexing") is extremely useful and powerful, and is one of the appeals of Python
    #...this is plot(x, y).
    #x is the distance north up the grid.
    #y is the elevation along all the rows, but only the 50th column (more slicing!), i.e., halfway along the grid and N-S

    print 'Completed loop ', i
 
print 'Completed the simulation. Plotting...'

#Finalize and plot:
#put a title on figure 1
pylab.figure(1)
pylab.title('N-S cross_section, nonlinear diffusion')
pylab.xlabel('Distance')
pylab.ylabel('Elevation')

#figure 2 is the map of the final elevations
pylab.figure(2)
im_nl = imshow_node_grid(mg, 'topographic_elevation')  # display a colored image

pylab.figure(3)
elev_r = mg.node_vector_to_raster(mg['node']['topographic_elevation']) #turn the 1-D array of elevation values into a spatially accurate 2-D gridded format, for plotting
im = pylab.plot(mg.dx*np.arange(nrows), elev_r[:,int(ncols//2)])

print('Done.')

#now do the linear diffusion:

##Reset the elevation field in the grid:
mg['node'][ 'topographic_elevation'] = z + np.random.rand(len(z))/100000.

# Display a message
print 'Running ...' 

##instantiate the components:
#diffuse = PerronNLDiffuse(mg, input_file)
#lin_diffuse = DiffusionComponent(grid=mg)
#lin_diffuse.initialize(input_file)

for i in xrange(nt):
    #This line performs the actual functionality of the component:
    #***NB: the nonlinear diffuser contains an "automatic" element of uplift. If you instead use the linear diffuser, you need to add the uplift manually...
    mg['node']['topographic_elevation'][uplifted_nodes] += uplift_per_step 
    mg = lin_diffuse.diffuse(mg, internal_uplift=False) #linear diffusion

    pylab.figure(4)
    elev_r = mg.node_vector_to_raster(mg['node']['topographic_elevation'])
    im = pylab.plot(mg.dx*np.arange(nrows), elev_r[:,int(ncols//2)])

    print 'Completed loop ', i

print 'Completed the simulation. Plotting...'

#Finalize and plot:
#put a title on figure 4
pylab.figure(4)
pylab.title('N-S cross_section, linear diffusion')
pylab.xlabel('Distance')
pylab.ylabel('Elevation')

#figure 5 is the map of the final elevations
pylab.figure(5)
im = imshow_node_grid(mg, 'topographic_elevation')

#superpose this final form onto figure 3:
pylab.figure(3)
elev_r = mg.node_vector_to_raster(mg['node']['topographic_elevation']) #turn the 1-D array of elevation values into a spatially accurate 2-D gridded format, for plotting
im = pylab.plot(mg.dx*np.arange(nrows), elev_r[:,int(ncols//2)])
pylab.xlabel('Distance')
pylab.ylabel('Elevation')

pylab.show() #this line displays all of the figures you've issued plot commands for, since you last called show()
