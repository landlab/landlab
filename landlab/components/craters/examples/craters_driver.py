from landlab.components.craters.dig_craters import impactor
from landlab import ModelParameterDictionary

from landlab import RasterModelGrid
import numpy as np
import time
import pylab

#get the needed properties to build the grid:
input_file = './craters_params_init.txt'
inputs = ModelParameterDictionary(input_file)
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
leftmost_elev = inputs.read_float('leftmost_elevation')
initial_slope = inputs.read_float('initial_slope')
nt = inputs.read_int('number_of_craters_per_loop')
loops = inputs.read_int('number_of_loops')

mg = RasterModelGrid(nrows, ncols, dx)
mg.set_looped_boundaries(True, True)

#create the fields in the grid
mg.create_node_array_zeros('planet_surface__elevation')
z = mg.create_node_array_zeros() + leftmost_elev
z += initial_slope*np.amax(mg.node_y) - initial_slope*mg.node_y
mg['node'][ 'planet_surface__elevation'] = z #+ np.random.rand(len(z))/10000.

# Display a message
print( 'Running ...' )
start_time = time.time()

#instantiate the component:
craters_component = impactor(mg, input_file)

#perform the loops:
x = np.empty(nt)
y = np.empty(nt)
r = np.empty(nt)
slope = np.empty(nt)
angle = np.empty(nt)
az = np.empty(nt)
mass_balance = np.empty(nt)
for i in xrange(loops):
    for j in xrange(nt):
        if j == 1:
            craters_component._ycoord = 0.45*(mg.get_grid_xdimension()-mg.dx)
            craters_component._radius /= 10.
        mg = craters_component.excavate_a_crater_furbish(mg)
        x[j] = craters_component.impact_xy_location[0]
        y[j] = craters_component.impact_xy_location[1]
        r[j] = craters_component.crater_radius
        slope[j] = craters_component.surface_slope_beneath_crater
        angle[j] = craters_component.impact_angle_to_normal
        az[j] = craters_component.impactor_travel_azimuth
        mass_balance[j] = craters_component.mass_balance
        print 'Completed loop ', j
    mystring = 'initial_topo'  #'craterssave'+str((i+1)*nt)
    np.save(mystring,mg['node']['planet_surface__elevation'])
    #Save the properties
    np.save(('x_'+str((i+1)*nt)),x)
    np.save(('y_'+str((i+1)*nt)),y)
    np.save(('r_'+str((i+1)*nt)),r)
    np.save(('slope_'+str((i+1)*nt)),slope)
    np.save(('angle_'+str((i+1)*nt)),angle)
    np.save(('az_'+str((i+1)*nt)),az)
    np.save(('mass_balance_'+str((i+1)*nt)),mass_balance)

#Finalize and plot
elev = mg['node']['planet_surface__elevation']
elev_r = mg.node_vector_to_raster(elev)
# Clear previous plots
pylab.figure(1)
pylab.close()
pylab.figure(1)
im = pylab.imshow(elev_r, cmap=pylab.cm.RdBu)  # display a colored image
pylab.colorbar(im)
pylab.title('Topography')

print('Done.')
print('Total run time = '+str(time.time()-start_time)+' seconds.')

pylab.show()

im=pylab.plot(elev_r[:,ncols//2])
pylab.show()