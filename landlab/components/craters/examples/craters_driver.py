from landlab.components.craters.dig_craters import impactor
from landlab import ModelParameterDictionary

from landlab import RasterModelGrid
import numpy as np
import time
import pylab

#get the needed properties to build the grid:
input_file = './craters_params.txt'
inputs = ModelParameterDictionary(input_file)
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
leftmost_elev = inputs.read_float('leftmost_elevation')
initial_slope = inputs.read_float('initial_slope')
nt = inputs.read_int('number_of_craters_per_loop')
loops = inputs.read_int('number_of_loops')

mg = RasterModelGrid(nrows, ncols, dx)
mg.set_inactive_boundaries(False, False, False, False)

#create the fields in the grid
mg.create_node_array_zeros('planet_surface__elevation')
z = mg.create_node_array_zeros() + leftmost_elev
z += initial_slope*np.amax(mg.node_y) - initial_slope*mg.node_y
mg['node'][ 'planet_surface__elevation'] = z + np.random.rand(len(z))/10000.

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
        mg = craters_component.excavate_a_crater(mg)
        x[j] = craters_component.impact_property_dict['x']
        y[j] = craters_component.impact_property_dict['y']
        r[j] = craters_component.impact_property_dict['r']
        slope[j] = craters_component.impact_property_dict['surface_slope']
        angle[j] = craters_component.impact_property_dict['normal_angle']
        az[j] = craters_component.impact_property_dict['impact_az']
        mass_balance[j] = craters_component.impact_property_dict['mass_balance']
        print 'Completed loop ', j
    mystring = 'craterssave'+str(i*nt)
    np.save(mystring,mg['node']['planet_surface__elevation'])
    #Save the properties
    np.save(('x_'+str(i*nt)),x)
    np.save(('y_'+str(i*nt)),y)
    np.save(('r_'+str(i*nt)),r)
    np.save(('slope_'+str(i*nt)),slope)
    np.save(('angle_'+str(i*nt)),angle)
    np.save(('az_'+str(i*nt)),az)
    np.save(('mass_balance_'+str(i*nt)),mass_balance)

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
