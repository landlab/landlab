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
mg['node'][ 'planet_surface__elevation'] = np.load('init.npy')

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
redug_crater = np.empty(nt)
for i in xrange(loops):
    for j in xrange(nt):
        mg = craters_component.excavate_a_crater_noangle(mg)
        x[j] = craters_component.impact_property_dict['x']
        y[j] = craters_component.impact_property_dict['y']
        r[j] = craters_component.impact_property_dict['r']
        slope[j] = craters_component.impact_property_dict['surface_slope']
        angle[j] = craters_component.impact_property_dict['normal_angle']
        az[j] = craters_component.impact_property_dict['impact_az']
        mass_balance[j] = craters_component.impact_property_dict['mass_balance']
        redug_crater[j] = craters_component.impact_property_dict['redug_crater']
        print 'Completed loop ', j
    mystring = 'craterssave'+str((i+1)*nt)
    np.save(mystring,mg['node']['planet_surface__elevation'])
    #Save the properties
    np.save(('x_'+str((i+1)*nt)),x)
    np.save(('y_'+str((i+1)*nt)),y)
    np.save(('r_'+str((i+1)*nt)),r)
    np.save(('slope_'+str((i+1)*nt)),slope)
    np.save(('angle_'+str((i+1)*nt)),angle)
    np.save(('az_'+str((i+1)*nt)),az)
    np.save(('mass_balance_'+str((i+1)*nt)),mass_balance)
    np.save(('redug_crater_'+str((i+1)*nt)),redug_crater)

#Finalize and plot
elev = mg['node']['planet_surface__elevation']
elev_r = mg.node_vector_to_raster(elev)
# Clear previous plots
#pylab.figure(1)
#pylab.close()
#pylab.figure(1)
#im = pylab.imshow(elev_r, cmap=pylab.cm.RdBu)  # display a colored image
#pylab.colorbar(im)
#pylab.title('Topography')

print('Done.')
print('Total run time = '+str(time.time()-start_time)+' seconds.')

#pylab.show()
