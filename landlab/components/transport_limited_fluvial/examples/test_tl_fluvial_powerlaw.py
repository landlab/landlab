from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.transport_limited_fluvial.tl_fluvial_monodirectional_v3 import TransportLimitedEroder
#from landlab.components.transport_limited_fluvial.tl_fluvial_polydirectional import TransportLimitedEroder
from landlab import ModelParameterDictionary
from landlab.plot import imshow
from landlab.plot.video_out import VideoPlotter
from landlab.plot import channel_profile as prf

from landlab import RasterModelGrid
import numpy as np
import pylab
from pylab import show
from copy import copy, deepcopy

from time import time

#get the needed properties to build the grid:
input_file = './stream_power_params_powerlaw.txt'
inputs = ModelParameterDictionary(input_file)
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
leftmost_elev = inputs.read_float('leftmost_elevation')
initial_slope = inputs.read_float('initial_slope')
uplift_rate = inputs.read_float('uplift_rate')

runtime = inputs.read_float('total_time')
dt = inputs.read_float('dt')

nt = int(runtime//dt)
uplift_per_step = uplift_rate * dt
print 'uplift per step: ', uplift_per_step

#instantiate the grid object
mg = RasterModelGrid(nrows, ncols, dx)

##create the elevation field in the grid:
#create the field
mg.create_node_array_zeros('topographic_elevation')
z = mg.create_node_array_zeros() + leftmost_elev
z += initial_slope*np.amax(mg.node_y) - initial_slope*mg.node_y
#put these values plus roughness into that field
mg['node'][ 'topographic_elevation'] = z + np.random.rand(len(z))/100000.

#set up grid's boundary conditions (bottom, left, top, right is inactive)
mg.set_inactive_boundaries(False, True, False, True)
mg.set_fixed_value_boundaries_at_grid_edges(True, False, True, False, value_of='topographic_elevation')
print 'fixed vals in grid: ', mg.fixed_value_node_properties['values']

# Display a message
print 'Running ...' 

#instantiate the components:
fr = FlowRouter(mg)
tl = TransportLimitedEroder(mg, input_file)
tl = TransportLimitedEroder(mg, input_file)
vid = VideoPlotter(mg, data_centering='node')

time_on = time()
#perform the loops:
for i in xrange(nt):
    #print 'loop ', i
    mg.at_node['topographic_elevation'][mg.core_nodes] += uplift_per_step
    mg = fr.route_flow(grid=mg)
    #mg.calculate_gradient_across_cell_faces(mg.at_node['topographic_elevation'])
    #neighbor_slopes = mg.calculate_gradient_along_node_links(mg.at_node['topographic_elevation'])
    #mean_slope = np.mean(np.fabs(neighbor_slopes),axis=1)
    #max_slope = np.max(np.fabs(neighbor_slopes),axis=1)
    #mg,_,capacity_out = tl.erode(mg,dt,slopes_at_nodes='steepest_slope')
    #mg,_,capacity_out = tl.erode(mg,dt,slopes_at_nodes=max_slope)
    mg_copy = deepcopy(mg)
    mg,_ = tl.erode(mg,dt,stability_condition='loose')
    if i%20 == 0:
        print 'loop ', i
        print 'subdivisions of dt used: ', tl.iterations_in_dt
        print 'max_slope', np.amax(mg.at_node['steepest_slope'][mg.core_nodes])
        pylab.figure("long_profiles")
        profile_IDs = prf.channel_nodes(mg, mg.at_node['steepest_slope'],
                        mg.at_node['drainage_area'], mg.at_node['upstream_ID_order'],
                        mg.at_node['flow_receiver'])
        dists_upstr = prf.get_distances_upstream(mg, len(mg.at_node['steepest_slope']),
                        profile_IDs, mg.at_node['links_to_flow_receiver'])
        prf.plot_profiles(dists_upstr, profile_IDs, mg.at_node['topographic_elevation'])
    #mg.update_boundary_nodes()
    #vid.add_frame(mg, 'topographic_elevation')
    
 
print 'Completed the simulation. Plotting...'

time_off = time()

#Finalize and plot

elev = mg['node']['topographic_elevation']
#imshow.imshow_node_grid(mg, elev)

print('Done.')
print 'Time: ', time_off-time_on

#pylab.show()

#vid.produce_video()
