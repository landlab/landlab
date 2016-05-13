# -*- coding: utf-8 -*-
from __future__ import print_function

from six.moves import range

from landlab.components.flow_routing import FlowRouter
from landlab.components.stream_power import SedDepEroder
from landlab import ModelParameterDictionary
from landlab.plot import imshow
from landlab.plot.video_out import VideoPlotter
from landlab.plot import channel_profile as prf
from landlab.plot.imshow import imshow_node_grid
from pylab import colorbar, show, plot, loglog, figure, savefig, close, ylim

from landlab import RasterModelGrid
import numpy as np
import pylab
from copy import copy, deepcopy

from time import time

#get the needed properties to build the grid:
input_file = './sed_dep_NMGparams2.txt'
#####remember to change the fixed y-axis dimension in the plots!!
y_max = 200
make_output_plots=True
out_interval=15 #was 15
inputs = ModelParameterDictionary(input_file)
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
uplift_rate = inputs.read_float('uplift_rate')

runtime = inputs.read_float('total_time')
dt = inputs.read_float('dt')

nt = int(runtime//dt)
uplift_per_step = uplift_rate * dt
print('uplift per step: ', uplift_per_step)

#check we have a plaubible grid
#mg = RasterModelGrid(nrows,ncols,dx)
assert mg.number_of_nodes == nrows*ncols
assert mg.node_spacing == dx

# Display a message
print('Running ...')

# instantiate the components:
fr = FlowRouter(mg)
sde = SedDepEroder(mg, input_file)
# don't allow overwriting of these, just in case
try:
    x_profiles
except NameError:
    x_profiles = []
    z_profiles = []
    S_profiles = []
    A_profiles = []

# plot init conds
if make_output_plots:
    mg = fr.route_flow(grid=mg)
    pylab.figure('long_profile_anim')
    ylim([0, y_max])
    prf.analyze_channel_network_and_plot(mg)
    savefig('0profile_anim_init.png')
    close('long_profile_anim')

(profile_IDs, dists_upstr) = prf.analyze_channel_network_and_plot(mg)
start_node = [profile_IDs[0]]

time_on = time()
#perform the loops:
for i in range(nt):
    #print 'loop ', i
    mg.at_node['topographic__elevation'][mg.core_nodes] += uplift_per_step
    mg = fr.route_flow()
    #mg.calc_grad_across_cell_faces(mg.at_node['topographic__elevation'])
    #neighbor_slopes = mg.calc_grad_along_node_links(mg.at_node['topographic__elevation'])
    #mean_slope = np.mean(np.fabs(neighbor_slopes),axis=1)
    #max_slope = np.max(np.fabs(neighbor_slopes),axis=1)
    #mg,_,capacity_out = tl.erode(mg,dt,slopes_at_nodes='topographic__steepest_slope')
    #mg,_,capacity_out = tl.erode(mg,dt,slopes_at_nodes=max_slope)
    mg_copy = deepcopy(mg)
    mg,_ = sde.erode(mg,dt)
    #print sde.iterations_in_dt
    #print 'capacity ', np.amax(capacity_out[mg.core_nodes])
    #print 'rel sed ', np.nanmax(sed_in[mg.core_nodes]/capacity_out[mg.core_nodes])
    if i%out_interval == 0:
        print('loop ', i)
        print('max_slope', np.amax(mg.at_node['topographic__steepest_slope'][mg.core_nodes]))
        pylab.figure("long_profiles")
        profile_IDs = prf.channel_nodes(mg, mg.at_node['topographic__steepest_slope'],
                                        mg.at_node['drainage_area'], mg.at_node['flow__receiver_node'])
        dists_upstr = prf.get_distances_upstream(mg, len(mg.at_node['topographic__steepest_slope']),
                                        profile_IDs, mg.at_node['flow__link_to_receiver_node'])
        prf.plot_profiles(dists_upstr, profile_IDs, mg.at_node['topographic__elevation'])
    if i%out_interval == 0:
        x_profiles.append(dists_upstr)
        z_profiles.append(mg.at_node['topographic__elevation'][profile_IDs])
        S_profiles.append(mg.at_node['topographic__steepest_slope'][profile_IDs])
        A_profiles.append(mg.at_node['drainage_area'][profile_IDs])
        if make_output_plots:
            pylab.figure('long_profile_anim')
            #prf.plot_profiles(dists_upstr, profile_IDs, mg.at_node['topographic_elevation'])
            plot(dists_upstr,mg.at_node['topographic_elevation'][profile_IDs])
            ylim([0,y_max])
            if i==0:
                savefig('profile_anim_000'+str(i)+'.png')
            elif i<100:
                savefig('profile_anim_00'+str(i)+'.png')
            elif i<1000:
                savefig('profile_anim_0'+str(i)+'.png')
            else:
                savefig('profile_anim_'+str(i)+'.png')
            close('long_profile_anim')
#vid.add_frame(mg, 'topographic__elevation')


print('Completed the simulation. Plotting...')

time_off = time()

#Finalize and plot

elev = mg['node']['topographic__elevation']
#imshow.imshow_node_grid(mg, elev)

print('Done.')
print('Time: ', time_off-time_on)

#pylab.show()

#vid.produce_video()
