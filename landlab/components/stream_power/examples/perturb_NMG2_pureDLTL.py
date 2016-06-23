# -*- coding: utf-8 -*-
#this driver runs 0.001->0.01 perturbation of NMG2 for both simple SP & simple transport limited response

from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.stream_power.stream_power import StreamPowerEroder
from landlab.components.transport_limited_fluvial.tl_fluvial_monodirectional import TransportLimitedEroder
from landlab import CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY
from landlab import ModelParameterDictionary
from landlab.plot import imshow
from landlab.plot.video_out import VideoPlotter
from landlab.plot import channel_profile as prf
from landlab.plot.imshow import imshow_node_grid
from pylab import colorbar, show, plot, loglog, figure, savefig, close, ylim, xlim, gca

from landlab import RasterModelGrid
import numpy as np
import pylab
from copy import copy, deepcopy

from time import time

show_figs_in_run = True #disable to run straight through
DL_or_TL = 'DL'


if DL_or_TL=='TL':
    init_interval = 20
else:
    init_interval = 100

#get the needed properties to build the grid:
input_file = './sed_dep_NMGparams2.txt'
#####remember to change the fixed y-axis dimension in the plots!!
y_max = 3000
make_output_plots=True
out_interval=15 #was 15
inputs = ModelParameterDictionary(input_file)
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
leftmost_elev = inputs.read_float('leftmost_elevation')
initial_slope = inputs.read_float('initial_slope')
uplift_rate = inputs.read_float('uplift_rate')

runtime = inputs.read_float('total_time')
dt = inputs.read_float('dt')

#check we have a plaubible grid
mg = RasterModelGrid(nrows,ncols,dx)
assert mg.number_of_nodes == nrows*ncols
assert mg.node_spacing == dx

# Display a message
print 'Running ...'
print 'Run to stable'
nt = int(runtime//dt)
uplift_per_step = uplift_rate * dt
print 'uplift per step: ', uplift_per_step

#create the field
z = mg.add_zeros('topographic__elevation', at='node')
z += leftmost_elev
z += initial_slope*np.amax(mg.node_y) - initial_slope*mg.node_y
#put these values plus roughness into that field
z += np.random.rand(len(z))/100000.

mg.status_at_node[mg.nodes_at_left_edge] = CLOSED_BOUNDARY
mg.status_at_node[mg.nodes_at_right_edge] = CLOSED_BOUNDARY

fr = FlowRouter(mg)
if DL_or_TL == 'TL':
    tle = TransportLimitedEroder(mg, input_file)
else:
    spe = StreamPowerEroder(mg, input_file)

for i in xrange(nt):
    # print 'loop ', i
    mg.at_node['topographic__elevation'][mg.core_nodes] += uplift_per_step
    mg = fr.route_flow(grid=mg)
    if DL_or_TL == 'TL':
        mg, _ = tle.erode(mg, dt)
    else:
        mg, _, _ = spe.erode(mg, dt=dt)
    if i % init_interval == 0:
        print 'loop ', i
        print 'max_slope', np.amax(mg.at_node['topographic__steepest_slope'][
            mg.core_nodes])
        pylab.figure("long_profiles_init")
        profile_IDs = prf.channel_nodes(
            mg, mg.at_node['topographic__steepest_slope'],
            mg.at_node['drainage_area'],
            mg.at_node['flow__receiver_node'])
        dists_upstr = prf.get_distances_upstream(mg, len(mg.at_node['topographic__steepest_slope']),
                                                profile_IDs, mg.at_node['flow__link_to_receiver_node'])
        prf.plot_profiles(dists_upstr, profile_IDs, mg.at_node['topographic__elevation'])

print 'completed run to steady state...'
if show_figs_in_run:
    show() #will cause a hang

#save a copy of the init conditions:
mg_init = deepcopy(mg)

#REinstantiate the components:
fr = FlowRouter(mg)
tle = TransportLimitedEroder(mg, input_file)
uplift_rate *= 10. #accelerate tenfold
runtime = 200000.
nt = int(runtime//dt)
uplift_per_step = uplift_rate * dt
print 'uplift per step: ', uplift_per_step

x_profiles = []
z_profiles = []
S_profiles = []
A_profiles = []

#plot init conds
if make_output_plots:
    mg = fr.route_flow(grid=mg)
    pylab.figure('long_profile_anim')
    ylim([0,y_max])
    prf.analyze_channel_network_and_plot(mg)
    savefig('0profile_anim_init.png')
    close('long_profile_anim')

(profile_IDs, dists_upstr) = prf.analyze_channel_network_and_plot(mg)
start_node = [profile_IDs[0]]

time_on = time()
#perform the loops:
for i in xrange(nt):
    #print 'loop ', i
    mg.at_node['topographic__elevation'][mg.core_nodes] += uplift_per_step
    mg = fr.route_flow(grid=mg)
    if DL_or_TL == 'TL':
        mg,_ = tle.erode(mg,dt)
    else:
        mg,_,_ = spe.erode(mg,dt=dt)
    if i%out_interval == 0:
        print 'loop ', i
        print 'max_slope', np.amax(mg.at_node['topographic__steepest_slope'][mg.core_nodes])
        pylab.figure("slope_area")
        loglog(mg.at_node['drainage_area'][profile_IDs[-1]], mg.at_node['topographic__steepest_slope'][profile_IDs[-1]], '-x')
    if i%out_interval == 0:
        x_profiles.append(dists_upstr[-1][:])
        z_profiles.append(mg.at_node['topographic__elevation'][profile_IDs[-1]])
        S_profiles.append(mg.at_node['topographic__steepest_slope'][profile_IDs[-1]])
        A_profiles.append(mg.at_node['drainage_area'][profile_IDs[-1]])
        if make_output_plots:
            pylab.figure('long_profile_anim')
            #prf.plot_profiles(dists_upstr, profile_IDs, mg.at_node['topographic_elevation'])
            plot(dists_upstr[-1],mg.at_node['topographic__elevation'][profile_IDs[-1]])
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

mg_perturbed = deepcopy(mg)

print 'Completed the simulation. Plotting...'

time_off = time()

#Finalize and plot

elev = mg['node']['topographic__elevation']
#imshow.imshow_node_grid(mg, elev)

print('Done.')
print 'Time: ', time_off-time_on

#pylab.show()

#vid.produce_video()

if True:
    pylab.figure('long_profile_anim_init')
    loglog(mg_init.at_node['drainage_area'][profile_IDs[-1]], mg_init.at_node['topographic__steepest_slope'][profile_IDs[-1]])
    ylim([1.e-3,1.e0])
    xlim(gca().get_xlim()[::-1]) #reverse the x axis for comparison with long profiles
    savefig('0profile_anim_init.png')
    close('long_profile_anim_init')
    for j in xrange(len(x_profiles)):
        i = j*15
        pylab.figure('long_profile_anim')
        #prf.plot_profiles(dists_upstr, profile_IDs, mg.at_node['topographic_elevation'])
        loglog(A_profiles[j], S_profiles[j])
        #ylim([0,y_max])
        ylim([1.e-3,1.e0])
        xlim(gca().get_xlim()[::-1]) #reverse the x axis for comparison with long profiles
        if i==0:
            savefig('profile_anim_000'+str(i)+'.png')
        elif i<100:
            savefig('profile_anim_00'+str(i)+'.png')
        elif i<1000:
            savefig('profile_anim_0'+str(i)+'.png')
        else:
            savefig('profile_anim_'+str(i)+'.png')
        close('long_profile_anim')
