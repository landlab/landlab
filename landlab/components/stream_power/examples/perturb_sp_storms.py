"""
A driver implementing Braun-Willett flow routing and then a
(non-fastscape) stream power component.
This version runs the model to something approximating steady
state, then perturbs the uplift rate to produce a propagating
wave, then stores the propagation as a gif.
"""
# DEJH, 09/15/14
from __future__ import print_function

from landlab.components.flow_routing import FlowRouter
from landlab.components.stream_power import StreamPowerEroder, FastscapeEroder
from landlab.components.uniform_precip import PrecipitationDistribution
from landlab.plot import channel_profile as prf
from landlab.plot import imshow as llplot
from landlab.plot.imshow import imshow_node_grid
from pylab import figure, plot, xlabel, ylabel, title, show

import numpy
from landlab import RasterModelGrid
from landlab import ModelParameterDictionary
import pylab
import time
import copy

input_file_string = './drive_sp_params_storms.txt'
inputs = ModelParameterDictionary(input_file_string)
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
dt = inputs.read_float('dt')
time_to_run = inputs.read_float('run_time')
#nt needs defining
uplift = inputs.read_float('uplift_rate')
init_elev = inputs.read_float('init_elev')

mg = RasterModelGrid(nrows, ncols, dx)

#create the fields in the grid
mg.add_zeros('topographic__elevation', at='node')
z = mg.zeros(at='node') + init_elev
mg['node'][ 'topographic__elevation'] = z + numpy.random.rand(len(z))/1000.
mg.add_zeros('node', 'water__unit_flux_in')

#make some K values in a field to test
#mg.at_node['K_values'] = 0.1+numpy.random.rand(nrows*ncols)/10.
mg.at_node['K_values'] = numpy.empty(nrows*ncols, dtype=float)
#mg.at_node['K_values'].fill(0.1+numpy.random.rand()/10.)
mg.at_node['K_values'].fill(0.001)

print( 'Running ...' )

#instantiate the components:
fr = FlowRouter(mg)
sp = StreamPowerEroder(mg, input_file_string)
#fsp = FastscapeEroder(mg, input_file_string)
precip = PrecipitationDistribution(input_file=input_file_string)

#load the Fastscape module too, to allow direct comparison
fsp = FastscapeEroder(mg, input_file_string)

try:
    #raise NameError
    mg = copy.deepcopy(mg_mature)
except NameError:
    print('building a new grid...')
    out_interval = 50000.
    last_trunc = time_to_run #we use this to trigger taking an output plot
    #run to a steady state:
    #We're going to cheat by running Fastscape SP for the first part of the solution
    for (interval_duration, rainfall_rate) in precip.yield_storm_interstorm_duration_intensity():
        if rainfall_rate != 0.:
            mg.at_node['water__unit_flux_in'].fill(rainfall_rate)
            mg = fr.route_flow()
            #print 'Area: ', numpy.max(mg.at_node['drainage_area'])
            mg,_,_ = sp.erode(mg, interval_duration, Q_if_used='surface_water__discharge', K_if_used='K_values')
        #add uplift
        mg.at_node['topographic__elevation'][mg.core_nodes] += uplift*interval_duration
        this_trunc = precip.elapsed_time//out_interval
        if this_trunc != last_trunc: #a new loop round
            print('made it to loop ', out_interval*this_trunc)
            last_trunc=this_trunc

    mg_mature = copy.deepcopy(mg)

else:
    #reinstantiate the components with the new grid
    fr = FlowRouter(mg)
    sp = StreamPowerEroder(mg, input_file_string)

x_profiles = []
z_profiles = []
S_profiles = []
A_profiles = []

if True:
    #perturb:
    time_to_run = 300000.
    precip_perturb = PrecipitationDistribution(input_file=input_file_string, mean_storm=10., mean_interstorm=40., total_t=time_to_run)
    out_interval = 5000.
    last_trunc = time_to_run #we use this to trigger taking an output plot
    for (interval_duration, rainfall_rate) in precip_perturb.yield_storm_interstorm_duration_intensity():
        if rainfall_rate != 0.:
            mg.at_node['water__unit_flux_in'].fill(rainfall_rate)
            mg = fr.route_flow() #the runoff_rate should pick up automatically
            #print 'Area: ', numpy.max(mg.at_node['drainage_area'])
            mg,_,_ = sp.erode(mg, interval_duration, Q_if_used='surface_water__discharge', K_if_used='K_values')

        #plot long profiles along channels
        this_trunc = precip_perturb.elapsed_time//out_interval
        if this_trunc != last_trunc: #a new loop round
            print('saving a plot at perturbed loop ', out_interval*this_trunc)
            pylab.figure("long_profiles")
            profile_IDs = prf.channel_nodes(mg, mg.at_node['topographic__steepest_slope'],
                            mg.at_node['drainage_area'], mg.at_node['flow__receiver_node'])
            dists_upstr = prf.get_distances_upstream(mg, len(mg.at_node['topographic__steepest_slope']),
                            profile_IDs, mg.at_node['flow__link_to_receiver_node'])
            prf.plot_profiles(dists_upstr, profile_IDs, mg.at_node['topographic__elevation'])
            last_trunc=this_trunc
            x_profiles.append(dists_upstr[:])
            z_profiles.append(mg.at_node['topographic_elevation'][profile_IDs])
            S_profiles.append(mg.at_node['steepest_slope'][profile_IDs])
            A_profiles.append(mg.at_node['drainage_area'][profile_IDs])
    
        #add uplift
        mg.at_node['topographic__elevation'][mg.core_nodes] += 5.*uplift*interval_duration

#Finalize and plot
elev = mg['node']['topographic__elevation']
elev_r = mg.node_vector_to_raster(elev)

# Clear previous plots
pylab.figure("topo")
pylab.close()

# Plot topography
pylab.figure("topo")
#im = pylab.imshow(elev_r, cmap=pylab.cm.RdBu)  # display a colored image
im = llplot.imshow_node_grid(mg, 'topographic__elevation')
#print elev_r
#pylab.colorbar(im)
#pylab.title('Topography')

pylab.figure("Xsec")
im = pylab.plot(dx*numpy.arange(nrows), elev_r[:,int(ncols//2)])  # display a colored image
pylab.title('Vertical cross section')

pylab.figure("Slope-Area")
im = pylab.loglog(mg.at_node['drainage_area'], mg.at_node['topographic__steepest_slope'],'.')
pylab.title('Slope-Area')

pylab.show()

print('Done.')


