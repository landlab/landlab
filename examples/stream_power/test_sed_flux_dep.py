# -*- coding: utf-8 -*-
from __future__ import print_function

from copy import deepcopy
from time import time

import numpy as np
import pylab
from six.moves import range

from landlab import ModelParameterDictionary, RasterModelGrid
from landlab.components import FlowAccumulator, SedDepEroder
from landlab.plot import channel_profile as prf
from landlab.plot.video_out import VideoPlotter

# get the needed properties to build the grid:
input_file = "./sed_dep_params_reliable.txt"  # './sed_dep_NMGparams2.txt'
inputs = ModelParameterDictionary(input_file)
nrows = inputs.read_int("nrows")
ncols = inputs.read_int("ncols")
dx = inputs.read_float("dx")
leftmost_elev = inputs.read_float("leftmost_elevation")
initial_slope = inputs.read_float("initial_slope")
uplift_rate = inputs.read_float("uplift_rate")

runtime = inputs.read_float("total_time")
dt = inputs.read_float("dt")

nt = int(runtime // dt)
uplift_per_step = uplift_rate * dt
print("uplift per step: ", uplift_per_step)

# instantiate the grid object
mg = RasterModelGrid((nrows, ncols), xy_spacing=dx)

# create the elevation field in the grid:
# create the field
mg.add_zeros("topographic__elevation", at="node")
z = mg.zeros(at="node") + leftmost_elev
z += initial_slope * np.amax(mg.node_y) - initial_slope * mg.node_y
# put these values plus roughness into that field
mg["node"]["topographic__elevation"] = z + np.random.rand(len(z)) / 100000.

# set up grid's boundary conditions (bottom, left, top, right is inactive)
mg.set_closed_boundaries_at_grid_edges(True, False, True, False)
mg.set_fixed_value_boundaries_at_grid_edges(
    False, True, False, True, value_of="topographic__elevation"
)
print("fixed vals in grid: ", mg.fixed_value_node_properties["values"])

# Display a message
print("Running ...")

# instantiate the components:
fr = FlowAccumulator(mg, flow_director="D8")
sde = SedDepEroder(mg, input_file)
vid = VideoPlotter(mg, data_centering="node")

time_on = time()
# perform the loops:
for i in range(nt):
    # print 'loop ', i
    mg.at_node["topographic__elevation"][mg.core_nodes] += uplift_per_step
    mg = fr.run_one_step()
    # mg.calc_grad_across_cell_faces(mg.at_node['topographic__elevation'])
    # neighbor_slopes = mg.calc_grad_along_node_links(mg.at_node['topographic__elevation'])
    # mean_slope = np.mean(np.fabs(neighbor_slopes),axis=1)
    # max_slope = np.max(np.fabs(neighbor_slopes),axis=1)
    # mg,_,capacity_out = tl.erode(mg,dt,slopes_at_nodes='topographic__steepest_slope')
    # mg,_,capacity_out = tl.erode(mg,dt,slopes_at_nodes=max_slope)
    mg_copy = deepcopy(mg)
    mg, _ = sde.erode(mg, dt)
    # print sde.iterations_in_dt
    # print 'capacity ', np.amax(capacity_out[mg.core_nodes])
    # print 'rel sed ',
    # np.nanmax(sed_in[mg.core_nodes]/capacity_out[mg.core_nodes])
    if i % 100 == 0:
        print("loop ", i)
        print(
            "max_slope",
            np.amax(mg.at_node["topographic__steepest_slope"][mg.core_nodes]),
        )
        pylab.figure("long_profiles")
        profile_IDs = prf.channel_nodes(
            mg,
            mg.at_node["topographic__steepest_slope"],
            mg.at_node["drainage_area"],
            mg.at_node["flow__receiver_node"],
        )
        dists_upstr = prf.get_distances_upstream(
            mg,
            len(mg.at_node["topographic__steepest_slope"]),
            profile_IDs,
            mg.at_node["flow__link_to_receiver_node"],
        )
        prf.plot_profiles(
            dists_upstr, profile_IDs, mg.at_node["topographic__elevation"]
        )
    # vid.add_frame(mg, 'topographic__elevation')


print("Completed the simulation. Plotting...")

time_off = time()

# Finalize and plot

elev = mg["node"]["topographic__elevation"]
# imshow.imshow_grid(mg, elev)

print("Done.")
print("Time: ", time_off - time_on)

# pylab.show()

# vid.produce_video()
