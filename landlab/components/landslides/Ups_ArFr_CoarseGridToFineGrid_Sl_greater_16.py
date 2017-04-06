
###########################################################################
## 26 Jan 2016 - SN & EI
## Authors: Sai Nudurupati & Erkan Istanbulluoglu
##
## Derived from Recharge_estimation.py
## - convert a coarser grid to finer grid
## - calculate fractions of upstream drainage area from each
##   Vic Grid cell
## - These fractions will be used to estimate recharge at
##   each grid cell of the finer grid
##
## This code is exactly the same as Upstream_AreaFraction_CoarserGridToFinerGrid.py
## - Only difference is that we will consider slopes > 16 degrees using a
##   mask
##
## Latest Edit - SN 26Jan16
###########################################################################

# %% Import required libraries
import numpy as np
#from landlab import RasterModelGrid as rmg
#from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.io.esri_ascii import read_esri_ascii
#from landlab.plot.imshow import imshow_grid, imshow_field
import copy
from collections import Counter
from datetime import datetime
import cPickle as pickle

# Start the clock
startTime = datetime.now()

# %% 
#Start of Script
grid,z = \
  read_esri_ascii('./Input_files/elevation.txt',\
    name='topographic__elevation')
grid.set_nodata_nodes_to_closed(grid['node']['topographic__elevation'], -9999.)

grid, flow_dir_arc = \
 read_esri_ascii('./Input_files/flow_direction.txt',\
   name='flow_dir',grid=grid)

grid,vic_ids = \
  read_esri_ascii('./Input_files/vic_idsnoca.txt',\
    name='vic_id', grid=grid)

grid,slp_g16 = \
  read_esri_ascii('./Input_files/slp_g16msk.txt',\
    name='slp_g16', grid=grid)

vic_ids = vic_ids.astype(int)

grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
grid.set_nodata_nodes_to_closed(grid['node']['flow_dir'], -9999.)
grid.set_nodata_nodes_to_closed(grid['node']['slp_g16'], -9999.)

r_arc_raw = np.log2(flow_dir_arc) # gives receiver of each node starting at
                                  # right and going clockwise
r_arc_raw = r_arc_raw.astype('int')
# %% Convert Arc flow_directions file to be represented in node ids

neigh_ = grid.neighbors_at_node
diag_ = grid.diagonals_at_node
neigh_ = np.fliplr(neigh_)
diag_ = np.fliplr(diag_)
a_n = np.hsplit(neigh_,4)
a_d = np.hsplit(diag_,4)
neighbors = np.hstack((a_n[-1],a_d[0],a_n[0],a_d[1],a_n[1],a_d[2],a_n[2],  \
                            a_d[3]))
# Now neighbors has node ids of neighboring nodes in cw order starting at
# right, hence the order of neighbors = [r, br, b, bl, l, tl, t, tr]
r = np.zeros(grid.number_of_nodes,dtype=int)
r[grid.core_nodes] = np.choose(r_arc_raw[grid.core_nodes], \
        np.transpose(neighbors[grid.core_nodes]))
#r = np.choose(r_arc_raw, np.transpose(neighbors))
# %% Source Routing Algorithm
# Note 1: This algorithm works on core nodes only because core nodes
# have neighbors that are real values and not -1s.
# Note 2: Nodes in the following comments in this section refer to core nodes.
core_nodes = grid.core_nodes
core_elev = z[core_nodes]
# Sort all nodes in the descending order of elevation
sor_z = core_nodes[np.argsort(core_elev)[::-1]]
# Create a list to record all nodes that have been visited
alr_counted = [] # To store nodes that have already been counted
r_core=r[core_nodes]
flow_accum = np.zeros(grid.number_of_nodes,dtype=int)
vic_upstr = {}
# Loop through all nodes
for i in sor_z:
    #print 'i= ', i
    # Check 1: Check if this node has been visited earlier. If yes,
    # then skip to next node
    if i in alr_counted:
        continue
    # Check 2: If the visited node is a sink
    if r[i]==i:
        vic_upstr.update({i:[vic_ids[i]]})
        flow_accum[i] += 1.
        alr_counted.append(i)
        continue
    # Check 3: Now, if the node is not a sink and hasn't been visited, it
    # belongs to a stream segment. Hence, all the nodes in the stream will
    # have to betraversed.
    # stream_buffer is a list that will hold the upstream contributing
    # node information for that particular segment until reaching the outlet.
    stream_buffer = []
    j = i
    switch_i = True
    a = 0.
    # Loop below will traverse the segment of the stream until an outlet
    # is reached.
    while True:
        # Following if loop is to execute the contents once the first node
        # in the segment is visited.
        if not switch_i:
            j = r[j]
            if j not in core_nodes:
                break
        #print 'wh_j = ',j
        # If the node being visited hasn't been visited before,
        # this if statement will executed.
        if flow_accum[j] == 0.:
            a += 1.
            alr_counted.append(j)
            stream_buffer.append(vic_ids[j])
        # Update number of upstream nodes.
        flow_accum[j] += a
        # If the node is being visited for the first time, the dictionary
        # 'vic_upstr' is being updated.
        if j in vic_upstr.keys():
            #print 'j=',j
            #print 'before ',vic_upstr[j]
            vic_upstr[j] += copy.copy(stream_buffer)
            #print 'after ',vic_upstr[j]
        # If the node has been already visited, then the upstream segment 
        # that was not accounted for in the main stem, would be added to
        # all the downstream nodes, one by one, until the outlet is reached.
        else:
            #print 'j=',j
            vic_upstr.update({j:copy.copy(stream_buffer)})
            #print 'after ',vic_upstr[j]
        # If the outlet is reached, the 'while' loop will be exited.
        if r[j]==j:
            break
        # This will be executed only for the first node of the stream segment.
        if switch_i:
            switch_i = False

OrderingTime = datetime.now() - startTime

# %% Algorithm to calculate coefficients of each upstream Vic ID
# vic_upstr = {1 : [21,22,23,24,25], 2: [23, 24, 25, 26, 25, 24]} #for testing
B = {}  # Holds unique upstream vic ids
C = {}  # Holds corresponding total numbers
coeff = {}  # Holds corresponding coefficients of contribution
for ke in vic_upstr.keys():
    cnt = Counter()
    for num in vic_upstr[ke]:
        cnt[num] += 1
    B.update({ke:cnt.keys()})
    buf = []
    for k in cnt.keys():
        buf.append(cnt[k])
    C.update({ke:buf})
    e = [s/float(sum(buf)) for s in buf]
    coeff.update({ke: e})

TotalTime = datetime.now() - startTime

# %% Plot or/and Save files
sim = '04Jan16_'
pickle.dump(B,open("dict_uniq_ids.p","wb"))    # Saving dict{30m id: uniq ids}
pickle.dump(coeff, open("dict_coeff.p","wb"))  # Saving dict{30m id: coeff}
np.save(sim + 'OrderingTime', OrderingTime)
np.save(sim + 'TotalTime', TotalTime)
