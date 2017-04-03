
###########################################################################
## 26 Jan 2016 - SN & EI
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
# convert Arc flow_directions file to be represented in node ids
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

core_nodes = grid.core_nodes
core_elev = z[core_nodes]
sor_z = core_nodes[np.argsort(core_elev)[::-1]]
alr_counted = [] # To store nodes that have already been counted
r_core=r[core_nodes]

flow_accum = np.zeros(grid.number_of_nodes,dtype=int)
vic_upstr = {}
for i in sor_z:
    #print 'i= ', i
    if i in alr_counted:
        continue
    if r[i]==i:
        vic_upstr.update({i:[vic_ids[i]]})
        flow_accum[i] += 1.
        alr_counted.append(i)
        continue
    carry_ = []
    j = i
    switch_i = True
    a = 0.
    while True:
        if not switch_i:
            j = r[j]
            if j not in core_nodes:
                break
        #print 'wh_j = ',j
        if flow_accum[j] == 0.:
            a += 1.
            alr_counted.append(j)
            carry_.append(vic_ids[j])

        flow_accum[j] += a

        if j in vic_upstr.keys():
            #print 'j=',j
            #print 'before ',vic_upstr[j]
            vic_upstr[j] += copy.copy(carry_)
            #print 'after ',vic_upstr[j]
        else:
            #print 'j=',j
            vic_upstr.update({j:copy.copy(carry_)})
            #print 'after ',vic_upstr[j]

        if r[j]==j:
            break
        if switch_i:
            switch_i = False

OrderingTime = datetime.now() - startTime

## Algorithm to calculate coefficients of each upstream Vic ID
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

sim = '04Jan16_'
pickle.dump(B,open("dict_uniq_ids.p","wb"))    # Saving dict{30m id: uniq ids}
pickle.dump(coeff, open("dict_coeff.p","wb"))  # Saving dict{30m id: coeff}
np.save(sim + 'OrderingTime', OrderingTime)
np.save(sim + 'TotalTime', TotalTime)
