#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 10:08:51 2020

@author: gtucker
"""

import numpy as np
from landlab import RasterModelGrid
from landlab.components import LakeMapperBarnes, FlowAccumulator
from landlab.components import FlowDirectorSteepest
from collections import deque


mg = RasterModelGrid((5, 6), xy_spacing=2.)
for edge in ('left', 'top', 'bottom'):
    mg.status_at_node[mg.nodes_at_edge(edge)] = mg.BC_NODE_IS_CLOSED
z = mg.add_zeros("topographic__elevation", at="node", dtype=float)
z.reshape(mg.shape)[2, 1:-1] = [2., 1., 0.5, 1.5]
z.reshape(mg.shape)[1, 1:-1] = [2.1, 1.1, 0.6, 1.6]
z.reshape(mg.shape)[3, 1:-1] = [2.2, 1.2, 0.7, 1.7]
z_init = z.copy()
fa = FlowAccumulator(mg)
lmb = LakeMapperBarnes(mg, method='Steepest', surface=z_init, 
                       fill_flat=False,
                       redirect_flow_steepest_descent=False,
                       track_lakes=False)

lmb.run_one_step()
z_out = np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,
                   0. ,  2.1,  1.5,  1.5,  1.6,  0. ,
                   0. ,  2. ,  1.5,  1.5,  1.5,  0. ,
                   0. ,  2.2,  1.5,  1.5,  1.7,  0. ,
                   0. ,  0. ,  0. ,  0. ,  0. ,  0. ])
assert(np.allclose(z, z_out))
#        True
assert(not np.all(np.equal(z, z_out)))  # those 1.5's are actually a bit > 1.5
#         False
try:
    lmb.lake_map  # not created, as we aren't tracking
except ValueError:
    print('ValueError was raised: ' +
          'Enable tracking to access information about lakes')
#         ValueError was raised: Enable tracking to access information about lakes
assert(not lmb.was_there_overfill)  # everything fine with slope adding
#         False

fd = FlowDirectorSteepest(mg)
fa = FlowAccumulator(mg)  # routing will work fine now
fd.run_one_step()
fa.run_one_step()
assert(np.all(mg.at_node['flow__sink_flag'][mg.core_nodes] == 0))
#         True
drainage_area = np.array([  0.,   0.,   0.,   0.,   0.,   0.,
              0.,   4.,   8.,  12.,   4.,   4.,
              0.,   4.,  16.,  36.,  40.,  40.,
              0.,   4.,   8.,   4.,   4.,   4.,
              0.,   0.,   0.,   0.,   0.,   0.])
assert(np.allclose(mg.at_node['drainage_area'], drainage_area))
#         True

#         Test two pits:

z[:] = mg.node_x.max() - mg.node_x
z[23] = 1.3
z[15] = 0.3
z[10] = 1.3  # raise "guard" exit nodes
z[7] = 2.  # is a lake on its own, if D8
z[9] = 0.5
z[14] = 0.6  # [9, 14, 15] is a lake in both methods
z[16] = 1.2
z[22] = 0.9  # a non-contiguous lake node also draining to 16 if D8
z_init = z.copy()
fa = FlowAccumulator(mg)
lmb = LakeMapperBarnes(mg, method='D8', fill_flat=True,
        track_lakes=True)
lmb.run_one_step()  # note the D8 routing now
assert(lmb.lake_dict == {22: deque([15, 9, 14])})
#         True
assert(lmb.number_of_lakes==1)
#         1
try:
    lmb.lake_depths  # z was both surface and 'fill_surface'
except ValueError:
    print('ValueError was raised: ' +
          'surface and fill_surface must be different fields ' +
          'or arrays to enable the property fill_depth!')
#         ValueError was raised: surface and fill_surface must be different fields or arrays to enable the property fill_depth!

z[:] = z_init
lmb = LakeMapperBarnes(mg, method='Steepest',
        fill_flat=False, track_lakes=True)
lmb.run_one_step()  # compare to the method='D8' lakes, above...
assert(lmb.lake_dict == {8: deque([7]), 16: deque([15, 9, 14, 22])})
#         True
assert(lmb.number_of_lakes == 2)
#         2
assert(np.allclose(lmb.lake_areas, np.array([ 16.,  4.])))
#         True
try:
    lmb.run_one_step()  # z already filled, so...
except ValueError:
    print('ValueError was raised: ' +
           'Pit is overfilled due to creation of two outlets as ')
#    +
#          'the minimum gradient gets applied. Suppress this ' +
#          'Error with the ignore_overfill flag at component ' +
#          'instantiation.')
#         ValueError was raised: Pit is overfilled due to creation of two outlets as the minimum gradient gets applied. Suppress this Error with the ignore_overfill flag at component instantiation.

#         Suppress this behaviour with ignore_overfill:

z[:] = z_init
lmb = LakeMapperBarnes(mg, method='Steepest',
        fill_flat=False, track_lakes=True,
        ignore_overfill=True)
lmb.run_one_step()
assert(lmb.lake_dict == {8: deque([7]), 16: deque([15, 9, 14, 22])})
#         True
lmb.run_one_step()
assert(np.allclose(lmb.lake_areas, np.array([ 16.,  4.])))  # found them!
#         True

#         The component can redirect flow to account for the fills that have
#         been carried out (all necessary fields get updated):

z[:] = z_init
fd.run_one_step()
init_flowdirs = mg.at_node['flow__receiver_node'].copy()
fa.run_one_step()
init_areas = mg.at_node['drainage_area'].copy()
init_qw = mg.at_node['surface_water__discharge'].copy()
fa = FlowAccumulator(mg)
lmb = LakeMapperBarnes(mg, method='Steepest',
        fill_flat=False, track_lakes=True,
        redirect_flow_steepest_descent=False,
        ignore_overfill=True)
lmb.run_one_step()
assert(np.all(mg.at_node['flow__receiver_node'] == init_flowdirs))
#         True

lmb = LakeMapperBarnes(mg, method='Steepest',
        fill_flat=False, track_lakes=True,
        redirect_flow_steepest_descent=True,
        ignore_overfill=True)
lmb.run_one_step()
assert(not np.all(mg.at_node['flow__receiver_node'] == init_flowdirs))
#         False

#         However, note that unless the reaccumulate_flow argument is also
#         set, the 'drainage_area' and 'surface_water__discharge' fields
#         *won't* also get updated:

assert(np.all(mg.at_node['drainage_area'] == init_areas))
#         True
assert(np.all(mg.at_node['surface_water__discharge'] == init_qw))
#         True

fa = FlowAccumulator(mg)
lmb = LakeMapperBarnes(mg, method='Steepest',
        fill_flat=False, track_lakes=True,
        redirect_flow_steepest_descent=True,
        reaccumulate_flow=True,
        ignore_overfill=True)
lmb.run_one_step()
assert(not np.all(mg.at_node['drainage_area'] == init_areas))
#         False
assert(not np.all(mg.at_node['surface_water__discharge'] == init_qw))
#         False

#         Be sure to set both redirect_flow_steepest_descent and
#         reaccumulate_flow to True if you want to reaccumulate flow...

try:
    lmb = LakeMapperBarnes(mg, method='Steepest',
            fill_flat=False, track_lakes=True,
            redirect_flow_steepest_descent=False,
            reaccumulate_flow=True,
            ignore_overfill=True)
except ValueError:
    print('Oops!')
#         Oops!

#         The component is completely happy with irregular grids:

from landlab import HexModelGrid, FieldError
hmg = HexModelGrid((5, 4), spacing=2.)
z_hex = hmg.add_zeros("topographic__elevation", at="node")
z_hex[:] = hmg.node_x
z_hex[11] = -3.
z_hex[12] = -1.
z_hex_init = z_hex.copy()
print(z_hex)
#         array([   2.,   4.,   6.,   8.,
#                1.,   3.,   5.,   7.,   9.,
#             0.,   2.,   -3.,  -1.,   8.,  10.,
#                1.,   3.,   5.,   7.,   9.,
#                   2.,   4.,  6.,   8.])

#         As you can see, nodes 11 and 12 are now a pit. If they were to fill
#         they would fill to the level of 2, the lowest downstream value.

fa = FlowAccumulator(hmg)
lmb = LakeMapperBarnes(
       hmg,
       method='Steepest',
       fill_flat=True,
       track_lakes=False)
lmb.run_one_step()
assert(np.allclose(z_hex[10:13], 2.))
#         True

hmg = HexModelGrid((5, 4), spacing=2.0)
z_hex = hmg.add_zeros("topographic__elevation", at="node")
z_hex[:] = z_hex_init
try:
    lmb = LakeMapperBarnes(hmg, method='Steepest',
            fill_flat=False,
            surface=z_hex_init,
            redirect_flow_steepest_descent=True,
            track_lakes=True)
except FieldError:
    print("Oops2!")  # flowdir field must already exist!
#         Oops2!
fd = FlowDirectorSteepest(hmg)
fa = FlowAccumulator(hmg)
lmb = LakeMapperBarnes(hmg, method='Steepest',
        fill_flat=False, surface=z_hex_init,
        redirect_flow_steepest_descent=True,
        track_lakes=True)
fd.run_one_step()
lmb.run_one_step()
assert(np.allclose(z_hex[10:13], 2.))
#         True
assert(z_hex[11] > z_hex[10])
#         True
assert(z_hex[12] > z_hex[11])
#         True
assert(np.allclose(lmb.lake_depths[10:14], np.array([ 0.,  5.,  3.,  0.])))
#         True
np.testing.assert_array_almost_equal(
       lmb.lake_volumes,
       27.712,
       decimal=3)

#         Together, all this means that we can now run a topographic growth
#         model that permits flooding as it runs:

import numpy as np
from landlab import RasterModelGrid
from landlab.components import LakeMapperBarnes, FlowAccumulator
from landlab.components import FlowDirectorSteepest
from landlab.components import FastscapeEroder
mg = RasterModelGrid((6, 8))
for edge in ('right', 'top', 'bottom'):
    mg.status_at_node[mg.nodes_at_edge(edge)] = mg.BC_NODE_IS_CLOSED

#         Because it is what we want the FastscapeEroder to see and work on,
#         it's actually the water surface that needs to go in as
#         'topographic__elevation'. We'll also need to keep track of the bed
#         elevation though, since the LakeMapper will need it. We start them
#         equal (i.e., topo starts dry).

z_water = mg.add_zeros("topographic__elevation", at="node", dtype=float)
z_water[:] = mg.node_x
z_water[11] = 1.5
z_water[19] = 0.5
z_water[34] = 1.1
z_bed = mg.add_zeros("bedrock__elevation", at="node", dtype=float)
z_bed[:] = z_water  # topo starts dry

#         Let's just take a look:

assert(np.all(np.equal(
    np.round(z_water, 2),
    np.array([0. , 1. , 2. , 3. , 4. , 5. , 6. , 7. ,
              0. , 1. , 2. , 1.5, 4. , 5. , 6. , 7. ,
              0. , 1. , 2. , 0.5, 4. , 5. , 6. , 7. ,
              0. , 1. , 2. , 3. , 4. , 5. , 6. , 7. ,
              0. , 1. , 1.1, 3. , 4. , 5. , 6. , 7. ,
              0. , 1. , 2. , 3. , 4. , 5. , 6. , 7. ]))))
#         True

fd = FlowDirectorSteepest(mg)
fa = FlowAccumulator(mg)
lmb = LakeMapperBarnes(mg, method='D8', fill_flat=True,
        surface='bedrock__elevation',
        fill_surface='topographic__elevation',
        redirect_flow_steepest_descent=True,
        reaccumulate_flow=True,
        track_lakes=True)
sp = FastscapeEroder(mg, K_sp=1., m_sp=0., n_sp=1.)
fd.run_one_step()
fa.run_one_step()  # node 18 is draining into the pit...
assert(np.isclose(mg.at_node['topographic__steepest_slope'][18], 1.5))
#         True
assert(np.allclose(mg.at_node['drainage_area'],
            np.array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        2.,  2.,  1.,  4.,  3.,  2.,  1.,  0.,
        1.,  1.,  1., 13.,  3.,  2.,  1.,  0.,
        2.,  2.,  1.,  4.,  3.,  2.,  1.,  0.,
        6.,  6.,  5.,  4.,  3.,  2.,  1.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])))
#        True
lmb.run_one_step()  # now node 18 drains correctly, outward ->  ***THIS IS THE LINE ***
assert(np.isclose(mg.at_node['topographic__steepest_slope'][18], 1.))
#         True
assert(np.allclose(mg.at_node['drainage_area'],
            np.array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
                              13., 13., 12.,  4.,  3.,  2.,  1.,  0.,
        2.,  2.,  1.,  7.,  3.,  2.,  1.,  0.,
        2.,  2.,  1.,  1.,  3.,  2.,  1.,  0.,
        7.,  7.,  6.,  4.,  3.,  2.,  1.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])))
#         True
assert(np.all(np.equal(
    np.round(z_water, 2),
    np.array([0. , 1. , 2. , 3. , 4. , 5. , 6. , 7. ,
              0. , 1. , 2. , 2. , 4. , 5. , 6. , 7. ,
              0. , 1. , 2. , 2. , 4. , 5. , 6. , 7. ,
              0. , 1. , 2. , 3. , 4. , 5. , 6. , 7. ,
              0. , 1. , 1.1, 3. , 4. , 5. , 6. , 7. ,
              0. , 1. , 2. , 3. , 4. , 5. , 6. , 7. ]))))
#         True

sp.run_one_step(0.05)  # note m=0 to illustrate effect of slopes
assert(np.all(np.equal(
    np.round(z_water, 2),
    np.array([0.  , 1.  , 2.  , 3.  , 4.  , 5.  , 6.  , 7.  ,
              0.  , 0.95, 1.95, 2.  , 3.9 , 4.95, 5.95, 7.  ,
              0.  , 0.95, 1.95, 2.  , 3.9 , 4.95, 5.95, 7.  ,
              0.  , 0.95, 1.95, 2.93, 3.93, 4.95, 5.95, 7.  ,
              0.  , 0.95, 1.09, 2.91, 3.95, 4.95, 5.95, 7.  ,
              0.  , 1.  , 2.  , 3.  , 4.  , 5.  , 6.  , 7.  ]))))
#         True

#         If we want to keep this going honouring the depths of the lakes try
#         this next in your loop:

z_bed[:] = np.minimum(z_water, z_bed)
assert(np.all(np.equal(
    np.round(z_bed, 2),
    np.array([0.  , 1.  , 2.  , 3.  , 4.  , 5.  , 6.  , 7.  ,
              0.  , 0.95, 1.95, 1.5 , 3.9 , 4.95, 5.95, 7.  ,
              0.  , 0.95, 1.95, 0.5 , 3.9 , 4.95, 5.95, 7.  ,
              0.  , 0.95, 1.95, 2.93, 3.93, 4.95, 5.95, 7.  ,
              0.  , 0.95, 1.09, 2.91, 3.95, 4.95, 5.95, 7.  ,
              0.  , 1.  , 2.  , 3.  , 4.  , 5.  , 6.  , 7.  ]))))
#         True
fd.run_one_step()
fa.run_one_step()
lmb.run_one_step()

#         Lake node depths are now updated in lmb:

np.testing.assert_array_equal(np.round(
    [lmb.lake_depths[lake] for lake in lmb.lake_dict.values()], 2),
        np.array([[ 0.45,  1.45]]))

#         ...and the "topography" (i.e., water surface) at the flooded nodes
#         has lowered itself as the lip of the outlet was eroded in the last
#         step:

assert(np.all(np.equal(
    np.round(z_water, 2),
    np.array([0.  , 1.  , 2.  , 3.  , 4.  , 5.  , 6.  , 7.  ,
              0.  , 0.95, 1.95, 1.95 , 3.9 , 4.95, 5.95, 7.  ,
              0.  , 0.95, 1.95, 1.95 , 3.9 , 4.95, 5.95, 7.  ,
              0.  , 0.95, 1.95, 2.93, 3.93, 4.95, 5.95, 7.  ,
              0.  , 0.95, 1.09, 2.91, 3.95, 4.95, 5.95, 7.  ,
              0.  , 1.  , 2.  , 3.  , 4.  , 5.  , 6.  , 7.  ]))))
#         True

# sp.run_one_step(0.05)
