# -*- coding: utf-8 -*-
"""
test_celllab_cts.py: Unit-test function for CellLabCTS.

Created on Thu Jul  9 08:20:06 2015

@author: gtucker
"""

from landlab import RasterModelGrid
from landlab.components.cellular_automata.celllab_cts import Transition
from landlab.components.cellular_automata.raster_cts import RasterCTS
from heapq import heappush
from heapq import heappop


def callback_function(ca, node1, node2, time_now):
    """
    This function is passed as an argument to a transition event, and then
    called automatically by CellLabCTSModel's do_event() method.
    """
    elapsed_time = time_now - ca.last_update_time[node1]
    if ca.node_state[node1]==1:
        ca.prop_data[ca.propid[node1]]+=100*elapsed_time
    elapsed_time = time_now - ca.last_update_time[node2]
    if ca.node_state[node2]==1:
        ca.prop_data[ca.propid[node2]]+=100*elapsed_time
    

def test_raster_cts():
    """
    Tests instantiation of a RasterCTS and implementation of one transition,
    with a callback function.
    """
    
    # Set up a small grid with no events scheduled
    mg = RasterModelGrid(4, 4, 1.0)
    mg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    node_state_grid = mg.add_ones('node', 'node_state_map', dtype=int)
    node_state_grid[6] = 0
    ns_dict = { 0 : 'black', 1 : 'white' }
    xn_list = []
    xn_list.append( Transition((1,0,0), (0,1,0), 0.1, '', True, callback_function))
    pd = mg.add_zeros('node', 'property_data', dtype=int)
    pd[5] = 50
    ca = RasterCTS(mg, ns_dict, xn_list, node_state_grid, prop_data=pd)

    # Test the data structures
    assert (ca.xn_to.size==4), 'wrong size for xn_to'
    assert (ca.xn_to.shape==(4, 1)), 'wrong size for xn_to'
    assert (ca.xn_to[2][0]==1), 'wrong value in xn_to'
    assert (len(ca.event_queue)==1), 'event queue has wrong size'
    assert (ca.num_link_states==4), 'wrong number of link states'
    assert (ca.prop_data[5]==50), 'error in property data'
    assert (ca.xn_rate[2][0]==0.1), 'error in transition rate array'
    assert (ca.node_active_links[1][6]==16), 'error in active link array'
    assert (ca.num_node_states==2), 'error in num_node_states'
    assert (ca.link_orientation[-1]==0), 'error in link orientation array'
    assert (ca.link_state_dict[(1, 0, 0)]==2), 'error in link state dict'
    assert (ca.n_xn[2]==1), 'error in n_xn'
    assert (ca.cell_pair[1]==(0, 1, 0)), 'error in cell_pair list'
    
    # Manipulate the data in the event queue for testing:
    
    # pop the scheduled event off the queue
    ev = heappop(ca.event_queue)
    assert (ca.event_queue==[]), 'event queue should now be empty but is not'
    
    # engineer an event
    ev.time = 1.0
    ev.link = 16
    ev.xn_to = 1
    ev.propswap = True
    ev.prop_update_fn = callback_function
    ca.next_update[16] = 1.0
    
    # push it onto the event queue
    heappush(ca.event_queue, ev)
    
    # run the CA
    ca.run(2.0)
    
    # some more tests. 
    # Is current time advancing correctly? (should only go to 1.0, not 2.0)
    # Did the two nodes (5 and 6) correctly exchange states?
    # Did the property ID and data arrays get updated? Note that the "propswap"
    # should switch propids between nodes 5 and 6, and the callback function
    # should increase the value of prop_data in the "swap" node from 50 to 150.
    assert (ca.current_time==1.0), 'current time incorrect'
    assert (ca.node_state[5]==0), 'error in node state 5'
    assert (ca.node_state[6]==1), 'error in node state 6'
    assert (ca.prop_data[ca.propid[6]]==150), 'error in prop swap'
    

if __name__=='__main__':
    test_raster_cts()
    