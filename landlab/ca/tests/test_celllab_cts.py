# -*- coding: utf-8 -*-
"""
test_celllab_cts.py: Unit-test function for CellLabCTS.

Created on Thu Jul  9 08:20:06 2015

@author: gtucker
"""

from nose.tools import assert_equal
from numpy.testing import assert_array_equal
from landlab import RasterModelGrid, HexModelGrid
from landlab.ca.celllab_cts import Transition, Event
from landlab.ca.raster_cts import RasterCTS
from landlab.ca.oriented_raster_cts import OrientedRasterCTS
from landlab.ca.hex_cts import HexCTS
from landlab.ca.oriented_hex_cts import OrientedHexCTS
from heapq import heappush
from heapq import heappop


def callback_function(ca, node1, node2, time_now):
    """
    This function is passed as an argument to a transition event, and then
    called automatically by CellLabCTSModel's do_event() method.
    """
    pass
#    elapsed_time = time_now - ca.last_update_time[node1]
#    if ca.node_state[node1]==1:
#        ca.prop_data[ca.propid[node1]]+=100*elapsed_time
#    elapsed_time = time_now - ca.last_update_time[node2]
#    if ca.node_state[node2]==1:
#        ca.prop_data[ca.propid[node2]]+=100*elapsed_time
#



def test_transition():
    """Test instantiation of Transition() object."""
    t = Transition((0, 0, 0), (1, 1, 0), 1.0, name='test',
                   swap_properties=False, prop_update_fn=None)
    assert_equal(t.from_state, (0,0,0))
    assert_equal(t.to_state, (1,1,0))
    assert_equal(t.rate, 1.0)
    assert_equal(t.name, 'test')
    assert_equal(t.swap_properties, False)
    assert_equal(t.prop_update_fn, None)


def test_event_init():
    """Test instantiation of Event() object"""
    e = Event(2.0, 3, 4, propswap=False, prop_update_fn=None)
    assert_equal(e.time, 2.0)
    assert_equal(e.link, 3)
    assert_equal(e.xn_to, 4)
    assert_equal(e.propswap, False)
    assert_equal(e.prop_update_fn, None)


def test_event_lt():
    """Test Event.__lt__()"""
    e1 = Event(2.0, 3, 4)
    e2 = Event(5.0, 2, 3)
    assert_equal(e1<e2, True)


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
    ca = RasterCTS(mg, ns_dict, xn_list, node_state_grid, prop_data=pd,
                   prop_reset_value=0)

    # Test the data structures
    assert (ca.xn_to.size==4), 'wrong size for xn_to'
    assert (ca.xn_to.shape==(4, 1)), 'wrong size for xn_to'
    assert (ca.xn_to[2][0]==1), 'wrong value in xn_to'
    assert (len(ca.event_queue)==1), 'event queue has wrong size'
    assert (ca.num_link_states==4), 'wrong number of link states'
    assert (ca.prop_data[5]==50), 'error in property data'
    assert (ca.xn_rate[2][0]==0.1), 'error in transition rate array'
    assert (ca._active_links_at_node[1][6]==8), 'error in active link array'
    assert (ca.num_node_states==2), 'error in num_node_states'
    assert (ca.link_orientation[-1]==0), 'error in link orientation array'
    assert (ca.link_state_dict[(1, 0, 0)]==2), 'error in link state dict'
    assert (ca.n_xn[2]==1), 'error in n_xn'
    assert (ca.node_pair[1]==(0, 1, 0)), 'error in cell_pair list'

    # Manipulate the data in the event queue for testing:

    # pop the scheduled event off the queue
    ev = heappop(ca.event_queue)
    assert (ca.event_queue==[]), 'event queue should now be empty but is not'

    # engineer an event
    ev.time = 1.0
    ev.link = 8
    ev.xn_to = 1
    ev.propswap = True
    ev.prop_update_fn = callback_function
    ca.next_update[8] = 1.0

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
    #assert (ca.prop_data[ca.propid[6]]==150), 'error in prop swap'


def test_oriented_raster_cts():
    """Tests instantiation of an OrientedRasterCTS() object"""
    mg = RasterModelGrid(3, 3, 1.0)
    nsd = {0 : 'oui', 1 : 'non'}
    xnlist = []
    xnlist.append(Transition((0,1,0), (1,1,0), 1.0, 'hopping'))
    nsg = mg.add_zeros('node', 'node_state_grid')
    orcts = OrientedRasterCTS(mg, nsd, xnlist, nsg)

    assert_equal(orcts.num_link_states, 8)
    #assert_array_equal(orcts.link_orientation, [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0])
    assert_array_equal(orcts.link_orientation, [0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0])


def test_hex_cts():
    """Tests instantiation of a HexCTS() object"""
    mg = HexModelGrid(3, 2, 1.0, orientation='vertical', reorient_links=True)
    nsd = {0 : 'zero', 1 : 'one'}
    xnlist = []
    xnlist.append(Transition((0,1,0), (1,1,0), 1.0, 'transitioning'))
    nsg = mg.add_zeros('node', 'node_state_grid')
    hcts = HexCTS(mg, nsd, xnlist, nsg)

    assert_equal(hcts.num_link_states, 4)
    assert_array_equal(hcts.link_orientation, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])


def test_oriented_hex_cts():
    """Tests instantiation of an OrientedHexCTS() object"""
    mg = HexModelGrid(3, 2, 1.0, orientation='vertical', reorient_links=True)
    nsd = {0 : 'zero', 1 : 'one'}
    xnlist = []
    xnlist.append(Transition((0,1,0), (1,1,0), 1.0, 'transitioning'))
    nsg = mg.add_zeros('node', 'node_state_grid')
    ohcts = OrientedHexCTS(mg, nsd, xnlist, nsg)
    
    assert_equal(ohcts.num_link_states, 12)
    #assert_array_equal(ohcts.link_orientation, [2, 1, 0, 0, 0, 2, 1, 0, 2, 1, 0])
    assert_array_equal(ohcts.link_orientation, [2, 0, 1, 0, 2, 0, 1, 0, 2, 0, 1])

if __name__ == '__main__':
    test_raster_cts()