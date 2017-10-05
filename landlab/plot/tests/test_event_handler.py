#! /usr/bin/env python
"""
Created on Thu Jun 29 11:02:28 2017

@author: njlyons
"""

from landlab import imshow_grid, RasterModelGrid
from landlab.plot.event_handler import query_grid_on_button_press
from matplotlib.pyplot import gcf
from matplotlib.backend_bases import Event
from numpy.testing import assert_equal


def test_query_grid_on_button_press():
    
    rmg = RasterModelGrid((5, 5))
    imshow_grid(rmg, rmg.nodes, cmap='RdYlBu')
    
    # Programmatically create an event near the grid center.
    event = Event('simulated_event', gcf().canvas)
    event.xdata = int(rmg.number_of_node_columns * 0.5)
    event.ydata = int(rmg.number_of_node_rows * 0.5)
    
    results = query_grid_on_button_press(event, rmg)
    x_coord = results['grid location']['x_coord']
    y_coord = results['grid location']['x_coord']
    
    msg = 'Items: Simulated matplotlib event and query results.'
    assert_equal(x_coord, event.xdata, msg)
    assert_equal(y_coord, event.ydata, msg)
    
    msg = 'Items: Node ID and grid coordinates of simulated matplotlib event.'
    node = rmg.grid_coords_to_node_id(event.xdata, event.ydata)
    assert_equal(results['node']['ID'], node, msg)
