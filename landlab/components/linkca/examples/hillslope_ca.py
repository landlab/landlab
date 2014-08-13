#!/usr/env/python

"""
hillslope_ca.py

Example of a pair-based cellular automaton model, which simulates the evolution
of a hillslope by disturbance-driven soil creep.

GT, August 2014
"""

_DEBUG = False

import time
import numpy
from landlab import RasterModelGrid
from landlab.components.linkca.link_ca import LinkCellularAutomaton, Transition, CAPlotter
from landlab.components.fracture_grid.fracture_grid import make_frac_grid


def setup_transition_list():
    """
    Creates and returns a list of Transition() objects to represent state
    transitions for a weathering model.
    
    The states and transitions are as follows (note: X-X means "horizontal"
    pair, X/X means "vertical" pair with first item beneath second):
    
    Node states
    -----------
    0 air
    1 moving left
    2 moving right
    3 immobile
        
    Pair state      Transition to       Process     Rate
    ----------      -------------       -------     ----
    0 (0-0)         
    1 (0-1)
    2 (0-2)
    3 (0-3)         4 (1-0)             left ejection
    4 (1-0)
    5 (1-1)      
    6 (1-2)
    7 (1-3)
    8 (2-0)
    9 (2-1)
    10 (2-2)
    11 (2-3)        15 (3-3)            demobilization (right wall)
    12 (3-0)        2 (0-2)             right ejection
    13 (3-1)        15 (3-3)            demobilization (left wall)
    14 (3-2)
    15 (3-3)
    16 (0/0)      
    17 (0/1) 
    18 (0/2)             -           
    19 (0/3)        20 (1/0)            downward ejection, left
                    24 (2/0)            downward ejection, right
    20 (1/0)
    21 (1/1)          
    22 (1/2)
    23 (1/3)
    24 (2/0)
    25 (2/1)
    26 (2/2)
    27 (2/3)
    28 (3/0)        17 (0/1)            upward ejection, left
                    18 (0/2)            upward ejection, right
    29 (3/1)
    30 (3/2)
    31 (3/3)
    
    """
    xn_list = []
    
    xn_list.append( Transition(3, 4, 1., 'left ejection') )
    xn_list.append( Transition(12, 2, 1., 'right ejection') )
    xn_list.append( Transition(19, 20, 1., 'downward ejection, left') )
    xn_list.append( Transition(19, 24, 1., 'downward ejection, right') )
    xn_list.append( Transition(28, 17, 1., 'upward ejection, left') )
    xn_list.append( Transition(28, 18, 1., 'upward ejection, right') )
    xn_list.append( Transition(11, 15, 10., 'demobilization (right wall)') )
    xn_list.append( Transition(13, 15, 10., 'demobilization (left wall)') )
    # still to add: motion l, r, u, d; demob for vertical pairs
        
    if _DEBUG:
        print
        print 'setup_transition_list(): list has',len(xn_list),'transitions:'
        for t in xn_list:
            print '  From state',t.from_state,'to state',t.to_state,'at rate',t.rate,'called',t.name
        
    return xn_list
    
    
def main():
    
    # INITIALIZE

    # User-defined parameters
    nr = 128
    nc = 128
    plot_interval = 0.25
    run_duration = 8.0
    report_interval = 5.0  # report interval, in real-time seconds
    
    # Initialize real time
    current_real_time = time.time()
    next_report = current_real_time + report_interval

    # Create grid and set up boundaries
    mg = RasterModelGrid(nr, nc, 1.0)
    
    # Transition data here represent a body of fractured rock, with rock 
    # represented by nodes with state 0, and saprolite (weathered rock)
    # represented by nodes with state 1. Node pairs (links) with 0-1 or 1-0
    # can undergo a transition to 1-1, representing chemical weathering of the
    # rock.
    ns_dict = { 0 : 'air', 1 : 'mobile left', 2 : 'mobile right', 3 : 'immobile' }
    xn_list = setup_transition_list()

    # Create the node-state map and attach it to the grid
    node_state_grid = mg.add_zeros('node', 'node_state_map', dtype=int)
    (lower_half,) = numpy.where(mg.node_y<nr/2)
    node_state_grid[lower_half] = 3
    
    # Create the CA model
    ca = LinkCellularAutomaton(mg, ns_dict, xn_list, node_state_grid)
    
    # Debug output if needed    
    if _DEBUG:
        n = ca.grid.number_of_nodes
        for r in range(ca.grid.number_of_node_rows):
            for c in range(ca.grid.number_of_node_columns):
                n -= 1
                print '{0:.0f}'.format(ca.node_state[n]),
            print

    ca_plotter = CAPlotter(ca)
    
    # RUN
    current_time = 0.0
    while current_time < run_duration:
        
        # Once in a while, print out simulation and real time to let the user
        # know that the sim is running ok
        current_real_time = time.time()
        if current_real_time >= next_report:
            print 'Current sim time',current_time,'(',100*current_time/run_duration,'%)'
            next_report = current_real_time + report_interval
        
        # Run the model forward in time until the next output step
        ca.run(current_time+plot_interval, ca.node_state, 
               plot_each_transition=False) #, plotter=ca_plotter)
        current_time += plot_interval
        
        # Plot the current grid
        ca_plotter.update_plot()

        # for debugging        
        if _DEBUG:
            n = ca.grid.number_of_nodes
            for r in range(ca.grid.number_of_node_rows):
                for c in range(ca.grid.number_of_node_columns):
                    n -= 1
                    print '{0:.0f}'.format(ca.node_state[n]),
                print
        
        
    # FINALIZE
    
    # Plot
    ca_plotter.finalize()
        

if __name__ == "__main__":
    main()
