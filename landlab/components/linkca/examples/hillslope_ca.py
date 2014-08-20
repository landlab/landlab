#!/usr/env/python

"""
hillslope_ca.py

Example of a pair-based cellular automaton model, which simulates the evolution
of a hillslope by disturbance-driven soil creep.

GT, August 2014
"""

_DEBUG = True

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
    1 (0-1)         4 (1-0)             leftward motion
    2 (0-2)
    3 (0-3)         4 (1-0)             left ejection
    4 (1-0)
    5 (1-1)      
    6 (1-2)
    7 (1-3)
    8 (2-0)         2 (0-2)             rightward motion
    9 (2-1)
    10 (2-2)
    11 (2-3)        15 (3-3)            demobilization (right wall)
    12 (3-0)        2 (0-2)             right ejection
    13 (3-1)        15 (3-3)            demobilization (left wall)
    14 (3-2)
    15 (3-3)
    16 (0/0)      
    17 (0/1)        20 (1/0)            downward motion
    18 (0/2)        24 (2/0)            downward motion         
    19 (0/3)        20 (1/0)            downward ejection, left
                    24 (2/0)            downward ejection, right
    20 (1/0)        17 (0/1)            upward motion
    21 (1/1)          
    22 (1/2)
    23 (1/3)
    24 (2/0)        18 (0/2)            upward motion
    25 (2/1)
    26 (2/2)
    27 (2/3)
    28 (3/0)        17 (0/1)            upward ejection, left
                    18 (0/2)            upward ejection, right
    29 (3/1)        31 (3/3)            demobilization (friction)
    30 (3/2)        31 (3/3)            demobilization (friction)
    31 (3/3)
    
    """
    xn_list = []
    
    xn_list.append( Transition(3, 4, 0.01, 'left ejection') )
    xn_list.append( Transition(12, 2, 0.01, 'right ejection') )
    xn_list.append( Transition(19, 20, 0.01, 'downward ejection, left') )
    xn_list.append( Transition(19, 24, 0.01, 'downward ejection, right') )
    xn_list.append( Transition(28, 17, 0.01, 'upward ejection, left') )
    xn_list.append( Transition(28, 18, 0.01, 'upward ejection, right') )
    xn_list.append( Transition(11, 15, 100., 'demobilization (right wall)') )
    xn_list.append( Transition(13, 15, 100., 'demobilization (left wall)') )
    xn_list.append( Transition(29, 31, 0.1, 'demobilization (friction)') )
    xn_list.append( Transition(30, 31, 0.1, 'demobilization (friction)') )
    xn_list.append( Transition(1, 4, 1.0, 'leftward motion') )
    xn_list.append( Transition(8, 2, 1.0, 'rightward motion') )
    xn_list.append( Transition(20, 17, 0.1, 'upward motion') )
    xn_list.append( Transition(24, 18, 0.1, 'upward motion') )
    xn_list.append( Transition(18, 24, 10.0, 'downward motion') )
    xn_list.append( Transition(17, 20, 10.0, 'downward motion') )
        
    if _DEBUG:
        print
        print 'setup_transition_list(): list has',len(xn_list),'transitions:'
        for t in xn_list:
            print '  From state',t.from_state,'to state',t.to_state,'at rate',t.rate,'called',t.name
        
    return xn_list
    
    
def main():
    
    # INITIALIZE

    # User-defined parameters
    nr = 9
    nc = 9
    plot_interval = 1.0
    run_duration = 100.0
    report_interval = 5.0  # report interval, in real-time seconds
    uplift_interval = 1
    
    # Initialize real time
    current_real_time = time.time()
    next_report = current_real_time + report_interval
    
    # Initialize next uplift time
    next_uplift = uplift_interval

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
    node_state_grid = mg.add_zeros('node', 'node_state', dtype=int)
    (lower_half,) = numpy.where(mg.node_y<nr/2)
    node_state_grid[lower_half] = 3
    
    # Set the left and right boundary conditions
    node_state_grid[mg.left_edge_node_ids()] = 0
    node_state_grid[mg.right_edge_node_ids()] = 0
    
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
        
        # Uplift
        if current_time >= next_uplift:
            print 'doing uplift at time', current_time
            next_uplift += uplift_interval
            ca.grid.roll_nodes_ud('node_state', 1, interior_only=True)
            #print 'after roll:',ca.node_state

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
