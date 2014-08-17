#!/usr/env/python

"""
hillslope_ca.py

Example of a pair-based cellular automaton model, which simulates the evolution
of a hillslope by disturbance-driven soil creep.
This version of the code incorporates all three of rock, regolith, and air.
To conduct an experiment just with regolith, no rock, set 
node_state_grid[lower_half] to 3, not 4.

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
    1 (0-1)         5 (1-0)             leftward motion
    2 (0-2)
    3 (0-3)         5 (1-0)             left ejection
    4 (0-4)         3 (0-3)             rock to regolith rightwards (bare surface)
    5 (1-0)
    6 (1-1)      
    7 (1-2)
    8 (1-3)
    9 (1-4)         8 (1-3)             rock to regolith rightwards (covered, mobile surface)
    10 (2-0)        2 (0-2)             rightward motion
    11 (2-1)
    12 (2-2)
    13 (2-3)        18 (3-3)            demobilization (right wall)
    14 (2-4)        13 (2-3)            rock to regolith rightwards (covered, mobile surface)
    15 (3-0)        2 (0-2)             right ejection
    16 (3-1)        18 (3-3)            demobilization (left wall)
    17 (3-2)
    18 (3-3)
    19 (3-4)        18 (3-3)            rock to regolith rightwards (covered surface)
    20 (4-0)        15 (3-0)            rock to regolith leftwards (bare surface)
    21 (4-1)        16 (3-1)            rock to regolith leftwards (covered, mobile surface)
    22 (4-2)        17 (3-2)            rock to regolith leftwards (covered, mobile surface)
    23 (4-3)        18 (3-3)            rock to regolith leftwards (covered surface)
    24 (4-4)
    25 (0/0)      
    26 (0/1)        30 (1/0)            downward motion
    27 (0/2)        35 (2/0)            downward motion         
    28 (0/3)        30 (1/0)            downward ejection, left
                    35 (2/0)            downward ejection, right
    29 (0/4)        28 (0/3)            rock to regolith upwards (bare surface)
    30 (1/0)        26 (0/1)            upward motion
    31 (1/1)          
    32 (1/2)
    33 (1/3)
    34 (1/4)        33 (1/3)            rock to regolith upwards (covered, mobile surface)
    35 (2/0)        27 (0/2)            upward motion
    36 (2/1)
    37 (2/2)
    38 (2/3)
    39 (2/4)        38 (2/3)            rock to regolith upwards (covered, mobile surface)
    40 (3/0)        26 (0/1)            upward ejection, left
                    27 (0/2)            upward ejection, right
    41 (3/1)        43 (3/3)            demobilization (friction)
    42 (3/2)        43 (3/3)            demobilization (friction)
    43 (3/3)
    44 (3/4)        43 (3/3)            rock to regolith upwards (covered surface)
    45 (4/0)        40 (3/0)            rock to regolith downwards (bare surface)
    46 (4/1)        41 (3/1)            rock to regolith downwards (covered, mobile surface)
    47 (4/2)        42 (3/2)            rock to regolith downwards (covered, mobile surface)
    48 (4/3)        43 (3/3)            rock to regolith downwards (covered surface)
    49 (4/4)
    
    """
    xn_list = []
    
    xn_list.append( Transition(3, 5, 0.01, 'left ejection') ) #0.01
    xn_list.append( Transition(15, 2, 0.01, 'right ejection') ) #0.01
    xn_list.append( Transition(28, 30, 0.05, 'downward ejection, left') ) #0.01
    xn_list.append( Transition(28, 35, 0.05, 'downward ejection, right') ) #0.01
    xn_list.append( Transition(40, 26, 0.005, 'upward ejection, left') ) #0.01
    xn_list.append( Transition(40, 27, 0.005, 'upward ejection, right') ) #0.01
    xn_list.append( Transition(13, 18, 100., 'demobilization (right wall)') )
    xn_list.append( Transition(16, 18, 100., 'demobilization (left wall)') )
    xn_list.append( Transition(41, 43, 0.1, 'demobilization (friction)') )
    xn_list.append( Transition(42, 43, 0.1, 'demobilization (friction)') )
    xn_list.append( Transition(1, 5, 1.0, 'leftward motion') )
    xn_list.append( Transition(10, 2, 1.0, 'rightward motion') )
    xn_list.append( Transition(30, 26, 0.1, 'upward motion') )
    xn_list.append( Transition(35, 27, 0.1, 'upward motion') )
    xn_list.append( Transition(26, 30, 10.0, 'downward motion') )
    xn_list.append( Transition(27, 35, 10.0, 'downward motion') )
    #DEJH adds rock-regolith properties:
    P_weath_bare = 0.002
    P_weath_mobile = 0.0005
    P_weath_cover = 0.0005
    P_weath_down = 0.05 #to simulate "undermining"
    xn_list.append( Transition(4, 3, P_weath_bare, 'weathering front, right, bare') )
    xn_list.append( Transition(9, 8, P_weath_mobile, 'weathering front, right, covered mobile') )
    xn_list.append( Transition(14, 13, P_weath_mobile, 'weathering front, right, covered mobile') )
    xn_list.append( Transition(19, 18, P_weath_cover, 'weathering front, right, covered') )
    xn_list.append( Transition(20, 15, P_weath_bare, 'weathering front, left, bare') )
    xn_list.append( Transition(21, 16, P_weath_mobile, 'weathering front, left, covered mobile') )
    xn_list.append( Transition(22, 17, P_weath_mobile, 'weathering front, left, covered mobile') )
    xn_list.append( Transition(23, 18, P_weath_cover, 'weathering front, left, covered') )
    xn_list.append( Transition(29, 28, P_weath_down, 'weathering front, up, bare') ) #here's our "undermined" transition
    xn_list.append( Transition(34, 33, P_weath_mobile, 'weathering front, up, covered mobile') )
    xn_list.append( Transition(39, 38, P_weath_mobile, 'weathering front, up, covered mobile') )
    xn_list.append( Transition(44, 43, P_weath_cover, 'weathering front, up, covered') )
    xn_list.append( Transition(45, 40, P_weath_bare, 'weathering front, down, bare') )
    xn_list.append( Transition(46, 41, P_weath_mobile, 'weathering front, down, covered mobile') )
    xn_list.append( Transition(47, 42, P_weath_mobile, 'weathering front, down, covered mobile') )
    xn_list.append( Transition(48, 43, P_weath_cover, 'weathering front, down, covered') )
        
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
    plot_interval = 5.0
    run_duration = 5000.0
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
    #ns_dict = { 0 : 'air', 1 : 'mobile left', 2 : 'mobile right', 3 : 'immobile' }
    ns_dict = { 0 : 'air', 1 : 'mobile left', 2 : 'mobile right', 3 : 'immobile', 4 : 'rock' }
    xn_list = setup_transition_list()

    # Create the node-state map and attach it to the grid
    node_state_grid = mg.add_zeros('node', 'node_state_map', dtype=int)
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
