#!/usr/env/python

"""
hillslope_ca.py

Example of a pair-based cellular automaton model, which simulates the evolution
of a hillslope by disturbance-driven soil creep.

GT, August 2014
"""
from six import print_ as print

_DEBUG = False

import time
import numpy
from landlab import RasterModelGrid
from landlab.components.linkca.link_ca import LinkCellularAutomaton, Transition, CAPlotter
#from landlab.components.fracture_grid.fracture_grid import make_frac_grid


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
    
    xn_list.append( Transition(3, 4, 2., 'left ejection') )
    xn_list.append( Transition(12, 2, 2., 'right ejection') )
    xn_list.append( Transition(19, 20, 2.e8, 'downward ejection, left') )
    xn_list.append( Transition(19, 24, 2.e8, 'downward ejection, right') )
    xn_list.append( Transition(28, 17, 1., 'upward ejection, left') )
    xn_list.append( Transition(28, 18, 1., 'upward ejection, right') )
    xn_list.append( Transition(11, 15, 3.0e7, 'demobilization (right wall)') )
    xn_list.append( Transition(13, 15, 3.0e7, 'demobilization (left wall)') )
    xn_list.append( Transition(29, 31, 2.0e6, 'demobilization (friction)') )
    xn_list.append( Transition(30, 31, 2.0e6, 'demobilization (friction)') )
    xn_list.append( Transition(1, 4, 3.0e8, 'leftward motion') )
    xn_list.append( Transition(8, 2, 3.0e8, 'rightward motion') )
    xn_list.append( Transition(20, 17, 2.0e6, 'upward motion') )
    xn_list.append( Transition(24, 18, 2.0e6, 'upward motion') )
    xn_list.append( Transition(18, 24, 2.0e8, 'downward motion') )
    xn_list.append( Transition(17, 20, 2.0e8, 'downward motion') )
        
    if _DEBUG:
        print()
        print('setup_transition_list(): list has',len(xn_list),'transitions:')
        for t in xn_list:
            print('  From state',t.from_state,'to state',t.to_state,'at rate',t.rate,'called',t.name)
        
    return xn_list
    
    
def extract_hillslope_profile(node_matrix):
    """
    Extracts the hillslope profile by finding, for each column, the highest row
    that contains a non-air node.
    """
    ncols = numpy.size(node_matrix, 1)
    z = numpy.zeros(ncols)
    for col in range(ncols):
        dirt = numpy.where(node_matrix[:,col]!=0)[0]
        if len(dirt)>0:
            z[col] = numpy.amax(dirt)
    return z
    
    
def main():
    """
    The basic algorithm for a LinkCA model is as follows:
        
        Initialize
            Set grid size, run duration, and plotting/output parameters
            Set up node states
            Set up link transitions
            Create and initialize the starting grid
            Create the CA model
            Create a CA plotter
            Set up (custom) boundary-control parameters
            Set up any other custom things
            
        Run
            Optionally report output to user in real time
            Figure out when the next interruption will be
            Run until that interruption
            If it's time to plot, do so
            If it's time to control boundaries, do so
            If it's time for file output, do so
            
        Finalize
            Do the last plot, if applicable
            Do any custom post-processing
    """
    
    # INITIALIZE

    # User-defined parameters: general
    nr = 129
    nc = 129
    opt_screen_display = True
    plot_interval = 1
    opt_file_output = False
    file_output_interval = 1
    run_duration = 127
    report_interval = 5.0  # report interval, in real-time seconds
    
    # User-defined parameters: boundary control
    boundary_update_interval = 0.1  # interval for baselevel lowering
    initial_hill_height = nr-3  # must be < nr-1
    
    # Initialize real time and intervals for reporting, plotting, and file output
    current_real_time = time.time()
    next_report = current_real_time + report_interval
    never = run_duration+1  # anything that scheduled for later than run_duration never happens!
    if opt_screen_display:
        next_plot = plot_interval
    else:
        next_plot = never
    if opt_file_output:
        next_file_write = file_output_interval
    else:
        next_file_write = never    
    
    # Custom: initialize information for baselevel lowering
    next_boundary_update = boundary_update_interval
    left_edge_alinks = (nc-2)*(nr-1)+(nc-1)*numpy.arange(nr-2)
    right_edge_alinks = left_edge_alinks + (nc-2)
    bl_height = initial_hill_height
    bl_left_alink = left_edge_alinks[bl_height-1]
    bl_right_alink = right_edge_alinks[bl_height-1]

    # Create grid and set up boundaries
    mg = RasterModelGrid(nr, nc, 1.0)
    
    # Set up node states and pair transitions
    ns_dict = { 0 : 'air', 1 : 'mobile left', 2 : 'mobile right', 3 : 'immobile' }
    xn_list = setup_transition_list()

    # Create the node-state map and attach it to the grid
    node_state_grid = mg.add_zeros('node', 'node_state', dtype=int)
    
    # Custom: set initial values in the node-state grid
    (initial_hill,) = numpy.where(mg.node_y<=initial_hill_height)
    node_state_grid[initial_hill] = 3
    
    # Custom: Set the left and right boundary conditions
    left_side_ids = mg.left_edge_node_ids()
    right_side_ids = mg.right_edge_node_ids()

    # Create the CA model
    ca = LinkCellularAutomaton(mg, ns_dict, xn_list, node_state_grid)
    
    # Debug output if needed    
    if _DEBUG:
        n = ca.grid.number_of_nodes
        for r in range(ca.grid.number_of_node_rows):
            for c in range(ca.grid.number_of_node_columns):
                n -= 1
                print('{0:.0f}'.format(ca.node_state[n]), end=' ')
            print()

    # Create a CAPlotter object to handle screen display
    if opt_screen_display:
        ca_plotter = CAPlotter(ca)
    
    
    # RUN
    current_time = 0.0
    while current_time < run_duration:
        
        # Once in a while, print out simulation and real time to let the user
        # know that the sim is running ok
        current_real_time = time.time()
        if current_real_time >= next_report:
            print('Current sim time',current_time,'(',100*current_time/run_duration,'%)')
            next_report = current_real_time + report_interval
            
        # Find out when the next interruption (for plotting, boundary control,
        # file I/O, etc.) will be
        next_interrupt = min(next_plot, next_file_write, next_boundary_update, run_duration)
        
        # Run the model forward in time until the next output step
        ca.run(next_interrupt, ca.node_state, 
               plot_each_transition=False) #, plotter=ca_plotter)
               
        # Update the current time
        current_time = next_interrupt
       
        # Plot the current grid
        if current_time >= next_plot:
            ca_plotter.update_plot()
            next_plot += plot_interval
            
        # Write the current grid to file
        if current_time >= next_file_write:
            print('FILE OUTPUT NOT YET IMPLEMENTED')
            next_file_write += file_output_interval
            
        # Custom boundary control: baselevel lowering
        if current_time >= next_boundary_update:
            #print 'doing bl drop at time', current_time
            next_boundary_update += boundary_update_interval
            if bl_height>0:
                ca.node_state[left_side_ids[bl_height]] = 0  # turn the former baselevel to air
                ca.node_state[right_side_ids[bl_height]] = 0  # turn the former baselevel to air
                ca.update_link_state(bl_left_alink, 0, current_time)
                ca.update_link_state(bl_right_alink, 0, current_time)
                bl_height -= 1
                bl_left_alink = left_edge_alinks[bl_height-1]
                bl_right_alink = right_edge_alinks[bl_height-1]
            
        # for debugging        
        if _DEBUG:
            n = ca.grid.number_of_nodes
            for r in range(ca.grid.number_of_node_rows):
                for c in range(ca.grid.number_of_node_columns):
                    n -= 1
                    print('{0:.0f}'.format(ca.node_state[n]), end=' ')
                print()
        
        
    # FINALIZE
    
    # Plot
    if opt_screen_display:
        ca_plotter.finalize()
        
    # Custom post-processing
    z = extract_hillslope_profile(mg.node_vector_to_raster(ca.node_state))
    #numpy.savetxt('h0822-01-u125.txt',z)
    print('PEAK ELEV = ',numpy.amax(z))
    

if __name__ == "__main__":
    #import profile
    #profile.run('main()')
    main()
