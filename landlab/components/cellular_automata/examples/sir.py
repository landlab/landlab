#!/usr/env/python
"""
sir.py: example of a Susceptible-Infectious-Recovered epidemiological 
cellular automaton model implemented on a hexagonal grid using stochastic
pair-transition rules.

GT Sep 2014
"""

_DEBUG = False

import time
from landlab import HexModelGrid
from numpy import where, logical_and
from landlab.components.cellular_automata.landlab_ca import Transition, CAPlotter
from landlab.components.cellular_automata.hex_lca import HexLCA


def setup_transition_list():
    """
    Creates and returns a list of Transition() objects to represent state
    transitions for the SIR model.
    
    Parameters
    ----------
    (none)
    
    Returns
    -------
    xn_list : list of Transition objects
        List of objects that encode information about the link-state transitions.
    
    Notes
    -----
    The states and transitions are as follows:

    Pair state      Transition to       Process
    ==========      =============       =======
    0 (0-0)         (none)              -       
    1 (0-1)         4 (1-1)             infection
                    2 (0-2)             recovery
    2 (0-2)         (none)              -
    3 (1-0)         4 (1-1)             infection
                    6 (2-0)             recovery
    4 (1-1)         5 (1-2)             recovery
                    6 (2-1)             recovery
    5 (1-2)         8 (2-2)             recovery
    6 (2-0)         (none)              -
    7 (2-1)         8 (2-2)             recovery
    8 (2-2)         (none)              -
    
    """
    xn_list = []
    
    xn_list.append( Transition((0,1), (1,1), 10., 'infection') )
    xn_list.append( Transition((0,1), (0,2), 1., 'recovery') )
    xn_list.append( Transition((1,0), (1,1), 10., 'infection') )
    xn_list.append( Transition((1,0), (2,0), 1., 'recovery') )
    xn_list.append( Transition((1,1), (1,2), 1., 'recovery') )
    xn_list.append( Transition((1,1), (2,1), 1., 'recovery') )
    xn_list.append( Transition((1,2), (2,2), 1., 'recovery') )
    xn_list.append( Transition((2,1), (2,2), 1., 'recovery') )
        
    if _DEBUG:
        print
        print 'setup_transition_list(): list has',len(xn_list),'transitions:'
        for t in xn_list:
            print '  From state',t.from_state,'to state',t.to_state,'at rate',t.rate,'called',t.name
        
    return xn_list
    
    
def main():
    
    # INITIALIZE
    
    # User-defined parameters
    nr = 10
    nc = 10
    plot_interval = 0.25
    run_duration = 4.0
    report_interval = 5.0  # report interval, in real-time seconds
    
    # Remember the clock time, and calculate when we next want to report
    # progress.
    current_real_time = time.time()
    next_report = current_real_time + report_interval

    # Create a grid
    hmg = HexModelGrid(nr, nc, 1.0)
    
    # Set up the states and pair transitions.
    # Transition data here represent the disease status of a population.
    ns_dict = { 0 : 'susceptible', 1 : 'infectious', 2: 'recovered' }
    xn_list = setup_transition_list()

    # Create data and initialize values
    node_state_grid = hmg.add_zeros('node', 'node_state_grid')
    is_middle_rows = logical_and(hmg.node_y>=0.4*nr, hmg.node_y<=0.6*nr)
    is_middle_cols = logical_and(hmg.node_x>=0.4*nc, hmg.node_x<=0.6*nc)
    print is_middle_rows
    print is_middle_cols
    middle_area = where(logical_and(is_middle_rows, is_middle_cols))[0]
    node_state_grid[middle_area] = 1
    print node_state_grid
    
    # Create the CA model
    ca = HexLCA(hmg, ns_dict, xn_list, node_state_grid)
    
    # Create a CAPlotter object for handling screen display
    ca_plotter = CAPlotter(ca)
    
    # Plot the initial grid
    ca_plotter.update_plot()

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
               plot_each_transition=True, plotter=ca_plotter)
        current_time += plot_interval
        
        # Plot the current grid
        ca_plotter.update_plot()


    # FINALIZE

    # Plot
    ca_plotter.finalize()


if __name__=='__main__':
    main()
