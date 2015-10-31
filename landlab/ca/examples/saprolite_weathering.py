#!/usr/env/python

"""
saprolite_weathering.py

Example of a continuous-time, stochastic, pair-based cellular automaton model,
which simulates weathering of rock into saprolite.

GT, August 2014 (adapted to new Landlab cellular automata framework Sep 2014)
"""
from __future__ import print_function

_DEBUG = False

import time
from landlab import RasterModelGrid
from landlab.ca.landlab_ca import Transition, CAPlotter
from landlab.ca.raster_lca import RasterLCA
from landlab.components.fracture_grid.fracture_grid import make_frac_grid


def setup_transition_list():
    """
    Creates and returns a list of Transition() objects to represent state
    transitions for a weathering model.

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

    Pair state      Transition to       Process     Rate
    ==========      =============       =======     ====
    0 (0-0)         (none)              -           -
    1 (0-1)         3 (1-1)             weathering  1.0
    2 (1-0)         3 (1-1)             weathering  1.0
    3 (1-1)         (none)              -           -

    """
    xn_list = []

    xn_list.append( Transition(1, 3, 1., 'weathering') ) # rock-sap to sap-sap
    xn_list.append( Transition(2, 3, 1., 'weathering') ) # sap-rock to sap-sap

    if _DEBUG:
        print()
        print('setup_transition_list(): list has',len(xn_list),'transitions:')
        for t in xn_list:
            print('  From state',t.from_state,'to state',t.to_state,'at rate',t.rate,'called',t.name)

    return xn_list


def main():

    # INITIALIZE

    # User-defined parameters
    nr = 128
    nc = 128
    fracture_spacing = 10  # fracture spacing, cell widths
    plot_interval = 0.25
    run_duration = 4.0
    report_interval = 5.0  # report interval, in real-time seconds

    # Remember the clock time, and calculate when we next want to report
    # progress.
    current_real_time = time.time()
    next_report = current_real_time + report_interval

    # Create grid
    mg = RasterModelGrid(nr, nc, 1.0)

    # Set up the states and pair transitions.
    # Transition data here represent a body of fractured rock, with rock
    # represented by nodes with state 0, and saprolite (weathered rock)
    # represented by nodes with state 1. Node pairs (links) with 0-1 or 1-0
    # can undergo a transition to 1-1, representing chemical weathering of the
    # rock.
    ns_dict = { 0 : 'rock', 1 : 'saprolite' }
    xn_list = setup_transition_list()

    # Create the node-state array and attach it to the grid
    node_state_grid = mg.add_zeros('node', 'node_state_map', dtype=int)

    # Initialize the node-state array as a "fracture grid" in which randomly
    # oriented fractures are represented as lines of saprolite embedded in
    # bedrock.
    node_state_grid[:] = make_frac_grid(fracture_spacing, model_grid=mg)

    # Create the CA model
    ca = RasterLCA(mg, ns_dict, xn_list, node_state_grid)

    # Debug output if needed
    if _DEBUG:
        n = ca.grid.number_of_nodes
        for r in range(ca.grid.number_of_node_rows):
            for c in range(ca.grid.number_of_node_columns):
                n -= 1
                print('{0:.0f}'.format(ca.node_state[n]), end=' ')
            print()

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
            print('Current sim time',current_time,'(',100*current_time/run_duration,'%)')
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
                    print('{0:.0f}'.format(ca.node_state[n]), end=' ')
                print()


    # FINALIZE

    # Plot
    ca_plotter.finalize()


if __name__ == "__main__":
    main()
