#!/usr/env/python

"""
diffusion_in_gravity.py

Example of a continuous-time, stochastic, pair-based cellular automaton model,
which simulates diffusion by random particle motion in a gravitational field.
The purpose of the example is to demonstrate the use of an OrientedRasterLCA.

GT, September 2014
"""
from __future__ import print_function

_DEBUG = False

import time
from numpy import where, bitwise_and
from landlab import RasterModelGrid
from landlab.ca.celllab_cts import Transition, CAPlotter
from landlab.ca.oriented_raster_cts import OrientedRasterCTS


def setup_transition_list():
    """
    Creates and returns a list of Transition() objects to represent state
    transitions for a biased random walk, in which the rate of downward
    motion is greater than the rate in the other three directions.

    Parameters
    ----------
    (none)

    Returns
    -------
    xn_list : list of Transition objects
        List of objects that encode information about the link-state transitions.

    Notes
    -----
    State 0 represents fluid and state 1 represents a particle (such as a
    sediment grain or dissolved heavy particle).

    The states and transitions are as follows:

    Pair state      Transition to       Process         Rate
    ==========      =============       =======         ====
    0 (0-0)         (none)              -               -
    1 (0-1)         2 (1-0)             left motion     1.0
    2 (1-0)         1 (0-1)             right motion    1.0
    3 (1-1)         (none)              -               -
    4 (0/0)         (none)              -               -
    5 (0/1)         6 (1/0)             down motion     1.1
    6 (1/0)         5 (0/1)             up motion       0.9
    7 (1/1)         (none)              -               -

    """
    xn_list = []

    xn_list.append( Transition((0,1,0), (1,0,0), 1., 'left motion') )
    xn_list.append( Transition((1,0,0), (0,1,0), 1., 'right motion') )
    xn_list.append( Transition((0,1,1), (1,0,1), 1.1, 'down motion') )
    xn_list.append( Transition((1,0,1), (0,1,1), 0.9, 'up motion') )

    if _DEBUG:
        print()
        print('setup_transition_list(): list has',len(xn_list),'transitions:')
        for t in xn_list:
            print('  From state',t.from_state,'to state',t.to_state,'at rate',t.rate,'called',t.name)

    return xn_list


def main():

    # INITIALIZE

    # User-defined parameters
    nr = 80
    nc = 80
    plot_interval = 2
    run_duration = 200
    report_interval = 5.0  # report interval, in real-time seconds

    # Remember the clock time, and calculate when we next want to report
    # progress.
    current_real_time = time.time()
    next_report = current_real_time + report_interval

    # Create grid
    mg = RasterModelGrid(nr, nc, 1.0)

    # Make the boundaries be walls
    mg.set_closed_boundaries_at_grid_edges(True, True, True, True)

    # Set up the states and pair transitions.
    ns_dict = { 0 : 'fluid', 1 : 'particle' }
    xn_list = setup_transition_list()

    # Create the node-state array and attach it to the grid
    node_state_grid = mg.add_zeros('node', 'node_state_map', dtype=int)

    # Initialize the node-state array
    middle_rows = where(bitwise_and(mg.node_y>0.45*nr, mg.node_y<0.55*nr))[0]
    node_state_grid[middle_rows] = 1

    # Create the CA model
    ca = OrientedRasterCTS(mg, ns_dict, xn_list, node_state_grid)

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
