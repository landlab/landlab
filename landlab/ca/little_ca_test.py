#!/usr/env/python
"""Test the creation and execution of a CellLab-CTS model.

Tests the creation and execution of a CellLab-CTS model, by creating a
simple two-state CA on a small grid.

Created by Greg Tucker, May 2015
"""
from __future__ import print_function

import time
from landlab import RasterModelGrid
from landlab.ca.celllab_cts import Transition, CAPlotter
from landlab.ca.raster_cts import RasterCTS


def setup_transition_list():
    """
    Creates and returns a list of Transition() objects to represent state
    transitions for a biased random walk, in which the rate of downward
    motion is greater than the rate in the other three directions.

    Returns
    -------
    xn_list : list of Transition objects
        List of objects that encode information about the link-state
        transitions.

    Notes
    -----
    This doesn't represent any particular process, but rather is simply used
    to test the CA code. The transition rules have 0-0 pairs transitioning to
    0-1 or 1-0 pairs (50/50 chance) and thence to 1-1 pairs, at which point
    there are no further transitions.

    The states and transitions are as follows::

        Pair state      Transition to       Process             Rate (cells/s)
        ==========      =============       =======             ==============
        0 (0-0)         1 (0-1)                                 0.5
                        2 (1-0)                                 0.5
        1 (0-1)         3 (1-1)                                 1.0
        2 (1-0)         3 (1-1)                                 1.0
        3 (1-1)         (none)                                  -

    """
    # Create an empty transition list
    xn_list = []

    # Append two transitions to the list.
    # Note that the arguments to the Transition() object constructor are:
    #  - Tuple representing starting pair state
    #    (left/bottom cell, right/top cell, orientation)
    #  - Tuple representing new pair state
    #    (left/bottom cell, right/top cell, orientation)
    #  - Transition rate (cells per time step, in this case 1 sec)
    #  - Name for transition
    xn_list.append(Transition((0, 0, 0), (0, 1, 0), 0.5, ''))
    xn_list.append(Transition((0, 0, 0), (1, 0, 0), 0.5, ''))
    xn_list.append(Transition((0, 1, 0), (1, 1, 0), 1., ''))
    xn_list.append(Transition((1, 0, 0), (1, 1, 0), 1., ''))

    return xn_list


def main():

    # INITIALIZE

    # User-defined parameters
    nr = 5  # number of rows in grid
    nc = 5  # number of columns in grid
    plot_interval = 10.0   # time interval for plotting, sec
    run_duration = 10.0   # duration of run, sec
    report_interval = 10.0  # report interval, in real-time seconds

    # Remember the clock time, and calculate when we next want to report
    # progress.
    current_real_time = time.time()
    next_report = current_real_time + report_interval

    # Create grid
    mg = RasterModelGrid(nr, nc, 1.0)

    # Make the boundaries be walls
    mg.set_closed_boundaries_at_grid_edges(True, True, True, True)

    # Set up the states and pair transitions.
    ns_dict = {0: 'black', 1: 'white'}
    xn_list = setup_transition_list()

    # Create the node-state array and attach it to the grid
    node_state_grid = mg.add_zeros('node', 'node_state_map', dtype=int)

    # For visual display purposes, set all boundary nodes to fluid
    node_state_grid[mg.closed_boundary_nodes] = 0

    # Create the CA model
    ca = RasterCTS(mg, ns_dict, xn_list, node_state_grid)

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
            print('Current sim time', current_time, '(',
                  100 * current_time / run_duration, '%)')
            next_report = current_real_time + report_interval

        # Run the model forward in time until the next output step
        ca.run(current_time + plot_interval, ca.node_state,
               plot_each_transition=True, plotter=ca_plotter)
        current_time += plot_interval

        # Plot the current grid
        ca_plotter.update_plot()

    # FINALIZE

    # Plot
    ca_plotter.finalize()

    print('ok, here are the keys')
    print(ca.__dict__.keys())


# If user runs this file, activate the main() function
if __name__ == "__main__":
    #import cProfile
    #fname = 'test_profiler_for_little_ca.txt'
    # cProfile.run('print main(); print') #, fname)
    main()
