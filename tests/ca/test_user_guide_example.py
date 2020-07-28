#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_user_guide_example.py

This test checks the example presented in the CellLab-CTS User Guide. If this
test fails as part of a compatibility-breaking major version update, the
corresponding code in the User Guide should also be updated. Some limits:
    - plotting is not tested (commented out here)
    - main is renamed here
    - run_duration is reduced to speed up test

Created on Thu Apr 30 11:43:50 2020

@author: gtucker
"""

import time

# import matplotlib (commented out for CI testing)
from numpy import where

from landlab import RasterModelGrid
from landlab.ca.celllab_cts import Transition  # , CAPlotter (# CI testing)
from landlab.ca.raster_cts import RasterCTS


def setup_transition_list():
    """
    Creates and returns a list of Transition() objects to represent state
    transitions for an unbiased random walk.

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
    sediment grain, tea leaf, or solute molecule).

    The states and transitions are as follows:

    Pair state      Transition to       Process             Rate (cells/s)
    ==========      =============       =======             ==============
    0 (0-0)         (none)              -                   -
    1 (0-1)         2 (1-0)             left/down motion    10.0
    2 (1-0)         1 (0-1)             right/up motion     10.0
    3 (1-1)         (none)              -                   -

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
    xn_list.append(Transition((0, 1, 0), (1, 0, 0), 10.0, "left/down motion"))
    xn_list.append(Transition((1, 0, 0), (0, 1, 0), 10.0, "right/up motion"))

    return xn_list


def test_user_guide_example():

    # INITIALIZE

    # User-defined parameters
    nr = 80  # number of rows in grid
    nc = 50  # number of columns in grid
    plot_interval = 0.5  # time interval for plotting, sec
    run_duration = 1.0  # duration of run, sec
    report_interval = 10.0  # report interval, in real-time seconds

    # Remember the clock time, and calculate when we next want to report
    # progress.
    current_real_time = time.time()
    next_report = current_real_time + report_interval

    # Create grid
    mg = RasterModelGrid(shape=(nr, nc), xy_spacing=1.0)

    # Make the boundaries be walls
    mg.set_closed_boundaries_at_grid_edges(True, True, True, True)

    # Create a node-state dictionary
    ns_dict = {0: "fluid", 1: "particle"}

    # Create the transition list
    xn_list = setup_transition_list()

    # Create an array containing the initial node-state values

    # Create the node-state array and attach it to the grid
    node_state_grid = mg.add_zeros("node", "node_state_map", dtype=int)

    # Initialize the node-state array: here, the initial condition is a pile of
    # resting grains at the bottom of a container.
    bottom_rows = where(mg.node_y < 0.1 * nr)[0]
    node_state_grid[bottom_rows] = 1

    # For visual display purposes, set all boundary nodes to fluid
    node_state_grid[mg.closed_boundary_nodes] = 0

    # Create the CA model
    ca = RasterCTS(mg, ns_dict, xn_list, node_state_grid)

    # Set up plotting
    # Set up colors for plotting
    # (commented out for CI testing)
    # grain = '#5F594D'
    # fluid = '#D0E4F2'
    # clist = [fluid,grain]
    # my_cmap = matplotlib.colors.ListedColormap(clist)

    # Create a CAPlotter object for handling screen display
    # (commented out for CI testing)
    # ca_plotter = CAPlotter(ca, cmap=my_cmap)

    # Plot the initial grid (commented out for CI testing)
    # ca_plotter.update_plot()

    # RUN
    current_time = 0.0
    while current_time < run_duration:

        # Once in a while, print out simulation real time to let the user
        # know that the sim is running ok
        current_real_time = time.time()
        if current_real_time >= next_report:
            print(
                "Current simulation time "
                + str(current_time)
                + "  \
                   ("
                + str(int(100 * current_time / run_duration))
                + "%)"
            )
            next_report = current_real_time + report_interval

        # Run the model forward in time until the next output step
        ca.run(current_time + plot_interval, ca.node_state, plot_each_transition=False)
        current_time += plot_interval

        # Plot the current grid (commented out for CI testing)
        # ca_plotter.update_plot()

    # Finalize plot (commented out for CI testing)
    # ca_plotter.finalize()
