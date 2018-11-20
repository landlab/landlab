"""
settle_test.py

Tests "external update" capability in landlab_ca.py.

GT Nov 2014
"""
from __future__ import print_function

import time

from numpy import where

from landlab import RasterModelGrid
from landlab.ca.celllab_cts import CAPlotter, Transition
from landlab.ca.oriented_raster_cts import OrientedRasterCTS

_DEBUG = False


def setup_transition_list():
    """
    Creates and returns a list of Transition() objects to represent state
    transitions for particles that simply settle under gravity.

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
    0 (0-0,0)
    1 (0-1,0)
    2 (1-0,0)
    3 (1-1,0)
    4 (0-0,1)
    5 (0-1,1)
    1 (0-1)         2 (1-0)             settling    1.0

    """
    xn_list = []

    xn_list.append(Transition(5, 6, 1., "settling"))

    if _DEBUG:
        print()
        print("setup_transition_list(): list has", len(xn_list), "transitions:")
        for t in xn_list:
            print(
                "  From state",
                t.from_state,
                "to state",
                t.to_state,
                "at rate",
                t.rate,
                "called",
                t.name,
            )

    return xn_list


def main():

    # INITIALIZE

    # User-defined parameters
    nr = 10
    nc = 10
    plot_interval = 0.25
    run_duration = 40.0
    report_interval = 5.0  # report interval, in real-time seconds

    # Remember the clock time, and calculate when we next want to report
    # progress.
    current_real_time = time.time()
    next_report = current_real_time + report_interval

    # Create grid
    mg = RasterModelGrid(nr, nc, 1.0)
    mg.set_closed_boundaries_at_grid_edges(True, True, True, True)

    # Set up the states and pair transitions.
    # Transition data here represent a body of fractured rock, with rock
    # represented by nodes with state 0, and saprolite (weathered rock)
    # represented by nodes with state 1. Node pairs (links) with 0-1 or 1-0
    # can undergo a transition to 1-1, representing chemical weathering of the
    # rock.
    ns_dict = {0: "air", 1: "particle"}
    xn_list = setup_transition_list()

    # Create the node-state array and attach it to the grid
    node_state_grid = mg.add_zeros("node", "node_state_map", dtype=int)
    node_state_grid[where(mg.node_y > nr - 3)[0]] = 1

    # Create the CA model
    ca = OrientedRasterCTS(mg, ns_dict, xn_list, node_state_grid)
    # ca = RasterCTS(mg, ns_dict, xn_list, node_state_grid)

    # Debug output if needed
    if _DEBUG:
        n = ca.grid.number_of_nodes
        for r in range(ca.grid.number_of_node_rows):
            for c in range(ca.grid.number_of_node_columns):
                n -= 1
                print("{0:.0f}".format(ca.node_state[n]), end=" ")
            print()

    # Create a CAPlotter object for handling screen display
    ca_plotter = CAPlotter(ca)

    # Plot the initial grid
    ca_plotter.update_plot()

    # RUN
    current_time = 0.0
    updated = False
    while current_time < run_duration:

        # Once in a while, print out simulation and real time to let the user
        # know that the sim is running ok
        current_real_time = time.time()
        if current_real_time >= next_report:
            print(
                "Current sim time",
                current_time,
                "(",
                100 * current_time / run_duration,
                "%)",
            )
            next_report = current_real_time + report_interval

        # Run the model forward in time until the next output step
        ca.run(
            current_time + plot_interval, ca.node_state, plot_each_transition=False
        )  # , plotter=ca_plotter)
        current_time += plot_interval

        # Add a bunch of particles
        if current_time > run_duration / 2. and not updated:
            print("updating...")
            node_state_grid[where(ca.grid.node_y > (nc / 2.0))[0]] = 1
            ca.update_link_states_and_transitions(current_time)
            updated = True

        # Plot the current grid
        ca_plotter.update_plot()

        # for debugging
        if _DEBUG:
            n = ca.grid.number_of_nodes
            for r in range(ca.grid.number_of_node_rows):
                for c in range(ca.grid.number_of_node_columns):
                    n -= 1
                    print("{0:.0f}".format(ca.node_state[n]), end=" ")
                print()

    # FINALIZE

    # Plot
    ca_plotter.finalize()


if __name__ == "__main__":
    main()
