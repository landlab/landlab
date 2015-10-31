#!/usr/env/python
"""
grains.py: simple cellular automaton model of granular movement under gravity.

This code provides an example of a Landlab cellular automaton that uses an
oriented hex grid. The model represents grains in fluid (such as air). The rules
are simple: a grain that lies either directly above, or above and to the side,
of a fluid-filled space, will fall to fill the space with a specified rate. This
model is not especially realistic; the reasons why include the lack of any
rebound/collision physics, and the fact that sideways motion of a grain is
independent of what's below it (fluid, or another grain). The point here is
simply to demonstrate the capability of an oriented hex grid.

GT Sep 2014
"""
from __future__ import print_function

_DEBUG = False

import time
from landlab import HexModelGrid
from numpy import where, logical_and, sqrt
from landlab.ca.celllab_cts import Transition, CAPlotter
from landlab.ca.oriented_hex_cts import OrientedHexCTS


def setup_transition_list():
    """
    Creates and returns a list of Transition() objects to represent state
    transitions for simple granular mechanics model.

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

    Pair state        Transition to       Process
    ==========        =============       =======
    0 (0-0-0)         (none)              -
    1 (0-1-0)         2 (1-0-0)           falling
    2 (1-0-0)         (none)              -
    3 (1-1-0)         (none)              -
    4 (0-0-1)         (none)              -
    5 (0-1-1)         6 (1-0-1)           falling
    6 (1-0-1)         (none)              -
    7 (1-1-1)         (none)              -
    8 (0-0-2)         (none)              -
    9 (0-1-2)         (none)              -
    10 (1-0-2)        9 (0-1-2)           falling
    11 (1-1-2)        (none)              -

    """
    xn_list = []

    xn_list.append( Transition((0,1,0), (1,0,0), 10., 'falling') )
    xn_list.append( Transition((0,1,1), (1,0,1), 1., 'falling') )
    xn_list.append( Transition((1,0,2), (0,1,2), 1., 'falling') )

    if _DEBUG:
        print()
        print('setup_transition_list(): list has',len(xn_list),'transitions:')
        for t in xn_list:
            print('  From state',t.from_state,'to state',t.to_state,'at rate',t.rate,'called',t.name)

    return xn_list


def main():

    # INITIALIZE

    # User-defined parameters
    nr = 21
    nc = 21
    plot_interval = 0.5
    run_duration = 25.0
    report_interval = 5.0  # report interval, in real-time seconds

    # Remember the clock time, and calculate when we next want to report
    # progress.
    current_real_time = time.time()
    next_report = current_real_time + report_interval

    # Create a grid
    hmg = HexModelGrid(nr, nc, 1.0, orientation='vertical', reorient_links=True)

    # Close the grid boundaries
    hmg.set_closed_nodes(hmg.open_boundary_nodes)

    # Set up the states and pair transitions.
    # Transition data here represent the disease status of a population.
    ns_dict = { 0 : 'fluid', 1 : 'grain' }
    xn_list = setup_transition_list()

    # Create data and initialize values. We start with the 3 middle columns full
    # of grains, and the others empty.
    node_state_grid = hmg.add_zeros('node', 'node_state_grid')
    middle = 0.25*(nc-1)*sqrt(3)
    is_middle_cols = logical_and(hmg.node_x<middle+1., hmg.node_x>middle-1.)
    node_state_grid[where(is_middle_cols)[0]] = 1

    # Create the CA model
    ca = OrientedHexCTS(hmg, ns_dict, xn_list, node_state_grid)

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
               plot_each_transition=False)
        current_time += plot_interval

        # Plot the current grid
        ca_plotter.update_plot()


    # FINALIZE

    # Plot
    ca_plotter.finalize()


if __name__=='__main__':
    main()
