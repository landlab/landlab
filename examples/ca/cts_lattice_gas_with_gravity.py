#!/usr/env/python
"""
cts_lattice_gas_with_gravity.py:
continuous-time stochastic version of a lattice-gas cellular
automaton model.

GT Sep 2014
"""
from __future__ import print_function

_DEBUG = False

import time
import random
from landlab import HexModelGrid
from landlab.ca.celllab_cts import Transition, CAPlotter
from landlab.ca.oriented_hex_cts import OrientedHexCTS


def setup_transition_list(g=1.0):
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
    """
    xn_list = []

    # Transitions for particle movement into an empty cell
    xn_list.append( Transition((1,0,0), (0,1,0), 1., 'motion') )
    xn_list.append( Transition((2,0,1), (0,2,1), 1., 'motion') )
    xn_list.append( Transition((3,0,2), (0,3,2), 1., 'motion') )
    xn_list.append( Transition((0,4,0), (4,0,0), 1., 'motion') )
    xn_list.append( Transition((0,5,1), (5,0,1), 1., 'motion') )
    xn_list.append( Transition((0,6,2), (6,0,2), 1., 'motion') )

    # Transitions for wall impact
    xn_list.append( Transition((1,8,0), (4,8,0), 1.0, 'wall rebound') )
    xn_list.append( Transition((2,8,1), (5,8,1), 1.0, 'wall rebound') )
    xn_list.append( Transition((3,8,2), (6,8,2), 1.0, 'wall rebound') )
    xn_list.append( Transition((8,4,0), (8,1,0), 1.0, 'wall rebound') )
    xn_list.append( Transition((8,5,1), (8,2,1), 1.0, 'wall rebound') )
    xn_list.append( Transition((8,6,2), (8,3,2), 1.0, 'wall rebound') )

    # Transitions for head-on collision
    xn_list.append( Transition((1,4,0), (3,6,0), 0.5, 'head-on collision') )
    xn_list.append( Transition((1,4,0), (5,2,0), 0.5, 'head-on collision') )
    xn_list.append( Transition((2,5,1), (4,1,1), 0.5, 'head-on collision') )
    xn_list.append( Transition((2,5,1), (6,3,1), 0.5, 'head-on collision') )
    xn_list.append( Transition((3,6,2), (1,4,2), 0.5, 'head-on collision') )
    xn_list.append( Transition((3,6,2), (5,2,2), 0.5, 'head-on collision') )

    # Transitions for glancing collision
    xn_list.append( Transition((1,3,0), (3,1,0), 1.0, 'glancing collision') )
    xn_list.append( Transition((1,5,0), (5,1,0), 1.0, 'glancing collision') )
    xn_list.append( Transition((2,4,0), (4,2,0), 1.0, 'glancing collision') )
    xn_list.append( Transition((6,4,0), (4,6,0), 1.0, 'glancing collision') )
    xn_list.append( Transition((2,4,1), (4,2,1), 1.0, 'glancing collision') )
    xn_list.append( Transition((2,6,1), (6,2,1), 1.0, 'glancing collision') )
    xn_list.append( Transition((1,5,1), (5,1,1), 1.0, 'glancing collision') )
    xn_list.append( Transition((3,5,1), (5,3,1), 1.0, 'glancing collision') )
    xn_list.append( Transition((3,1,2), (1,3,2), 1.0, 'glancing collision') )
    xn_list.append( Transition((3,5,2), (5,3,2), 1.0, 'glancing collision') )
    xn_list.append( Transition((2,6,2), (6,2,2), 1.0, 'glancing collision') )
    xn_list.append( Transition((4,6,2), (6,4,2), 1.0, 'glancing collision') )

    # Transitions for oblique-from-behind collisions
    xn_list.append( Transition((1,2,0), (2,1,0), 1.0, 'oblique') )
    xn_list.append( Transition((1,6,0), (6,1,0), 1.0, 'oblique') )
    xn_list.append( Transition((3,4,0), (4,3,0), 1.0, 'oblique') )
    xn_list.append( Transition((5,4,0), (4,5,0), 1.0, 'oblique') )
    xn_list.append( Transition((2,1,1), (1,2,1), 1.0, 'oblique') )
    xn_list.append( Transition((2,3,1), (3,2,1), 1.0, 'oblique') )
    xn_list.append( Transition((4,5,1), (5,4,1), 1.0, 'oblique') )
    xn_list.append( Transition((6,5,1), (5,6,1), 1.0, 'oblique') )
    xn_list.append( Transition((3,2,2), (2,3,2), 1.0, 'oblique') )
    xn_list.append( Transition((3,4,2), (4,3,2), 1.0, 'oblique') )
    xn_list.append( Transition((1,6,2), (6,1,2), 1.0, 'oblique') )
    xn_list.append( Transition((5,6,2), (6,5,2), 1.0, 'oblique') )

    # Transitions for direct-from-behind collisions
    xn_list.append( Transition((1,1,0), (2,6,0), 0.5, 'behind') )
    xn_list.append( Transition((1,1,0), (6,2,0), 0.5, 'behind') )
    xn_list.append( Transition((4,4,0), (3,5,0), 0.5, 'behind') )
    xn_list.append( Transition((4,4,0), (5,3,0), 0.5, 'behind') )
    xn_list.append( Transition((2,2,1), (1,3,1), 0.5, 'behind') )
    xn_list.append( Transition((2,2,1), (3,1,1), 0.5, 'behind') )
    xn_list.append( Transition((5,5,1), (4,6,1), 0.5, 'behind') )
    xn_list.append( Transition((5,5,1), (6,4,1), 0.5, 'behind') )
    xn_list.append( Transition((3,3,2), (2,4,2), 0.5, 'behind') )
    xn_list.append( Transition((3,3,2), (4,2,2), 0.5, 'behind') )
    xn_list.append( Transition((6,6,2), (1,5,2), 0.5, 'behind') )
    xn_list.append( Transition((6,6,2), (5,1,2), 0.5, 'behind') )

    # Transitions for collision with stationary (resting) particle
    xn_list.append( Transition((1,7,0), (7,2,0), 0.5, 'rest') )
    xn_list.append( Transition((1,7,0), (7,6,0), 0.5, 'rest') )
    xn_list.append( Transition((7,4,0), (3,7,0), 0.5, 'rest') )
    xn_list.append( Transition((7,4,0), (5,7,0), 0.5, 'rest') )
    xn_list.append( Transition((2,7,1), (7,1,1), 0.5, 'rest') )
    xn_list.append( Transition((2,7,1), (7,3,1), 0.5, 'rest') )
    xn_list.append( Transition((7,5,1), (4,7,1), 0.5, 'rest') )
    xn_list.append( Transition((7,5,1), (6,7,1), 0.5, 'rest') )
    xn_list.append( Transition((3,7,2), (7,2,2), 0.5, 'rest') )
    xn_list.append( Transition((3,7,2), (7,4,2), 0.5, 'rest') )
    xn_list.append( Transition((7,6,2), (1,7,2), 0.5, 'rest') )
    xn_list.append( Transition((7,6,2), (5,7,2), 0.5, 'rest') )

    # Gravity rules
    xn_list.append( Transition((1,0,0), (7,0,0), g, 'up to rest') )
    xn_list.append( Transition((1,1,0), (7,1,0), g, 'up to rest') )
    xn_list.append( Transition((1,2,0), (7,2,0), g, 'up to rest') )
    xn_list.append( Transition((1,3,0), (7,3,0), g, 'up to rest') )
    xn_list.append( Transition((1,4,0), (7,4,0), g, 'up to rest') )
    xn_list.append( Transition((1,5,0), (7,5,0), g, 'up to rest') )
    xn_list.append( Transition((1,6,0), (7,6,0), g, 'up to rest') )
    xn_list.append( Transition((1,7,0), (7,7,0), g, 'up to rest') )
    xn_list.append( Transition((0,1,0), (0,7,0), g, 'up to rest') )
    xn_list.append( Transition((1,1,0), (1,7,0), g, 'up to rest') )
    xn_list.append( Transition((2,1,0), (2,7,0), g, 'up to rest') )
    xn_list.append( Transition((3,1,0), (3,7,0), g, 'up to rest') )
    xn_list.append( Transition((4,1,0), (4,7,0), g, 'up to rest') )
    xn_list.append( Transition((5,1,0), (5,7,0), g, 'up to rest') )
    xn_list.append( Transition((6,1,0), (6,7,0), g, 'up to rest') )
    xn_list.append( Transition((7,1,0), (7,7,0), g, 'up to rest') )

    xn_list.append( Transition((7,0,0), (4,0,0), g, 'rest to down') )
    xn_list.append( Transition((7,1,0), (4,1,0), g, 'rest to down') )
    xn_list.append( Transition((7,2,0), (4,2,0), g, 'rest to down') )
    xn_list.append( Transition((7,3,0), (4,3,0), g, 'rest to down') )
    xn_list.append( Transition((7,4,0), (4,4,0), g, 'rest to down') )
    xn_list.append( Transition((7,5,0), (4,5,0), g, 'rest to down') )
    xn_list.append( Transition((7,6,0), (4,6,0), g, 'rest to down') )
    xn_list.append( Transition((7,7,0), (4,7,0), g, 'rest to down') )
    xn_list.append( Transition((0,7,0), (0,4,0), g, 'rest to down') )
    xn_list.append( Transition((1,7,0), (1,4,0), g, 'rest to down') )
    xn_list.append( Transition((2,7,0), (2,4,0), g, 'rest to down') )
    xn_list.append( Transition((3,7,0), (3,4,0), g, 'rest to down') )
    xn_list.append( Transition((4,7,0), (4,4,0), g, 'rest to down') )
    xn_list.append( Transition((5,7,0), (5,4,0), g, 'rest to down') )
    xn_list.append( Transition((6,7,0), (6,4,0), g, 'rest to down') )
    xn_list.append( Transition((7,7,0), (7,4,0), g, 'rest to down') )

    xn_list.append( Transition((2,0,1), (3,0,1), g, 'right up to right down') )
    xn_list.append( Transition((2,1,1), (3,1,1), g, 'right up to right down') )
    xn_list.append( Transition((2,2,1), (3,2,1), g, 'right up to right down') )
    xn_list.append( Transition((2,3,1), (3,3,1), g, 'right up to right down') )
    xn_list.append( Transition((2,4,1), (3,4,1), g, 'right up to right down') )
    xn_list.append( Transition((2,5,1), (3,5,1), g, 'right up to right down') )
    xn_list.append( Transition((2,6,1), (3,6,1), g, 'right up to right down') )
    xn_list.append( Transition((2,7,1), (3,7,1), g, 'right up to right down') )
    xn_list.append( Transition((0,2,1), (0,3,1), g, 'right up to right down') )
    xn_list.append( Transition((1,2,1), (1,3,1), g, 'right up to right down') )
    xn_list.append( Transition((2,2,1), (2,3,1), g, 'right up to right down') )
    xn_list.append( Transition((3,2,1), (3,3,1), g, 'right up to right down') )
    xn_list.append( Transition((4,2,1), (4,3,1), g, 'right up to right down') )
    xn_list.append( Transition((5,2,1), (5,3,1), g, 'right up to right down') )
    xn_list.append( Transition((6,2,1), (6,3,1), g, 'right up to right down') )
    xn_list.append( Transition((7,2,1), (7,3,1), g, 'right up to right down') )

    xn_list.append( Transition((6,0,2), (5,0,2), g, 'left up to left down') )
    xn_list.append( Transition((6,1,2), (5,1,2), g, 'left up to left down') )
    xn_list.append( Transition((6,2,2), (5,2,2), g, 'left up to left down') )
    xn_list.append( Transition((6,3,2), (5,3,2), g, 'left up to left down') )
    xn_list.append( Transition((6,4,2), (5,4,2), g, 'left up to left down') )
    xn_list.append( Transition((6,5,2), (5,5,2), g, 'left up to left down') )
    xn_list.append( Transition((6,6,2), (5,6,2), g, 'left up to left down') )
    xn_list.append( Transition((6,7,2), (5,7,2), g, 'left up to left down') )
    xn_list.append( Transition((0,6,2), (0,5,2), g, 'left up to left down') )
    xn_list.append( Transition((1,6,2), (1,5,2), g, 'left up to left down') )
    xn_list.append( Transition((2,6,2), (2,5,2), g, 'left up to left down') )
    xn_list.append( Transition((3,6,2), (3,5,2), g, 'left up to left down') )
    xn_list.append( Transition((4,6,2), (4,5,2), g, 'left up to left down') )
    xn_list.append( Transition((5,6,2), (5,5,2), g, 'left up to left down') )
    xn_list.append( Transition((6,6,2), (6,5,2), g, 'left up to left down') )
    xn_list.append( Transition((7,6,2), (7,5,2), g, 'left up to left down') )

    if _DEBUG:
        print()
        print('setup_transition_list(): list has',len(xn_list),'transitions:')
        for t in xn_list:
            print('  From state',t.from_state,'to state',t.to_state,'at rate',t.rate,'called',t.name)

    return xn_list


def main():

    # INITIALIZE

    # User-defined parameters
    nr = 41
    nc = 61
    g = 0.05
    plot_interval = 1.0
    run_duration = 100.0
    report_interval = 5.0  # report interval, in real-time seconds
    p_init = 0.1  # probability that a cell is occupied at start
    plot_every_transition = False

    # Remember the clock time, and calculate when we next want to report
    # progress.
    current_real_time = time.time()
    next_report = current_real_time + report_interval

    # Create a grid
    hmg = HexModelGrid(nr, nc, 1.0, orientation='vertical', reorient_links=True)

    # Close the grid boundaries
    #hmg.set_closed_nodes(hmg.open_boundary_nodes)

    # Set up the states and pair transitions.
    # Transition data here represent particles moving on a lattice: one state
    # per direction (for 6 directions), plus an empty state, a stationary
    # state, and a wall state.
    ns_dict = { 0 : 'empty',
                1 : 'moving up',
                2 : 'moving right and up',
                3 : 'moving right and down',
                4 : 'moving down',
                5 : 'moving left and down',
                6 : 'moving left and up',
                7 : 'rest',
                8 : 'wall'}
    xn_list = setup_transition_list(g)

    # Create data and initialize values.
    node_state_grid = hmg.add_zeros('node', 'node_state_grid')

    # Make the grid boundary all wall particles
    node_state_grid[hmg.boundary_nodes] = 8

    # Seed the grid interior with randomly oriented particles
    for i in hmg.core_nodes:
        if random.random()<p_init:
            node_state_grid[i] = random.randint(1, 7)

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
               plot_each_transition=plot_every_transition, plotter=ca_plotter)
        current_time += plot_interval

        # Plot the current grid
        ca_plotter.update_plot()


    # FINALIZE

    # Plot
    ca_plotter.finalize()


if __name__=='__main__':
    main()
