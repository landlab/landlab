#!/usr/env/python
"""
sir.py: 

Example of a Susceptible-Infectious-Recovered epidemiological 
cellular automaton model implemented on a hexagonal grid using stochastic
pair-transition rules.

GT Sep 2014; updated for Landlab 2.0beta, March 2020
"""

import sys
import time
from landlab import HexModelGrid
import numpy as np
from landlab.ca.celllab_cts import Transition, CAPlotter
from landlab.ca.hex_cts import HexCTS
import matplotlib.pyplot as plt


def setup_transition_list(infection_rate):
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
    
    xn_list.append( Transition((0,1,0), (1,1,0), infection_rate, 'infection') )
    xn_list.append( Transition((0,1,0), (0,2,0), 1., 'recovery') )
    xn_list.append( Transition((1,0,0), (1,1,0), infection_rate, 'infection') )
    xn_list.append( Transition((1,0,0), (2,0,0), 1., 'recovery') )
    xn_list.append( Transition((1,1,0), (1,2,0), 1., 'recovery') )
    xn_list.append( Transition((1,1,0), (2,1,0), 1., 'recovery') )
    xn_list.append( Transition((1,2,0), (2,2,0), 1., 'recovery') )
    xn_list.append( Transition((2,1,0), (2,2,0), 1., 'recovery') )
        
    return xn_list

def set_colormap():
    """Create and return colormap for plotting."""
    import matplotlib
    susceptible_color = (0.5, 0.5, 0.5)  # gray
    infectious_color = (0.5, 0.0, 0.0)  # dark red
    recovered_color = (0.0, 0.0, 1.0)  # blue
    clist = [susceptible_color, infectious_color, recovered_color]
    return matplotlib.colors.ListedColormap(clist)


class StochasticCellularSIRmodel(object):
    """Continuous-time stochastic cellular Susceptible-Infectious-Recovered
    model.
    """
    def __init__(self,
                 infection_rate=2.0,
                 number_of_node_rows=80,
                 number_of_base_node_columns=41,
                 run_duration=5.0,
                 report_interval=10.0,
                 output_interval=1.0,
                 output_file_base_name='sirmodel',
                 show_plots=False,
                 seed=0,
                 ):
        """Initialize."""

        # Create grid
        self.grid = HexModelGrid(shape=(number_of_node_rows,
                                        number_of_base_node_columns),
                                 spacing=1.0)

        # Create node-state field
        node_state_grid = self.grid.add_zeros('node', 'node_state_grid',
                                              dtype=np.int)
        node_state_grid[self.grid.number_of_nodes//2] = 1
        node_state_grid[0] = 2  # full color range: set node 0 to 'recovered'

        # Set up the states and pair transitions.
        # Transition data here represent the disease status of a population.
        ns_dict = { 0 : 'susceptible', 1 : 'infectious', 2: 'recovered' }
        xn_list = setup_transition_list(infection_rate)

        # Create the CA model
        self.ca = HexCTS(self.grid, ns_dict, xn_list, node_state_grid,
                         seed=seed)

        # Set up parameters and timing
        self.run_duration = run_duration
        self.report_interval = report_interval
        self.outfilename = output_file_base_name

        # Set up output and plotting
        self.output_interval = output_interval
        self.time_slice = 0
        my_cmap = set_colormap()
        self.ca_plotter = CAPlotter(self.ca, cmap=my_cmap)

        # First plot/save
        self.plot_and_save(0)

    def plot_and_save(self, time_slice_number):
        """Plot to screen and file."""
        self.ca_plotter.update_plot()
        plt.axis('off')
        savename = self.outfilename + str(time_slice_number)
        plt.savefig(savename+'.pdf', format='pdf')

    def run(self):
        """Run for complete duration."""
        current_time = 0.0
        time_slice = 1
        next_report = time.time() + self.report_interval
        while current_time < self.run_duration:

            # Once in a while, print out simulation and real time to let the 
            # user know that the sim is running ok
            current_real_time = time.time()
            if current_real_time >= next_report:
                print('Current simulation time: ', current_time, '(',
                      100 * current_time/self.run_duration, '%)')
                next_report = current_real_time + self.report_interval

            # Run the model forward in time until the next output step
            self.ca.run(current_time + self.output_interval, self.ca.node_state, 
                        plot_each_transition=False)
            current_time += self.output_interval

            # Plot the current grid
            self.plot_and_save(time_slice)
            time_slice += 1

    def finalize(self):
        """Finish up."""
        self.ca_plotter.finalize()
        plt.clf()

        # Report
        print('Done. Close figure window to finish.')



def run_sir(params=None):
    """Run the SIR model with given parameters or params in given input file"""

    # Handle input parameters
    if type(params) is not dict:
        if len(sys.argv) > 1:
            from landlab import load_params
            print('Reading parameters from file "' + sys.argv[1] + '"...',
                  end='')
            params = load_params(sys.argv[1])
            print('done.')
        else:
            params = {}

    # Initial and run the model, then clean up
    sir_model = StochasticCellularSIRmodel(**params)
    sir_model.run()
    sir_model.finalize()


if __name__=='__main__':
    run_sir()
