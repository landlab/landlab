#!/usr/env/python

"""
rock_weathering.py

CellLab-CTS model that simulates the weathering of rock to saprolite around
a network of fractures.

Created (and translated from earlier code by) by Greg Tucker, Jul 2015
"""
from __future__ import print_function

import time
import numpy as np
from landlab import RasterModelGrid
from landlab.ca.celllab_cts import Transition, CAPlotter
from landlab.ca.raster_cts import RasterCTS
from landlab.components.fracture_grid.fracture_grid import make_frac_grid
import matplotlib
from landlab.io.netcdf import write_netcdf


def setup_transition_list():
    """
    Creates and returns a list of Transition() objects to represent the
    grain-by-grain transformation of bedrock to saprolite.
    
    Returns
    -------
    xn_list : list of Transition objects
        List of objects that encode information about the link-state transitions.
    
    Notes
    -----
    Weathering here is treated very simply: a bedrock particle adjacent to a
    saprolite particle has a specified probability (rate) of weathering to
    saprolite; in other words, a rock-saprolite pair can turn into a
    saprolite-saprolite pair.
    
    The states and transitions are as follows:

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
    xn_list.append( Transition((0,1,0), (1,1,0), 1., 'weathering') )
    xn_list.append( Transition((1,0,0), (1,1,0), 1., 'weathering') )
    
    return xn_list
    
    
def main():
    
    # INITIALIZE

    # User-defined parameters
    nr = 200  # number of rows in grid
    nc = 200  # number of columns in grid
    plot_interval = 0.05   # time interval for plotting (unscaled)
    run_duration = 5.0   # duration of run (unscaled)
    report_interval = 10.0  # report interval, in real-time seconds
    frac_spacing = 10  # average fracture spacing, nodes
    outfilename = 'wx' # name for netCDF files
    
    # Remember the clock time, and calculate when we next want to report
    # progress.
    current_real_time = time.time()
    next_report = current_real_time + report_interval
    
    # Counter for output files
    time_slice = 0

    # Create grid
    mg = RasterModelGrid(nr, nc, 1.0)
    
    # Make the boundaries be walls
    mg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    
    # Set up the states and pair transitions.
    ns_dict = { 0 : 'rock', 1 : 'saprolite' }
    xn_list = setup_transition_list()

    # Create the node-state array and attach it to the grid.
    # (Note use of numpy's uint8 data type. This saves memory AND allows us
    # to write output to a netCDF3 file; netCDF3 does not handle the default
    # 64-bit integer type)
    node_state_grid = mg.add_zeros('node', 'node_state_map', dtype=np.uint8)
    
    node_state_grid[:] = make_frac_grid(frac_spacing, model_grid=mg)    
    
    # Create the CA model
    ca = RasterCTS(mg, ns_dict, xn_list, node_state_grid)

    # Set up the color map
    rock_color = (0.8, 0.8, 0.8)
    sap_color = (0.4, 0.2, 0)
    clist = [rock_color, sap_color]
    my_cmap = matplotlib.colors.ListedColormap(clist)
    
    # Create a CAPlotter object for handling screen display
    ca_plotter = CAPlotter(ca, cmap=my_cmap)
    
    # Plot the initial grid
    ca_plotter.update_plot()
    
    # Output the initial grid to file
    write_netcdf((outfilename+str(time_slice)+'.nc'), mg, 
                 #format='NETCDF3_64BIT',
                 names='node_state_map')

    # RUN
    current_time = 0.0
    while current_time < run_duration:
        
        # Once in a while, print out simulation and real time to let the user
        # know that the sim is running ok
        current_real_time = time.time()
        if current_real_time >= next_report:
            print('Current sim time', current_time, '(',
                  100 * current_time/run_duration, '%)')
            next_report = current_real_time + report_interval
        
        # Run the model forward in time until the next output step
        ca.run(current_time+plot_interval, ca.node_state, 
               plot_each_transition=False)
        current_time += plot_interval
        
        # Plot the current grid
        ca_plotter.update_plot()
        
        # Output the current grid to a netCDF file
        time_slice += 1
        write_netcdf((outfilename+str(time_slice)+'.nc'), mg, 
                     #format='NETCDF3_64BIT',
                     names='node_state_map')        
        

    # FINALIZE

    # Plot
    ca_plotter.finalize()



# If user runs this file, activate the main() function
if __name__ == "__main__":
    main()
