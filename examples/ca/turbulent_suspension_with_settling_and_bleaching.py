#!/usr/env/python

"""
isotropic_turbulent_suspension_with_settling_and_bleaching.py

Example of a continuous-time, stochastic, pair-based cellular automaton model,
which simulates the diffusion of suspended particles in a turbulent fluid.
Particles start with an accumulated luminescence signal L = 1, and are bleached
by exposure to light at a rate that depends on distance below the upper surface.

Written by Greg Tucker, July 2015
"""

from __future__ import print_function  # for both python 2 and 3 compability
import time
import matplotlib
from pylab import figure, show, clf
from numpy import where, exp, amin
from landlab import RasterModelGrid, ModelParameterDictionary
from landlab.plot.imshow import imshow_node_grid
from landlab.components.cellular_automata.celllab_cts import Transition, CAPlotter
from landlab.components.cellular_automata.oriented_raster_cts import OrientedRasterCTS




class TurbulentSuspensionAndBleachingModel(OrientedRasterCTS):
    """
    Examples
    ---------
    >>> from six import StringIO
    >>> p = StringIO('''
    ... model_grid_row__count: number of rows in grid
    ... 4
    ... model_grid_column__count: number of columns in grid
    ... 4
    ... plot_interval: interval for plotting to display, s
    ... 2.0
    ... model__run_time: duration of model run, s
    ... 1.0
    ... model__report_interval: time interval for reporting progress, real-time seconds
    ... 1.0e6
    ... surface_bleaching_time_scale: time scale for OSL bleaching, s
    ... 2.42
    ... light_attenuation_length: length scale for light attenuation, cells (1 cell = 1 mm)
    ... 2.0
    ... ''')
    >>> tsbm = TurbulentSuspensionAndBleachingModel(p)
    >>> tsbm.node_state
    array([1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
    >>> tsbm.grid.at_node['osl']
    array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.,  0.,  0.,
            0.,  0.,  0.])
    >>> tsbm.n_xn
    array([0, 1, 1, 0, 0, 1, 1, 0])
    >>> tsbm.fluid_surface_height
    3.5
    """

    def __init__(self, input_stream):
        """
        Reads in parameters and initializes the model.

        Examples
        --------
        >>> from six import StringIO
        >>> p = StringIO('''
        ... model_grid_row__count: number of rows in grid
        ... 4
        ... model_grid_column__count: number of columns in grid
        ... 4
        ... plot_interval: interval for plotting to display, s
        ... 2.0
        ... model__run_time: duration of model run, s
        ... 1.0
        ... model__report_interval: time interval for reporting progress, real-time seconds
        ... 1.0e6
        ... surface_bleaching_time_scale: time scale for OSL bleaching, s
        ... 2.42
        ... light_attenuation_length: length scale for light attenuation, cells (1 cell = 1 mm)
        ... 2.0
        ... ''')
        >>> tsbm = TurbulentSuspensionAndBleachingModel(p)
        >>> tsbm.node_state
        array([1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
        >>> tsbm.grid.at_node['osl']
        array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.,  0.,  0.,
                0.,  0.,  0.])
        >>> tsbm.n_xn
        array([0, 1, 1, 0, 0, 1, 1, 0])
        >>> tsbm.fluid_surface_height
        3.5
        """
        # Get a source for input parameters.
        params = ModelParameterDictionary(input_stream)

        # Read user-defined parameters
        nr = params.read_int('model_grid_row__count')    # number of rows (CSDMS Standard Name [CSN])
        nc = params.read_int('model_grid_column__count') # number of cols (CSN)
        self.plot_interval = params.read_float('plot_interval')  # interval for plotting output, s
        self.run_duration = params.read_float('model__run_time')   # duration of run, sec (CSN)
        self.report_interval = params.read_float('model__report_interval')  # report interval, in real-time seconds
        self.bleach_T0 = params.read_float('surface_bleaching_time_scale')  # time scale for bleaching at fluid surface, s
        self.zstar = params.read_float('light_attenuation_length')  # length scale for light attenuation in fluid, CELLS

        # Derived parameters
        self.fluid_surface_height = nr-0.5

        # Calculate when we next want to report progress.
        self.next_report = time.time() + self.report_interval

        # Create grid
        mg = RasterModelGrid(nr, nc, 1.0)

        # Make the boundaries be walls
        mg.set_closed_boundaries_at_grid_edges(True, True, True, True)

        # Set up the states and pair transitions.
        ns_dict = { 0 : 'fluid', 1 : 'particle' }
        xn_list = self.setup_transition_list()

        # Create the node-state array and attach it to the grid
        node_state_grid = mg.add_zeros('node', 'node_state_map', dtype=int)

        # For visual display purposes, set all boundary nodes to fluid
        node_state_grid[mg.closed_boundary_nodes] = 0

        # Initialize the node-state array: here, the initial condition is a pile of
        # resting grains at the bottom of a container.
        bottom_rows = where(mg.node_y<0.4*nr)[0]
        node_state_grid[bottom_rows] = 1

        # Create a data array for bleaching.
        # Here, osl=optically stimulated luminescence, normalized to the original
        # signal (hence, initially all unity). Over time this signal may get
        # bleached out due to exposure to light.
        self.osl = mg.add_zeros('node', 'osl')
        self.osl[bottom_rows] = 1.0
        self.osl_display = mg.add_zeros('node', 'osl_display')
        self.osl_display[bottom_rows] = 1.0

        # We'll need an array to track the last time any given node was
        # updated, so we can figure out the duration of light exposure between
        # update events
        self.last_update_time = mg.add_zeros('node','last_update_time')

        # Call the  base class (RasterCTS) init method
        super(TurbulentSuspensionAndBleachingModel, \
              self).__init__(mg, ns_dict, xn_list, node_state_grid, prop_data=self.osl)

        # Set up plotting (if plotting desired)
        if self.plot_interval <= self.run_duration:
            self.initialize_plotting()


    def initialize_plotting(self):
        """
        Creates a CA plotter object, sets its colormap, and plots the initial
        model state.
        """
        # Set up some plotting information
        grain = '#5F594D'
        bleached_grain = '#CC0000'
        fluid = '#D0E4F2'
        clist = [fluid,bleached_grain,grain]
        my_cmap = matplotlib.colors.ListedColormap(clist)

        # Create a CAPlotter object for handling screen display
        self.ca_plotter = CAPlotter(self, cmap=my_cmap)

        # Plot the initial grid
        self.ca_plotter.update_plot()

        # Make a colormap for use in showing the bleaching of each grain
        clist = [(0.0, (1.0, 1.0, 1.0)), (0.49, (0.8, 0.8, 0.8)), (1.0, (0.0, 0.0, 0.0))]
        self.cmap_for_osl = matplotlib.colors.LinearSegmentedColormap.from_list('osl_cmap', clist)


    def setup_transition_list(self):
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
        sediment grain, tea leaf, or dissolved heavy particle).

        The states and transitions are as follows:

        Pair state      Transition to       Process             Rate (cells/s)
        ==========      =============       =======             ==============
        0 (0-0)         (none)              -                   -
        1 (0-1)         2 (1-0)             left motion         10.0
        2 (1-0)         1 (0-1)             right motion        10.0
        3 (1-1)         (none)              -                   -
        4 (0-0)         (none)              -                   -
        5 (0-1)         2 (1-0)             down motion         10.55
        6 (1-0)         1 (0-1)             up motion            9.45
        7 (1-1)         (none)              -                   -

        """

        # Create an empty transition list
        xn_list = []

        # Append four transitions to the list.
        # Note that the arguments to the Transition() object constructor are:
        #  - Tuple representing starting pair state
        #    (left cell, right cell, orientation [0=horizontal])
        #  - Tuple representing new pair state
        #    (bottom cell, top cell, orientation [1=vertical])
        #  - Transition rate (cells per time step, in this case 1 sec)
        #  - Name for transition
        #  - Flag indicating that the transition involves an exchange of properties
        #  - Function to be called after each transition, to update a property
        #    (in this case, to simulate bleaching of the luminescence signal)
        xn_list.append( Transition((0,1,0), (1,0,0), 10., 'left motion', True, self.update_bleaching) )
        xn_list.append( Transition((1,0,0), (0,1,0), 10., 'right motion', True, self.update_bleaching) )
        xn_list.append( Transition((0,1,1), (1,0,1), 10.55, 'down motion', True, self.update_bleaching) )
        xn_list.append( Transition((1,0,1), (0,1,1), 9.45, 'up motion', True, self.update_bleaching) )

        return xn_list


    def bleach_grain(self, node, dt):
        """
        Updates the luminescence signal at node.

        Examples
        --------
        >>> from six import StringIO
        >>> p = StringIO('''
        ... model_grid_row__count: number of rows in grid
        ... 10
        ... model_grid_column__count: number of columns in grid
        ... 3
        ... plot_interval: interval for plotting to display, s
        ... 2.0
        ... model__run_time: duration of model run, s
        ... 1.0
        ... model__report_interval: time interval for reporting progress, real-time seconds
        ... 1.0e6
        ... surface_bleaching_time_scale: time scale for OSL bleaching, s
        ... 2.42
        ... light_attenuation_length: length scale for light attenuation, cells (1 cell = 1 mm)
        ... 6.5
        ... ''')
        >>> tsbm = TurbulentSuspensionAndBleachingModel(p)
        >>> tsbm.bleach_grain(10, 1.0)
        >>> int(tsbm.prop_data[tsbm.propid[10]]*1000)
        858
        """
        depth = self.fluid_surface_height - self.grid.node_y[node]
        T_bleach = self.bleach_T0*exp( depth/self.zstar)
        self.prop_data[self.propid[node]] *= exp( -dt/T_bleach )


    def update_bleaching(self, ca_unused, node1, node2, time_now):
        """
        Updates the luminescence signal at a pair of nodes that have just
        undergone a transition, if either or both nodes is a grain.

        Examples
        --------
        >>> from six import StringIO
        >>> p = StringIO('''
        ... model_grid_row__count: number of rows in grid
        ... 10
        ... model_grid_column__count: number of columns in grid
        ... 3
        ... plot_interval: interval for plotting to display, s
        ... 2.0
        ... model__run_time: duration of model run, s
        ... 1.0
        ... model__report_interval: time interval for reporting progress, real-time seconds
        ... 1.0e6
        ... surface_bleaching_time_scale: time scale for OSL bleaching, s
        ... 2.42
        ... light_attenuation_length: length scale for light attenuation, cells (1 cell = 1 mm)
        ... 6.5
        ... ''')
        >>> tsbm = TurbulentSuspensionAndBleachingModel(p)
        >>> tsbm.update_bleaching(tsbm, 10, 13, 1.0)
        >>> int(tsbm.prop_data[tsbm.propid[10]]*1000)
        858
        >>> tsbm.prop_data[tsbm.propid[13]]
        0.0
        """
        if self.node_state[node1]==1:
            dt = time_now - self.last_update_time[self.propid[node1]]
            self.bleach_grain(node1, dt)
            self.last_update_time[self.propid[node1]] = time_now
        if self.node_state[node2]==1:
            dt = time_now - self.last_update_time[self.propid[node2]]
            self.bleach_grain(node2, dt)
            self.last_update_time[self.propid[node2]] = time_now


    def synchronize_bleaching(self, sync_time):
        """
        Brings all nodes up to the same time, sync_time, by applying bleaching
        up to this time, and updating last_update_time.

        Notes
        -----
        In a CellLab-CTS model, the "time" is usually different for each node:
        some will have only just recently undergone a transition and had their
        properties (in this case, OSL bleaching) updated, while others will
        have last been updated a long time ago, and some may never have had a
        transition. If we want to plot the properties at a consistent time, we
        need to bring all node properties (again, in this case, OSL) up to
        date. This method does so.
            We multiply elapsed time (between last update and "sync time") by
        the node state, because we only want to update the solid particles---
        because the state of a particle is 1 and fluid 0, this multiplication
        masks out the fluid nodes.
            We don't call bleach_grain(), because we want to take advantage of
        numpy array operations rather than calling a method for each node.

        Examples
        --------
        >>> from six import StringIO
        >>> p = StringIO('''
        ... model_grid_row__count: number of rows in grid
        ... 10
        ... model_grid_column__count: number of columns in grid
        ... 3
        ... plot_interval: interval for plotting to display, s
        ... 2.0
        ... model__run_time: duration of model run, s
        ... 1.0
        ... model__report_interval: time interval for reporting progress, real-time seconds
        ... 1.0e6
        ... surface_bleaching_time_scale: time scale for OSL bleaching, s
        ... 2.42
        ... light_attenuation_length: length scale for light attenuation, cells (1 cell = 1 mm)
        ... 6.5
        ... ''')
        >>> tsbm = TurbulentSuspensionAndBleachingModel(p)
        >>> tsbm.synchronize_bleaching(1.0)
        >>> int(tsbm.osl[10]*100000)
        85897
        """
        dt = (sync_time - self.last_update_time[self.propid])*self.node_state
        assert (amin(dt)>=0.0), 'sync_time must be >= 0 everywhere'
        depth = self.fluid_surface_height - self.grid.node_y
        T_bleach = self.bleach_T0*exp( depth/self.zstar)
        self.prop_data[self.propid] *= exp( -dt/T_bleach )
        self.last_update_time[self.propid] = sync_time*self.node_state


    def go(self):
        """
        Runs the model.
        """
        # RUN
        while self.current_time < self.run_duration:

            # Once in a while, print out simulation and real time to let the user
            # know that the sim is running ok
            current_real_time = time.time()
            if current_real_time >= self.next_report:
                print('Current sim time',self.current_time,'(',100*self.current_time/self.run_duration,'%)')
                self.next_report = current_real_time + self.report_interval

            # Run the model forward in time until the next output step
            self.run(self.current_time+self.plot_interval, self.node_state,
                   plot_each_transition=False)
            self.current_time += self.plot_interval
            self.synchronize_bleaching(self.current_time)

            if self.plot_interval <= self.run_duration:

                # Plot the current grid
                self.ca_plotter.update_plot()

                # Display the OSL content of grains
                figure(3)
                clf()
                self.osl_display[:] = self.osl[self.propid]+self.node_state
                imshow_node_grid(self.grid, 'osl_display', limits=(0.0, 2.0),
                                 cmap=self.cmap_for_osl)
                show()
                figure(1)


    def finalize(self):

        # FINALIZE

        # Plot
        self.ca_plotter.finalize()



# If user runs this file, activate the main() function.
if __name__ == "__main__":

    # Parse command-line argument, if any
    import sys
    if len(sys.argv)>1:
        input_file_name = sys.argv[1]
    else:
        input_file_name = 'tsbm_inputs.txt'

    # Instantiate the model
    ca_model = TurbulentSuspensionAndBleachingModel(input_file_name)

    # Run the model
    ca_model.go()

    # Clean up
    ca_model.finalize()
