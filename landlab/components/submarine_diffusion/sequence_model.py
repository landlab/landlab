# -*- coding: utf-8 -*-
"""
erosion_model.py: generic base class for an erosion model.

Created on Thu Dec 24 12:28:31 2015

@author: gtucker
"""

from landlab import RasterModelGrid
from landlab import CLOSED_BOUNDARY
from landlab.core import load_params
from landlab.io.netcdf import write_netcdf
import sys

_VERBOSE = True

class SequenceModel(object):
    """
    INSERT DESCRIPTION HERE...
    """

    def __init__(self, input_file=None, params=None):
        """
        Handle inputs and set params.
        """
        
        # Make sure user has given us an input file or parameter dictionary
        # (but not both)
        if input_file is None and params is None:
            print('You must specify either an input_file or params dict')
            sys.exit(1)
        if input_file is not None and params is not None:
            print('ErosionModel constructor takes EITHER')
            print('input_file or params, but not both.')
            sys.exit(1)

        # If we have an input file, let's read it
        ### CHANGE TO USE LOAD_PARAMS AND YAML
        if input_file is None:
            self.params = params  # assume params is a dictionary
        else:
            self.params = load_params(input_file)                

        # create a grid
        nrows = self.params['number_of_node_rows']
        ncols = self.params['number_of_node_columns']
        self.grid = RasterModelGrid((nrows, ncols), self.params['dx'])
        self.z = self.grid.add_zeros('node', 'topographic__elevation')

        # A little reporting
        if _VERBOSE:
            print('Grid has ' + str(self.grid.number_of_node_rows) +
                  ' rows and ' + str(self.grid.number_of_node_columns) +
                  ' columns, with spacing of ' +str(self.grid.dx))

        # Set grid boundaries: top and bottom closed, left and right open
        # (which we assume is the default)
        self.grid.status_at_node[self.grid.nodes_at_bottom_edge] = CLOSED_BOUNDARY
        self.grid.status_at_node[self.grid.nodes_at_top_edge] = CLOSED_BOUNDARY

        ### CALL INITIAL PROFILE FN AROUND HERE

    def write_output(self, params, iteration, field_names=None, silent=False):
        """Write output to file (currently netCDF)."""
        if not silent:
            print('Iteration ' + str(iteration))
        filename = self.params['output_filename'] + str(iteration).zfill(4) \
                    + '.nc'
        write_netcdf(filename, self.grid, names=field_names)
        ### SOMEHOW FACTOR IN LAYER OUTPUT HERE

    def run_one_step(self, dt):
        """
        Run each component for one time step.
        
        This base-class method does nothing. Derived classes should override
        it to run each component in turn for a time period dt.
        """
        ### HERE CALCULATE TRANSPORT, ISOSTASY, COMPACTION, AND ITERATE
        pass


    def run_for(self, dt, runtime):
        """
        Run model without interruption for a specified time period.
        """
        elapsed_time = 0.
        keep_running = True
        while keep_running:
            if elapsed_time+dt >= runtime:  # make sure we don't exceed runtime
                dt = runtime-elapsed_time
                keep_running = False
            self.run_one_step(dt)
            elapsed_time += dt


    def run(self, output_fields=None):
        """
        Run the model until complete.
        """
        total_run_duration = self.params['run_duration']
        output_interval = self.params['output_interval']
        iteration = 1
        time_now = 0.0
        while time_now < total_run_duration:
            next_run_pause = min(time_now + output_interval, total_run_duration)
            self.run_for(self.params['dt'], next_run_pause - time_now)
            time_now = next_run_pause
            self.write_output(self.params, iteration, field_names=output_fields)
            iteration += 1


    def finalize(self, params):
        """
        Finalize the model.
        """
        pass


def main():
    """Executes model."""
    try:
        infile = sys.argv[1]
    except IndexError:
        print('Must include input file name on command line')
        sys.exit(1)
        
    erosion_model = SequenceModel(input_file=infile)
    erosion_model.run()
    
    # MIGHT WANT TO CHANGE PARAMS in the middle of a run, like so ...
    # sequence_model.run_for(...)
    # sequence_model.params['rock_uplift__rate'] = 0.004
    # sequence_model.run_for(...)


if __name__ == '__main__':
    main()