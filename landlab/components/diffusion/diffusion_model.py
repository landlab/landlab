#! /usr/env/python
"""

Landlab model of 2D diffusion, using DiffusionComponent.

Created July 2013 GT
Last updated August 2013 GT

"""

import sys                                    # for command-line arguments
from landlab import ModelParameterDictionary  # for input file
from landlab import create_and_initialize_grid
from landlab import RasterModelGrid
from landlab.plot.imshow import imshow_node_grid
import diffusion

class DiffusionModel(object):
    
    def initialize(self, input_stream=None):
    
        # If no input file/stream specified, use default (or command line?)
        # (or default values?)
        if input_stream is None:
            input_stream = str(raw_input('Enter name of input file: '))
        
        # Open input file
        inputs = ModelParameterDictionary(input_stream)
    
        # Create grid
        self.grid = create_and_initialize_grid(inputs)
        
        # Create a diffusion component
        self.diffusion_component = diffusion.DiffusionComponent(self.grid)
        self.diffusion_component.initialize(input_stream)
        
        # Read parameters
        self.run_duration = inputs.get('RUN_DURATION', ptype=float)
        self.opt_netcdf_output = inputs.get('OPT_FILE_OUTPUT', ptype='bool')
        self.opt_display_output = inputs.get('OPT_DISPLAY_OUTPUT', ptype='bool')
        
        self.setup_output_timing(inputs)
        
        
    def setup_output_timing(self, inputs):
        
        # Setup output timing
        if self.opt_netcdf_output:
            self.netcdf_output_interval = inputs.get('FILE_OUTPUT_INTERVAL',
                                                     ptype=float)
            self.next_file_output = self.netcdf_output_interval
        else:
            self.next_file_output = self.run_duration+1.0
            
        if self.opt_display_output:
            self.display_output_interval = inputs.get('DISPLAY_OUTPUT_INTERVAL',
                                                      ptype=float)
            self.next_display_output = self.display_output_interval
        else:
            self.next_display_output = self.run_duration+1.0
            
        self.find_next_stop_time()
                    
        # Time
        self.current_time = 0.0
        
        
    def handle_output(self):
        
        if self.current_time >= self.next_display_output:
            self.do_display_output()
            self.next_display_output += self.display_output_interval
            
        if self.current_time >= self.next_file_output:
            self.do_file_output()
            self.next_file_output += self.file_output_interval
            
            
    def find_next_stop_time(self):
        
        self.next_stop_time = min(self.run_duration, self.next_file_output,
                                  self.next_display_output)
            
            
    def do_display_output(self):
        
        # ok, here one problem is that we don't know what grid type we have
        if type(self.grid) is RasterModelGrid:
            print 'This is a raster grid'
        else:
            print 'non-raster grid'
            
        imshow_node_grid(self.grid, self.diffusion_component.z)
        
        
    def do_file_output(self):
        
        print 'File output goes here'
        
        
    def update(self):
        
        self.diffusion_component.run_one_step()
        
        
    def run(self):
        
        while self.current_time < self.run_duration:
            self.diffusion_component.run_until(self.next_stop_time)
            self.current_time = self.next_stop_time
            print self.current_time
            self.handle_output()
            self.find_next_stop_time()
        
        
    def finalize(self):
        
        pass
        
        
def main():
    
    # Instantiate the model
    difmod = DiffusionModel()

    # Handle command-line arguments (if any)
    if len(sys.argv)==1:
        input_file_name = None
    else:
        input_file_name = sys.argv[1]
    
    # Initialize the model
    difmod.initialize(input_file_name)

    # Run the model
    difmod.run()
    
    # Finalize
    difmod.finalize()
    


if __name__ == "__main__":
    main()
    
    
# Issues arising (for group discussion):
#  - how to get argument for input source to model: command-line argument?
#    command-line prompt? window??? is there a way to toggle between file and
#    command-line input? How best to be compatible with CMF/CMT?
#  - naming conventions for components, models, etc.
#  - How to design for maximum simplicity?
#  - How to get parameters/data passed between functions? Use class with data
#    members? Use big Data object?
#  - We need a way to handle defaults in MPD
#  - We need a way to create and return a grid whose subclass is determined by
#    the input file specs, not hardwired
#  - Better error handling for MPD
#  - Exception handling generally!

# Sketch of a class-based model:
