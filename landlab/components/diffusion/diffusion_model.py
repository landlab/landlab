#! /usr/env/python
"""

Landlab model of 2D diffusion, using DiffusionComponent.

Created July 2013 GT
Last updated August 2013 GT

"""

import sys                                    # for command-line arguments
from landlab import ModelParameterDictionary  # for input file
from landlab import create_and_initialize_grid
import diffusion

class DiffusionModel():
    
    def initialize(self, input_stream=None):
    
        # If no input file/stream specified, use default (or command line?)
        # (or default values?)
        if input_stream==None:
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
        
        # Create state variables
        self.z = self.grid.create_node_dvector()
        
    def update(self):
        
        self.diffusion_component.run_one_step(self.z)
        
    def run(self):
        
        self.diffusion_component.run_until(self.run_duration, self.z)
        
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
