#! /usr/env/python
"""

Component that models 2D diffusion using an explicit finite-volume method.

Created July 2013 GT
Last updated July 2013 GT

"""

from landlab import ModelParameterDictionary

class DiffusionComponent():
    
    def __init__(self, grid=None):
        
        self.grid = grid
        
    def initialize(self, input_stream):
        
        # Create a ModelParameterDictionary for the inputs
        inputs = ModelParameterDictionary(input_stream)
        
        # Read input/configuration parameters
        self.kd = inputs.get('KD', ptype=float)
        
        # Create grid if one doesn't already exist
        if self.grid==None:
            pass