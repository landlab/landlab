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
        self.kd = inputs.get('DIFMOD_KD', ptype=float)
        self.uplift_rate = inputs.get('DIFMOD_UPLIFT_RATE', 0., ptype=float)
        
        # Create grid if one doesn't already exist
        if self.grid==None:
            pass   # to be implemented