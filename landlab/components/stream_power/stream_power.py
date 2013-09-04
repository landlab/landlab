import numpy as np

class StreamPower(object):
    """
    UNDER DEVELOPMENT. DO NOT ATTEMPT TO USE YET. DEJH Sept 2013.
    """
    
    def __init__(self, grid, data):
        self.m = 0.5
        self.n = 1.
        self.K = 1.e-5
        #This will us the MPD once finalized
        
    def self.initialize(self, grid, data):
        self.__init__(grid, data)
        
    def stream_power_erosion(self, grid, data):
        grads_at_links = grid.calculate_gradients_at_active_links(data.elev) #Make this!