import numpy as np

class StreamPower(object):
    """
    UNDER DEVELOPMENT. DO NOT ATTEMPT TO USE YET. DEJH Sept 2013.
    """
    
    def __init__(self, grid, data, tstep):
        self.initialize(grid, data, tstep)
            
#This draws attention to a potential problem. It will be easy to have modules update z, but because noone "owns" the data, to forget to also update dz/dx...
#How about a built in grid utility that updates "derived" data (i.e., using only grid fns, e.g., slope, curvature) at the end of any given tstep loop?
#Or an explicit flagging system for all variables in the modelfield indicating if they have been updated this timestep. (Currently implemented)
#Or wipe the existance of any derived grid data at the end of a timestep entirely, so modules find they don't have it next timestep.
        
    def initialize(self, grid, data, tstep):
        self.m = 0.5
        self.n = 1.
        self.K = 1.e-3 * tstep
        self.unit_stream_power = grid.zeros(centering='node')
        #Flags for self-building of derived data:
        self.made_max_gradients = 0
        self.made_link_gradients = 0
        #This will us the MPD once finalized
        #Now perform checks for existance of needed data items:
        try:
            data.node_max_gradients
        except:
            try:
                data.node_max_gradients, dummy = grid.calculate_steepest_descent_on_nodes(data.elev, data.link_gradients)
                self.made_max_gradients = 1
            except:
                data.link_gradients = grid.calculate_gradients_at_active_links(data.elev)
                data.node_max_gradients, dummy = grid.calculate_steepest_descent_on_nodes(data.elev, data.link_gradients)
                self.made_max_gradients = 1
                self.made_link_gradients = 1
        try:
            data.flowacc.shape
        except:
            print 'Flow accumulation data is not available! Run will crash out.'
        
    def stream_power_erosion(self, grid, data):
        """
        """
        #Perform checks on self-build flags:
        if self.made_link_gradients:
            data.link_gradients = grid.calculate_gradients_at_active_links(data.elev)
        if self.made_max_gradients:
            data.node_max_gradients, dummy = grid.calculate_steepest_descent_on_nodes(data.elev, data.link_gradients)
        
        #Operate the main function:
        active_nodes = grid.get_active_cell_node_ids()
        node_max_gradients = data.node_max_gradients * np.where(data.node_max_gradients>0., 1., 0.) #Don't permit erosion on upward-pointing or flat surfaces
        unit_stream_power_active_nodes = data.flowacc[active_nodes]**self.m * node_max_gradients[active_nodes]**self.n
        self.unit_stream_power[active_nodes] = unit_stream_power_active_nodes
        
        data.elev -= (self.K * self.unit_stream_power)
        
        return data.elev
