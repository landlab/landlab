
from landlab.components.flow_accum import FlowAccumulator 
from landlab.components.flow_director import FlowDirectorD4 as FlowDirector
#from landlab.components.flow_accum import flow_accum_bw 

class FlowAccumulatorD4(FlowAccumulator):
    """ 
    Info here
    """
    
    _name = 'FlowAccumulator_D4'
    
    # of _name, _input_var_names, _output_var_names, _var_units, _var_mapping, 
    # and _var_doc , only _name needs to change. 
    
    def __init__(self, grid, surface='topographic__elevation', depression_finder=None):
        super(FlowAccumulatorD4, self).__init__(grid, surface)

        # save method as attribute
        self.method = 'D4'
        
        # save 
        self.flow_director=FlowDirector()
        self.depression_finder=depression_finder
        
        
    def run_one_step(self):
        
        # step 1. Find flow directions by specified method
        self.FlowDirector.run_one_step()
        
        # step 2. Get r (and potentially p) array(s)        
        #delta = flow__data_structure_delta=flow_accum._get
        #D = delf.FlowDirection.r        
        
        
        # step 2. Stack, D, delta construction
        
       #p ut theese in grid so that d finder can use it.         
        # store the generated data in the grid
        
        
        
        self._grid['node']['drainage_area'][:] = a
        self._grid['node']['flow__receiver_node'][:] = receiver
        self._grid['node']['topographic__steepest_slope'][:] = steepest_slope
        self._grid['node']['surface_water__discharge'][:] = q
        self._grid['node']['flow__upstream_node_order'][:] = s
        self._grid['node']['flow__link_to_receiver_node'][:] = recvr_link
        
        

        # step 3. Depression finder/router if called for
        
        
        # step 4. Accumulate (to one or to N depending on direction method. )
        
if __name__ == '__main__':
    import doctest
    doctest.testmod()    