
from landlab.components.flow_accum import FlowAccumulator 
from landlab.components.flow_director import FlowDirector_D4 as FlowDirector
#from landlab.components.flow_accum import flow_accum_bw 

class FlowAccumulator_D4(FlowAccumulator):
    """ 
    Info here
    """
    
    _name = 'FlowAccumulator_D4'
    
    # of _name, _input_var_names, _output_var_names, _var_units, _var_mapping, 
    # and _var_doc , only _name needs to change. 
    
    def __init__(self, grid, depression_finder=None, **kwds):
        super(FlowAccumulator_D4, self).__init__(**kwds)

        # save method as attribute
        self.method = 'D4'
        
        # need to pass surface to FlowDirection
        self.FlowDirector=FlowDirector()
        self.DepressionFinder=depression_finder
        
        
    def run_one_step(self):
        
        # step 1. Find flow directions by specified method
        self.FlowDirector.run_one_step()
        
        # step 2. Get r (and potentially p) array(s)        
        #delta = flow__data_structure_delta=flow_accum._get
        #D = delf.FlowDirection.r        
        
        
        # step 2. Stack, D, delta construction
        
       #p ut theese in grid so that d finder can use it.         
        # store the generated data in the grid
        
        elevs = self._grid['node']['topographic__elevation']

        node_cell_area = self._grid.cell_area_at_node.copy()
        node_cell_area[self._grid.closed_boundary_nodes] = 0.
        
        
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