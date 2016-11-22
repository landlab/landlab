from landlab.components.flow_accum.flow_accumulator import FlowAccumulator 
from landlab.components.flow_director import FlowDirectorD4 as FlowDirector
from landlab.components.flow_accum import flow_accum_bw 

class FlowAccumulatorD4(FlowAccumulator):
    """ 
    Info here
    """
    
    _name = 'FlowAccumulatorD4'
    
    # of _name, _input_var_names, _output_var_names, _var_units, _var_mapping, 
    # and _var_doc , only _name needs to change. 
    
    def __init__(self, grid, surface='topographic__elevation', depression_finder=None):
        super(FlowAccumulatorD4, self).__init__(grid, surface)

        # save method as attribute
        self.method = 'D4'
        
        # save 
        self.fd = FlowDirector(self._grid, self.elevs)
        self.df = depression_finder
        
    def run_one_step(self):
        
        # step 0. Check and update BCs
        if self._bc_set_code != self.grid.bc_set_code:
            self.updated_boundary_conditions()
            self._bc_set_code = self.grid.bc_set_code        
        
        # step 1. Find flow directions by specified method
        self.fd.run_one_step()
        
        # step 2. Get r (and potentially p) array(s)        
        r = self._grid['node']['flow__receiver_node']
        s = self.fd.baselevel_nodes       
        
        # step 2. Stack, D, delta construction
        nd = flow_accum_bw._make_number_of_donors_array(r)
        delta = flow_accum_bw._make_delta_array(nd)
        D = flow_accum_bw._make_array_of_donors(r, delta)
        s = flow_accum_bw.make_ordered_node_array(r, s)
        
        #put theese in grid so that depression finder can use it.         
        # store the generated data in the grid
        self._grid['node']['flow__data_structure_delta'][:] = delta
        self._grid['node']['flow__data_structure_D'][:] = D
        self._grid['node']['flow__upstream_node_order'][:] = s
        
        # step 3. Initialize and Run depression finder if passed 
        if self.df:
            df=self.df(self.grid)
            df.run_one_step()
        
        # step 4. Accumulate (to one or to N depending on direction method. )
        a, q = flow_accum_bw.find_drainage_area_and_discharge(self._grid['node']['flow__upstream_node_order'], 
                                                              self._grid['node']['flow__receiver_node'], 
                                                              node_cell_area=self.node_cell_area,
                                                              runoff_rate=self._grid.at_node['water__unit_flux_in'])
        
        self._grid['node']['drainage_area'][:] = a
        self._grid['node']['surface_water__discharge'][:] = q
        

    
if __name__ == '__main__':
    import doctest
    doctest.testmod()    