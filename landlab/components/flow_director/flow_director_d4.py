from landlab.components.flow_director import FlowDirectorToOne
from landlab.components.flow_routing import flow_direction_DN
from landlab import FIXED_VALUE_BOUNDARY, FIXED_GRADIENT_BOUNDARY
import numpy

class FlowDirectorD4(FlowDirectorToOne):
    """
    """

    _name = 'FlowDirectorD4'

    def __init__(self, grid, surface='topographic_elevation'):
        super(FlowDirectorD4, self).__init__(grid, surface)
        self.method = 'D4'
       
    def run_one_step(self):   
        
        # step 0. Check and update BCs
        if self._bc_set_code != self.grid.bc_set_code:
            self.updated_boundary_conditions()
            self._bc_set_code = self.grid.bc_set_code           
        
        # step 1. Calculate link slopes. 
        link_slope = - self._grid.calc_grad_of_active_link(
                self.elevs)
                
        # Step 2. Find and save base level nodes. 
        (self.baselevel_nodes, ) = numpy.where(
            numpy.logical_or(self._grid.status_at_node == FIXED_VALUE_BOUNDARY,
                             self._grid.status_at_node == FIXED_GRADIENT_BOUNDARY))
                             
        # Calculate flow directions
        num_d4_active = self._grid.number_of_active_links  
        receiver, steepest_slope, sink, recvr_link = \
        flow_direction_DN.flow_directions(self.elevs, self._active_links,
                                         self._activelink_tail[:num_d4_active],
                                         self._activelink_head[:num_d4_active],
                                         link_slope,
                                         grid=self._grid,
                                         baselevel_nodes=self.baselevel_nodes)
                                         
       
       # Save the four ouputs of this component.                                  
        self._grid['node']['flow__receiver_node'][:] = receiver
        self._grid['node']['topographic__steepest_slope'][:] = steepest_slope
        self._grid['node']['flow__link_to_receiver_node'][:] = recvr_link
        self._grid['node']['flow__sink_flag'][:] = numpy.zeros_like(receiver,
                                                                    dtype=bool)
        self._grid['node']['flow__sink_flag'][sink] = True
    
        
if __name__ == '__main__':
    import doctest
    doctest.testmod()