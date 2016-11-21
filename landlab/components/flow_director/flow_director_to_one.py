from __future__ import print_function

from landlab import FieldError
from landlab.components.flow_director import FlowDirector
from landlab.utils.decorators import use_field_name_or_array
import numpy

class FlowDirectorToOne(FlowDirector):
    """
    """
    # of _name, _input_var_names, _output_var_names, _var_units, _var_mapping, 
    # and _var_doc , all need to change. 
    
    _name = 'FlowDirectorToOne'

    _input_var_names = ('topographic__elevation',
                        )

    _output_var_names = ('flow__receiver_node',
                         'topographic__steepest_slope',
                         'flow__upstream_node_order',
                         'flow__sink_flag',
                         )

    _var_units = {'topographic__elevation': 'm',
                  'flow__receiver_node': '-',
                  'topographic__steepest_slope': '-',
                  'flow__upstream_node_order': '-',
                  'flow__link_to_receiver_node': '-',
                  'flow__sink_flag': '-',
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'flow__receiver_node': 'node',
                    'topographic__steepest_slope': 'node',
                    'flow__link_to_receiver_node': 'node',
                    'flow__sink_flag': 'node',
                    }

    _var_doc = {
        'topographic__elevation': 'Land surface topographic elevation',
        'flow__receiver_node':
            'Node array of receivers (node that receives flow from current '
            'node)',
        'topographic__steepest_slope':
            'Node array of steepest *downhill* slopes',
        'surface_water__discharge': 'Discharge of water through each node',
        'flow__link_to_receiver_node':
            'ID of link downstream of each node, which carries the discharge',
        'flow__sink_flag': 'Boolean array, True at local lows',
    }
    

    @use_field_name_or_array
    def __init__(self, grid, surface):
        super(FlowDirectorToOne, self).__init__(grid, surface)
        
        # initialize new fields
        
        try:
            self.receiver = grid.add_zeros('flow__receiver_node', at='node',
                                           dtype=int)
        except FieldError:
            self.receiver = grid.at_node['flow__receiver_node']
        
        try:
            self.steepest_slope = grid.add_zeros(
                'topographic__steepest_slope', at='node', dtype=float)
        except FieldError:
            self.steepest_slope = grid.at_node['topographic__steepest_slope']
        
        try:
            self.links_to_receiver = grid.add_zeros(
                'flow__link_to_receiver_node', at='node', dtype=int)
        except FieldError:
            self.links_to_receiver = grid.at_node[
                'flow__link_to_receiver_node']
        
        grid.add_zeros('flow__sink_flag', at='node', dtype=numpy.int8,
                       noclobber=False)
       
    
    def run_one_step(self):
        raise NotImplementedError('run_one_step()')
        
if __name__ == '__main__':
    import doctest
    doctest.testmod()