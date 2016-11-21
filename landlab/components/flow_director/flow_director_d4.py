from __future__ import print_function

from landlab import FieldError
from landlab.components.flow_director import FlowDirectorToOne
from landlab.utils.decorators import use_field_name_or_array
import numpy

class FlowDirectorD4(FlowDirectorToOne):
    """
    """
    # of _name, _input_var_names, _output_var_names, _var_units, _var_mapping, 
    # and _var_doc , all need to change. 
    
    _name = 'FlowDirectorD4'



    @use_field_name_or_array
    def __init__(self, grid, surface='topographic_elevation'):
        super(FlowDirectorD4, self).__init__(grid, surface)
        
        # load 
       
    
    def run_one_step(self):
        test=1
        
if __name__ == '__main__':
    import doctest
    doctest.testmod()