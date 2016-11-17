

from __future__ import print_function

from landlab.components.flow_director import FlowDirector
from landlab.utils.decorators import use_file_name_or_kwds

class FlowDirector_D4(FlowDirector):
    """
    """

    _name = 'FlowDirector_D4'

    _input_var_names = ('topographic__elevation')

    _output_var_names = ('flow__sink_flag'
                         )

    _var_units = {'flow__sink_flag': '-',
                  }

    _var_mapping = {'flow__sink_flag': 'node',
                    }

    _var_doc = {'flow__sink_flag': 'Boolean array, True at local lows',
    }
    
    # of _name, _input_var_names, _output_var_names, _var_units, _var_mapping, 
    # and _var_doc , all need to change. 
    @use_file_name_or_kwds
    def __init__(self, grid, surface='topographic_elevation', **kwds):
        super(FlowDirector_D4, self).__init__(grid, surface, **kwds)

        # save method as attribute
        self.method = 'D4'
        
        # load correct flow direction module 
    
    def run_one_step(self):
        print('testing yay')
        #self._grid['node']['flow__sink_flag'][:] = numpy.zeros_like(receiver,
        #                                                          dtype=bool)
        #self._grid['node']['flow__sink_flag'][sink] = True

        #return self._grid
        
if __name__ == '__main__':
    import doctest
    doctest.testmod()