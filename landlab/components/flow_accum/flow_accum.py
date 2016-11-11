#!/usr/env/python

from landlab import (ModelParameterDictionary, Component, FieldError,
                     FIXED_VALUE_BOUNDARY, CLOSED_BOUNDARY, CORE_NODE)
from six.moves import range
from .cfuncs import _add_to_stack
import numpy

class FlowAccumulation_ToOne(Component):

    """

    """
    _name = 'FlowAccumulation'

    _input_var_names = ('flow__receiver_node',
                        'runoff_rate',
                        )

    _output_var_names = ('drainage_area',
                         'surface_water__discharge',
                         'flow__upstream_node_order'
                         )

    _var_units = {'flow__receiver_node': 'm',
                  'runoff_rate': 'm/s',
                  'drainage_area': 'm**2',
                  'surface_water__discharge': 'm**3/s',
                  'flow__upstream_node_order':'-'
                  }

    _var_mapping = {'flow__receiver_node': 'node',
                    'runoff_rate': 'node',
                    'drainage_area': 'node',
                    'surface_water__discharge': 'node',
                    'flow__upstream_node_order': 'node'
                    }
    _var_doc = {
        'flow__receiver_node':
            'Node array of receivers (node that receives flow from current '
            'node)',
        'drainage_area':
            "Upstream accumulated surface area contributing to the node's "
            "discharge",
        'surface_water__discharge': 'Discharge of water through each node',
        'flow__upstream_node_order':
            'Node array containing downstream-to-upstream ordered list of '
            'node IDs',
        'runoff_rate':
            'Local runoff rate at each cell.'
            }
    def __init__(self, grid):
        """

        """
        self._grid = grid
        self._bc_set_code = self.grid.bc_set_code


        # Put all the other init checking things here.










        def run_one_step(self, **kwds):
        """Accumulate surface-water flow over a landscape.

        Examples
        --------

        """
        self.flow_accumulation(**kwds)

        @property
        def node_drainage_area(self):
        return self._grid['node']['drainage_area']

        @property
        def node_water_discharge(self):
        return self._grid['node']['surface_water__discharge']

        @property
        def node_order_upstream(self):
        return self._grid['node']['flow__upstream_node_order']


    if __name__ == '__main__':
        import doctest
        doctest.testmod()
