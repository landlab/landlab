# -*- coding: utf-8 -*-
"""
Created on Sat May 30 14:01:10 2015

@author: gtucker
"""

from landlab import ModelParameterDictionary, Component, FieldError, FIXED_VALUE_BOUNDARY
from landlab.core.model_parameter_dictionary import MissingKeyError


class DepressionFinderAndRouter(Component):
    
    _name = 'DepressionFinderAndRouter'
    
    _input_var_names = set(['topographic__elevation',
                            ])
    
    _output_var_names = set(['depression__depth',  # depth below spill point
                             'depression__outlet_node_ID',
                             ])
                             
    _var_units = {'depression__depth' : 'm',
                  'depression__outlet_node_ID' : '-'
                  }
    
    _var_mapping = {'depression__depth' : 'node',
                    'depression__outlet_node_ID' : 'node'
                    }
    
    _var_defs = {'topographic__elevation' : 'Surface topographic elevation',
                 'depression__depth' : 'Depth of depression below its spillway point'
                  }
    
    def __init__(self, grid, input_stream=None, current_time=0.):
        self._grid = grid
        self.current_time = current_time
        self.initialize(input_stream)

    def initialize(self, input_stream=None):
        
        # Create a ModelParameterDictionary for the inputs
        if input_stream is None:
            inputs = None
        elif type(input_stream)==ModelParameterDictionary:
            inputs = input_stream
        else:
            inputs = ModelParameterDictionary(input_stream)
        
        # Make sure the grid includes elevation data. This means either:
        #  1. The grid has a node field called 'topographic__elevation', or
        #  2. The input file has an item called 'ELEVATION_FIELD_NAME' *and*
        #     a field by this name exists in the grid.
        try:
            self._elev = self._grid.at_node['topographic__elevation']
        except FieldError:
            try:
                print 'type(inputs)=',type(inputs)
                topo_field_name = inputs.read_string('ELEVATION_FIELD_NAME')
            except AttributeError:
                print 'Error: Because your grid does not have a node field called'
                print '"topographic__elevation", you need to pass the name of'
                print 'a text input file or ModelParameterDictionary, and this'
                print 'file or dictionary needs to include the name of another'
                print 'field in your grid that contains your elevation data.'
                raise AttributeError
            except MissingKeyError:
                print 'Error: Because your grid does not have a node field called'
                print '"topographic__elevation", your input file (or'
                print 'ModelParameterDictionary) must include an entry with the'
                print 'key "ELEVATION_FIELD_NAME", which gives the name of a'
                print 'field in your grid that contains your elevation data.'
                raise MissingKeyError('ELEVATION_FIELD_NAME')
            try:
                self._elev = self._grid.at_node[topo_field_name]
            except AttributeError:
                print 'Your grid does not seem to include a node field called',topo_field_name
                
                
    def find_pits(self):
        
        self.is_pit = self._grid.add_ones('node', 'is_pit', dtype=bool)
        self.is_pit[self._grid.boundary_nodes] = False
        # For each core node in the grid
        for link in self._grid.active_links:
            h = self._grid.node_index_at_link_head[link]
            t = self._grid.node_index_at_link_tail[link]
            if self._elev[h] > self._elev[t]:
                self.is_pit[h] = False
            elif self._elev[t] > self._elev[h]:
                self.is_pit[t] = False
            elif self._elev[h] == self._elev[t]:
                if self._grid.node_boundary_status[h]==FIXED_VALUE_BOUNDARY:
                    self.is_pit[t] = False
                elif self._grid.node_boundary_status[t]==FIXED_VALUE_BOUNDARY:
                    self.is_pit[h] = False
        print self.is_pit
  
                
    def map_depressions(self):
        
        self.find_pits()
    

def main():
    """
    temporary: sketch out the basic algorithm; later we'll divide this component-style.
    """
    print 'howdy'
    from landlab import RasterModelGrid
    from numpy.random import rand
    grid = RasterModelGrid(4, 5, 1.0)
    z = grid.add_zeros('node', 'topographic__elevation')
    z[:] = rand(len(z))*100
    print z
    dep_finder = DepressionFinderAndRouter(grid)
    dep_finder.initialize('/Users/gtucker/Dev/Landlab/gt_tests/test_inputs_for_depression_mapper.txt')
    dep_finder.map_depressions()
    
    n=0
    for r in range(grid.number_of_node_rows):
        for c in range(grid.number_of_node_columns):
            print int(z[n]),'(',dep_finder.is_pit[n],')',
            n+=1
        print
    
    
if __name__=='__main__':
    main()

    