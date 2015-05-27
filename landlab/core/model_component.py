#! /usr/bin/env python
from landlab import ModelParameterDictionary as mpd
from landlab import MissingKeyError

class Component(object):
    _input_var_names = set()
    _output_var_names = set()
    _var_units = dict()

    def __init__(self, grid, map_vars=None, **kwds):
        #ISSUES:
        #params must be named 'input_params', or use a dict containing all poss names
        #It's not obvious the var_names exist to the novice coder
        #...?
        map_vars = map_vars or {}
        self._grid = grid

        for (location, vars) in map_vars.items():
            for (dest, src) in vars.items():
                grid.add_field(location, dest,
                               grid.field_values(location, src))

        self.initialize(grid, **kwds)
        self._var_names = {}
        _var_set = self._input_var_names.union(self.output_var_names)
        
        #attempt to build a dict translating user defined field names
        #to param determined field names:
        try:
            MPD_in = mpd(kwds['input_params'])
        except KeyError:
            for i in _var_set:
                self._var_names[i] = i
        else:
            for i in _var_set:
                try:
                    new_key = MPD_in.read_string(i)
                except MissingKeyError:
                    self._var_names[i] = i
                else:
                    self._var_names[new_key] = i
                    
        
        
        

    @property
    def input_var_names(self):
        return self._input_var_names

    @property
    def output_var_names(self):
        return self._output_var_names
        
    @property
    def var_names(self):
        return self._var_names

    @property
    def name(self):
        return self._name

    @property
    def units(self):
        return self._var_units
    
    @property
    def var_units(self):
        return self._var_units
        
    @property
    def var_definitions(self):
        return self._var_defs
    
    @property
    def var_mapping(self):
        """var_mapping
        This is 'node', 'cell', 'active_link', etc.
        """
        return self._var_mapping

    @property
    def shape(self):
        return self.grid._shape

    @property
    def grid(self):
        return self._grid

    @property
    def coords(self):
        return (self.grid.node_x, self.grid.node_y)

    def imshow(self, name, **kwds):
        self._grid.imshow(name, **kwds)
