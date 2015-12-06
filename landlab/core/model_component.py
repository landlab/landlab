#! /usr/bin/env python


class Component(object):
    _input_var_names = set()
    _output_var_names = set()
    _var_units = dict()

    def __init__(self, grid, map_vars=None):
        map_vars = map_vars or {}
        self._grid = grid

        for (location, vars) in map_vars.items():
            for (dest, src) in vars.items():
                grid.add_field(location, dest,
                               grid.field_values(location, src))

    @property
    def input_var_names(self):
        return self._input_var_names

    @property
    def output_var_names(self):
        return self._output_var_names

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
        return self._var_doc

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
