#! /usr/bin/env python

class Component(object):
    _input_var_names = set()
    _output_var_names = set()
    _var_units = dict()

    def __init__(self, grid, map_vars={}):
        self._grid = grid

        for (dest, src) in map_vars.items():
            grid[dest] = grid[src]

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
    def shape(self):
        return self.grid._shape

    @property
    def grid(self):
        return self._grid

    @property
    def coords(self):
        return (self.grid.get_cell_x_coords(), self.grid.get_cell_y_coords())

    def imshow(self, name, **kwds):
        self._grid.imshow(name, **kwds)
