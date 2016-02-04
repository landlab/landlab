

def setup_var_grid(var_map):
    return dict(enumerate(set(var_map.values())))


class BmiMixIn(object):
    def __init__(self, *args, **kwds):

        at = set(self._var_mapping.values())
        for name, at in self._var_mapping.items():
            self._var_grid[name] = 

    def initialize(self, conf):
        pass

    def update(self):
        pass

    def update_until(self, stop):
        pass

    def finalize(self):
        pass

    def get_component_name(self):
        return self._name

    def get_input_item_count(self):
        return len(self._input_var_names)

    def get_output_item_count(self):
        return len(self._input_var_names)

    def get_input_var_names(self):
        return tuple(self._input_var_names)

    def get_output_var_names(self):
        return tuple(self._output_var_names)

    def get_var_type(self, name):
        return self.get_value(name).dtype

    def get_var_units(self, name):
        self._var_units[name]

    def get_var_nbytes(self, name):
        return self.get_value_ptr(name).nbytes

    def get_var_itemsize(self, name):
        return self.get_value_ptr(name).itemsize

    def get_var_grid(self, name):
        return self._var_grid[name]

    def get_grid_type(self, grid):
        pass

    def get_grid_rank(self, grid):
        return 2

    def get_grid_size(self, grid):
        pass

    def get_grid_x(self, grid):
        pass

    def get_grid_y(self, grid):
        pass

    def get_grid_z(self, grid):
        pass

    def get_value(self, name, dest):
        dest[:] = self.get_value_ptr(name)[:]
        return dest

    def get_value_ptr(self, name):
        return self._grid[self._var_mapping[name]][name]

    def set_value(self, name, src):
        dest = self.get_value_ptr(name)
        dest[:] = src[:]
