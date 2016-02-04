import numpy as np
import yaml

from ..core.model_component import Component
from ..grid import RasterModelGrid


class TimeStepper(object):
    def __init__(self, start=0., stop=None, step=1.):
        self._start = start
        self._stop = stop
        self._step = step

        self._time = start

    @property
    def time(self):
        return self._time

    @property
    def start(self):
        return self._start

    @property
    def stop(self):
        return self._stop

    @property
    def step(self):
        return self._step

    @step.setter
    def step(self, new_val):
        self._step = new_val

    def advance(self):
        self._time += self.step
        if self._time > self._stop:
            raise StopIteration()


def wrap_as_bmi(cls):
    if not issubclass(cls, Component):
        raise TypeError('class must inherit from Component')

    class BmiWrapper(object):
        __doc__ = cls.__doc__
        _cls = cls
        def __init__(self):
            self._base = None
            self._clock = None
            super(BmiWrapper, self).__init__()

        def get_component_name(self):
            return self._cls.name

        def get_input_var_names(self):
            return self._cls.input_var_names

        def get_output_var_names(self):
            return self._cls.output_var_names

        def get_current_time(self):
            return self._clock.time

        def get_end_time(self):
            return self._clock.stop

        def get_start_time(self):
            return self._clock.start

        def get_time_step(self):
            return self._clock.step

        def get_time_units(self):
            raise NotImplementedError('get_time_units not implemented')

        def initialize(self, fname):
            with open(fname, 'r') as fp:
                params = yaml.load(fp)

            grid_params = params.pop('grid')
            grid = RasterModelGrid(grid_params['shape'],
                                   spacing=grid_params.get('spacing',
                                                           (1., 1.)),
                                   origin=grid_params.get('origin', (0., 0.)))
            clock_params = params.pop('clock')
            self._clock = TimeStepper(**clock_params)

            self._base = self._cls(grid, **params)

        def update(self):
            self._clock.advance()

        def update_frac(self, frac):
            time_step = self.get_time_step()
            self._clock.step = time_step * frac
            self.update()
            self._clock.step = time_step

        def update_until(self, then):
            n_steps = (then - self.get_current_time()) / self.get_time_step()
            for _ in xrange(int(n_steps)):
                self.update()
            self.update_frac(n_steps - int(n_steps))

        def finalize(self):
            pass

        def get_var_grid(self, name):
            return 0

        def get_var_itemsize(self, name):
            return np.dtype('float').itemsize

        def get_var_nbytes(self, name):
            return self.get_itemsize(name) * self._base.grid.number_of_nodes

        def get_var_type(self, name):
            return str(np.dtype('float'))

        def get_var_units(self, name):
            return self._cls.var_units(name)

        def get_value_ref(self, name):
            return self._base.grid.at_node[name]

        def get_value(self, name):
            return self._base.grid.at_node[name].copy()

        def set_value(self, name, vals):
            if name in self.get_input_var_names():
                if name in self._base.grid.at_node:
                    self._base.grid.at_node[name][:] = vals.flat
                else:
                    self._base.grid.at_node[name] = vals
            else:
                raise KeyError('{name} is not an input item'.format(name=name))

        def get_grid_origin(self, gid):
            return (self._base.grid.node_y[0], self._base.grid.node_x[0])

        def get_grid_rank(self, gid):
            return 2

        def get_grid_shape(self, gid):
            return (self._base.grid.number_of_node_rows,
                    self._base.grid.number_of_node_columns)

        def get_grid_spacing(self, gid):
            return (self._base.grid.dy, self._base.grid.dx)

        def get_grid_type(self, gid):
            return 'uniform_rectilinear'


    BmiWrapper.__name__ = cls.__name__
    return BmiWrapper
