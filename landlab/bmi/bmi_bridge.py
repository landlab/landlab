"""
========================================================================================
Wrap landlab component with the Basic Modeling Interface (:mod:`landlab.bmi.bmi_bridge`)
========================================================================================

.. sectionauthor:: Eric Hutton

Function reference
------------------

The `wrap_as_bmi` function wraps a landlab component class so that it
exposes a Basic Modelling Interface.

"""
import os

import numpy as np
import yaml

from ..core.model_component import Component
from ..grid import RasterModelGrid

__all__ = ['TimeStepper', 'wrap_as_bmi']


class TimeStepper(object):

    """Step through time.

    Parameters
    ----------
    start : float, optional
        Clock start time.
    stop : float, optional
        Stop time.
    step : float, optional
        Time step.

    Examples
    --------
    >>> from landlab.bmi import TimeStepper
    >>> time_stepper = TimeStepper()
    >>> time_stepper.start
    0.0
    >>> time_stepper.stop is None
    True
    >>> time_stepper.step
    1.0
    >>> time_stepper.time
    0.0
    >>> for _ in range(10): time_stepper.advance()
    >>> time_stepper.time
    10.0
    >>> time_stepper = TimeStepper(1., 13., 2.)
    >>> [time for time in time_stepper]
    [1.0, 3.0, 5.0, 7.0, 9.0, 11.0]
    """

    def __init__(self, start=0., stop=None, step=1.):
        self._start = start
        self._stop = stop
        self._step = step

        self._time = start

    def __iter__(self):
        if self.stop is None:
            while 1:
                yield self._time
                self._time += self._step
        else:
            while self._time < self._stop:
                yield self._time
                self._time += self._step
        raise StopIteration()

    @property
    def time(self):
        """Current time."""
        return self._time

    @property
    def start(self):
        """Start time."""
        return self._start

    @property
    def stop(self):
        """Stop time."""
        return self._stop

    @property
    def step(self):
        """Time Step."""
        return self._step

    @step.setter
    def step(self, new_val):
        """Change the time step."""
        self._step = new_val

    def advance(self):
        """Advance the time stepper by one time step."""
        self._time += self.step
        if self._stop is not None and self._time > self._stop:
            raise StopIteration()


def wrap_as_bmi(cls):
    """Wrap a landlab class so it exposes a BMI.

    Parameters
    ----------
    cls : class
        A landlab class that inherits from `Component`.

    Returns
    -------
    class
        A wrapped class that exposes a BMI.

    Examples
    --------
    >>> from landlab.bmi import wrap_as_bmi
    >>> from landlab.components.flexure import Flexure

    >>> BmiFlexure = wrap_as_bmi(Flexure)
    >>> flexure = BmiFlexure()

    >>> config = \"\"\"
    ... eet: 10e3
    ... method: flexure
    ... clock:
    ...     start: 0.
    ...     stop: 10.
    ...     step: 2.
    ... grid:
    ...     type: raster
    ...     shape: [20, 40]
    ...     spacing: [1000., 2000.]
    ... \"\"\"
    >>> flexure.initialize(config)
    >>> flexure.get_output_var_names()
    ('lithosphere__elevation_increment',)
    >>> flexure.get_var_grid('lithosphere__elevation_increment')
    0
    >>> flexure.get_grid_shape(0)
    (20, 40)
    >>> dz = flexure.get_value('lithosphere__elevation_increment')
    >>> dz.shape
    (800,)
    >>> np.all(dz == 0.)
    True
    >>> flexure.get_current_time()
    0.0

    >>> flexure.get_input_var_names()
    ('lithosphere__overlying_pressure_increment',)
    >>> load = np.zeros((20, 40), dtype=float)
    >>> load[0, 0] = 1.
    >>> flexure.set_value('lithosphere__overlying_pressure_increment', load)
    >>> flexure.update()
    >>> flexure.get_current_time()
    2.0
    >>> dz = flexure.get_value('lithosphere__elevation_increment')
    >>> np.all(dz == 0.)
    False
    """
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
            if os.path.isfile(fname):
                with open(fname, 'r') as fp:
                    params = yaml.load(fp)
            else:
                params = yaml.load(fname)

            grid_params = params.pop('grid')
            gtype = grid_params.pop('type')
            if gtype == 'raster':
                cls = RasterModelGrid
            else:
                raise ValueError(
                    'unrecognized grid type {gtype}'.format(gtype=gtype))

            grid = cls.from_dict(grid_params)

            clock_params = params.pop('clock')
            self._clock = TimeStepper(**clock_params)

            self._base = self._cls(grid, **params)

        def update(self):
            if hasattr(self._base, 'update'):
                self._base.update()
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
