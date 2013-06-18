#import numpy as np

from landlab import Implements, ImplementsOrRaise
from landlab.framework.interfaces import BmiBase, BmiNoGrid


#@ImplementsOrRaise(BmiBase)
@Implements(BmiBase)
class Sample1(object):
    """
    A sample component.
    """
    __implements__ = (BmiBase, BmiNoGrid, )

    _input_var_names = [
        'air__temperature',
        'surface__elevation',
    ]
    _output_var_names = [
        'deposition__rate',
    ]

    model_name = 'Sample 1'
    author_name = 'Eric Hutton'
    version = '0.1'
    time_units = 's'
    time_step_type = 'fixed'
    step_method = 'explicit'
    grid_type = 'none'

    _vars = {
        'deposition__rate': [1.]
    }

    def initialize(self, name):
        pass

    def update(self):
        pass

    def finalize(self):
        pass

    def get_input_var_names(self):
        return self._input_var_names

    def get_output_var_names(self):
        return self._output_var_names

    def get_var_rank(self, name):
        return 0

    def get_start_time(self):
        return 0.

    def get_current_time(self):
        return 0.

    def get_end_time(self):
        return 100.

    def get_time_step(self):
        return 1.

    def get_var_type(self, name):
        return 'float64'

    def get_var_units(self, name):
        return 'm'

    def set_value(self, name, value):
        pass

    def get_value(self, name):
        return self._vars[name]


@ImplementsOrRaise(BmiBase)
class Sample2(object):
    """
    A sample component.
    """
    __implements__ = (BmiBase, BmiNoGrid, )

    _input_var_names = [
        'deposition__rate',
    ]
    _output_var_names = [
        'air__temperature',
        'surface__elevation',
    ]

    model_name = 'Sample 2'
    author_name = 'Eric Hutton'
    version = '0.1'
    time_units = 's'
    time_step_type = 'fixed'
    step_method = 'explicit'
    grid_type = 'none'

    _vars = {
        'air__temperature': [1.],
        'surface__elevation': [1.],
    }

    def initialize(self, name):
        pass

    def update(self):
        pass

    def finalize(self):
        pass

    def get_input_var_names(self):
        return self._input_var_names

    def get_output_var_names(self):
        return self._output_var_names

    def get_var_rank(self, name):
        return 0

    def get_start_time(self):
        return 0.

    def get_current_time(self):
        return 0.

    def get_end_time(self):
        return 100.

    def get_time_step(self):
        return 1.

    def get_var_type(self, name):
        return 'float64'

    def get_var_units(self, name):
        return 'm'

    def get_value(self, name):
        return self._vars[name]

    def set_value(self, name, value):
        pass
