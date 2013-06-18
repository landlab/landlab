#! /usr/bin/env python

import numpy as np
from scipy import ndimage

from landlab import Implements
from landlab.framework.interfaces import (BmiBase, BmiUniformRectilinear,
                                          BadVarNameError, GRID_TYPE_UNIFORM)


#@ImplementsOrRaise(BmiBase, BmiUniformRectilinear)
@Implements(BmiBase, BmiUniformRectilinear)
class Diff2D(object): #pylint:disable-msg=R0902,R0904
    """
    Basic Modeling Interface implementation for a 2d diffusion process.
    """

    __implements__ = (BmiBase, BmiUniformRectilinear, )

    var_units = {'surface__elevation': 'meters'}
    name = '2D Diffusion'

    _input_var_names = [
        'surface__elevation',
    ]

    _output_var_names = [
        'surface__elevation',
    ]

    def __init__(self):
        self._time_step = 0.
        self._shape = (0, 0)
        self._spacing = (0., 0.)
        self._origin = (0., 0.)
        self._start_time = 0.
        self._end_time = np.finfo ('d').max
        self._time = self._start_time

        self._elevation = None
        self._temp_elevation = None

        self._var_names = set(self.output_vars + self.input_vars)
        self._value = {}

        self._stencil = np.array((3, 3), dtype=np.float64)

    @staticmethod
    def _set_bc(data):
        """
        Set boundary conditions for 2D data.

        Values on left, right, and top edge are 0. Values on the bottom edge
        are a parabola centered on the edge with roots on the left and right
        boundaries.

        :data: 2D numpy array
        """
        data[:, 0] = 0.
        data[:, -1] = 0.
        data[0, :] = 0.

        top_column = data.shape[1] - 1
        data[-1, :] = .25*top_column**2 - (np.arange(top_column + 1) -
                                           .5*top_column)**2

    def initialize(self, file_name=None):
        """
        Initialve the model from a file, or use default values. The
        initialization file is one line with values,

        ${TIME_STEP}, ${NUMBER_OF_ROWS}, ${NUMBER_OF_COLUMNS}

        :keyword file_name: Name of initialization file, or None for defaults
        """
        if file_name is not None:
            with open(file_name, 'r') as init_file:
                vals = init_file.readline().split(',')
            self._time_step, row_count, column_count = (
                float(vals[0]),
                int(vals[1]),
                int(vals[2]),
            )
            self._shape = (row_count, column_count)
        else:
            self._time_step = 1.
            self._shape = (20, 10)

        self._spacing = (1., 1.)
        self._origin = (0., 0.)

        self._elevation = np.zeros(self._shape)
        self._set_bc(self._elevation)
        self._temp_elevation = self._elevation.copy()

        self._value['surface__elevation'] = self._elevation

        row_spacing, column_spacing = self._spacing
        self._stencil = np.array(
            [[0., column_spacing**2, 0.],
             [row_spacing**2, 0., row_spacing**2],
             [0., column_spacing**2, 0.]]
        ) / (2. * (row_spacing**2 + column_spacing**2)) * self._time_step

    def update(self):
        """
        Update the model for a single time step.
        """
        ndimage.convolve(self._elevation, self._stencil,
                         output=self._temp_elevation)
        self._set_bc(self._temp_elevation)
        self._elevation[:] = self._temp_elevation

        self._time += self._time_step

    def update_until(self, stop_time):
        """
        Update the model until a given time.

        :stop_time: Update model until this time
        """
        while 1:
            next_time = self.get_current_time() + self._time_step
            if next_time > stop_time:
                break
            self.update()

        partial_step = stop_time - self.get_current_time()
        if partial_step > 0:
            self._stencil *= partial_step / self._time_step
            self.update()
            self._stencil *= self._time_step / partial_step

    def finalize(self):
        """
        Reset model values, clean up resources used.
        """
        self._time_step = 0.
        self._time = 0.

        self._elevation = np.array([])
        self._temp_elevation = np.array([])

    def get_var_type(self, long_var_name):
        """
        Type of a model variable.

        :long_var_name: Standard name string

        :return: Type as a numpy.dtype string
        """
        return str (self.get_value (long_var_name).dtype)

    def get_var_units(self, long_var_name):
        """
        Units of a model variable.

        :long_var_name: Standard name string

        :returns: Units as a udunits string
        """
        return self.var_units[long_var_name]

    def get_var_rank(self, long_var_name):
        """
        Rank of a model variable.

        :long_var_name: Standard name string

        :returns: Number of dimension
        """
        if long_var_name in self._var_names:
            return len(self._spacing)
        else:
            raise BadVarNameError(long_var_name)

    def get_value(self, long_var_name):
        """
        Values of a model variable.

        :long_var_name: Standard name string

        :returns: Numpy ND array of values
        """
        return self._value[long_var_name]

    def get_value_at_indices(self, long_var_name, indices):
        """
        Values of a model variable.

        :long_var_name: Standard name string
        :indices: Indices of (flattened) value array

        :returns: Numpy ND array of values
        """
        return self.get_value(long_var_name)[indices]

    def set_value(self, long_var_name, src):
        """
        Set values of a model variable.

        :long_var_name: Standard name string
        :src: Array of values to set
        """
        val = self.get_value(long_var_name)
        val[:] = src

    def set_value_at_indices(self, long_var_name, indices, src):
        """
        Set values of a model variable at given indices.

        :long_var_name: Standard name string
        :indices: Indices of (flattened) value array
        :src: Array of values to set
        """
        val = self.get_value(long_var_name)
        val[indices] = src

    def get_component_name(self):
        """
        Name of the component

        :returns: Name of the component as a string
        """
        return self.name

    def get_input_var_names(self):
        """
        List of input var names as standard name strings

        :returns: List of standard names
        """
        return self.input_vars

    def get_output_var_names(self):
        """
        List of output var names as standard name strings

        :returns: List of standard names
        """
        return self.output_vars

    def get_grid_shape(self, long_var_name):
        """
        Shape of the grid on which long_var_name is defined.

        :long_var_name: Standard name string

        :returns: Shape of the variable as a tuple
        """
        return self.get_value(long_var_name).shape

    def get_grid_spacing(self, long_var_name):
        """
        Row and column spacing of the grid on which long_var_name is defined.

        :long_var_name: Standard name string

        :returns: Spacing of the variable as a tuple
        """
        if long_var_name in self._var_names:
            return self._spacing
        else:
            raise BadVarNameError(long_var_name)

    def get_grid_origin(self, long_var_name):
        """
        Origin of the grid on which long_var_name is defined.

        :long_var_name: Standard name string

        :returns: Origin coordinates of the variable as a tuple
        """
        if long_var_name in self._var_names:
            return self._origin
        else:
            raise BadVarNameError(long_var_name)

    def get_grid_type(self, long_var_name):
        """
        Grid-type on which long_var_name is defined.

        :long_var_name: Standard name string

        :returns: Grid-type as a BmiGridType
        """
        if long_var_name in self._var_names:
            return GRID_TYPE_UNIFORM
        else:
            raise BadVarNameError(long_var_name)

    def get_time_step(self):
        """
        Model time step.

        :returns: Time step
        """
        return self._time_step

    def get_start_time(self):
        """
        Model start time.

        :returns: Start time
        """
        return self._start_time

    def get_end_time(self):
        """
        Model end time.

        :returns: End time
        """
        return self._end_time

    def get_current_time(self):
        """
        Current time of the model.

        :returns: Current model time
        """
        return self._time

def main ():
    """
    Run some unit tests
    """
    import doctest
    doctest.testmod ()

if __name__ == '__main__':
    main ()
