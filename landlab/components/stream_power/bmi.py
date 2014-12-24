#! /usr/bin/env python
import numpy as np

from ..flow_routing.route_flow_dn import FlowRouter
from .stream_power import StreamPowerEroder

from ...core.model_parameter_dictionary import ModelParameterDictionary
from ...grid.raster import RasterModelGrid


class BmiStreamPower(object):
    """BMI implementation of the landlab stream power component.

    Examples
    --------
    Create an instance of the model. Note that although an instance has been
    created, the component has not yet been initialized. This means that not
    all model data has been set. For instance, the model does not yet have a
    time step. The model is initialized through an input file by calling
    the `initialize` method.

    >>> from landlab.components.stream_power.bmi import BmiStreamPower
    >>> sp = BmiStreamPower()
    >>> sp.get_time_step() is None
    True
    >>> sp.get_input_var_names() # As CSDMS Standard Names
    ('land_surface__elevation', 'bedrock__uplift_rate')
    >>> sp.get_output_var_names()
    ('land_surface__elevation',)

    >>> from StringIO import StringIO
    >>> config = StringIO('''
    ... nrows:
    ... 100
    ... ncols:
    ... 100
    ... dx:
    ... 1.
    ... dt:
    ... 1.
    ... uplift_rate: #per time, not per timestep
    ... 0.001
    ... run_time:
    ... 100.
    ... init_elev:
    ... 0.
    ... rock_density:
    ... 2.7
    ... sed_density:
    ... 2.7
    ... K_sp:
    ... array
    ... m_sp:
    ... 0.5
    ... n_sp:
    ... 1.
    ... threshold_sp:
    ... 0.
    ... ''')

    After initialization, the model instance now has all of its data set.

    >>> sp.initialize('examples/drive_sp_params.txt')
    >>> sp.get_current_time(), sp.get_time_units()
    (0.0, 'd')
    >>> sp.get_time_step()
    365.0

    To advance the model in time, use the `update` method to advance by a
    single time step. Use the `update_until` to update to a particular time.

    >>> for _ in xrange(100):
    ...     sp.update()
    >>> sp.get_current_time()
    36500.0
    >>> sp.update_until(36501.5)
    >>> sp.get_current_time()
    36501.5
    """
    _name = 'Stream Power'
    _input_var_names = ('land_surface__elevation', 'bedrock__uplift_rate')
    _output_var_names = ('land_surface__elevation', )
    _var_units = {
        'land_surface__elevation': 'm',
        'bedrock__uplift_rate': 'm/d'
    }

    def __init__(self):
        self._stream_power = None
        self._flow_router = None
        self._grid = None
        self._values = {}
        self._time_step = None
        self._time = None

    def initialize(self, filename):
        """Initialize the component from a file.

        Parameters
        ----------
        filename : str
            Name of input file.
        """
        params = ModelParameterDictionary(filename)

        self._time_step = float(params['dt']) * 365.
        self._time = 0.

        self._grid = RasterModelGrid(int(params['nrows']),
                                     int(params['ncols']),
                                     float(params['dx']))

        z = self._grid.add_zeros('node', 'topographic_elevation')
        z[:] = (z + float(params['init_elev']) +
                np.random.rand(z.size) / 1000.)

        k = self._grid.add_zeros('node', 'K_values')
        k[:] = k + 0.1 + np.random.rand(k.size) / 10.

        uplift_rate = np.ones_like(z) * float(params['uplift_rate']) / 365.

        self._stream_power = StreamPowerEroder(self._grid, filename)
        self._flow_router = FlowRouter(self._grid)

        self._values = {
            'land_surface__elevation': z,
            'bedrock__uplift_rate': uplift_rate,
        }

    def update(self):
        """Update the model by one time step.
        """
        self.update_frac(1.)

    def update_until(self, then):
        """Update model until a time.

        Parameters
        ----------
        then : float
            Time to run model to.

        Raises
        ------
        ValueError : *then* is less than now.
        """
        n_steps = (then - self.get_current_time()) / self.get_time_step()
        if n_steps < 0:
            raise ValueError('negative time step')

        for _ in xrange(int(n_steps)):
            self.update()
        self.update_frac(n_steps - int(n_steps))

    def update_frac(self, dt_frac):
        """Run model for a fraction of a time step.

        Parameters
        ----------
        dt_frac : float
            Fraction of time step to run model for.

        Raises
        ------
        ValueError : *dt_frac* is outside of (0., 1.]
        """
        if dt_frac <= 0. or dt_frac > 1.:
            raise ValueError('dt_frac out of bounds')

        dt = self.get_time_step() * dt_frac

        self._flow_router.route_flow(grid=self._grid)
        self._stream_power.erode(self._grid, dt / 365.,
                                 node_drainage_areas='drainage_area',
                                 slopes_at_nodes='steepest_slope',
                                 K_if_used='K_values')

        z = self._values['land_surface__elevation']
        uplift_rate = self._values['bedrock__uplift_rate']
        z[self._grid.core_nodes] += uplift_rate[self._grid.core_nodes] * dt

        self._time += dt

    def get_component_name(self):
        """The name of the component.

        Returns
        -------
        name : str
            Name of the model.
        """
        return self._name

    def get_input_var_names(self):
        """Standard names for input variables.

        Returns
        -------
        names : tuple
            Names of input variables as CSDMS Standard Names.
        """
        return self._input_var_names

    def get_output_var_names(self):
        """Standard names for output variables.

        Returns
        -------
        names : tuple
            Names of output variables as CSDMS Standard Names.
        """
        return self._output_var_names

    def get_grid_shape(self, name):
        """Shape of the grid for a variable.

        Parameters
        ----------
        name : str
            Variable name

        Returns
        -------
        shape : tuple of ints
            Grid shape.
        """
        return (self._grid.number_of_node_rows,
                self._grid.number_of_node_columns)

    def get_grid_spacing(self, name):
        """Spacing of grid rows and columns for a variable.

        Parameters
        ----------
        name : str
            Variable name

        Returns
        -------
        spacing : tuple of floats
            Grid spacing.
        """
        return (self._grid.node_spacing, self._grid.node_spacing)

    def get_grid_origin(self, name):
        """Coordinates of lower-left corner of grid variable.

        Parameters
        ----------
        name : str
            Variable name

        Returns
        -------
        origin : tuple of floats
            Grid origin.
        """
        return (0., 0.)

    def get_value(self, name):
        """Copy of grid variable data.

        Parameters
        ----------
        name : str
            Variable name

        Returns
        -------
        data : ndarray
            Variable data.
        """
        return self.get_value_ptr(name).copy()

    def get_value_ptr(self, name):
        """Grid variable data.

        Parameters
        ----------
        name : str
            Variable name

        Returns
        -------
        data : ndarray
            Variable data.
        """
        if name in self._output_var_names:
            return self._values[name]
        else:
            raise KeyError(name)

    def set_value(self, name, src):
        """Set the values for a grid variable.

        Parameters
        ----------
        name : str
            Variable name
        src : ndarray
            New values for the variable.
        """
        if name in self._input_var_names:
            val = self._values[name]
            val[:] = src.flat
        else:
            raise KeyError(name)

    def get_start_time(self):
        """Model start time.

        Returns
        -------
        time : float
            Start time.
        """
        return 0.

    def get_end_time(self):
        """Model end time.

        Returns
        -------
        time : float
            End time.
        """
        return np.finfo('d').max

    def get_current_time(self):
        """Current time of model.

        Returns
        -------
        time : float
            Current time.
        """
        return self._time

    def get_time_step(self):
        """Model time step.

        Returns
        -------
        dt : float
            Time step.
        """
        return self._time_step

    def get_time_units(self):
        """Units that all time is reported in.

        Returns
        -------
        units : str
            Units for time values.
        """
        return 'd'
