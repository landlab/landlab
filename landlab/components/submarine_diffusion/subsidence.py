#! /usr/bin/env python
import numpy as np
from scipy.interpolate import interp1d
from landlab import RasterModelGrid, Component


class SubsidenceTimeSeries(Component):

    _name = 'Subsider'

    _time_units = 'y'

    _input_var_names = ()

    _output_var_names = (
        'bedrock_surface__increment_of_elevation',
        'bedrock_surface__elevation',
    )

    _var_units = {
        'bedrock_surface__increment_of_elevation': 'm',
        'bedrock_surface__elevation': 'm',
    }

    _var_mapping = {
        'bedrock_surface__increment_of_elevation': 'node',
        'bedrock_surface__elevation': 'node',
    }

    _var_doc = {
        'bedrock_surface__increment_of_elevation': 'Increment of elevation',
        'bedrock_surface__elevation': 'Surface elevation',
    }

    def __init__(self, grid, filepath=None, kind='linear', **kwds):
        """Generate subsidence rates.

        Parameters
        ----------
        grid: RasterModelGrid
            A landlab grid.
        filepath: str
            Name of csv-formatted subsidence file.
        kind: str, optional
            Kind of interpolation as a string (one of 'linear',
            'nearest', 'zero', 'slinear', 'quadratic', 'cubic').
            Default is 'linear'.
        """
        super(SubsidenceTimeSeries, self).__init__(grid, **kwds)

        data = np.loadtxt(filepath, delimiter=',')
        subsidence = interp1d(data[:, 0], data[:, 1], kind=kind,
                              copy=True, assume_sorted=True,
                              bounds_error=True)
        inc = self.grid.add_empty('bedrock_surface__increment_of_elevation',
                                  at='node').reshape(self.grid.shape)
        inc[:] = subsidence(self.grid.x_of_node[self.grid.nodes_at_bottom_edge])

        self._dz = inc.copy()
        self._time = 0.

    @property
    def time(self):
        return self._time

    def run_one_step(self, dt):
        dz = self.grid.at_node['bedrock_surface__increment_of_elevation']
        z = self.grid.at_node['bedrock_surface__elevation']

        dz = dz.reshape(self.grid.shape)
        z = z.reshape(self.grid.shape)

        dz[:] = self._dz * dt
        z[:] += dz

        self._time += dt
