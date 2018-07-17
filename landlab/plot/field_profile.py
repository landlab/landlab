#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 10:26:17 2017

@author: njlyons
"""

from landlab import FieldError
from landlab.plot.imshow import imshow_grid
import matplotlib.pyplot as plt
import numpy as np


class FieldProfiler:

    def __init__(self, grid, field, start, end):
        """Create profiles of model grid field values.

        The profile extends from *start* to *end*. The field is sampled at
        approximately the resolution of *grid*. Profile traces that are
        horizontal or vertical will

        Parameters
        ----------
        grid : RasterModelGrid
            A landlab RasterModelGrid.
        field : string
            The field in which to draw a profile.
        start : tuple or integer
            The start of the profile as coordinates (x, y) or node.
        end : tuple or integer
            The end of the profile as coordinates (x, y) or node.
        """

        if field not in grid.at_node.keys():
            raise FieldError('the field, {} must be a field of the input grid '
                             'to create a profile'.format(field))

        # Store inputs.
        self._grid = grid
        self._field = field

        x_in = np.empty(2)
        y_in = np.empty(2)

        # Handle profile endpoint inputs that can be integers or tuples.
        x_in[0], y_in[0] = self._get_node_x_y(start)
        x_in[1], y_in[1] = self._get_node_x_y(end)

        # Get sample coordinates.

        dx = grid.dx

        total_distance = np.hypot(x_in[1] - x_in[0], y_in[1] - y_in[0])
        n = np.floor(total_distance / dx)

        n_samples = np.floor(np.hypot(x_in[1] - x_in[0], y_in[1] - y_in[0]) / dx)
        print(np.hypot(x_in[1] - x_in[0], y_in[1] - y_in[0]), n_samples, n)
        x_samples = np.linspace(x_in[0], x_in[1], n, endpoint=False)
        y_samples = np.linspace(y_in[0], y_in[1], n, endpoint=False)
        self.coordinates = list(zip(x_samples, y_samples))

        self.distance = np.hypot(x_samples - x_in[0], y_samples - y_in[0])

        # Get the nodes nearest to the sample coordinates in order to get the
        # field value at these nodes.

        self.nodes = self._grid.find_nearest_node((x_samples, y_samples))

        self.field_value = self._grid.at_node[self._field][self.nodes]

    def _get_node_x_y(self, point):
        if isinstance(point, (float, int)):
            x = self._grid.x_of_node[point]
            y = self._grid.y_of_node[point]
        elif isinstance(point, tuple):
            x = point[0]
            y = point[1]

        return x, y

    def plot(self, distance_unit_label, field_y_axis_label, **kwds):
        """Create a figure of the grid field profile.

        The figure has two subplots:
        1.  The trace of the profile plotted on the grid shaded by the field
            values.
        2.  The profile illustrating the profile samples (distance vs the
            field value).

        Parameters
        ----------
        distance_unit_label : string
            The distance unit label for the profile figure x-axis.
        field_y_axis_label : string
            The field value label for the profile figure y-axis.
        **kwds : tuple or integer
            The same parameters in :func:`imshow_grid_at_node` can be used
            in this method.
        """

        mg = self._grid

        plt.figure()

        # Plot field values as grid.
        plt.subplot(2, 1, 1)
        z = mg.at_node[self._field]
        z_min = min(z[z != -9999])
        z_max = max(z[z != -9999])

        # Set colorbar label to the profile y axis label unless specified in
        # kwds.
        if 'colorbar_label' not in kwds:
            kwds['colorbar_label'] = field_y_axis_label

        imshow_grid(self._grid, self._field, limits=(z_min, z_max), **kwds)

        # Plot profile trace on grid.
        x0 = self.coordinates[0][0]
        y0 = self.coordinates[0][1]
        xn = self.coordinates[-1][0]
        yn = self.coordinates[-1][1]
        plt.plot((x0, xn), (y0, yn), 'co-')

        # Plot profile.
        plt.subplot(2 ,1, 2)
        plt.plot(self.distance, self.field_value, 'k')
        plt.xlabel('Distance ({})'.format(distance_unit_label))
        plt.ylabel(field_y_axis_label)

        plt.tight_layout()
