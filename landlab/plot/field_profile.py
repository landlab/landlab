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
        """
        Create profiles of model grid field values.

        The profile extends from start to end. The grid field is sampled at the
        resolution of the grid. Profile endpoints and samples are snapped to
        node centers.

        Parameters
        ----------
        grid : RasterModelGrid
            A landlab RasterModelGrid.
        field : string
            The field in which to draw a profile.
        start : tuple or integer
            The start of the profile.
        end : tuple or integer
            The end of the profile
        """

        if field not in grid.at_node.keys():
            raise FieldError('A grid field, {} is required to '
                             'update a stream species range.'.format(field))

        self._grid = grid
        self._field = field

        # Handle inputs that can be integers or tuples for the profile
        # endpoints.
        start_node, x0, y0 = self._get_node_x_y(start)
        end_node, x1, y1  = self._get_node_x_y(end)

        # Get sample coordinates.
        n_samples = np.ceil(np.hypot(x1 - x0, y1 - y0) / grid.dx) + 1
        x_samples = np.linspace(x0, x1, n_samples)
        y_samples = np.linspace(y0, y1, n_samples)

        # Set sample parameters.
        self.nodes = []
        self.distance = []
        self.field_value = []

        for xi, yi in zip(x_samples, y_samples):
            node = self._grid.find_nearest_node((xi, yi))

            if node not in self.nodes and not grid.node_is_boundary(node):

                if len(self.nodes) == 0:
                    di = 0
                else:
                    prior_node = self.nodes[-1]
                    xp = grid.x_of_node[prior_node]
                    yp = grid.y_of_node[prior_node]
                    di = np.hypot(xi - xp, yi - yp) + self.distance[-1]
     
                zi = self._grid.at_node[self._field][node]
                self.nodes.append(node)
                self.distance.append(di)
                self.field_value.append(zi)

    def _get_node_x_y(self, point): 
        if isinstance(point, (float, int)):
            node = int(point)
            x = self._grid.x_of_node[point]
            y = self._grid.y_of_node[point]
        elif isinstance(point, tuple):
            node = self._grid.find_nearest_node((point[0], point[1]))
            x = point[0]
            y = point[1]

        return node, x, y

    def plot(self, distance_unit_label, field_y_axis_label, **kwds):
        """
        Create a figure of the grid field profile.

        The figure has two subplots:
            1. The trace of the profile plotted on the grid shaded by the
               field values.
            2. The profile illustrating the profile samples (distance vs the
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

        self.fig = plt.subplots(nrows=2, ncols=1)

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
        x1 = mg.x_of_node[self.nodes[0]]
        y1 = mg.y_of_node[self.nodes[0]]
        x2 = mg.x_of_node[self.nodes[-1]]
        y2 = mg.y_of_node[self.nodes[-1]]
        plt.plot((x1, x2), (y1, y2), 'co-')

        # Plot profile.
        plt.subplot(2 ,1, 2)
        plt.plot(self.distance, self.field_value, 'k')
        plt.xlabel('Distance ({})'.format(distance_unit_label))
        plt.ylabel(field_y_axis_label)

        plt.tight_layout()
