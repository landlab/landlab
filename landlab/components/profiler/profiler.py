#! /usr/bin/env python
"""
"""


from six.moves import range

import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError:
    import warnings
    warnings.warn('matplotlib not found', ImportWarning)

from landlab.plot import imshow_grid

from landlab.utils.return_array import return_array_at_node

from landlab import RasterModelGrid

from landlab import Component

class Profiler(Component):
    """
    """
    def __init__(self, grid, threshold):
        super(Profiler, self).__init__(grid)
        self._grid = grid

        if 'drainage_area' in grid.at_node:
            self.drainage_area = grid.at_node['drainage_area']
        else:
            msg = ''
            raise ValueError(msg)

        if 'flow__receiver_node' in grid.at_node:
            self.flow_receiver = grid.at_node['flow__receiver_node']
        else:
            msg = ''
            raise ValueError(msg)

        if 'flow__link_to_receiver_nodes':
            self.link_to_flow_receiver = grid.at_node['flow__link_to_receiver_node']
        else:
            msg = ''
            raise ValueError(msg)

    def plot_profiles(field='topographic_elevation'):
        """
        Plot distance-upstream vs arbitrary quantity, default when calling through
        analyze_channel_network_and_plot is topographic_elevation.

        Parameters
        ----------
        field : nnode array, required
            Array of  the at-node-field to plot against distance upstream.
        """
        quantity = return_array(field)

        # for each stream network
        for i in range(len(self._profile_structure)):
            network_nodes = self._profile_structure[i]
            network_distance = self._distances_upstream[i]

            # for each stream segment in the network
            for j in range(len(network_nodes)):

                # identify the nodes and distances upstream for this channel segment
                the_nodes = network_nodes[j]
                the_distances = network_distance[j]
                plt.plot(the_distances, quantity[the_nodes])


    def plot_profiles_in_map_view(grid, field='topographic__elevation',  **kwargs):
        """
        Plot profile locations in map view on a frame.

        Parameters
        ----------
        field, name or nnode long array to plot with imshow_grid
        **kwargs: additional parameters to pass to imshow_grid
        """
        # make imshow_grid background
        imshow_grid(grid, field, **kwargs)

        # for each stream network
        for i in range(len(self._profile_structure)):
            network_nodes = self._profile_structure[i]

            # for each stream segment in the network
            for j in range(len(network_nodes)):

                # identify the nodes and distances upstream for this channel segment
                the_nodes = network_nodes[j]
                plt.plot(grid.x_of_node[the_nodes], grid.y_of_node[the_nodes])
