# coding: utf8
#! /usr/env/python
"""
"""

from six.moves import range

try:
    import matplotlib.pyplot as plt
except ImportError:
    import warnings
    warnings.warn('matplotlib not found', ImportWarning)

from landlab.plot import imshow_grid

from landlab.utils.return_array import return_array_at_node

from landlab import Component

class _Profiler(Component):
    """
    """
    def __init__(self, grid):
        super(Profiler, self).__init__(grid)
        self._grid = grid

        if 'drainage_area' in grid.at_node:
            self._drainage_area = grid.at_node['drainage_area']
        else:
            msg = ''
            raise ValueError(msg)

        if 'flow__receiver_node' in grid.at_node:
            self._flow_receiver = grid.at_node['flow__receiver_node']
        else:
            msg = ''
            raise ValueError(msg)

        if 'flow__link_to_receiver_nodes':
            self._link_to_flow_receiver = grid.at_node['flow__link_to_receiver_node']
        else:
            msg = ''
            raise ValueError(msg)

    def plot_profiles(self, field='topographic__elevation'):
        """
        Plot distance-upstream vs arbitrary quantity, default when calling through
        analyze_channel_network_and_plot is topographic__elevation.

        Parameters
        ----------
        field : nnode array, required
            Array of  the at-node-field to plot against distance upstream.
        """
        quantity = return_array_at_node(self._grid, field)

        # for each stream network
        for i in range(len(self.profile_structure)):
            network_nodes = self.profile_structure[i]
            network_distance = self.distances_upstream[i]

            # for each stream segment in the network
            for j in range(len(network_nodes)):

                # identify the nodes and distances upstream for this channel segment
                the_nodes = network_nodes[j]
                the_distances = network_distance[j]
                plt.plot(the_distances, quantity[the_nodes])


    def plot_profiles_in_map_view(self, field='topographic__elevation',  **kwargs):
        """
        Plot profile locations in map view on a frame.

        Parameters
        ----------
        field, name or nnode long array to plot with imshow_grid
        **kwargs: additional parameters to pass to imshow_grid
        """
        # make imshow_grid background
        imshow_grid(self._grid, field, **kwargs)

        # for each stream network
        for i in range(len(self.profile_structure)):
            network_nodes = self.profile_structure[i]

            # for each stream segment in the network
            for j in range(len(network_nodes)):

                # identify the nodes and distances upstream for this channel segment
                the_nodes = network_nodes[j]
                plt.plot(self._grid.x_of_node[the_nodes], self._grid.y_of_node[the_nodes])
