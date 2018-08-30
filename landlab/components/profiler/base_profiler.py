# coding: utf8
#! /usr/env/python
"""
"""

from six.moves import range

from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt

from landlab.plot import imshow_grid
from landlab.utils.return_array import return_array_at_node
from landlab import Component

class _NetworkProfiler(Component):
    """
    """
    def __init__(self, grid):
        super(_NetworkProfiler, self).__init__(grid)
        self._grid = grid

        if 'drainage_area' in grid.at_node:
            self._drainage_area = grid.at_node['drainage_area']
        else:
            msg = 'drainage_area is a required field to run a _NetworkProfiler.'
            raise ValueError(msg)

        if 'flow__receiver_node' in grid.at_node:
            self._flow_receiver = grid.at_node['flow__receiver_node']
        else:
            msg = 'flow__receiver_node is a required field to run a _NetworkProfiler.'
            raise ValueError(msg)

        if 'flow__link_to_receiver_node' in grid.at_node:
            self._link_to_flow_receiver = grid.at_node['flow__link_to_receiver_node']
        else:
            msg = 'flow__link_to_receiver_node is a required field to run a _NetworkProfiler.'
            raise ValueError(msg)

    def plot_profiles(self,
                      field='topographic__elevation',
                      colors=None,
                      xlabel='Distance Upstream',
                      ylabel='Plotted Quantity',
                      title='Channel Long Profile'):
        """
        Plot distance-upstream vs at at-node or size (nnodes,) quantity.

        Parameters
        ----------
        field : field name or nnode array
            Array of  the at-node-field to plot against distance upstream.
            Default value is the at-node field 'topographic__elevation'.
        colors : sequence of RGB tuples, optional
            Sequence of RGB tuples to use with each stream segment.
        xlabel : str, optional
            X-axis label, default is "Distance Upstream".
        ylabel : str, optional
            Y-axis label, default value is "Plotted Quantity".
        title : str, optional
            Plot title, default value is "Channel Long Profile".
        """
        quantity = return_array_at_node(self._grid, field)

        # for each stream network
        segments = []
        for i in range(len(self._profile_structure)):
            network_nodes = self._profile_structure[i]
            network_distance = self._distances_upstream[i]

            # for each stream segment in the network
            for j in range(len(network_nodes)):

                # identify the nodes and distances upstream for this channel segment
                the_nodes = network_nodes[j]
                the_distances = network_distance[j]
                segments.append(list(zip(the_distances, quantity[the_nodes])))

        # We need to set the plot limits.
        fig, ax = plt.subplots()
        ax.set_xlim(min(min(min(self._distances_upstream))),
                    max(max(max(self._distances_upstream))))
        ax.set_ylim(quantity.min(), quantity.max())

        line_segments = LineCollection(segments, colors=colors)
        ax.add_collection(line_segments)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

    def plot_profiles_in_map_view(self,
                                  field='topographic__elevation',
                                  colors=None,
                                  **kwargs):
        """
        Plot profile locations in map view.

        Parameters
        ----------
        field : field name or nnode array
            Array of  the at-node-field to plot as the 2D map values.
            Default value is the at-node field 'topographic__elevation'.
        colors : f..
        **kwargs : keyword arguments, arguments
            Additional parameters to pass to imshow_grid
        """
        # make imshow_grid background
        imshow_grid(self._grid, field, **kwargs)
        ax = plt.gca()

        segments = []

        # for each stream network
        for i in range(len(self._profile_structure)):
            network_nodes = self._profile_structure[i]

            # for each stream segment in the network
            for j in range(len(network_nodes)):

                # identify the nodes and distances upstream for this channel segment
                the_nodes = network_nodes[j]
                segments.append(list(zip(self._grid.x_of_node[the_nodes], self._grid.y_of_node[the_nodes])))

        line_segments = LineCollection(segments, colors=colors)
        ax.add_collection(line_segments)
