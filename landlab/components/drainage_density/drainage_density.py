# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 10:00:25 2016

@author: Charlie Shobe

for West Valley Project: Landlab component to calculate drainage density

"""
from landlab import Component
import numpy as np


class DrainageDensity(Component):

    """Calculate drainage density over a DEM.

    calc_drainage_density function returns drainage density for the model
    domain.

    Landlab component that implements the distance to channel algorithm of
    Tucker et al., 2001.

    Written by C. Shobe on 7/11/2016

    Construction::

        DrainageDensity(grid, channel_network_name='string')

    Parameters
    ----------
    grid : ModelGrid
    channel_network_name : string naming a boolean grid field with 'True'
        in channels and 'False' elsewhere.

    Examples
    ---------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.flow_routing import FlowRouter
    >>> from landlab.components import FastscapeEroder

    >>> mg = RasterModelGrid((10, 10), 1.0)
    >>> _ = mg.add_zeros('node', 'topographic__elevation')

    >>> np.random.seed(50)
    >>> noise = np.random.rand(100)
    >>> mg.at_node['topographic__elevation'] += noise
    >>> mg.at_node['topographic__elevation'] # doctest: +NORMALIZE_WHITESPACE
    array([ 0.49460165,  0.2280831 ,  0.25547392,  0.39632991,  0.3773151 ,
        0.99657423,  0.4081972 ,  0.77189399,  0.76053669,  0.31000935,
        0.3465412 ,  0.35176482,  0.14546686,  0.97266468,  0.90917844,
        0.5599571 ,  0.31359075,  0.88820004,  0.67457307,  0.39108745,
        0.50718412,  0.5241035 ,  0.92800093,  0.57137307,  0.66833757,
        0.05225869,  0.3270573 ,  0.05640164,  0.17982769,  0.92593317,
        0.93801522,  0.71409271,  0.73268761,  0.46174768,  0.93132927,
        0.40642024,  0.68320577,  0.64991587,  0.59876518,  0.22203939,
        0.68235717,  0.8780563 ,  0.79671726,  0.43200225,  0.91787822,
        0.78183368,  0.72575028,  0.12485469,  0.91630845,  0.38771099,
        0.29492955,  0.61673141,  0.46784623,  0.25533891,  0.83899589,
        0.1786192 ,  0.22711417,  0.65987645,  0.47911625,  0.07344734,
        0.13896007,  0.11230718,  0.47778497,  0.54029623,  0.95807105,
        0.58379231,  0.52666409,  0.92226269,  0.91925702,  0.25200886,
        0.68263261,  0.96427612,  0.22696165,  0.7160172 ,  0.79776011,
        0.9367512 ,  0.8537225 ,  0.42154581,  0.00543987,  0.03486533,
        0.01390537,  0.58890993,  0.3829931 ,  0.11481895,  0.86445401,
        0.82165703,  0.73749168,  0.84034417,  0.4015291 ,  0.74862   ,
        0.55962945,  0.61323757,  0.29810165,  0.60237917,  0.42567684,
        0.53854438,  0.48672986,  0.49989164,  0.91745948,  0.26287702])

    >>> fr = FlowRouter(mg)
    >>> fsc = FastscapeEroder(mg, K_sp=.01, m_sp=.5, n_sp=1)
    >>> for x in range(100):
    ...     fr.run_one_step()
    ...     fsc.run_one_step(dt = 10.0)
    ...     mg.at_node['topographic__elevation'][mg.core_nodes] += .01

    >>> channels = mg.at_node['drainage_area'] > 5
    >>> _ = mg.add_field('channel_network', channels, at='node')

    >>> dd = DrainageDensity(mg, channel_network_name='channel_network')
    >>> mean_drainage_density = dd.calc_drainage_density()
    >>> np.isclose(mean_drainage_density, 0.3831100571)
    True
    """

    _name = 'DrainageDensity'

    _input_var_names = (
        'flow__receiver_node',
        'flow__link_to_receiver_node',
        'channel_network',
    )

    _output_var_names = (
        'distance_to_channel',
    )

    _var_units = {
        'flow__receiver_node': '-',
        'flow__link_to_receiver_node': '-',
        'channel_network': '-',
        'distance_to_channel': 'm',
    }

    _var_mapping = {
        'flow__receiver_node': 'node',
        'flow__link_to_receiver_node': 'node',
        'channel_network': 'node',
        'distance_to_channel': 'node',
    }

    _var_doc = {
        'flow__receiver_node':
            'Node array of receivers (node that receives flow from current '
            'node)',
        'flow__link_to_receiver_node':
            'ID of link downstream of each node, which carries the discharge',
        'channel_network':
            'Logical map of at which grid nodes channels are present',
        'distance_to_channel':
            'Distance from each node to the nearest channel',
    }

    def __init__(self, grid, channel_network_name='channel_network', **kwds):
        """Initialize the DrainageDensity component.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        channel_network_name : string, optional (defaults to 'channel_network')
            Name of Landlab grid field that holds 1's where channels exist
            and 0's elsewhere
        """

        # Store grid
        self._grid = grid

        # Create fields... Input fields should raise an error if not present,
        # rather than silently creating a nonsensical blank field.
        # Channel network
        self.channel_network_name = channel_network_name
        if channel_network_name in grid.at_node:
            self.channel_network = grid.at_node[channel_network_name].astype(
                int)
            # for this component to work with Cython acceleration,
            # the channel_network must be int, not bool...
        else:
            raise FieldError(
                'DrainageDensity needs the field channel_network!')

        # Flow receivers
        if 'flow__receiver_node' in grid.at_node:
            self.flow_receivers = grid.at_node['flow__receiver_node']
        else:
            raise FieldError(
                'DrainageDensity needs the field flow__receiver_node!')

        # Links to receiver nodes
        if 'flow__link_to_receiver_node' in grid.at_node:
            # ^flow__link_to_receiver_node
            self.stack_links = grid.at_node['flow__link_to_receiver_node']
        else:
            raise FieldError(
                'DrainageDensity needs the field flow__link_to_receiver_node!')

        #   Distance to channel
        if 'distance_to_channel' in grid.at_node:
            self.distance_to_channel = grid.at_node['distance_to_channel']
        else:
            self.distance_to_channel = grid.add_zeros(
                'node', 'distance_to_channel', dtype=float)

    def calc_drainage_density(self, **kwds):
        # ^there is no 'run_one_step' methid b/c this is a tool, not a model.
        """Calculate distance to channel and drainage density, after
        Tucker et al., 2001.

        Returns
        -------
        landscape_drainage_density : float (1/m)
            Drainage density over the model domain.
        """
        from .cfuncs import _calc_dists_to_channel
        _calc_dists_to_channel(self.channel_network,
                               self.flow_receivers,
                               self.grid._length_of_link_with_diagonals,
                               self.stack_links,
                               self.distance_to_channel,
                               self.grid.number_of_nodes)
        landscape_drainage_density = 1. / (2.0 * np.mean(self.grid.at_node[
            'distance_to_channel'][self.grid.core_nodes]))
        # self.distance_to_channel))  # this is THE drainage density
        return landscape_drainage_density
