# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 10:00:25 2016

@author: Charlie Shobe

Landlab component to calculate drainage density

"""
from warnings import warn

import numpy as np

from landlab import Component, FieldError


class DrainageDensity(Component):

    r"""Calculate drainage density over a DEM.

    Landlab component that implements the distance to channel algorithm of
    Tucker et al., 2001.

    calc_drainage_density function returns drainage density for the model
    domain.

    calc_drainage_density calculates the distance from every node to the
    nearest channel node :math:`L` along the flow line of steepest descent
    (assuming D8 routing). The drainage density is then (after Tucker et al.,
    2001):

    .. math::

        D_d=\frac{1}{2\overline{L}}

    where :math:`\overline{L}` is the mean L for the model domain.

    This component requires EITHER a channel__mask array with 1's
    where channels exist and 0's elsewhere, OR a set of coefficients
    and exponents for a slope-area relationship and a
    channelization threshold to compare against that relationship.

    Written by C. Shobe on 7/11/2016, modified 10/7/2016.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator, FastscapeEroder
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
    >>> fr = FlowAccumulator(mg, flow_director='D8')
    >>> fsc = FastscapeEroder(mg, K_sp=.01, m_sp=.5, n_sp=1)
    >>> for x in range(100):
    ...     fr.run_one_step()
    ...     fsc.run_one_step(dt = 10.0)
    ...     mg.at_node['topographic__elevation'][mg.core_nodes] += .01
    >>> channels = mg.at_node['drainage_area'] > 5
    >>> dd = DrainageDensity(mg, channel__mask=channels)
    >>> mean_drainage_density = dd.calc_drainage_density()
    >>> np.isclose(mean_drainage_density, 0.3831100571)
    True
    """

    _name = "DrainageDensity"

    _input_var_names = (
        "flow__receiver_node",
        "flow__link_to_receiver_node",
        "topographic__steepest_slope",
        "channel__mask",
        "area_coefficient",
        "slope_coefficient",
        "area_exponent",
        "slope_exponent",
        "channelization_threshold",
    )

    _output_var_names = ("surface_to_channel__minimum_distance",)

    _var_units = {
        "flow__receiver_node": "-",
        "flow__link_to_receiver_node": "-",
        "topographic__steepest_slope": "-",
        "channel__mask": "-",
        "surface_to_channel__minimum_distance": "m",
    }

    _var_mapping = {
        "flow__receiver_node": "node",
        "flow__link_to_receiver_node": "node",
        "topographic__steepest_slope": "node",
        "channel__mask": "node",
        "surface_to_channel__minimum_distance": "node",
    }

    _var_doc = {
        "flow__receiver_node": "Node array of receivers (node that receives flow from current "
        "node)",
        "flow__link_to_receiver_node": "ID of link downstream of each node, which carries the discharge",
        "topographic__steepest_slope": "Topographic slope at each node",
        "channel__mask": "Logical map of at which grid nodes channels are present",
        "surface_to_channel__minimum_distance": "Distance from each node to the nearest channel",
    }

    def __init__(
        self,
        grid,
        channel__mask=None,
        area_coefficient=None,
        slope_coefficient=None,
        area_exponent=None,
        slope_exponent=None,
        channelization_threshold=None,
        **kwds
    ):
        """Initialize the DrainageDensity component.

        Parameters
        ----------
        grid : ModelGrid
        channel__mask : Array that holds 1's where
            channels exist and 0's elsewhere
        area_coefficient : coefficient to multiply drainage area by,
            for calculating channelization threshold
        slope_coefficient : coefficient to multiply slope by,
            for calculating channelization threshold
        area_exponent : exponent to raise drainage area to,
            for calculating channelization threshold
        slope_exponent : exponent to raise slope to,
            for calculating channelization threshold
        channelization_threshold : threshold value above
            which channels exist
        """
        if grid.at_node["flow__receiver_node"].size != grid.size("node"):
            msg = (
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that DrainageDensity is compatible with "
                "route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )
            raise NotImplementedError(msg)

        if channel__mask is not None:
            if area_coefficient is not None:
                warn(
                    "Channel mask and area "
                    "coefficient supplied. Defaulting "
                    "to channel mask, ignoring area "
                    "coefficient."
                )
            if slope_coefficient is not None:
                warn(
                    "Channel mask and slope "
                    "coefficient supplied. Defaulting "
                    "to channel mask, ignoring slope "
                    "coefficient."
                )
            if area_exponent is not None:
                warn(
                    "Channel mask and area "
                    "exponent supplied. Defaulting "
                    "to channel mask, ignoring area "
                    "exponent."
                )
            if slope_exponent is not None:
                warn(
                    "Channel mask and slope "
                    "exponent supplied. Defaulting "
                    "to channel mask, ignoring slope "
                    "exponent."
                )
            if channelization_threshold is not None:
                warn(
                    "Channel mask and channelization "
                    "threshold supplied. Defaulting "
                    "to channel mask, ignoring "
                    "threshold."
                )
            if grid.number_of_nodes != len(channel__mask):
                raise ValueError(
                    "Length of channel mask is not equal to " "number of grid nodes"
                )
            if "channel__mask" in grid.at_node:
                warn("Existing channel__mask grid field was overwritten.")

            grid.at_node["channel__mask"] = channel__mask

        if channel__mask is None:
            if area_coefficient is None:
                raise FieldError(
                    "No channel mask and no area "
                    "coefficient supplied. Either "
                    "a channel mask or all 5 threshold "
                    "parameters are needed."
                )
            if slope_coefficient is None:
                raise FieldError(
                    "No channel mask and no slope "
                    "coefficient supplied. Either "
                    "a channel mask or all 5 threshold "
                    "parameters are needed."
                )
            if area_exponent is None:
                raise FieldError(
                    "No channel mask and no area "
                    "exponent supplied. Either "
                    "a channel mask or all 5 threshold "
                    "parameters are needed."
                )
            if slope_exponent is None:
                raise FieldError(
                    "No channel mask and no slope "
                    "exponent supplied. Either "
                    "a channel mask or all 5 threshold "
                    "parameters are needed."
                )
            if channelization_threshold is None:
                raise FieldError(
                    "No channel mask and no channelization "
                    "threshold supplied. Either "
                    "a channel mask or all 5 threshold "
                    "parameters are needed."
                )
            channel__mask = (
                area_coefficient
                * np.power(grid.at_node["drainage_area"], area_exponent)
                * slope_coefficient
                * np.power(grid.at_node["topographic__steepest_slope"], slope_exponent)
            ) > channelization_threshold
            grid.at_node["channel__mask"] = channel__mask

        required = (
            "flow__receiver_node",
            "flow__link_to_receiver_node",
            "topographic__steepest_slope",
        )
        for name in required:
            if name not in grid.at_node:
                raise FieldError("{name}: missing required field".format(name=name))

        # Store grid
        self._grid = grid

        # for this component to work with Cython acceleration,
        # the channel_network must be uint8, not bool...
        self.channel_network = grid.at_node["channel__mask"].view(dtype=np.uint8)

        # Flow receivers
        self.flow_receivers = grid.at_node["flow__receiver_node"]

        # Links to receiver nodes
        self.stack_links = grid.at_node["flow__link_to_receiver_node"]

        # Distance to channel
        try:
            self.distance_to_channel = grid.at_node[
                "surface_to_channel__minimum_distance"
            ]
        except KeyError:
            self.distance_to_channel = grid.add_zeros(
                "surface_to_channel__minimum_distance", at="node", dtype=float
            )

    def calc_drainage_density(self):
        """Calculate drainage density.

        Returns
        -------
        float
            The drainage density.
        """

        # ^there is no 'run_one_step' method b/c this is a tool, not a model.
        """Calculate distance to channel and drainage density, after
        Tucker et al., 2001.

        Returns
        -------
        landscape_drainage_density : float (1/m)
            Drainage density over the model domain.
        """
        from .cfuncs import _calc_dists_to_channel

        _calc_dists_to_channel(
            self.channel_network,
            self.flow_receivers,
            self.grid.length_of_d8,
            self.stack_links,
            self.distance_to_channel,
            self.grid.number_of_nodes,
        )
        landscape_drainage_density = 1. / (
            2.0
            * np.mean(
                self.grid.at_node["surface_to_channel__minimum_distance"][
                    self.grid.core_nodes
                ]
            )
        )
        # self.distance_to_channel))  # this is THE drainage density
        return landscape_drainage_density
