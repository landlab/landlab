# -*- coding: utf-8 -*-
r"""
stream_power_smooth_threshold.py: Defines the StreamPowerSmoothThresholdEroder,
which is a version of the FastscapeEroder (derived from it).

StreamPowerSmoothThresholdEroder uses a mathematically smooth threshold
formulation, rather than one with a singularity. The erosion rate is defined as
follows:

$\omega = K A^m S$

$E = \omega - \omega_c \left[ 1 - \exp ( -\omega / \omega_c ) \right]$

Created on Sat Nov 26 08:36:49 2016

@author: gtucker
"""

import numpy as np

from landlab import BAD_INDEX_VALUE

from .cfuncs import smooth_stream_power_eroder_solver
from .fastscape_stream_power import FastscapeEroder

UNDEFINED_INDEX = BAD_INDEX_VALUE


class StreamPowerSmoothThresholdEroder(FastscapeEroder):
    """Stream erosion component with smooth threshold function.

    Parameters
    ----------
    grid : ModelGrid
        A grid.
    K_sp : float, array, or field name
        K in the stream power equation (units vary with other parameters).
    m_sp : float, optional
        m in the stream power equation (power on drainage area).
    n_sp : float, optional, ~ 0.5<n_sp<4.
        n in the stream power equation (power on slope). NOTE: NOT PRESENTLY
        HONORED BY StreamPowerSmoothThresholdEroder (TODO)
    threshold_sp : float (TODO: array, or field name)
        The threshold stream power.
    rainfall_intensity : float; optional
        NOT PRESENTLY HONORED (TODO)

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> rg = RasterModelGrid((3, 4))
    >>> rg.set_closed_boundaries_at_grid_edges(False, True, True, True)
    >>> z = rg.add_zeros('node', 'topographic__elevation')
    >>> z[5] = 2.0
    >>> z[6] = 1.0
    >>> from landlab.components import FlowAccumulator
    >>> fr = FlowAccumulator(rg, flow_director='D4')
    >>> fr.run_one_step()
    >>> from landlab.components import StreamPowerSmoothThresholdEroder
    >>> sp = StreamPowerSmoothThresholdEroder(rg, K_sp=1.0)
    >>> sp.thresholds
    1.0
    >>> sp.run_one_step(dt=1.0)
    >>> import numpy as np
    >>> np.round(z[5:7], 3)
    array([ 1.646,  0.667])
    >>> z[5] = 2.0
    >>> z[6] = 1.0
    >>> import numpy as np
    >>> q = np.zeros(rg.number_of_nodes) + 0.25
    >>> q[6] = 100.0
    >>> sp = StreamPowerSmoothThresholdEroder(rg, K_sp=1.0, use_Q=q)
    >>> sp.run_one_step(dt=1.0)
    >>> np.round(z[5:7], 3)
    array([ 1.754,  0.164])
    """

    def __init__(
        self,
        grid,
        K_sp=None,
        m_sp=0.5,
        n_sp=1.0,
        threshold_sp=1.0,
        rainfall_intensity=1.0,
        use_Q=None,
        **kwargs
    ):
        """Initialize StreamPowerSmoothThresholdEroder."""
        if "flow__receiver_node" in grid.at_node:
            if grid.at_node["flow__receiver_node"].size != grid.size("node"):
                msg = (
                    "A route-to-multiple flow director has been "
                    "run on this grid. The landlab development team has not "
                    "verified that StreamPowerSmoothThresholdEroder is compatible "
                    "with route-to-multiple methods. Please open a GitHub Issue "
                    "to start this process."
                )
                raise NotImplementedError(msg)

        if n_sp != 1.0:
            raise ValueError(
                ("StreamPowerSmoothThresholdEroder currently only " "supports n_sp = 1")
            )

        # Call base-class init
        super(StreamPowerSmoothThresholdEroder, self).__init__(
            grid,
            K_sp=K_sp,
            m_sp=m_sp,
            n_sp=n_sp,
            threshold_sp=threshold_sp,
            rainfall_intensity=rainfall_intensity,
            **kwargs
        )

        # Handle "use_Q" option (ideally should be done by base class, but
        # FastscapeEroder, unlike StreamPowerEroder, lacks this option)
        if use_Q is not None:
            if type(use_Q) is str:  # if str, assume it's a field name
                self.area_or_discharge = self._grid.at_node[use_Q]
            elif type(use_Q) is np.ndarray:  # if array, use it
                self.area_or_discharge = use_Q
            else:
                print("Warning: use_Q must be field name or array")
                self.area_or_discharge = self._grid.at_node["drainage_area"]
        else:
            self.area_or_discharge = self._grid.at_node["drainage_area"]

        # Arrays with parameters for use in implicit solver
        self.gamma = grid.empty(at="node")
        self.delta = grid.empty(at="node")

    def run_one_step(self, dt, flooded_nodes=None, runoff_rate=None, **kwds):
        """Run one forward iteration of duration dt.

        Parameters
        ----------
        dt : float
            Time step size
        flooded_nodes : ndarray of int (optional)
            Indices of nodes in lakes/depressions
        runoff_rate : (not used yet)
            (to be added later)

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> rg = RasterModelGrid((3, 3))
        >>> rg.set_closed_boundaries_at_grid_edges(False, True, True, True)
        >>> z = rg.add_zeros('node', 'topographic__elevation')
        >>> z[4] = 1.0
        >>> from landlab.components import FlowAccumulator
        >>> fr = FlowAccumulator(rg, flow_director='D4')
        >>> fr.run_one_step()
        >>> from landlab.components import StreamPowerSmoothThresholdEroder
        >>> sp = StreamPowerSmoothThresholdEroder(rg, K_sp=1.0)
        >>> sp.run_one_step(dt=1.0)
        >>> sp.alpha
        array([ 0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.])
        >>> sp.gamma
        array([ 0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.])
        >>> sp.delta
        array([ 0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.])
        """
        if self._grid.at_node["flow__receiver_node"].size != self._grid.size("node"):
            msg = (
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that StreamPowerSmoothThresholdEroder is compatible "
                "with route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )
            raise NotImplementedError(msg)
        # Set up needed arrays
        #
        # Get shorthand for elevation field ("z"), and for up-to-downstream
        # ordered array of node IDs ("upstream_order_IDs")
        node_id = np.arange(self._grid.number_of_nodes)
        upstream_order_IDs = self._grid["node"]["flow__upstream_node_order"]
        z = self._grid["node"]["topographic__elevation"]
        flow_receivers = self._grid["node"]["flow__receiver_node"]

        # Get an array of flow-link length for each node that has a defined
        # receiver (i.e., that drains to another node).
        defined_flow_receivers = np.not_equal(
            self._grid["node"]["flow__link_to_receiver_node"], UNDEFINED_INDEX
        )
        defined_flow_receivers[flow_receivers == node_id] = False
        if flooded_nodes is not None:
            defined_flow_receivers[flooded_nodes] = False
        flow_link_lengths = self._grid.length_of_d8[
            self._grid.at_node["flow__link_to_receiver_node"][defined_flow_receivers]
        ]

        # (Later on, add functionality for a runoff rate, or discharge, or
        # variable K)

        # Handle possibility of spatially varying K
        if type(self.K) is np.ndarray:
            K = self.K[defined_flow_receivers]
        else:
            K = self.K

        # Handle possibility of spatially varying threshold
        if isinstance(self.thresholds, np.ndarray):
            thresh = self.thresholds[defined_flow_receivers]
        else:
            thresh = self.thresholds

        # Set up alpha, beta, delta arrays
        #
        #   First, compute drainage area or discharge raised to the power m.
        np.power(self.area_or_discharge, self.m, out=self.A_to_the_m)

        #   Alpha
        self.alpha[~defined_flow_receivers] = 0.0
        self.alpha[defined_flow_receivers] = (
            K * dt * self.A_to_the_m[defined_flow_receivers] / flow_link_lengths
        )

        #   Gamma
        self.gamma[~defined_flow_receivers] = 0.0
        self.gamma[defined_flow_receivers] = dt * thresh

        #   Delta
        self.delta[~defined_flow_receivers] = 0.0
        if isinstance(self.thresholds, np.ndarray):
            self.delta[defined_flow_receivers] = (
                K * self.A_to_the_m[defined_flow_receivers]
            ) / (thresh * flow_link_lengths)

            self.delta[defined_flow_receivers][thresh == 0.0] = 0.0
        else:
            if thresh == 0:
                self.delta[defined_flow_receivers] = 0.0
            else:
                self.delta[defined_flow_receivers] = (
                    K * self.A_to_the_m[defined_flow_receivers]
                ) / (thresh * flow_link_lengths)

        # Iterate over nodes from downstream to upstream, using scipy's
        # 'newton' function to find new elevation at each node in turn.
        smooth_stream_power_eroder_solver(
            upstream_order_IDs, flow_receivers, z, self.alpha, self.gamma, self.delta
        )
