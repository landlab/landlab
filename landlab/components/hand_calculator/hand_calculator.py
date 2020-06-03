# -*- coding: utf-8 -*-
"""Landlab component to calculate height above nearest drainage.

@author: D Litwin
"""
from warnings import warn

import numpy as np

from landlab import Component
from landlab.utils import return_array_at_node


class HeightAboveDrainageCalculator(Component):
    """
    Calculate the elevation difference between each node and its nearest
    drainage node in a DEM.

    This component implements the method described by Nobre et al (2011). A
    single direction flow director (D8 or steepest descent) must be run prior
    to HeightAboveDrainageCalculator to supply the flow directions. This component does
    not fill depressions in a DEM, but rather it treats them as drainage nodes.
    For best results, please run one of the available pit filling components
    prior to HeightAboveDrainageCalculator.

    Examples
    --------

    >>> import numpy as np
    >>> from numpy.testing import assert_equal

    >>> from landlab import RasterModelGrid
    >>> from landlab.components import HeightAboveDrainageCalculator, FlowAccumulator

    >>> mg = RasterModelGrid((4, 5))
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> mg.set_status_at_node_on_edges(right=mg.BC_NODE_IS_CLOSED, bottom=mg.BC_NODE_IS_FIXED_VALUE, \
                                  left=mg.BC_NODE_IS_CLOSED, top=mg.BC_NODE_IS_CLOSED)
    >>> elev = np.array([[2,1,0,1,2],[3,2,1,2,3],[4,3,2,3,4],[5,4,4,4,5]])
    >>> z[:] = elev.reshape(len(z))
    >>> elev
    array([[2, 1, 0, 1, 2],
       [3, 2, 1, 2, 3],
       [4, 3, 2, 3, 4],
       [5, 4, 4, 4, 5]])

    >>> fa = FlowAccumulator(mg, flow_director="D8")
    >>> fa.run_one_step()

    >>> channel__mask = mg.zeros(at="node")
    >>> channel__mask[[2,7]] = 1
    >>> channel__mask.reshape(elev.shape)
    array([[ 0.,  0.,  1.,  0.,  0.],
       [ 0.,  0.,  1.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.]])

    >>> hd = HeightAboveDrainageCalculator(mg, channel_mask=channel__mask)
    >>> hd.run_one_step()

    >>> mg.at_node["height_above_drainage__elevation"].reshape(elev.shape)  # doctest: +NORMALIZE_WHITESPACE
    array([[ 2.,  0.,  0.,  0.,  0.],
           [ 3.,  2.,  0.,  2.,  3.],
           [ 4.,  2.,  1.,  2.,  4.],
           [ 5.,  4.,  4.,  4.,  5.]])

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    Nobre, A. D., Cuartas, L. A., Hodnett, M., Rennó, C. D., Rodrigues, G.,
    Silveira, A., et al. (2011). Height Above the Nearest Drainage – a
    hydrologically relevant new terrain model. Journal of Hydrology, 404(1),
    13–29. https://doi.org/10.1016/j.jhydrol.2011.03.051

    """

    _name = "HeightAboveDrainageCalculator"

    _unit_agnostic = True

    _info = {
        "channel__mask": {
            "dtype": np.uint8,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "Logical map of at which grid nodes channels are present",
        },
        "flow__receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        "flow__upstream_node_order": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array containing downstream-to-upstream ordered list of node IDs",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "height_above_drainage__elevation": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Elevation above the nearest channel node",
        },
    }

    def __init__(self, grid, channel_mask="channel__mask"):
        """
        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        channel_mask : field name, array of uint8
            Logical map of nodes where drainage is present
        """
        super().__init__(grid)

        if grid.at_node["flow__receiver_node"].size != grid.size("node"):
            msg = (
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that HeightAboveDrainageCalculator is compatible with "
                "route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )
            raise NotImplementedError(msg)

        self._grid = grid
        self._channel_mask = return_array_at_node(self._grid, channel_mask)
        self._elev = grid.at_node["topographic__elevation"]
        self._receivers = grid.at_node["flow__receiver_node"]
        self._node_order = grid.at_node["flow__upstream_node_order"]

        # height above nearest drainage
        if "height_above_drainage__elevation" in grid.at_node:
            self._hand = grid.at_node["height_above_drainage__elevation"]
        else:
            self._hand = grid.add_zeros(
                "height_above_drainage__elevation", at="node", dtype=float
            )

    @property
    def channel_mask(self):
        return self._channel_mask

    @channel_mask.setter
    def channel_mask(self, new_val):
        self._channel_mask = return_array_at_node(self._grid, new_val)

    def run_one_step(self):

        is_drainage_node = self._channel_mask
        is_drainage_node[self._grid.open_boundary_nodes] = 1

        # check for pits
        self_draining_nodes = np.where(
            self._receivers == np.arange(self._grid.number_of_nodes)
        )
        pits = np.setxor1d(self_draining_nodes, self._grid.boundary_nodes)
        if pits.any():
            warn(
                "Pits detected in the flow directions supplied. "
                "Pits will be treated as drainage nodes."
            )
            is_drainage_node[pits] = 1

        # iterate downstream through stack to find nearest drainage elevation
        nearest_drainage_elev = np.empty(self._elev.shape)
        for n in self._node_order:
            r = self._receivers[n]
            # if not drainage node set drainage elevation to downstream.
            if not is_drainage_node[n]:
                nearest_drainage_elev[n] = nearest_drainage_elev[r]
            else:  # set elevation of drainage to self.
                nearest_drainage_elev[n] = self._elev[n]

        self._hand[:] = self._elev - nearest_drainage_elev
