# -*- coding: utf-8 -*-
"""Landlab component to calculate height above nearest drainage."""
from warnings import warn

import numpy as np

from landlab import Component


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
    >>> elev
    >>> z[:] = elev.reshape(len(z))

    >>> fa = FlowAccumulator(mg, flow_director="D8")
    >>> fa.run_one_step()

    >>> channel__mask = mg.zeros(at="node")
    >>> channel__mask[[2,7]] = 1
    >>> channel__mask.reshape(elev.shape)

    >>> hd = HeightAboveDrainage(mg, channel__mask)
    >>> hd.run_one_step()

    >>> mg.at_node["height_above_drainage__elevation"].reshape(elev.shape)



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

    _name = "HeightAboveDrainage"

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
        "downstream_drainage__node": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "node array of nearst drainage node ID",
        },
    }

    def __init__(self, grid, channel__mask):

        super().__init__(grid)

        if grid.at_node["flow__receiver_node"].size != grid.size("node"):
            msg = (
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that HeightAboveDrainage is compatible with "
                "route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )
            raise NotImplementedError(msg)

        self._grid = grid
        self._channel_mask = channel__mask
        self._elev = grid.at_node["topographic__elevation"]
        self._receivers = grid.at_node["flow__receiver_node"]

        # Downstream drainage node
        if "downstream_drainage__node" in grid.at_node:
            self._downstream_drainage_id = grid.at_node["downstream_drainage__node"]
        else:
            self._downstream_drainage_id = grid.add_zeros(
                "downstream_drainage__node", at="node", dtype=int
            )

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
        self._channel_mask = new_val

    def run_one_step(self):

        self._downstream_drainage_id[:] = 0
        is_drainage_node = self._channel_mask
        is_drainage_node[self._grid.open_boundary_nodes] = 1

        # check for pits
        self_draining_nodes = np.where(self._receivers == np.arange(self._grid.number_of_nodes))
        pits = np.setxor1d(self_draining_nodes,self._grid.boundary_nodes)
        if pits.any():
            warn(
                "Pits detected in the flow directions supplied. "
                "Pits will be treated as drainage nodes."
            )
            is_drainage_node[pits] = 1

        # find drainage nodes
        for i in range(self._grid.number_of_nodes):

            if (
                i == self._receivers[i] or is_drainage_node[i]
            ):  # started on a boundary, depression, or in channel
                self._downstream_drainage_id[i] = i

            else:
                cur_node = i
                reached_drainage = False
                while not reached_drainage:
                    downstream_id = self._receivers[cur_node]
                    if is_drainage_node[downstream_id]:
                        self._downstream_drainage_id[i] = downstream_id
                        reached_drainage = True
                    else:
                        cur_node = downstream_id

        # calculate height above drainage node
        nearest_drainage_elev = self._elev[self._downstream_drainage_id]
        self._hand[:] = self._elev - nearest_drainage_elev
