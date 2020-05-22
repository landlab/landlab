# -*- coding: utf-8 -*-
"""Landlab component to calculate height above nearest drainage."""
from warnings import warn

import numpy as np

from landlab import Component


class HeightAboveDrainage(Component):
    """
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
        "flow__link_to_receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "ID of link downstream of each node, which carries the discharge",
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
        "downstream_drainage__node": {
            "dtype": float,
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

        if channel__mask is None:
            raise ValueError(
                "No channel mask supplied. "
                "A channel mask is needed "
                "to determine nearest drainage."
            )

        self._grid = grid
        self._channel_mask = channel__mask
        self._elev = grid.at_node["topographic__elevation"]
        self._receivers = grid.at_node["flow__receiver_node"]
        self._downstream_drainage_id = grid.at_node["downstream_drainage__node"]
        self._hand = grid.at_node["height_above_drainage__elevation"]


        @property
        def channel_mask(self):
            return self._channel_mask

        @channel_mask.setter
        def channel_mask(self,new_val):
            self._channel_mask = new_val


        def run_one_step(self):

            self_draining_nodes = sum(self._receivers == np.arange(self._grid.number_of_nodes))
            if self_draining_nodes != len(self._grid.boundary_nodes):
                warn(
                "Pits detected in the flow directions supplied. "
                "Pits will be treated as drainage nodes."
                )

            self._downstream_drainage_id[:] = 0
            is_drainage_node = self._channel_mask
            is_drainage_node[self._grid.open_boundary_nodes] = 1

            for i in range(self._grid.number_of_nodes):

                if i == receivers[i] or is_drainage_node[i]: #started on a boundary, depression, or in channel
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

            nearest_drainage_elev = self._elev[self._downstream_drainage_id]
            self._hand = self._elev-nearest_drainage_elev
