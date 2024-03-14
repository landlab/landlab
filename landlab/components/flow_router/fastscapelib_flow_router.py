"""fastscapelib_flow_router.py: Component to compute flow routes, accumulate
flow and calculate drainage area.

It provides the FastscapelibFlowRouter components which exposes the flow routing
capabilities available in the Fastscapelib library
(https://fastscapelib.readthedocs.io/). This includes:

- computing flow routes on different kinds of grids and using simple to advanced
  strategies (single vs. multiple direction flow, closed depression resolving)
- accumulate flow
- compute drainage area
"""

import numpy as np
import fastscapelib as fs

from landlab import (
    Component,
    RasterModelGrid,
)
from landlab.utils import return_array_at_node


def _get_nodes_status(grid):
    """Return a dictionary of Fastscapelib NodeStatus from a landlab grid."""

    node_status = {}
    for i in grid.fixed_value_boundary_nodes:
        node_status[i] = fs.NodeStatus.FIXED_VALUE
    for i in grid.fixed_gradient_boundary_nodes:
        node_status[i] = fs.NodeStatus.FIXED_GRADIENT

    return node_status


class FastscapelibFlowRouter(Component):

    _name = "FastscapelibFlowRouter"

    _unit_agnostic = True

    _info = {
        # input grid fields
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        # output grid fields
        "flow__link_to_receiver_node": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "ID of link downstream of each node, which carries the discharge",
        },
        "flow__receiver_node": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        "flow__upstream_node_order": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array containing downstream-to-upstream ordered list of node IDs",
        },
        "flood_status_code": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Map of flood status (_PIT, _CURRENT_LAKE, _UNFLOODED, or _FLOODED).",
        },
        "topographic__steepest_slope": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "The steepest *downhill* slope",
        },
        "drainage_area": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**2",
            "mapping": "node",
            "doc": "Upstream accumulated surface area contributing to the node's discharge",
        },
        "surface_water__discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": "Volumetric discharge of surface water",
        },
    }

    def __init__(
        self,
        grid,
        operators,
        surface="topographic__elevation",
    ):
        super().__init__(grid)

        node_status = _get_nodes_status(grid)

        if isinstance(grid, RasterModelGrid):
            self._fs_grid = fs.RasterGrid(
                list(grid.shape), grid.spacing, fs.NodeStatus.CORE
            )
        else:
            raise TypeError("other grids than raster are not yet supported")

        self._fs_flow_graph = fs.FlowGraph(self._fs_grid, operators)
        base_levels = [i for i, status in node_status.items() if status == fs.NodeStatus.FIXED_VALUE]
        self._fs_flow_graph.base_levels = base_levels

        # do not include closed nodes in the flow graph
        mask = np.zeros(grid.shape, dtype=bool)
        mask.ravel()[grid.closed_boundary_nodes] = True
        self._fs_flow_graph.mask = mask

        self._surface = surface
        self._surface_values = return_array_at_node(grid, surface)

        self.initialize_output_fields()

        self._drainage_area = grid.at_node["drainage_area"]

    def _maybe_reshape(self, arr) -> np.ndarray:
        if isinstance(self._grid, RasterModelGrid):
            return arr.reshape(self._grid.shape)
        else:
            return arr

    @property
    def surface_values(self):
        """Values of the surface over which flow is routed and accumulated."""
        return self._surface_values

    @property
    def node_drainage_area(self):
        """Return the drainage area."""
        return self._grid["node"]["drainage_area"]

    @property
    def node_water_discharge(self):
        """Return the surface water discharge."""
        return self._grid["node"]["surface_water__discharge"]

    def _set_link_to_receiver_node(self, receivers):
        # TODO: more grid-agnostic and efficient way to do this?
        ncols = self._grid.shape[1]
        rel_idx = receivers - np.arange(receivers.size)

        # translate relative receiver index to landlab relative link index
        map_idx = [1, ncols, -1, -ncols, ncols + 1, ncols - 1, -ncols - 1, -ncols + 1]

        temp = rel_idx.copy()

        for i, idx in enumerate(map_idx):
            temp = np.where(rel_idx == idx, i, temp)

        rcvr_links = np.take_along_axis(self._grid.d8s_at_node, temp[:, None], axis=1).squeeze()
        self._grid.at_node["flow__link_to_receiver_node"][:] = np.where(rel_idx == 0, -1, rcvr_links)

    def update_routes(self):
        filled_elevation = self._fs_flow_graph.update_routes(self._maybe_reshape(self.surface_values))
        filled_elevation_flat = filled_elevation.ravel()

        elevation = self.surface_values
        receivers = self._fs_flow_graph.impl().receivers.astype("int").squeeze()
        receivers_distance = self._fs_flow_graph.impl().receivers_distance.squeeze()
        node_order = self._fs_flow_graph.impl().nodes_indices_bottomup.astype("int")
        slope = (filled_elevation_flat - filled_elevation_flat[receivers]) / np.where(receivers_distance == 0, 1, receivers_distance)
        flood_code = np.where(elevation < filled_elevation_flat, 3, 0)   # 0: unflooded, 3: flooded
        flood_code[self._fs_flow_graph.base_levels] = 0

        self._grid.at_node["flow__receiver_node"][:] = receivers
        self._grid.at_node["flow__upstream_node_order"][:] = node_order
        self._grid.at_node["topographic__steepest_slope"][:] = slope
        self._grid.at_node["flood_status_code"] = flood_code
        self._set_link_to_receiver_node(receivers)

    def accumulate_flow(self):
        self._fs_flow_graph.accumulate(self._maybe_reshape(self.node_drainage_area), 1)
        self._fs_flow_graph.accumulate(self._maybe_reshape(self.node_water_discharge), 1)

    def run_one_step(self):
        self.update_routes()
        self.accumulate_flow()
