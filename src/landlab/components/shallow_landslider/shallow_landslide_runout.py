from __future__ import annotations

import heapq
import logging

import numpy as np

from landlab.core.model_component import Component

try:
    from tqdm import tqdm
except Exception:
    tqdm = None

logger = logging.getLogger("landslider")


class ShallowLandslideRunout(Component):
    """
    Subcomponent that routes shallow-landslide material downslope
    using flow-receiver topology and a distance-limited runout law.
    """

    _name = "ShallowLandslideRunout"

    _info = {
        "soil__depth": {
            "dtype": float,
            "intent": "inout",
            "mapping": "node",
            "units": "m",
            "optional": False,
            "doc": "Soil thickness at nodes.",
        },
        "hill_flow__receiver_node": {
            "dtype": int,
            "intent": "in",
            "mapping": "node",
            "units": "-",
            "optional": False,
            "doc": "Node IDs receiving routed hillslope flow from each node.",
        },
        "hill_flow__receiver_proportions": {
            "dtype": float,
            "intent": "in",
            "mapping": "node",
            "units": "-",
            "optional": False,
            "doc": "Proportions of hillslope flow routed to receiver nodes.",
        },
        "topographic__elevation": {
            "dtype": float,
            "units": "m",
            "intent": "in",
            "mapping": "node",
            "optional": False,
            "doc": "Land surface topographic elevation.",
        },
    }

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------

    def __init__(self, grid):
        super().__init__(grid)

        required = (
            "hill_flow__receiver_node",
            "hill_flow__receiver_proportions",
        )

        missing = [f for f in required if f not in grid.at_node]

        if missing:
            raise RuntimeError(
                "ShallowLandslideRunout requires flow routing fields, "
                f"but missing: {missing}. "
                "Run a FlowAccumulator or router before enabling runout."
            )

    def run_one_step(
        self,
        failed_nodes: np.ndarray,
        runout_distance: np.ndarray,
        return_paths: bool = False,
    ):
        """
        Route failed shallow-landslide material downslope and update
        soil thickness on the grid.
        """

        if failed_nodes.size == 0:
            return None

        paths, proportions, details = self._trace_paths_landslides(
            failed_nodes, runout_distance
        )

        soil, erosion, deposition = self._update_soil_depth(
            paths,
            proportions,
            self._grid.at_node["soil__depth"],
        )

        self._grid.at_node["soil__depth"][:] = soil

        # cache diagnostics (optional but useful)
        self._last_erosion = erosion
        self._last_deposition = deposition

        if return_paths:
            return details

    # %% Helpers
    def _calculate_node_distance(self, node1, node2):
        grid = self._grid
        x1, y1 = grid.node_x[node1], grid.node_y[node1]
        x2, y2 = grid.node_x[node2], grid.node_y[node2]
        return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

    def _trace_paths_landslides(
        self,
        starting_nodes,
        newmark_distances,
    ):
        grid = self._grid

        receiver_nodes = grid.at_node["hill_flow__receiver_node"]
        receiver_props = grid.at_node["hill_flow__receiver_proportions"]
        boundary_nodes = set(grid.boundary_nodes)

        final_paths = []
        final_proportions = []
        path_details = {}

        logger.info("Tracing landslide runout paths...")

        for node in starting_nodes:
            max_dist = newmark_distances[node]
            stack = [(0.0, node, 1.0, [node])]
            path_details[node] = []

            while stack:
                dist, current, prop, path = heapq.heappop(stack)

                if dist >= max_dist:
                    final_paths.append(tuple(path))
                    final_proportions.append(prop)
                    path_details[node].append((path, prop))
                    continue

                recs = receiver_nodes[current]
                props = receiver_props[current]

                # pit cell
                if recs[0] == current:
                    final_paths.append(tuple(path))
                    final_proportions.append(prop)
                    path_details[node].append((path, prop))
                    continue

                for r, p in zip(recs, props):
                    if r == -1:
                        continue

                    new_dist = dist + self._calculate_node_distance(current, r)
                    new_prop = prop * p
                    new_path = path + [r]

                    if r in boundary_nodes or np.isnan(
                        grid.at_node["topographic__elevation"][r]
                    ):
                        final_paths.append(tuple(path))
                        final_proportions.append(prop)
                        path_details[node].append((path, prop))
                    else:
                        heapq.heappush(stack, (new_dist, r, new_prop, new_path))

        return final_paths, final_proportions, path_details

    def _update_soil_depth(
        self,
        paths,
        proportions,
        soil_depth,
    ):
        soil = soil_depth.copy()
        erosion = np.zeros_like(soil)
        deposition = np.zeros_like(soil)

        logger.info("Updating soil depth...")
        for path, prop in zip(paths, proportions):
            for i in range(len(path) - 1):
                src = path[i]
                dst = path[i + 1]

                moved = soil_depth[src] * prop
                moved = min(moved, soil[src])

                erosion[src] += moved
                deposition[dst] += moved

        soil -= erosion
        soil += deposition

        return soil, erosion, deposition
