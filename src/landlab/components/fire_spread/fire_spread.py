"""
FireSpread – Rothermel fire propagation on Landlab RasterModelGrid
"""

from __future__ import annotations

import numpy as np
from landlab import Component, RasterModelGrid
from .fuel_models import ANDERSON_13
from .utils import byram_flame_length, rothermel_rate_of_spread, reaction_intensity


class FireSpread(Component):
    """
    Rothermel (1972) surface fire spread using level-set fast marching.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((40, 60), xy_spacing=30.0)
    >>> grid.add_zeros("fuel__model", at="cell", dtype=int)[:] = 1
    >>> grid.add_zeros("fuel__moisture", at="cell")[:] = 0.08
    >>> fs = FireSpread(grid, ignition_row=20, ignition_col=30)
    >>> for i in range(300):
    ...     fs.run_one_step(dt=60)
    """

    _name = "FireSpread"
    _unit_agnostic = False
    _time_units = "s"
    _cite_as = """@techreport{rothermel1972mathematical,
      author = {Rothermel, Richard C},
      title = {A mathematical model for predicting fire spread in wildland fuels},
      institution = {USDA Forest Service},
      year = {1972}
    }"""

    _info = {
        "fuel__model": {"dtype": int, "intent": "in", "mapping": "cell"},
        "fuel__moisture": {"dtype": float, "intent": "in", "mapping": "cell"},
        "wind__speed": {"dtype": float, "intent": "in", "mapping": "cell"},
        "wind__direction": {"dtype": float, "intent": "in", "mapping": "cell"},
        "topographic__slope_steepness": {"dtype": float, "intent": "in", "mapping": "cell"},
        "fire__arrival_time": {"dtype": float, "intent": "out", "mapping": "cell"},
        "fire__flame_length": {"dtype": float, "intent": "out", "mapping": "cell"},
        "fire__reaction_intensity": {"dtype": float, "intent": "out", "mapping": "cell"},
    }

    def __init__(
        self,
        grid: RasterModelGrid,
        ignition_row: int | None = None,
        ignition_col: int | None = None,
        ignition_cells: list[int] | None = None,
        ignition_time: float = 0.0,
        **kwds,
    ):
        if not isinstance(grid, RasterModelGrid):
            raise TypeError("FireSpread requires a RasterModelGrid")

        super().__init__(grid, **kwds)

        self._time = ignition_time
        self._dx = grid.dx

        # Output fields
        self.arrival = grid.add_field("fire__arrival_time", -9999.0, at="cell", dtype=float)
        self.flame = grid.add_zeros("fire__flame_length", at="cell", dtype=float)
        self.intensity = grid.add_zeros("fire__reaction_intensity", at="cell", dtype=float)

        # Set ignition
        if ignition_cells:
            self.arrival[np.array(ignition_cells)] = ignition_time
        elif ignition_row is not None and ignition_col is not None:
            cell = grid.grid_coords_to_node_id(ignition_row, ignition_col, "cell")
            self.arrival[cell] = ignition_time
        else:
            raise ValueError("Specify ignition_row/col or ignition_cells")

        # Neighbor offsets for 4-connectivity
        ny, nx = grid.shape
        self._offsets = np.array([[-1, 0], [1, 0], [0, -1], [0, 1]])

    def run_one_step(self, dt: float = 60.0):
        """Advance fire by dt seconds."""
        arrival = self.arrival
        burned = arrival >= 0

        if not np.any(burned):
            return

        # Compute ROS at all currently burning cells
        fuel = self.grid.at_cell["fuel__model"][burned].astype(int)
        Mf = self.grid.at_cell["fuel__moisture"][burned]
        U = self.grid.at_cell.get("wind__speed", np.zeros_like(burned, dtype=float))[burned]
        tan_phi = self.grid.at_cell.get("topographic__slope_steepness", np.zeros_like(burned))[burned]

        ros = np.zeros_like(Mf)
        for fm in np.unique(fuel):
            mask = fuel == fm
            ros[mask] = rothermel_rate_of_spread(fm, Mf[mask], U[mask], tan_phi[mask], ANDERSON_13)

        # Candidate frontier cells
        frontier = set()
        for dy, dx in self._offsets:
            ny, nx = self.grid.shape
            idx = np.where(burned)[0]
            rows = idx // nx + dy
            cols = idx % nx + dx
            valid = (rows >= 0) & (rows < ny) & (cols >= 0) & (cols < nx)
            frontier.update(rows[valid] * nx + cols[valid])

        for cell in frontier:
            if arrival[cell] >= 0:
                continue  # already burned

            # Find fastest neighbor
            nbrs = self.grid.neighbors_at_cell[cell]
            valid = (nbrs != -1) & (arrival[nbrs] >= 0)
            if not np.any(valid):
                continue

            nbr_times = arrival[nbrs[valid]]
            nbr_ros = ros[nbrs[valid] - np.where(burned)[0][0]]  # rough mapping
            best_idx = np.argmin(nbr_times)
            t_arrive = nbr_times[best_idx] + self._dx / nbr_ros[best_idx]

            if arrival[cell] < 0 or t_arrive < arrival[cell]:
                arrival[cell] = t_arrive
                Ir = reaction_intensity(fuel[0], Mf[0], None, ANDERSON_13)  # simplified
                self.flame[cell] = byram_flame_length(Ir)
                self.intensity[cell] = Ir

        self._time += dt

    @property
    def current_time(self):
        return self._time

    @property
    def burned_area_m2(self):
        return np.sum(self.arrival >= 0) * self._dx ** 2
