# landlab/components/fire_spread/fire_spread.py

from __future__ import annotations

import numpy as np

from landlab import Component
from landlab import RasterModelGrid

from .fuel_models import ANDERSON_13
from .utils import byram_flame_length
from .utils import reaction_intensity
from .utils import rothermel_rate_of_spread


class FireSpread(Component):
    """Rothermel (1972) surface fire spread on a RasterModelGrid."""

    _name = "FireSpread"
    _time_units = "s"

    _info = {
        "fuel__model": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "mapping": "cell",
            "doc": "Fuel model number (Anderson 13)",
        },
        "fuel__moisture": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "mapping": "cell",
            "doc": "Dead fuel moisture fraction",
        },
        "wind__speed": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "mapping": "cell",
            "doc": "Mid-flame wind speed (m/s)",
        },
        "wind__direction": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "mapping": "cell",
        },
        "topographic__slope_steepness": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "mapping": "cell",
            "doc": "tan(slope angle)",
        },
        "fire__arrival_time": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "mapping": "cell",
            "doc": "Time fire arrives at cell (s)",
        },
        "fire__flame_length": {
            "dtype": float,
            "intent": "out",
            "optional": True,
            "mapping": "cell",
            "doc": "Byram flame length (m)",
        },
        "fire__reaction_intensity": {
            "dtype": float,
            "intent": "out",
            "optional": True,
            "mapping": "cell",
            "doc": "Reaction intensity (kW/m²)",
        },
    }

    def __init__(
        self,
        grid: RasterModelGrid,
        ignition_row: int | None = None,
        ignition_col: int | None = None,
        ignition_cells: list[int] | None = None,
        ignition_time: float = 0.0,
    ):
        super().__init__(grid)

        if not isinstance(grid, RasterModelGrid):
            raise TypeError("FireSpread requires a RasterModelGrid")

        self._dx = grid.dx
        self._time = ignition_time

        # current simulation time

        # Output fields
        self.arrival = grid.add_field(
            "fire__arrival_time",
            -9999.0 * np.ones(grid.number_of_cells),
            at="cell",
            dtype=float,
        )
        self.flame = grid.add_zeros("fire__flame_length", at="cell")
        self.intensity = grid.add_zeros("fire__reaction_intensity", at="cell")

        # === Set ignition ===
        if ignition_cells is not None:
            self.arrival[np.asarray(ignition_cells)] = ignition_time
        elif ignition_row is not None and ignition_col is not None:
            ny, nx = grid.shape
            node_id = ignition_row * nx + ignition_col
            cell_id = grid.cell_at_node[node_id]
            if cell_id == -1:
                raise ValueError("Ignition node is on the boundary (no cell there)")
            self.arrival[cell_id] = ignition_time
        else:
            raise ValueError("Specify ignition_row/col or ignition_cells")

        # Pre-compute 4-connected cell neighbours
        ny, nx = grid.shape
        cell_grid = np.arange(grid.number_of_cells).reshape(ny - 2, nx - 2)

        north = np.full_like(cell_grid, -1)
        south = np.full_like(cell_grid, -1)
        west = np.full_like(cell_grid, -1)
        east = np.full_like(cell_grid, -1)

        north[:-1, :] = cell_grid[1:, :]
        south[1:, :] = cell_grid[:-1, :]
        west[:, :-1] = cell_grid[:, 1:]
        east[:, 1:] = cell_grid[:, :-1]

        self._neighbors = np.stack(
            [north.ravel(), south.ravel(), west.ravel(), east.ravel()]
        )

    def run_one_step(self, dt: float = 60.0):
        """Advance fire front by dt seconds."""
        grid = self.grid
        arrival = self.arrival

        # Input fields (with defaults if missing)
        fuel = grid.at_cell["fuel__model"].astype(int)
        Mf = grid.at_cell["fuel__moisture"]
        U = grid.at_cell.get("wind__speed", np.zeros(grid.number_of_cells))
        tanφ = grid.at_cell.get(
            "topographic__slope_steepness", np.zeros(grid.number_of_cells)
        )

        # Rate of spread for every cell
        ros = np.zeros(grid.number_of_cells)
        for fm in np.unique(fuel):
            m = fuel == fm
            ros[m] = rothermel_rate_of_spread(fm, Mf[m], U[m], tanφ[m], ANDERSON_13)

        # Find candidate cells (have at least one burning neighbour)
        burning = arrival >= 0
        candidates = np.unique(self._neighbors[:, burning])
        candidates = candidates[candidates != -1]

        updated = False
        for cell in candidates:
            if arrival[cell] >= 0:  # already burned
                continue

            # burning neighbours of this cell
            nbrs = self._neighbors[:, cell]
            burning_nbrs = nbrs[(nbrs != -1) & burning[nbrs]]

            if len(burning_nbrs) == 0:
                continue

            travel = self._dx / np.maximum(ros[burning_nbrs], 1e-12)
            new_time = (arrival[burning_nbrs] + travel).min()

            if arrival[cell] < 0 or new_time < arrival[cell]:
                arrival[cell] = new_time

                # fire behaviour for the newly ignited cell
                p = ANDERSON_13[fuel[cell]]
                Ir = reaction_intensity(
                    p["w0"] * 20.83,  # tons/acre → kg/m²
                    p["sigma"],
                    p["h"] * 2326,  # BTU/lb → kJ/kg
                    Mf[cell],
                    None,
                )
                self.intensity[cell] = Ir
                self.flame[cell] = byram_flame_length(Ir)
                updated = True

        if updated:
            self._time += dt

    @property
    def current_time(self):
        return self._time

    @property
    def burned_area_m2(self):
        return np.count_nonzero(self.arrival >= 0) * self._dx**2
