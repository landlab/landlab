#!/usr/bin/env python3
"""
Created on Fri May 29 08:41:50 2026

@author: matheuswanderleydealmeida
"""

import warnings

warnings.filterwarnings("ignore")


class ErosionChanger:

    def __init__(
        self,
        grid,
        vegetation="fuel_availability",
        rivers_unburn=None,
        K_sed0=1e-5,
        K_boost=1e2,
    ):
        self.grid = grid
        self.vegetation = grid.at_node[vegetation]
        self.rivers_unburn = rivers_unburn
        self.K_sed0 = K_sed0
        self.K_boost = K_boost

    def erodibility_update(self, ha):
        updated_K = self.K_sed0 * (self.K_boost ** (1 - self.vegetation))

        if self.rivers_unburn is not None:
            updated_K[self.rivers_unburn] = self.K_sed0

        ha.K_sed = updated_K

    def run_one_step(self, ha, dt):
        self.erodibility_update(ha)


class DiffusionChanger:

    def __init__(
        self, grid, vegetation="fuel_availability", diff_0=0.1, diffusion_boost=10
    ):
        self.grid = grid
        self.vegetation = grid.at_node[vegetation]
        self.diff_0 = diff_0
        self.diffusion_boost = diffusion_boost

    def diffusion_update(self, nld):
        veg_link = self.grid.map_value_at_max_node_to_link(
            "topographic__elevation", "fuel_availability"
        )
        nld._K = self.diff_0 * (self.diffusion_boost ** (1 - veg_link))

    def run_one_step(self, nld, dt):
        self.diffusion_update(nld)
