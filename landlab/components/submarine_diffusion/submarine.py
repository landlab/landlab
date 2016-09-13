#! /usr/bin/env python
import numpy as np

from .. import LinearDiffuser


class SubmarineDiffuser(LinearDiffuser):

    def __init__(self, grid, k0=1., sea_level=0., k_land=1., **kwds):
        self._k0 = k0
        self._sea_level = sea_level
        self._k_land = k_land
        kwds['linear_diffusivity'] = 'kd'

        grid.add_zeros('kd', at='node')

        super(SubmarineDiffuser, self).__init__(grid, **kwds)

    @property
    def k0(self):
        return self._k0

    @property
    def k_land(self):
        return self._k_land

    @property
    def sea_level(self):
        return self._sea_level

    @sea_level.setter
    def sea_level(self, sea_level):
        self._sea_level = sea_level

    def calc_diffusion_coef(self):
        water_depth = (self.sea_level -
                       self._grid.at_node['topographic__elevation'])
        under_water = water_depth > 0.

        k = self._grid.at_node['kd']

        x = self._grid.x_of_node[under_water]
        x_sh = 0.
        a = b = 1e-6
        big_z = 1.

        k[under_water] = self.k0 * (1. /
                                    (water_depth[under_water] + b)) * big_z
        k[~ under_water] = self.k_land

        return k

    def run_one_step(self, dt):
        self.calc_diffusion_coef()
        super(SubmarineDiffuser, self).run_one_step(dt)
