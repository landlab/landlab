#! /usr/bin/env python
import numpy as np

from .. import LinearDiffuser
from .shoreline import find_shoreline


class SubmarineDiffuser(LinearDiffuser):

    _name = 'Submarine Diffusion'

    _time_units = 'y'

    _input_var_names = (
        'topographic__elevation',
        'sea_level__elevation',
    )

    _output_var_names = (
        'topographic__elevation',
    )

    _var_units = {
        'topographic__elevation': 'm',
        'sea_level__elevation': 'm',
    }

    _var_mapping = {
        'topographic__elevation': 'node',
        'sea_level__elevation': 'grid',
    }

    _var_doc = {
        'topographic__elevation': 'land and ocean bottom elevation',
        'sea_level__elevation': 'Position of sea level',
    }

    def __init__(self, grid, sea_level=0., ksh=100., wave_base=60.,
                 shelf_depth=15., alpha=.0005, shelf_slope=.001, load=3.,
                 **kwds):
        """Diffuse the ocean bottom.

        Parameters
        ----------
        grid: RasterModelGrid
            A landlab grid.
        sea_level: float, optional
            The current sea level (m).
        ksh: float, optional
            Diffusion coefficient at the shoreline.
        wave_base: float, optional
            Wave base (m).
        shelf_depth: float, optional
            Water depth of the shelf/slope break (m).
        alpha: float, optional
            Some coefficient.
        shelf_slope: float, optional
            Slope of the shelf (m / m).
        load: float, optional
            Sediment load entering the profile.
        """
        self._ksh = float(ksh)
        self._wave_base = float(wave_base)
        self._shelf_depth = float(shelf_depth)
        self._alpha = float(alpha)
        self._shelf_slope = float(shelf_slope)
        self._load = float(load)
        self._sea_level = sea_level

        grid.add_zeros('kd', at='node')

        self._time = 0.

        kwds.setdefault('linear_diffusivity', 'kd')
        super(SubmarineDiffuser, self).__init__(grid, **kwds)

    @property
    def k0(self):
        return self._k0

    @property
    def k_land(self):
        return self._ksh

    @property
    def time(self):
        return self._time

    @property
    def sea_level(self):
        return self.grid.at_grid['sea_level__elevation']

    @sea_level.setter
    def sea_level(self, sea_level):
        self.grid.at_grid['sea_level__elevation'] = sea_level

    def calc_diffusion_coef(self, x_of_shore):
        sea_level = self.grid.at_grid['sea_level__elevation']
        water_depth = (sea_level -
                       self._grid.at_node['topographic__elevation'])

        under_water = water_depth > 0.
        deep_water = water_depth > self._wave_base
        land = ~ under_water

        k = self.grid.at_node['kd']

        x = self.grid.x_of_node
        b = (self._shelf_depth * self._alpha + self._shelf_slope) * self.grid.dx
        
        k[under_water] = self._load * (
            (x[under_water] - x_of_shore) + self.grid.dx) / (water_depth[under_water] + b) 

        k[deep_water] *= np.exp(- (water_depth[deep_water] - self._wave_base) /
                                self._wave_base)

        k[land] = self._load / .0008 # self._shelf_slope

        return k

    def run_one_step(self, dt):
        z_before = self.grid.at_node['topographic__elevation'].copy()

        shore = find_shoreline(self.grid.x_of_node[self.grid.node_at_cell], 
                               z_before[self.grid.node_at_cell], self.sea_level)
        self.calc_diffusion_coef(shore)

        super(SubmarineDiffuser, self).run_one_step(dt)

        z_after = self.grid.at_node['topographic__elevation']
        water_depth = self.grid.at_grid['sea_level__elevation'] - z_after
        dz = (z_after - z_before)[self.grid.node_at_cell]

        self.grid.layers.add(dz,
                             age=self._time,
                             water_depth=water_depth[self.grid.node_at_cell],
                             t0=dz.clip(0.))

        self._time += dt
