#! /usr/bin/env python
import numpy as np

from .. import LinearDiffuser
from .shoreline import find_shoreline


class SubmarineDiffuser(LinearDiffuser):

    def __init__(self, grid, shore, sl_array, dictionary, **kwds):
        self._ksh = dictionary['ksh']
        self._sl_array   = sl_array
        self._sea_level = sl_array[0]
        self._wavebase  = dictionary['wave_base']
        self._hgt       = dictionary['hgt']
        self._alpha     = dictionary['alpha']
        self._sl_sh     = dictionary['sl_sh']
        self._load      = dictionary['load']
        kwds['linear_diffusivity'] = 'kd'
        self._model_step = 0
        grid.add_zeros('kd', at='node')

        super(SubmarineDiffuser, self).__init__(grid, **kwds)

    @property
    def k0(self):
        return self._k0

    @property
    def k_land(self):
        return self._ksh

    @property
    def sea_level(self):
        return self._sea_level

    @sea_level.setter
    def sea_level(self, sea_level):
        self._sea_level = sea_level

    def calc_diffusion_coef(self,shore):
        self.sea_level = self._sl_array[self._model_step]
        water_depth = (self.sea_level -
                       self._grid.at_node['topographic__elevation'])
        wavebase = self._wavebase
        under_water = water_depth > 0.
        deep_water = water_depth > wavebase

        k = self._grid.at_node['kd']
    

        x = self._grid.x_of_node
        x_sh = shore
        a = dx = self.grid.dx
        b = (self._hgt*self._alpha + self._sl_sh)*dx
        
        k_land = self._load/.0008 #self._sl_sh
        k[under_water] = self._load * ((x[under_water]-x_sh) +a) / (
                                    water_depth[under_water] + b) 
        k[deep_water] = k[deep_water] * np.exp(-1*(water_depth[deep_water]-
                                                wavebase)/wavebase)
        k[~ under_water] = k_land

        return k

    def run_one_step(self, dt):
        z = self._grid.at_node['topographic__elevation'] - self.sea_level #+ subsidence_array
        shore = find_shoreline(self.grid.x_of_node[self.grid.core_nodes], 
                               z[self.grid.core_nodes], self.sea_level)
        self.calc_diffusion_coef(shore)
        super(SubmarineDiffuser, self).run_one_step(dt)
        self._model_step = self._model_step + 1
