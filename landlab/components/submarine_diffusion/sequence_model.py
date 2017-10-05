#! /usr/bin/env python
from .sea_level import SinusoidalSeaLevel, SeaLevelTimeSeries
from .subsidence import SubsidenceTimeSeries
from .submarine import SubmarineDiffuser
from ..flexure.sediment_flexure import SedimentFlexure
from .raster_model import RasterModel


class SequenceModel(RasterModel):

    DEFAULT_PARAMS = {
        'grid': {
            'shape': [3, 100],
            'spacing': 100.,
            'origin': 0.,
            'bc': {
                'top': 'closed',
                'bottom': 'closed',
            },
        },
        'clock': {
            'start': -20000.,
            'stop': 0.,
            'step': 100.,
        },
        'submarine_diffusion': {
            'ksh': 100.,
            'wave_base': 60.,
            'shelf_depth': 15.,
            'alpha': .0005,
            'shelf_slope': .001,
            'load': 3.,
        },
        'sea_level': {
            'amplitude': 10.,
            'wave_length': 1000.,
            'phase': 0.,
        },
        'subsidence': {
            'filepath': 'subsidence.csv',
        },
        'flexure': {
            'method': 'airy',
            'rho_mantle': 3300.,
        },
    }

    LONG_NAME = {
        'z': 'topographic__elevation',
        'z0': 'bedrock_surface__elevation',
    }

    def __init__(self, grid=None, clock=None, submarine_diffusion=None,
                 sea_level=None, subsidence=None, flexure=None):
        RasterModel.__init__(self, grid=grid, clock=clock)

        z0 = self.grid.add_empty('bedrock_surface__elevation', at='node')
        z = self.grid.add_empty('topographic__elevation', at='node')

        z0[:] = - .01 * self.grid.x_of_node + 10.
        z[:] = z0

        self.grid.layers.add(0.,
                             age=self.clock.start,
                             water_depth=-z[self.grid.core_nodes],
                             t0=0.)

        self._sea_level = SinusoidalSeaLevel(self.grid, start=clock['start'],
                                             **sea_level)

        self._subsidence = SubsidenceTimeSeries(self.grid, **subsidence)

        self._submarine_diffusion = SubmarineDiffuser(self.grid,
                                                      **submarine_diffusion)
        self._flexure = SedimentFlexure(self.grid, **flexure)

        self._components = (
            self._sea_level,
            self._subsidence,
            self._submarine_diffusion,
            self._flexure,
        )

    def advance_components(self, dt):
        for component in self._components:
            component.run_one_step(dt)

        dz = self.grid.at_node['sediment_deposit__thickness']
        water_depth = (self.grid.at_grid['sea_level__elevation'] -
                       self.grid.at_node['topographic__elevation'])

        self.grid.layers.add(dz[self.grid.node_at_cell],
                             age=self.clock.time,
                             water_depth=water_depth[self.grid.node_at_cell],
                             t0=dz[self.grid.node_at_cell].clip(0.))


if __name__ == '__main__':
    main()
