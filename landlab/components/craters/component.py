#!/usr/bin/env python

import numpy as np

from landlab.grid.raster import RasterModelGrid
from landlab import Component
from . import craters


class CratersComponent(Component):
    _name = 'Craters'

    _input_var_names = [
        'topographic__elevation',
    ]
    _output_var_names = [
        'topographic__elevation',
        'topographic__elevation_increment',
    ]

    _var_units = {
        'topographic__elevation': 'm',
        'topographic__elevation_increment': 'm',
    }

    def __init__(self, grid, **kwds):
        self._impactor_mean_size = kwds.pop('mean_size', .75)
        self._impactor_size_std_dev = kwds.pop('std_dev', .1)
        seed = kwds.pop('seed', None)

        super(CratersComponent, self).__init__(grid, **kwds)

        self._vectors = craters.data(self.grid)

        self._grid.add_field('node', 'topographic__elevation',
                             self._vectors.elev, units='m')

        self._grid.add_field('node', 'topographic__elevation_increment',
                             np.zeros(self._grid.shape, dtype=np.float),
                             units='m')

        if seed is not None:
            np.random.seed(seed)

    def update(self):
        self.dig_one_crater(self.next_impact_position(),
                            self.next_impactor_radius())

    def dig_one_crater(self, impact_loc, radius):
        """
        Dig a crater that impacts the surface at *impact_loc*, which is a
        tuple of an x, y position (in km) on the planet surface. *radius* is
        the size of the impactor in km.
        """
        cr = craters.impactor()

        cr._radius = radius
        cr.set_depth_from_size()
        cr.set_crater_volume()
        cr._xcoord, cr._ycoord = (impact_loc[1], impact_loc[0])

        z0 = self._grid.at_node['topographic__elevation'].copy()

        vertices_array = self.grid.get_nodes_around_point(impact_loc[1],
                                                          impact_loc[0])

        distances_to_vertices = []
        for node_id in vertices_array:
            try:
                distances_to_vertices.append(
                    np.sqrt((impact_loc[1] - self.grid.node_x[node_id]) ** 2. +
                            (impact_loc[0] - self.grid.node_y[node_id]) ** 2.))
            except IndexError:
                return

        cr.closest_node_index = vertices_array[np.argmin(distances_to_vertices)]
        cr.closest_node_elev = self._vectors.elev[cr.closest_node_index]
        cr._angle_to_horizontal = np.pi * 0.5 * 3. / 3.
        cr._azimuth_of_travel = np.pi * 1.5

        cr.set_crater_mean_slope_v2(self.grid, self._vectors)
        cr.set_elev_change_across_grid(self.grid, self._vectors)

        self._grid.at_node['topographic__elevation_increment'][:].flat = (
            self._grid.at_node['topographic__elevation'].flat - z0)

    def next_impactor(self):
        cr = craters.impactor()

        (impact_x, impact_y) = self.next_impact_position()
        cr._radius = self.next_impactor_radius()
        cr.set_depth_from_size()
        cr.set_crater_volume()
        cr._xcoord, cr._ycoord = (impact_x, impact_y)

        return cr

    def next_impact_position(self):
        return (np.random.uniform() * self.grid.get_grid_xdimension(),
                np.random.uniform() * self.grid.get_grid_ydimension())

    def next_impactor_radius(self, method='weibull'):
        assert(method in ['weibull', 'normal'])

        if method == 'weibull':
            self._weibull_a = 1.
            self._weibull_lambda = self._impactor_mean_size / 5
            radius = (np.random.weibull(self._weibull_a) *
                      self._weibull_lambda)
        else:
            radius = np.random.normal(self._impactor_mean_size,
                                      self._impactor_size_std_dev)
            while radius <= 0:
                radius = np.random.normal(self._impactor_mean_size,
                                          self._impactor_size_std_dev)

        return radius + .005


if __name__ == "__main__":
    import doctest
    doctest.testmod()
