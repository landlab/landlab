#! /usr/bin/env python

import numpy as np

from landlab.components.flexure import Flexure
from landlab import RasterModelGrid


SHAPE = (100, 100)
SPACING = (10e3, 10e3)
LOAD_LOCS = [
    (SHAPE[0] / 2, SHAPE[1] / 2),
    (SHAPE[0] / 4, SHAPE[1] * 3 / 4),
]


def put_two_point_loads_on_grid(grid):
    load = grid.field_values('node',
                             'lithosphere__overlying_pressure_increment')
    load = load.view()
    load.shape = grid.shape
    for loc in LOAD_LOCS:
        load[loc] = 10e6


def create_lithosphere_elevation_with_bulge(grid):
    grid.add_zeros('node', 'lithosphere_surface__elevation')

    z = grid.field_values('node', 'lithosphere_surface__elevation').view()
    z.shape = grid.shape

    (y, x) = np.meshgrid(np.linspace(0, np.pi * .5, grid.shape[0]),
                         np.linspace(0, np.pi * .5, grid.shape[1]))
    (x0, y0) = (np.pi / 3, np.pi / 8)
    np.sin((x - x0) ** 2 + (y - y0) ** 2, out=z)


def main():

    grid = RasterModelGrid(SHAPE[0], SHAPE[1], SPACING[0])

    create_lithosphere_elevation_with_bulge(grid)

    flex = Flexure(grid, method='flexure')

    put_two_point_loads_on_grid(grid)

    flex.update()

    grid.at_node['lithosphere_surface__elevation'] += grid.at_node['lithosphere_surface__elevation_increment']
    grid.imshow('node', 'lithosphere_surface__elevation',
                symmetric_cbar=False, show=True)


if __name__ == '__main__':
    main()
