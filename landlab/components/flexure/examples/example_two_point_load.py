#! /usr/bin/env python

import numpy as np

from landlab.components.flexure import FlexureComponent
from landlab import RasterModelGrid
from landlab.plot import imshow_field


SHAPE = (100, 100)
SPACING = (10e3, 10e3)
LOAD_LOCS = [
    (SHAPE[0] / 2, SHAPE[1] / 2),
    (SHAPE[0] / 4, SHAPE[1] * 3 / 4),
]


def put_two_point_loads_on_grid(grid):
    load = grid.field_values('node', 'lithosphere__overlying_pressure')
    load = load.view()
    load.shape = grid.shape
    for loc in LOAD_LOCS:
        load[loc] = 1e15


def create_lithosphere_elevation_with_bulge(grid):
    grid.add_zeros('node', 'lithosphere__elevation')

    z = grid.field_values('node', 'lithosphere__elevation').view()
    z.shape = grid.shape

    (y, x) = np.meshgrid(np.linspace(0, np.pi * .5, grid.shape[0]),
                         np.linspace(0, np.pi * .5, grid.shape[1]))
    (x0, y0) = (np.pi / 3, np.pi / 8)
    np.sin((x - x0) ** 2 + (y - y0) ** 2, out=z)


def main():

    grid = RasterModelGrid(SHAPE[0], SHAPE[1], SPACING[0])

    create_lithosphere_elevation_with_bulge(grid)

    flex = FlexureComponent(grid, method='flexure')

    put_two_point_loads_on_grid(grid)

    flex.update()

    grid.imshow('node', 'lithosphere__elevation', symmetric_cbar=False,
                show=True) 


if __name__ == '__main__':
    main()
