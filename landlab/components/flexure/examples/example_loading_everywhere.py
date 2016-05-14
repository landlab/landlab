#! /usr/bin/env python

import numpy as np

from landlab.components.flexure import Flexure
from landlab import RasterModelGrid


def get_random_load_magnitudes(n_loads):
    return np.random.normal(1e3, 10e6, n_loads)


def put_loads_on_grid(grid, load_sizes):
    load = grid.at_node['lithosphere__overlying_pressure_increment']
    load[:] = load_sizes

    #for (loc, size) in zip(load_locations, load_sizes):
    #    load.flat[loc] = size


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--shape', type=int, default=200,
                        help='Number rows and columns')
    parser.add_argument('--spacing', type=int, default=5e3,
                        help='Spading between rows and columns (m)')
    parser.add_argument('--n-procs', type=int, default=1,
                        help='Number of processors to use')
    parser.add_argument('--plot', action='store_true', default=False,
                        help='Plot an image of the total deflection')

    args = parser.parse_args()

    shape = (args.shape, args.shape)
    spacing = (args.spacing, args.spacing)

    load_sizes = get_random_load_magnitudes(args.shape * args.shape)

    grid = RasterModelGrid(shape[0], shape[1], spacing[0])

    flex = Flexure(grid, method='flexure')

    put_loads_on_grid(grid, load_sizes)

    flex.update(n_procs=args.n_procs)

    if args.plot:
        grid.imshow('node', 'lithosphere_surface__elevation_increment',
                    symmetric_cbar=False, cmap='spectral', show=True)

if __name__ == '__main__':
    main()
