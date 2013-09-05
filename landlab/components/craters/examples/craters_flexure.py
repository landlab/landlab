#! /usr/bin/env python

from landlab import RasterModelGrid
from landlab.components.craters.component import CratersComponent
import landlab.components.flexure as flexure


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--impact-count', type=int, default=128,
                        help='Number of impacts')
    parser.add_argument('--grid-size', type=int, default=128,
                        help='Number of grid rows and columns')
    parser.add_argument('--grid-spacing', type=float, default=.025,
                        help='Spacing between grid rows and columns')

    seed_parser = parser.add_mutually_exclusive_group()
    seed_parser.add_argument('--with-seed', type=int, default=1973,
                             help='Seed for random number generator')
    seed_parser.add_argument('--without-seed', action='store_true',
                             default=False,
                             help='Seed for random number generator')
    args = parser.parse_args()

    shape = (args.grid_size, args.grid_size)
    spacing = (args.grid_spacing, args.grid_spacing)
    impact_count = args.impact_count
    if args.without_seed:
        seed = None
    else:
        seed = args.with_seed

    grid = RasterModelGrid(shape[0], shape[1], spacing[0])

    craters_comp = CratersComponent(grid, seed=seed)
    flexure_comp = flexure.FlexureComponent(
        grid,
        map_vars={
            'lithosphere__elevation':
                'planet_surface__elevation',
            'planet_surface_sediment__deposition_increment':
                'planet_surface__elevation_increment'
        }
    )

    for _ in xrange(impact_count):
        craters_comp.update()
        flexure_comp.update()

    craters_comp.imshow('planet_surface__elevation', grid_units=('km', 'km'),
                        symmetric_cbar=True)


if __name__ == "__main__":
    main()
