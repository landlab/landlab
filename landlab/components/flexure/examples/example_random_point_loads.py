#! /usr/bin/env python

import pylab
import numpy as np

from landlab.components.flexure import Flexure

#_SHAPE = (200, 200)
#_SPACING = (5e3, 5e3)
#_N_LOADS = 128


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--n-loads', type=int, default=1,
                        help='Number of loads to apply')
    parser.add_argument('--shape', type=int, default=128,
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

    load_locs = np.random.random_integers(0, shape[0] * shape[1] - 1,
                                          args.n_loads)
    load_sizes = np.random.normal(1e9, 1e7, args.n_loads)

    flex = Flexure(shape, spacing, (0., 0.))
    load = flex['lithosphere__overlying_pressure']

    for (loc, size) in zip(load_locs, load_sizes):
        load.flat[loc] = size

    flex.update(n_procs=args.n_procs)

    if args.plot:
        pylab.imshow(flex['lithosphere__elevation'])
        pylab.colorbar()
        pylab.show()


if __name__ == '__main__':
    main()
