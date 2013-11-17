#! /usr/env/python

#-----------------------------------------------------------------------
#
# diffusion2Dtest.py: Tests ModelGrid() class by implementing a 2D
# diffusion code.
#
# GT, July 2010
#-----------------------------------------------------------------------

import numpy as np

from landlab import RasterModelGrid

import pylab
from pylab import plot, draw, show, contour


def main():
    # Initialize

    # User-defined parameters
    nr = 130     # number of rows
    nc = 3     # number of columns
    dx = 0.25   # cell spacing
    s = 0.001  # source
    k = 0.01   # flux coefficient
    dt = 0.25*dx**2.0/k   # time step
    print 'dt=', dt
    nt = 40000    # number of iterations

    # Setup
    mg = RasterModelGrid(nr, nc, dx)
    u = mg.zeros(centering='cell')
    interior_cells = mg.get_interior_cells()
    dudt = mg.zeros(centering='cell')

    # Boundaries
    mg.set_noflux_boundaries(False, True, False, True)

    # Process
    for i in range(0, nt):
        print i
        g = mg.calculate_face_gradients(u)
        q = -k*g
        dqds = mg.calculate_flux_divergences(q)
        for c in interior_cells:
            dudt[c] = s - dqds[c]
        u = u + dudt * dt
        u = mg.update_noflux_boundaries(u)

    # Finalize
    ur = mg.cell_vector_to_raster(u)
    print 'Max u = ', max(u)
    contour(ur)
    show()


if __name__ == '__main__':
    main()
