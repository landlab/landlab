#! /usr/env/python
"""
diffusion2Dtest.py: Tests ModelGrid() class by implementing a 2D
diffusion code.

GT, July 2010
"""
import numpy as np
from pylab import plot, draw, show, contour

from landlab import RasterModelGrid


def set_flux_coefficients(mg, dx):
    """
    Sets the spatial pattern of "K" (diffusion coefficient at cell
    faces)
    """
    # here we set it up such that K ranges from 0 at the edge of
    # cell 1 (left edge) to 1 at right edge of last cell (x=L)
    Kmax = 1.0
    Kexp = 1.0
    yf = mg.get_face_y_coords()
    xf = mg.get_face_x_coords()
    #L = self.dx * (self.nnodes-1)
    #self.K = Kmax * ((self.x[:self.nnodes-1] + self.dx/2) / L)**2.0

    # increases downward
    #self.K = Kmax * (((self.dx*(self.nr - 1)) - yf) / (self.dx*self.nr))**Kexp

    # increases upward
    #self.K = Kmax * ((yf-self.dx/2) / (self.dx*(self.nr-2)))**Kexp

    # increases to left
    #self.K = Kmax * ((self.dx*(self.nc - 2) - (xf-self.dx/2)) /
    #                 (self.dx*(self.nc-2)))**Kexp

    # increases to right
    #self.K = Kmax * ((xf-self.dx/2) / (self.dx*(self.nc-2)))**Kexp

    # Exponential decline from a point at center top
    decay_scale = 0.5
    xf = xf - dx/2.0
    yf = yf - dx/2.0
    x0 = 0.5 * (mg.number_of_node_columns - 2) * dx
    y0 = (mg.number_of_node_columns - 2) * dx
    dist = np.sqrt((xf - x0) ** 2.0 + (yf - y0) ** 2.0)
    K = Kmax * np.exp(-dist / decay_scale)

    if False:
        print 'dist:', dist
        print 'yf:', yf
        print 'xf:', xf
        print 'K:', K

    return K


def main():
    # Initialize

    # User-defined parameters
    nr = 12     # number of rows
    nc = 12     # number of columns
    dx = 0.5   # cell spacing
    s = 0.001  # source
    #k0 = 0.01   # flux coefficient

    # Setup
    mg = RasterModelGrid(nr, nc, dx)
    u = mg.zeros(centering='cell')
    interior_cells = mg.get_interior_cells()
    dudt = mg.zeros(centering='cell')
    #K = mg.zeros(centering='cell')
    fx = mg.get_face_x_coords()
    fy = mg.get_face_y_coords()
    #k = k0 * (fy / max(fy))
    #k = k0 * (fx / max(fx))
    #print 'k=',k

    opt_plot = True

    # A "channel"
    # k[:] = k0/10.0
    # for i in xrange(0, len(k)):
    #     if fx[i] > 14 and fx[i] < 16 and fy[i] > 14:
    #         k[i] = k0

    # Radially symmetric, exponential decay
    k = set_flux_coefficients(mg, dx)

    dt = 0.25*dx**2.0/max(k)   # time step
    #print 'dt=',dt
    run_time = 100.0
    nt = int(round(run_time / dt))    # number of iterations

    #ri = raw_input('frog')

    # Boundaries
    mg.set_noflux_boundaries(True, True, False, True)

    #for i in range(0, len(u)):
    #    u[i] = mg.x(i)
    #print 'u=',u

    # Process

    for i in range(0, nt):
        print i
        g = mg.calculate_face_gradients(u)
        q = -k*g
        #print 'g=',g
        #print 'q=',q
        dqds = mg.calculate_flux_divergences(q)
        #print 'dqds=',dqds
        for c in interior_cells:
            dudt[c] = s - dqds[c]
        u = u + dudt * dt
        u = mg.update_noflux_boundaries(u)
        #print 'new u:',u

        #plot(x, u)
        #draw()

    # Finalize
    print 'Max u = ', max(u)
    if opt_plot:
        ur = mg.cell_vector_to_raster(u)
        #x = 0.5*dx + dx*arange(0, nr)
        #x = 0.5*dx + dx*arange(0, nc)
        #plot(x, ur[1,:])
        contour(ur)
        show()


if __name__ == '__main__':
    main()
