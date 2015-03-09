#! /usr/env/python
"""
diffusion2Dtest.py: Tests ModelGrid() class by implementing a 2D
diffusion code.

version 3 same as 2 but uses source as a dvector

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

    # Exponential decline from a point at center top
    decay_scale = 0.5
    xf = xf - dx/2.0
    yf = yf - dx/2.0
    x0 = 0.5 * (mg.number_of_node_columns - 2) * dx
    y0 = (mg.number_of_node_columns - 2) * dx
    dist = np.sqrt((xf - x0) ** 2.0 + (yf - y0) ** 2.0)
    K = Kmax * np.exp(-dist / decay_scale)

    return K


def main():
    # Initialize

    # User-defined parameters
    nr = 12     # number of rows
    nc = 12     # number of columns
    dx = 0.5   # cell spacing
    s0 = 0.001  # source

    # Setup
    mg = RasterModelGrid(nr, nc, dx)  # Create the grid
    u = mg.zeros(centering='cell') # Dependent variable (temperature)
    interior_cells = mg.get_interior_cells()  # ID's of interior cells
    dudt = mg.zeros(centering='cell') # Rate of change of temperature
    s = mg.zeros(centering='cell') # Source term (e.g., radioactive decay)
    s[interior_cells] = s0  # Set source at interior cells
    opt_plot = True  # Option for plotting output

    # Radially symmetric, exponential decay
    k = set_flux_coefficients(mg, dx)       # Set up diffusion coefficients

    dt = 0.25*dx**2.0/max(k)                  # Set time step
    run_time = 100.0                          # Set run time
    nt = int(round(run_time / dt))        # number of iterations

    mg.set_noflux_boundaries(True, True, False, True)  # Boundary conditions

    #-------------------------------------------------------------------
    # Process

    for i in range(0, nt):

        print i
        g = mg.calculate_face_gradients(u)  # Thermal gradients
        q = -k*g  # Heat flux across faces
        dqds = mg.calculate_flux_divergences(q)  # Divergence of heat flux
        dudt = s - dqds  # Rate of change of temperature
        u = u + dudt * dt  # Update temperature field
        u = mg.update_noflux_boundaries(u)  # Update boundaries

    #-------------------------------------------------------------------
    # Finalize

    if opt_plot:
        ur = mg.cell_vector_to_raster(u)   # Get a raster temperature field
        contour(ur)
        show()


if __name__ == '__main__':
    main()
