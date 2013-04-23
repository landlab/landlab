#! /usr/env/python

#-----------------------------------------------------------------------
#
# diffusion2Dtest.py: Tests ModelGrid() class by implementing a 2D
# diffusion code.
#
# GT, July 2010
#-----------------------------------------------------------------------

import numpy
from numpy import *

import model_grid
from model_grid import *

import pylab
from pylab import plot, draw, show, contour

# Initialize

# User-defined parameters
nr = 130     # number of rows
nc = 3     # number of columns
dx = 0.25   # cell spacing
s = 0.001  # source
k = 0.01   # flux coefficient
dt = 0.25*dx**2.0/k   # time step
print 'dt=',dt
nt = 40000    # number of iterations

# Setup
mg = RasterModelGrid()
mg.initialize( nr, nc, dx )
u = mg.create_cell_dvector()
interior_cells = mg.get_interior_cells()
dudt = mg.create_cell_dvector()


# Boundaries
mg.set_noflux_boundaries( False, True, False, True )

#for i in range( 0, len( u ) ):
#	u[i] = mg.x(i)
#print 'u=',u

# Process

for i in range( 0, nt ):

	print i
	g = mg.calculate_face_gradients( u )
	q = -k*g
	#print 'g=',g
	#print 'q=',q
	dqds = mg.calculate_flux_divergences( q )
	#print 'dqds=',dqds
	for c in interior_cells:
		dudt[c] = s - dqds[c]
	u = u + dudt * dt
	u = mg.update_noflux_boundaries( u )
	#print 'new u:',u
	
	#plot( x, u )
	#draw()
	
	
# Finalize
ur = mg.cell_vector_to_raster( u )
#x = 0.5*dx + dx*arange( 0, nr )
#x = 0.5*dx + dx*arange( 0, nc )
#plot( x, ur[1,:] )
print 'Max u = ',max(u)
contour( ur )
show()


	