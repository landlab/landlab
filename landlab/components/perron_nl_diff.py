import numpy
from collections import deque
import itertools
import scipy.sparse as sparse
import landlab.components.craters as craters #Note this brings in the data class

import time

#these ones only so we can run this module ad-hoc:
#import pylab
from pylab import plot, draw, show, contour, imshow, colorbar
#from copy import copy

#Things to add: 1. Explicit stability check. 2. Implicit handling of scenarios where kappa*dt exceeds critical step - subdivide dt automatically. 3. Don't use loops.

class perron_nl_diff_faster(object):
    '''
    This module uses Taylor Perron's implicit (2011) method to solve the nonlinear hillslope diffusion equation across a rectangular grid for a single timestep. Note it works with the mass flux implicitly, and thus does not actually calculate it. Grid must be at least 5x5.
    Built DEJH early June 2013.
    '''
    def __init__(self, grid, data, dt):
        self._delta_t = dt
        self._uplift = 0.
        self._rock_density = 2.7
        self._sed_density = 2.7
        self._kappa = 5.e-5 # ==_a
        self._S_crit = 32.*numpy.pi/180.
        self._delta_x = grid.get_grid_spacing()
        self._delta_y = self._delta_x
        self._one_over_delta_x = 1./self._delta_x
        self._one_over_delta_y = 1./self._delta_y
        self._one_over_delta_x_sqd = self._one_over_delta_x**2.
        self._one_over_delta_y_sqd = self._one_over_delta_y**2.
        self._b = 1./self._S_crit**2.
        #self._grid = grid
        #self._data = data
        
        ncols = grid.get_count_of_cols()
        nrows = grid.get_count_of_rows()
        nnodes = grid.get_count_of_all_nodes()
        #Superceded by next text block (faster)
#        self._interior_cells = list(grid.get_interior_cells())
#        _core_cells = self._interior_cells[:] #we'll thin this down into the other lists
#        _interior_corners = []
#        _interior_edges = []
#        _interior_corners.append(_core_cells.pop())
#        _interior_edges.extend(_core_cells[(4-ncols):])
#        _core_cells[(4-ncols):] = []
#        _interior_corners.append(_core_cells.pop())
#        
#        _interior_corners.append(_core_cells.pop(0))
#        _interior_edges.extend(_core_cells[:(ncols-4)])
#        _core_cells[:(ncols-4)] = []
#        _interior_corners.append(_core_cells.pop(0))
#        
#        _interior_edges.append(_core_cells.pop(0))
#        i=1
#        while 1:
#            try:
#                _interior_edges.append(_core_cells.pop(i*(ncols-4)))
#            except:
#                break
#            try:
#                _interior_edges.append(_core_cells.pop(i*(ncols-4)))
#            except:
#                break
#            else:
#                i = i+1
    
        self._interior_corners = numpy.array(ncols+1,2*ncols-2,nnodes-2*ncols+1,nnodes-ncols-2)
        _left_list = range(2*ncols+1,nnodes-2*ncols,ncols)
        _right_list = range(3*ncols-2,nnodes-2*ncols,ncols)
        _left_right = list(itertools.chain.fromiterable(itertools.izip(_left_list, _right_list)))
        self._interior_edges = numpy.array(itertools.chain(range((ncols+2),(2*ncols-2)),_left_right,range(nnodes-2*ncols+2,nnodes-ncols-2)))
        self._core_cells = numpy.zeros(grid.get_count_of_interior_nodes()-2*ncols-4-2*(nrows-4),dtype=int)
        for i in range(nrows-4):
            _core_cells[i*(ncols-4):(i+1)*(ncols-4)] = range((2+i)*ncols+2,(3+i)*ncols-2) #fill IDs into core cell matrix
        
        #self._core_cells = numpy.array(_core_cells)
        #self._interior_edges = numpy.array(_interior_edges) #order is [ncols-4 of TOP]+[ncols-4 of BOTTOM]+[(nrows-4)*(LEFT,RIGHT) pairs]
        #self._interior_corners = numpy.array(_interior_corners) #order is topright,topleft,bottomleft,bottomright
        #print _interior_corners
        #print _interior_edges
        #print _core_cells
        #setup is giving the correct interior cell lists here


    def set_variables(self, grid, data):
        '''
        This function sets the variables needed for update().
        '''
        #self._grid = grid
        #self._data = data
        n_interior_cells = grid.get_count_of_interior_cells()
        _operating_matrix = sparse.lil_matrix((n_interior_cells, n_interior_cells), dtype=float)
        #_interior_elevs = [-1] * n_interior_cells

        #Initialize the local builder lists
        _mat_RHS = numpy.zeros(n_interior_cells)
        
        self.set_variables_for_core_cells(grid, data, _operating_matrix, _mat_RHS)
        self.set_variables_for_interior_edges(grid, data, _operating_matrix, _mat_RHS)
        self.set_variables_for_interior_corners(grid, data, _operating_matrix, _mat_RHS)
    
        self._operating_matrix = _operating_matrix.tocsr()
        self._mat_RHS = _mat_RHS


    def set_variables_for_core_cells(self, grid, data, _operating_matrix, _mat_RHS):
    
        elev = data.elev
        ncols = grid.get_count_of_cols()
        _delta_t = self._delta_t
        _one_over_delta_x = self._one_over_delta_x
        _one_over_delta_x_sqd = self._one_over_delta_x_sqd
        _one_over_delta_y = self._one_over_delta_y
        _one_over_delta_y_sqd = self._one_over_delta_y_sqd
        _kappa = self._kappa
        _b = self._b
        _S_crit = self._S_crit
        count = 0
        interior_grid_width = ncols-2
        core_cell_width = ncols-4
        for i in self._core_cells:
            n = (count//core_cell_width+1)*interior_grid_width + (count%core_cell_width) + 1 #This is the ID within the interior grid
            #n_test = self._interior_cells.index(i)
            #assert n == n_test
            cell_neighbors = grid.get_neighbor_list(i)
            cell_diagonals = grid.get_diagonal_list(i)
            _z_x = (data.elev[cell_neighbors[0]]-data.elev[cell_neighbors[2]])*0.5*_one_over_delta_x
            _z_y = (data.elev[cell_neighbors[1]]-data.elev[cell_neighbors[3]])*0.5*_one_over_delta_y
            _z_xx = (data.elev[cell_neighbors[0]]-2.*data.elev[i]+data.elev[cell_neighbors[2]])*_one_over_delta_x_sqd
            _z_yy = (data.elev[cell_neighbors[1]]-2.*data.elev[i]+data.elev[cell_neighbors[3]])*_one_over_delta_y_sqd
            _z_xy = (data.elev[cell_diagonals[0]] - data.elev[cell_diagonals[1]] - data.elev[cell_diagonals[3]] + data.elev[cell_diagonals[2]])*0.25*_one_over_delta_x*_one_over_delta_y
            _d = 1./(1.-_b*(_z_x*_z_x+_z_y*_z_y))
            
            #assert type(_z_x) is numpy.float64
            #assert type(_z_y) is numpy.float64
            #assert type(_z_xx) is numpy.float64
            #assert type(_z_yy) is numpy.float64
            #assert type(_z_xy) is numpy.float64
            
            _abd_sqd = _kappa*_b*_d*_d
            _F_ij = -2.*_kappa*_d*(_one_over_delta_x_sqd+_one_over_delta_y_sqd) - 4.*_abd_sqd*(_z_x*_z_x*_one_over_delta_x_sqd+_z_y*_z_y*_one_over_delta_y_sqd)
            _F_ijminus1 = _kappa*_d*_one_over_delta_x_sqd - _abd_sqd*_z_x*(_z_xx+_z_yy)*_one_over_delta_x - 4.*_abd_sqd*_b*_d*(_z_x*_z_x*_z_xx+_z_y*_z_y*_z_yy+2.*_z_x*_z_y*_z_xy)*_z_x*_one_over_delta_x - 2.*_abd_sqd*(_z_x*_z_xx*_one_over_delta_x-_z_x*_z_x*_one_over_delta_x_sqd+_z_y*_z_xy*_one_over_delta_x)
            _F_ijplus1 = _kappa*_d*_one_over_delta_x_sqd + _abd_sqd*_z_x*(_z_xx+_z_yy)*_one_over_delta_x + 4.*_abd_sqd*_b*_d*(_z_x*_z_x*_z_xx+_z_y*_z_y*_z_yy+2.*_z_x*_z_y*_z_xy)*_z_x*_one_over_delta_x + 2.*_abd_sqd*(_z_x*_z_xx*_one_over_delta_x+_z_x*_z_x*_one_over_delta_x_sqd+_z_y*_z_xy*_one_over_delta_x)
            _F_iminus1j = _kappa*_d*_one_over_delta_y_sqd - _abd_sqd*_z_y*(_z_xx+_z_yy)*_one_over_delta_y - 4.*_abd_sqd*_b*_d*(_z_x*_z_x*_z_xx+_z_y*_z_y*_z_yy+2.*_z_x*_z_y*_z_xy)*_z_y*_one_over_delta_y - 2.*_abd_sqd*(_z_y*_z_yy*_one_over_delta_y-_z_y*_z_y*_one_over_delta_y_sqd+_z_x*_z_xy*_one_over_delta_y)
            _F_iplus1j = _kappa*_d*_one_over_delta_y_sqd + _abd_sqd*_z_y*(_z_xx+_z_yy)*_one_over_delta_y + 4.*_abd_sqd*_b*_d*(_z_x*_z_x*_z_xx+_z_y*_z_y*_z_yy+2.*_z_x*_z_y*_z_xy)*_z_y*_one_over_delta_y + 2.*_abd_sqd*(_z_y*_z_yy*_one_over_delta_y+_z_y*_z_y*_one_over_delta_y_sqd+_z_x*_z_xy*_one_over_delta_y)
            _F_iplus1jplus1 = _abd_sqd*_z_x*_z_y*_one_over_delta_x*_one_over_delta_y
            _F_iminus1jminus1 = _F_iplus1jplus1
            _F_iplus1jminus1 = -_F_iplus1jplus1
            _F_iminus1jplus1 = _F_iplus1jminus1
            
            #assert type(_F_ij) is numpy.float64
            #assert type(_F_ijminus1) is numpy.float64
            #assert type(_F_ijplus1) is numpy.float64
            #assert type(_F_iminus1j) is numpy.float64
            #assert type(_F_iplus1j) is numpy.float64
            #assert type(_F_iplus1jplus1) is numpy.float64

            #RHS of equ 6 (see para [20])
            _func_on_z = self._rock_density/self._sed_density*self._uplift + _kappa*((_z_xx+_z_yy)/(1.-(_z_x*_z_x+_z_y*_z_y)/_S_crit*_S_crit) + 2.*(_z_x*_z_x*_z_xx+_z_y*_z_y*_z_yy+2.*_z_x*_z_y*_z_xy)/(_S_crit*_S_crit*(1.-(_z_x*_z_x+_z_y*_z_y)/_S_crit*_S_crit)**2.))
                        
            #assert type(_func_on_z) is numpy.float64

            _mat_RHS[n] = _mat_RHS[n] + data.elev[i] + _delta_t*(_func_on_z - (_F_ij*data.elev[i]+_F_ijminus1*data.elev[cell_neighbors[2]]+_F_ijplus1*data.elev[cell_neighbors[0]]+_F_iminus1j*data.elev[cell_neighbors[3]]+_F_iplus1j*data.elev[cell_neighbors[1]]+_F_iminus1jminus1*data.elev[cell_diagonals[2]]+_F_iplus1jplus1*data.elev[cell_diagonals[0]]+_F_iplus1jminus1*data.elev[cell_diagonals[1]]+_F_iminus1jplus1*data.elev[cell_diagonals[3]]))
            
            #build the operating matrix. No logic operations needed as these are all core cells
            _operating_matrix[n,(n-1):(n+2)] += [-_delta_t*_F_ijminus1, 1.-_delta_t*_F_ij, -_delta_t*_F_ijplus1]
            #_operating_matrix[n,n] = _operating_matrix[n,n]+1.-_delta_t*_F_ij
            #_operating_matrix[n,n-1] = _operating_matrix[n,n-1]-_delta_t*_F_ijminus1
            #_operating_matrix[n,n+1] = _operating_matrix[n,n+1]-_delta_t*_F_ijplus1
            _operating_matrix[n,(n-ncols+1):(n-ncols+4)] -= [_F_iminus1jminus1, _F_iminus1j, _F_iminus1jplus1]*_delta_t
            #_operating_matrix[n,n-ncols+2] = _operating_matrix[n,n-ncols+2]-_delta_t*_F_iminus1j
            #_operating_matrix[n,n-ncols+1] = _operating_matrix[n,n-ncols+1]-_delta_t*_F_iminus1jminus1
            #_operating_matrix[n,n-ncols+3] = _operating_matrix[n,n-ncols+3]-_delta_t*_F_iminus1jplus1
            _operating_matrix[n,(n+ncols-3):(n+ncols)] -= [_F_iplus1jminus1, _F_iplus1j, _F_iplus1jplus1]*_delta_t
            #_operating_matrix[n,n+ncols-2] = _operating_matrix[n,n+ncols-2]-_delta_t*_F_iplus1j
            #_operating_matrix[n,n+ncols-1] = _operating_matrix[n,n+ncols-1]-_delta_t*_F_iplus1jplus1
            #_operating_matrix[n,n+ncols-3] = _operating_matrix[n,n+ncols-3]-_delta_t*_F_iplus1jminus1
        
            count = count + 1


    def set_variables_for_interior_corners(self, grid, data, _operating_matrix, _mat_RHS):
        ncols = grid.get_count_of_cols()
        count = 0
        corner_ids = [-1,-(ncols-2),0,ncols-3] #topright,topleft,bottomleft,bottomright. This is in the ID reference frame of the INTERIOR GRID
        for i in self._interior_corners:
            #print i
            n = corner_ids[count]
            cell_neighbors = grid.get_neighbor_list(i)
            cell_diagonals = grid.get_diagonal_list(i)
            _z_x = (data.elev[cell_neighbors[0]]-data.elev[cell_neighbors[2]])*0.5*self._one_over_delta_x
            _z_y = (data.elev[cell_neighbors[1]]-data.elev[cell_neighbors[3]])*0.5*self._one_over_delta_y
            _z_xx = (data.elev[cell_neighbors[0]]-2.*data.elev[i]+data.elev[cell_neighbors[2]])*self._one_over_delta_x_sqd
            _z_yy = (data.elev[cell_neighbors[1]]-2.*data.elev[i]+data.elev[cell_neighbors[3]])*self._one_over_delta_y_sqd
            _z_xy = (data.elev[cell_diagonals[0]] - data.elev[cell_diagonals[1]] - data.elev[cell_diagonals[3]] + data.elev[cell_diagonals[2]])*0.25*self._one_over_delta_x*self._one_over_delta_y
            _d = 1./(1.-self._b*(_z_x**2.+_z_y**2.))
            
            #assert type(_z_x) is numpy.float64
            #assert type(_z_y) is numpy.float64
            #assert type(_z_xx) is numpy.float64
            #assert type(_z_yy) is numpy.float64
            #assert type(_z_xy) is numpy.float64
            
            _abd_sqd = self._kappa*self._b*_d**2.
            _F_ij = -2.*self._kappa*_d*(self._one_over_delta_x_sqd+self._one_over_delta_y_sqd) - 4.*_abd_sqd*(_z_x**2.*self._one_over_delta_x_sqd+_z_y**2.*self._one_over_delta_y_sqd)
            _F_ijminus1 = self._kappa*_d*self._one_over_delta_x_sqd - _abd_sqd*_z_x*(_z_xx+_z_yy)*self._one_over_delta_x - 4.*_abd_sqd*self._b*_d*(_z_x**2.*_z_xx+_z_y**2.*_z_yy+2.*_z_x*_z_y*_z_xy)*_z_x*self._one_over_delta_x - 2.*_abd_sqd*(_z_x*_z_xx*self._one_over_delta_x-_z_x**2.*self._one_over_delta_x_sqd+_z_y*_z_xy*self._one_over_delta_x)
            _F_ijplus1 = self._kappa*_d*self._one_over_delta_x_sqd + _abd_sqd*_z_x*(_z_xx+_z_yy)*self._one_over_delta_x + 4.*_abd_sqd*self._b*_d*(_z_x**2.*_z_xx+_z_y**2.*_z_yy+2.*_z_x*_z_y*_z_xy)*_z_x*self._one_over_delta_x + 2.*_abd_sqd*(_z_x*_z_xx*self._one_over_delta_x+_z_x**2.*self._one_over_delta_x_sqd+_z_y*_z_xy*self._one_over_delta_x)
            _F_iminus1j = self._kappa*_d*self._one_over_delta_y_sqd - _abd_sqd*_z_y*(_z_xx+_z_yy)*self._one_over_delta_y - 4.*_abd_sqd*self._b*_d*(_z_x**2.*_z_xx+_z_y**2.*_z_yy+2.*_z_x*_z_y*_z_xy)*_z_y*self._one_over_delta_y - 2.*_abd_sqd*(_z_y*_z_yy*self._one_over_delta_y-_z_y**2.*self._one_over_delta_y_sqd+_z_x*_z_xy*self._one_over_delta_y)
            _F_iplus1j = self._kappa*_d*self._one_over_delta_y_sqd + _abd_sqd*_z_y*(_z_xx+_z_yy)*self._one_over_delta_y + 4.*_abd_sqd*self._b*_d*(_z_x**2.*_z_xx+_z_y**2.*_z_yy+2.*_z_x*_z_y*_z_xy)*_z_y*self._one_over_delta_y + 2.*_abd_sqd*(_z_y*_z_yy*self._one_over_delta_y+_z_y**2.*self._one_over_delta_y_sqd+_z_x*_z_xy*self._one_over_delta_y)
            _F_iplus1jplus1 = _abd_sqd*_z_x*_z_y*self._one_over_delta_x*self._one_over_delta_y
            _F_iminus1jminus1 = _F_iplus1jplus1
            _F_iplus1jminus1 = -_F_iplus1jplus1
            _F_iminus1jplus1 = _F_iplus1jminus1
            
            #assert type(_F_ij) is numpy.float64
            #assert type(_F_ijminus1) is numpy.float64
            #assert type(_F_ijplus1) is numpy.float64
            #assert type(_F_iminus1j) is numpy.float64
            #assert type(_F_iplus1j) is numpy.float64
            #assert type(_F_iplus1jplus1) is numpy.float64

            #RHS of equ 6 (see para [20])
            _func_on_z = self._rock_density/self._sed_density*self._uplift + self._kappa*((_z_xx+_z_yy)/(1.-(_z_x**2.+_z_y**2.)/self._S_crit**2.) + 2.*(_z_x**2.*_z_xx+_z_y**2.*_z_yy+2.*_z_x*_z_y*_z_xy)/(self._S_crit**2.*(1.-(_z_x**2.+_z_y**2.)/self._S_crit**2.)**2.))
                        
            #assert type(_func_on_z) is numpy.float64

            _mat_RHS[n] = _mat_RHS[n] + data.elev[i] + self._delta_t*(_func_on_z - (_F_ij*data.elev[i]+_F_ijminus1*data.elev[cell_neighbors[2]]+_F_ijplus1*data.elev[cell_neighbors[0]]+_F_iminus1j*data.elev[cell_neighbors[3]]+_F_iplus1j*data.elev[cell_neighbors[1]]+_F_iminus1jminus1*data.elev[cell_diagonals[2]]+_F_iplus1jplus1*data.elev[cell_diagonals[0]]+_F_iplus1jminus1*data.elev[cell_diagonals[1]]+_F_iminus1jplus1*data.elev[cell_diagonals[3]]))
            
            #build the operating matrix
            _operating_matrix[n,n] = _operating_matrix[n,n]+1.-self._delta_t*_F_ij
            if count==0: #topright
                _operating_matrix[n,n-1] = _operating_matrix[n,n-1]-self._delta_t*_F_ijminus1 #left
                _operating_matrix[n,n-ncols+2] = _operating_matrix[n,n-ncols+2]-self._delta_t*_F_iminus1j #below
                _operating_matrix[n,n-ncols+1] = _operating_matrix[n,n-ncols+1]-self._delta_t*_F_iminus1jminus1 #leftbelow
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, +ncols-1, _F_iplus1jminus1) #aboveleft
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, ncols, _F_iplus1j) #above
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, +ncols+1, _F_iplus1jplus1) #aboveright
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, 1, _F_ijplus1) #right
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, -ncols+1, _F_iminus1jplus1) #belowright
            if count==1: #topleft
                _operating_matrix[n,n+1] = _operating_matrix[n,n+1]-self._delta_t*_F_ijplus1 #right
                _operating_matrix[n,n-ncols+2] = _operating_matrix[n,n-ncols+2]-self._delta_t*_F_iminus1j #below
                _operating_matrix[n,n-ncols+3] = _operating_matrix[n,n-ncols+3]-self._delta_t*_F_iminus1jplus1 #belowright
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, +ncols-1, _F_iplus1jminus1) #aboveleft
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, ncols, _F_iplus1j) #above
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, +ncols+1, _F_iplus1jplus1) #aboveright
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, -1, _F_ijminus1) #left
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, -ncols-1, _F_iminus1jminus1) #belowleft
            if count==2: #bottomleft
                _operating_matrix[n,n+1] = _operating_matrix[n,n+1]-self._delta_t*_F_ijplus1 #right
                _operating_matrix[n,n+ncols-2] = _operating_matrix[n,n+ncols-2]-self._delta_t*_F_iplus1j #above
                _operating_matrix[n,n+ncols-1] = _operating_matrix[n,n+ncols-1]-self._delta_t*_F_iplus1jplus1 #aboveright
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, +ncols-1, _F_iplus1jminus1) #aboveleft
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, -1, _F_ijminus1) #left
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, -ncols-1, _F_iminus1jminus1) #belowleft
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, -ncols, _F_iminus1j) #below
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, -ncols+1, _F_iminus1jplus1) #belowright
            if count==3: #bottomright
                _operating_matrix[n,n-1] = _operating_matrix[n,n-1]-self._delta_t*_F_ijminus1 #left
                _operating_matrix[n,n+ncols-2] = _operating_matrix[n,n+ncols-2]-self._delta_t*_F_iplus1j #above
                _operating_matrix[n,n+ncols-3] = _operating_matrix[n,n+ncols-3]-self._delta_t*_F_iplus1jminus1 #aboveleft
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, +ncols+1, _F_iplus1jplus1) #aboveright
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, 1, _F_ijplus1) #right
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, -ncols+1, _F_iminus1jplus1) #belowright
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, -ncols, _F_iminus1j) #below
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, -ncols-1, _F_iminus1jminus1) #belowleft
            count = count+1
                

    def set_variables_for_interior_edges(self, grid, data, _operating_matrix, _mat_RHS):
        count = 0
        side_IDs = []
        nrows = grid.get_count_of_rows()
        ncols = grid.get_count_of_cols()
        for i in range(nrows-4):
            side_IDs.extend([(1+i)*(ncols-2), (2+i)*(ncols-2)-1])
        edge_ids = range(-(ncols-3),-1)+range(1,(ncols-3))+side_IDs #order is [ncols-4 of TOP]+[ncols-4 of BOTTOM]+[(nrows-4)*(LEFT,RIGHT) pairs]; ID to the INTERIOR GRID
        #print edge_ids
        for i in self._interior_edges:
            #print '---'
            #print count
            #print i
            n = edge_ids[count]
            #print n
            cell_neighbors = grid.get_neighbor_list(i)
            cell_diagonals = grid.get_diagonal_list(i)
            _z_x = (data.elev[cell_neighbors[0]]-data.elev[cell_neighbors[2]])*0.5*self._one_over_delta_x
            _z_y = (data.elev[cell_neighbors[1]]-data.elev[cell_neighbors[3]])*0.5*self._one_over_delta_y
            _z_xx = (data.elev[cell_neighbors[0]]-2.*data.elev[i]+data.elev[cell_neighbors[2]])*self._one_over_delta_x_sqd
            _z_yy = (data.elev[cell_neighbors[1]]-2.*data.elev[i]+data.elev[cell_neighbors[3]])*self._one_over_delta_y_sqd
            _z_xy = (data.elev[cell_diagonals[0]] - data.elev[cell_diagonals[1]] - data.elev[cell_diagonals[3]] + data.elev[cell_diagonals[2]])*0.25*self._one_over_delta_x*self._one_over_delta_y
            _d = 1./(1.-self._b*(_z_x*_z_x+_z_y*_z_y))
            
            #assert type(_z_x) is numpy.float64
            #assert type(_z_y) is numpy.float64
            #assert type(_z_xx) is numpy.float64
            #assert type(_z_yy) is numpy.float64
            #assert type(_z_xy) is numpy.float64
            
            _abd_sqd = self._kappa*self._b*_d*_d
            _F_ij = -2.*self._kappa*_d*(self._one_over_delta_x_sqd+self._one_over_delta_y_sqd) - 4.*_abd_sqd*(_z_x*_z_x*self._one_over_delta_x_sqd+_z_y*_z_y*self._one_over_delta_y_sqd)
            _F_ijminus1 = self._kappa*_d*self._one_over_delta_x_sqd - _abd_sqd*_z_x*(_z_xx+_z_yy)*self._one_over_delta_x - 4.*_abd_sqd*self._b*_d*(_z_x*_z_x*_z_xx+_z_y*_z_y*_z_yy+2.*_z_x*_z_y*_z_xy)*_z_x*self._one_over_delta_x - 2.*_abd_sqd*(_z_x*_z_xx*self._one_over_delta_x-_z_x*_z_x*self._one_over_delta_x_sqd+_z_y*_z_xy*self._one_over_delta_x)
            _F_ijplus1 = self._kappa*_d*self._one_over_delta_x_sqd + _abd_sqd*_z_x*(_z_xx+_z_yy)*self._one_over_delta_x + 4.*_abd_sqd*self._b*_d*(_z_x*_z_x*_z_xx+_z_y*_z_y*_z_yy+2.*_z_x*_z_y*_z_xy)*_z_x*self._one_over_delta_x + 2.*_abd_sqd*(_z_x*_z_xx*self._one_over_delta_x+_z_x*_z_x*self._one_over_delta_x_sqd+_z_y*_z_xy*self._one_over_delta_x)
            _F_iminus1j = self._kappa*_d*self._one_over_delta_y_sqd - _abd_sqd*_z_y*(_z_xx+_z_yy)*self._one_over_delta_y - 4.*_abd_sqd*self._b*_d*(_z_x*_z_x*_z_xx+_z_y*_z_y*_z_yy+2.*_z_x*_z_y*_z_xy)*_z_y*self._one_over_delta_y - 2.*_abd_sqd*(_z_y*_z_yy*self._one_over_delta_y-_z_y*_z_y*self._one_over_delta_y_sqd+_z_x*_z_xy*self._one_over_delta_y)
            _F_iplus1j = self._kappa*_d*self._one_over_delta_y_sqd + _abd_sqd*_z_y*(_z_xx+_z_yy)*self._one_over_delta_y + 4.*_abd_sqd*self._b*_d*(_z_x*_z_x*_z_xx+_z_y*_z_y*_z_yy+2.*_z_x*_z_y*_z_xy)*_z_y*self._one_over_delta_y + 2.*_abd_sqd*(_z_y*_z_yy*self._one_over_delta_y+_z_y*_z_y*self._one_over_delta_y_sqd+_z_x*_z_xy*self._one_over_delta_y)
            _F_iplus1jplus1 = _abd_sqd*_z_x*_z_y*self._one_over_delta_x*self._one_over_delta_y
            _F_iminus1jminus1 = _F_iplus1jplus1
            _F_iplus1jminus1 = -_F_iplus1jplus1
            _F_iminus1jplus1 = _F_iplus1jminus1
            
            #assert type(_F_ij) is numpy.float64
            #assert type(_F_ijminus1) is numpy.float64
            #assert type(_F_ijplus1) is numpy.float64
            #assert type(_F_iminus1j) is numpy.float64
            #assert type(_F_iplus1j) is numpy.float64
            #assert type(_F_iplus1jplus1) is numpy.float64

            #RHS of equ 6 (see para [20])
            _func_on_z = self._rock_density/self._sed_density*self._uplift + self._kappa*((_z_xx+_z_yy)/(1.-(_z_x*_z_x+_z_y*_z_y)/self._S_crit*self._S_crit) + 2.*(_z_x*_z_x*_z_xx+_z_y*_z_y*_z_yy+2.*_z_x*_z_y*_z_xy)/(self._S_crit*self._S_crit*(1.-(_z_x*_z_x+_z_y*_z_y)/self._S_crit*self._S_crit)**2.))
                        
            #assert type(_func_on_z) is numpy.float64

            _mat_RHS[n] = _mat_RHS[n] + data.elev[i] + self._delta_t*(_func_on_z - (_F_ij*data.elev[i]+_F_ijminus1*data.elev[cell_neighbors[2]]+_F_ijplus1*data.elev[cell_neighbors[0]]+_F_iminus1j*data.elev[cell_neighbors[3]]+_F_iplus1j*data.elev[cell_neighbors[1]]+_F_iminus1jminus1*data.elev[cell_diagonals[2]]+_F_iplus1jplus1*data.elev[cell_diagonals[0]]+_F_iplus1jminus1*data.elev[cell_diagonals[1]]+_F_iminus1jplus1*data.elev[cell_diagonals[3]]))
            
            #build the operating matrix
            _operating_matrix[n,n] = _operating_matrix[n,n]+1.-self._delta_t*_F_ij
            if count<ncols-4: #top
                _operating_matrix[n,n-1] = _operating_matrix[n,n-1]-self._delta_t*_F_ijminus1 #left
                _operating_matrix[n,n+1] = _operating_matrix[n,n+1]-self._delta_t*_F_ijplus1 #right
                _operating_matrix[n,n-ncols+2] = _operating_matrix[n,n-ncols+2]-self._delta_t*_F_iminus1j #below
                _operating_matrix[n,n-ncols+1] = _operating_matrix[n,n-ncols+1]-self._delta_t*_F_iminus1jminus1 #leftbelow
                _operating_matrix[n,n-ncols+3] = _operating_matrix[n,n-ncols+3]-self._delta_t*_F_iminus1jplus1 #belowright
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, +ncols-1, _F_iplus1jminus1) #aboveleft
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, ncols, _F_iplus1j) #above
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, +ncols+1, _F_iplus1jplus1) #aboveright
            elif count<2*(ncols-4): #bottom
                _operating_matrix[n,n-1] = _operating_matrix[n,n-1]-self._delta_t*_F_ijminus1 #left
                _operating_matrix[n,n+1] = _operating_matrix[n,n+1]-self._delta_t*_F_ijplus1 #right
                _operating_matrix[n,n+ncols-2] = _operating_matrix[n,n+ncols-2]-self._delta_t*_F_iplus1j #above
                _operating_matrix[n,n+ncols-1] = _operating_matrix[n,n+ncols-1]-self._delta_t*_F_iplus1jplus1 #aboveright
                _operating_matrix[n,n+ncols-3] = _operating_matrix[n,n+ncols-3]-self._delta_t*_F_iplus1jminus1 #aboveleft
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, -ncols-1, _F_iminus1jminus1) #belowleft
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, -ncols, _F_iminus1j) #below
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, -ncols+1, _F_iminus1jplus1) #belowright
            elif (count-2*(ncols-4))%2: #right
                _operating_matrix[n,n-1] = _operating_matrix[n,n-1]-self._delta_t*_F_ijminus1 #left
                _operating_matrix[n,n+ncols-2] = _operating_matrix[n,n+ncols-2]-self._delta_t*_F_iplus1j #above
                _operating_matrix[n,n+ncols-3] = _operating_matrix[n,n+ncols-3]-self._delta_t*_F_iplus1jminus1 #aboveleft
                _operating_matrix[n,n-ncols+2] = _operating_matrix[n,n-ncols+2]-self._delta_t*_F_iminus1j #below
                _operating_matrix[n,n-ncols+1] = _operating_matrix[n,n-ncols+1]-self._delta_t*_F_iminus1jminus1 #leftbelow
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, +ncols+1, _F_iplus1jplus1) #aboveright
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, 1, _F_ijplus1) #right
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, -ncols+1, _F_iminus1jplus1) #belowright
            else: #left
                _operating_matrix[n,n+1] = _operating_matrix[n,n+1]-self._delta_t*_F_ijplus1 #right
                _operating_matrix[n,n+ncols-2] = _operating_matrix[n,n+ncols-2]-self._delta_t*_F_iplus1j #above
                _operating_matrix[n,n+ncols-1] = _operating_matrix[n,n+ncols-1]-self._delta_t*_F_iplus1jplus1 #aboveright
                _operating_matrix[n,n-ncols+2] = _operating_matrix[n,n-ncols+2]-self._delta_t*_F_iminus1j #below
                _operating_matrix[n,n-ncols+3] = _operating_matrix[n,n-ncols+3]-self._delta_t*_F_iminus1jplus1 #belowright
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, +ncols-1, _F_iplus1jminus1) #aboveleft
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, -1, _F_ijminus1) #left
                self.set_bc_cell(grid, data, _operating_matrix, _mat_RHS, n, -ncols-1, _F_iminus1jminus1) #belowleft

            count = count+1
        
    #This needs updating to reflect Greg's new BC handling in model_grid
    def set_bc_cell(self, grid, data, _operating_matrix, _mat_RHS, n, modulator, _F_value):
        bc_pos = grid.boundary_cells.tolist().index(self._interior_cells[n]+modulator)
        if grid.default_bc.boundary_code[bc_pos] == 3: #tracks cell
            #Find the cell it's linked to in *this* internal-cells-only matrix:
            l = self._interior_cells.tolist().index(grid.default_bc.tracks_cell[bc_pos])
            _operating_matrix[n,l] = _operating_matrix[n,l]-self._delta_t*_F_value
        elif grid.default_bc.boundary_code[bc_pos] == 2: #fixed gradient
            #Which cell is it linked to?
            l = self._interior_cells.tolist().index(grid.default_bc.tracks_cell[bc_pos])
            #One part goes in this cell...
            _operating_matrix[n,l] = _operating_matrix[n,l]-self._delta_t*_F_value
            #And one part goes on the RHS
            _mat_RHS[n] = _mat_RHS[n]+self._delta_t*_F_value*self._delta_x*grid.default_bc.boundary_gradient[bc_pos]
        elif grid.default_bc.boundary_code[bc_pos] == 1: #fixed value
            #Goes to the RHS
            _mat_RHS[n] = _mat_RHS[n]+self._delta_t*_F_value*data.elev[self._interior_cells[n]+modulator]
        else:
            try:
                raise IndexError
            except IndexError:
                print "Boundary code was not of an expected value!"
                raise


    def update(self, grid, data_in):
        #Initialize the variables for the step:
        self.set_variables(grid, data_in)
        #Solve interior of grid:
        _interior_elevs = sparse.linalg.spsolve(self._operating_matrix, self._mat_RHS)
        #print _interior_elevs.shape
        #this fn solves Ax=B for x
        
        #Load the new elevs into the actual elev data
        j=0 #counter for location in the matrix _interior_elevs
        ncells = grid.get_number_of_nodes()
        for i in range(ncells):
            if grid.is_interior(i):
                data_in.elev[i] = _interior_elevs[j]
                j = j+1
            else:
                bc_pos = grid.boundary_cells.tolist().index(i)
                if grid.default_bc.boundary_code[bc_pos] == 3:
                    data_in.elev[i] = _interior_elevs[0,self._interior_cells.tolist().index(grid.default_bc.tracks_cell[bc_pos])]
                elif grid.default_bc.boundary_code[bc_pos] == 2:
                    data_in.elev[i] = _interior_elevs[0,self._interior_cells.tolist().index(grid.default_bc.tracks_cell[bc_pos])] + grid.default_bc.gradient[bc_pos]*self._delta_x
                elif grid.default_bc.boundary_code[bc_pos] == 1:
                    continue #The old value of elev is still correct
                else:
                    try:
                        raise IndexError
                    except IndexError:
                        print "Boundary code was not of an expected value!"
                        raise


def main():

    dt = 1.
    nt = 500
    cr, mg, vectors = craters.dig_one_crater(120,120, 0.025, 0.5, 0.5, 1.)
    start_time = time.time()
    diff = perron_nl_diff_faster(mg, vectors, dt)
    for i in range (nt):
        diff.update(mg, vectors)
        print "Loop number ", i+1, " completed!"

    #Finalize
    elev_raster = mg.cell_vector_to_raster(vectors.elev)

    end_time = time.time()
    print('Elapsed time was %g seconds' % (end_time - start_time))

    #contour(elev_raster)
    flipped_elev_raster = numpy.empty_like(elev_raster)
    for i in range(0,mg.nrows):
        flipped_elev_raster[i,:] = elev_raster[(mg.nrows-i-1),:]
    imshow(flipped_elev_raster)
    colorbar()
    show()
    vectors.viewing_raster = flipped_elev_raster


if __name__=='__main__':
    main()