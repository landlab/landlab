import numpy
import scipy.sparse as sparse
import landlab.components.craters as craters #Note this brings in the data class

#these ones only so we can run this module ad-hoc:
#import pylab
from pylab import plot, draw, show, contour, imshow, colorbar
#from copy import copy


class perron_nl_diff(object):
    '''
    This module uses Taylor Perron's implicit (2011) method to solve the nonlinear hillslope diffusion equation across a rectangular grid for a single timestep. Note it works with the mass flux implicitly, and thus does not actually calculate it.
    '''
    def __init__(self, grid, data, dt):
        self._delta_t = dt
        self._uplift = 0.
        self._rock_density = 2.7
        self._sed_density = 2.7
        self._kappa = 1.e-4 # ==_a
        self._S_crit = 32.*numpy.pi/180.
        self._delta_x = grid.dx
        self._delta_y = grid.dx
        self._one_over_delta_x = 1./self._delta_x
        self._one_over_delta_y = 1./self._delta_y
        self._one_over_delta_x_sqd = self._one_over_delta_x**2.
        self._one_over_delta_y_sqd = self._one_over_delta_y**2.
        self._b = 1./self._S_crit**2.
        self._grid = grid
        self._data = data
        
        self._interior_cells = self._grid.get_interior_cells()
        self._func_on_z = [-1] * self._grid.ncells


    def make_fns(self):
        '''
        This method builds and returns the "F functions", and sets the related helper functions needed to run Perron's implicit method - excluding the actual ncellsxncells inverted matrix. It exists to minimise memory usage by the module. All excess 1xncells lists will get garbage collected when the method terminates.
        '''
        #Initialize the local builder lists
        _equ23_RHS = []        
        _F_ = [numpy.zeros((3,3))] * self._grid.ncells

        for i in self._interior_cells:
            cell_neighbors = self._grid.get_neighbor_list(i)
            cell_diagonals = self._grid.get_diagonal_list(i)
            _z_x = (self._data.elev[cell_neighbors[0]]-self._data.elev[cell_neighbors[2]])*0.5*self._one_over_delta_x
            _z_y = (self._data.elev[cell_neighbors[1]]-self._data.elev[cell_neighbors[3]])*0.5*self._one_over_delta_y
            _z_xx = (self._data.elev[cell_neighbors[0]]-2.*self._data.elev[i]+self._data.elev[cell_neighbors[2]])*self._one_over_delta_x_sqd
            _z_yy = (self._data.elev[cell_neighbors[1]]-2.*self._data.elev[i]+self._data.elev[cell_neighbors[3]])*self._one_over_delta_y_sqd
            _z_xy = (self._data.elev[cell_diagonals[0]] - self._data.elev[cell_diagonals[1]] - self._data.elev[cell_diagonals[3]] + self._data.elev[cell_diagonals[2]])*0.25*self._one_over_delta_x*self._one_over_delta_y
            _d = 1./(1.-self._b*(_z_x**2.+_z_y**2.))
            
            assert type(_z_x) is numpy.float64
            assert type(_z_y) is numpy.float64
            assert type(_z_xx) is numpy.float64
            assert type(_z_yy) is numpy.float64
            assert type(_z_xy) is numpy.float64
            
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
            
            assert type(_F_ij) is numpy.float64
            assert type(_F_ijminus1) is numpy.float64
            assert type(_F_ijplus1) is numpy.float64
            assert type(_F_iminus1j) is numpy.float64
            assert type(_F_iplus1j) is numpy.float64
            assert type(_F_iplus1jplus1) is numpy.float64

            #RHS of equ 6 (see para [20])
            self._func_on_z[i] = self._rock_density/self._sed_density*self._uplift + self._kappa*((_z_xx+_z_yy)/(1.-(_z_x**2.+_z_y**2.)/self._S_crit**2.) + 2.*(_z_x**2.*_z_xx+_z_y**2.*_z_yy+2.*_z_x*_z_y*_z_xy)/(self._S_crit**2.*(1.-(_z_x**2.+_z_y**2.)/self._S_crit**2.)**2.))
                        
            assert type(self._func_on_z[i]) is numpy.float64

            _equ23_RHS.append(self._data.elev[i] + self._delta_t*(self._func_on_z[i] - (_F_ij*self._data.elev[i]+_F_ijminus1*self._data.elev[cell_neighbors[2]]+_F_ijplus1*self._data.elev[cell_neighbors[0]]+_F_iminus1j*self._data.elev[cell_neighbors[3]]+_F_iplus1j*self._data.elev[cell_neighbors[1]]+_F_iminus1jminus1*self._data.elev[cell_diagonals[2]]+_F_iplus1jplus1*self._data.elev[cell_diagonals[0]]+_F_iplus1jminus1*self._data.elev[cell_diagonals[1]]+_F_iminus1jplus1*self._data.elev[cell_diagonals[3]])))
        
            _F_[i] = numpy.array([[_F_iminus1jminus1, _F_iminus1j, _F_iminus1jplus1], [_F_ijminus1, _F_ij, _F_ijplus1], [_F_iplus1jminus1, _F_iplus1j, _F_iplus1jplus1]])
        
        assert type(_equ23_RHS[0]) is numpy.float64
        assert len(_equ23_RHS) == len(self._interior_cells)

        #self._mat_RHS = numpy.matrix(_equ23_RHS) #Put this in set_variables
        

        return _F_, _equ23_RHS

    
    def set_variables(self, grid, data):
        self._grid = grid
        self._data = data
        self._operating_matrix = sparse.lil_matrix((self._grid.n_interior_cells, self._grid.n_interior_cells), dtype=float)
        self._interior_elevs = [-1] * self._grid.n_interior_cells
        
        _F_, _equ23_RHS = self.make_fns()

        #assemble the operation matrix
        #This is where the BoundaryConditions get complex...
        #print self._grid.n_interior_cells
        for n in range(self._grid.n_interior_cells): #this is for each cell in the grid_interior_cells[i]
        #1st row is all boundary
        #core rows are one boundary cell at either end
        #Note we have to add onto the existing value, as we might have already put a value in the cell before the iteration on i gets there
            for i in range(3):
                for j in range(3):
                    #print n
                    modulator = (i-1)*self._grid.ncols + (j-1)
                    #print modulator
                    if self._grid.is_interior(self._interior_cells[n]+modulator):
                        if modulator == 0:
                            self._operating_matrix[n,n] = self._operating_matrix[n,n]+1.-self._delta_t*_F_[self._interior_cells[n]][i,j]
                        else: #Need to move into "interior cell" indices
                            if modulator > 1:
                                self._operating_matrix[n,n+modulator-2] = self._operating_matrix[n,n+modulator-2]-self._delta_t*_F_[self._interior_cells[n]][i,j]
                            elif modulator < -1:
                                self._operating_matrix[n,n+modulator+2] = self._operating_matrix[n,n+modulator+2]-self._delta_t*_F_[self._interior_cells[n]][i,j]
                            else:
                                self._operating_matrix[n,n+modulator] = self._operating_matrix[n,n+modulator]-self._delta_t*_F_[self._interior_cells[n]][i,j]
                    else:
                        bc_pos = self._grid.boundary_cells.tolist().index(self._interior_cells[n]+modulator)
                        if self._grid.default_bc.boundary_code[bc_pos] == 3: #tracks cell
                            #Find the cell it's linked to in *this* internal-cells-only matrix:
                            l = self._interior_cells.tolist().index(self._grid.default_bc.tracks_cell[bc_pos])
                            #Put the value in that cell:
                            self._operating_matrix[n,l] = self._operating_matrix[n,l]-self._delta_t*_F_[self._interior_cells[n]][i,j]
                        elif self._grid.default_bc.boundary_code[bc_pos] == 2: #fixed gradient
                            #Which cell is it linked to?
                            l = self._interior_cells.tolist().index(self._grid.default_bc.tracks_cell[bc_pos])
                            #One part goes in this cell...
                            self._operating_matrix[n,l] = self._operating_matrix[n,l]-self._delta_t*_F_[self._interior_cells[n]][i,j]
                            #And one part goes on the RHS
                            _equ23_RHS[n] = _equ23_RHS[n]+self._delta_t*_F_[self._interior_cells[n]][i,j]*self._grid.dx*self._grid.default_bc.boundary_gradient[bc_pos]
                        elif self._grid.default_bc.boundary_code[bc_pos] == 1: #fixed value
                            #Goes to the RHS
                            _equ23_RHS[n] = _equ23_RHS[n]+self._delta_t*_F_[self._interior_cells[n]][i,j]*self._data.elev[self._interior_cells[n]+modulator]
                        else:
                            try:
                                raise IndexError
                            except IndexError:
                                print "Boundary code was not of an expected value!"
                                raise

        self._operating_matrix = self._operating_matrix.tocsr()
        self._mat_RHS = numpy.matrix(_equ23_RHS)


    def update(self, grid, data_in):
        #Initialize the variables for the step:
        self.set_variables(grid, data_in)
        #Solve interior of grid:
        _interior_elevs = sparse.linalg.spsolve(self._operating_matrix, self._mat_RHS)
        #this fn solves Ax=B for x
        
        #Load the new elevs into the actual elev data
        j=0 #counter for location in the matrix _interior_elevs
        for i in range(self._grid.ncells):
            if self._grid.is_interior(i):
                data_in.elev[i] = _interior_elevs[0,j]
                j = j+1
            else:
                bc_pos = self._grid.boundary_cells.tolist().index(i)
                if self._grid.default_bc.boundary_code[bc_pos] == 3:
                    data_in.elev[i] = _interior_elevs[0,self._interior_cells.tolist().index(self._grid.default_bc.tracks_cell[bc_pos])]
                elif self._grid.default_bc.boundary_code[bc_pos] == 2:
                    data_in.elev[i] = _interior_elevs[0,self._interior_cells.tolist().index(self._grid.default_bc.tracks_cell[bc_pos])] + self._grid.default_bc.gradient[bc_pos]*self._grid.dx
                elif self._grid.default_bc.boundary_code[bc_pos] == 1:
                    continue #The old value of elev is still correct
                else:
                    try:
                        raise IndexError
                    except IndexError:
                        print "Boundary code was not of an expected value!"
                        raise


def main():
    
    dt = 1.
    nt = 10
    cr, mg, vectors = craters.dig_one_crater(120,120, 0.025, 0.5, 0.5, 1.)
    diff = perron_nl_diff(mg, vectors, dt)
    for i in range (nt):
        diff.update(mg, vectors)
        print "Loop number ", i+1, " completed!"

    #Finalize
    elev_raster = mg.cell_vector_to_raster(vectors.elev)
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