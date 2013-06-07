import numpy
import craters

class perron_nl_diff(object):
    '''
    This module uses Taylor Perron's implicit (2011) method to solve the
    nonlinear hillslope diffusion equation across a rectangular grid for a
    single timestep. Note it works with the mass flux implicitly, and thus does
    not actually calculate it.
    '''
    def __init__(self):
#        self._current_elevs = []
#        self._elevs_next_tstep = []
#        self._operating_matrix = []
#        self._z_x = []
#        self._z_y = []
#        self._z_xx = []
#        self._z_yy = []
#        self._z_xy = []
#        self._b = numpy.nan
#        self._d = []
#        self._F_ij = []
#        self._F_ijminus1 = []
#        self._F_ijplus1 = []
#        self._F_iminus1j = []
#        self._F_iplus1j = []
#        self._F_iplus1jplus1 = []
        self._delta_t = numpy.nan
        self._uplift = 0.
        self._rock_density = 2.7
        self._sed_density = 2.7
        self._kappa = 1. # ==_a
        self._delta_x = numpy.nan
        self._delta_y = numpy.nan
        self._S_crit = 32.*numpy.pi/180.
            
    def initialize(self, grid, data):
        #self._current_elevs = numpy.matrix(data.elev)
        #self._elevs_next_tstep = numpy.zeros(grid.ncells)
        self._operating_matrix = numpy.zeros((grid.n_interior_cells, grid.n_interior_cells), dtype=float)
        self._delta_x = grid.dx
        self._delta_y = grid.dx
        one_over_delta_x = 1./self._delta_x
        one_over_delta_y = 1./self._delta_y
        one_over_delta_x_sqd = one_over_delta_x**2.
        one_over_delta_y_sqd = one_over_delta_y**2.
        self._delta_t = 100.
        self._b = 1./self._S_crit**2.
        self._z_x = [-1] * grid.ncells #ncells, or n_interior_cells?
        self._z_y = [-1] * grid.ncells
        self._z_xx = [-1] * grid.ncells
        self._z_yy = [-1] * grid.ncells
        self._z_xy = [-1] * grid.ncells
        self._d = [-1] * grid.ncells
        self._F_ij = [-1] * grid.ncells
        self._F_ijminus1 = [-1] * grid.ncells
        self._F_ijplus1 = [-1] * grid.ncells
        self._F_iminus1j = [-1] * grid.ncells
        self._F_iplus1j = [-1] * grid.ncells
        self._F_iplus1jplus1 = [-1] * grid.ncells
        self._F_iminus1jminus1 = [-1] * grid.ncells
        self._F_iplus1jminus1 = [-1] * grid.ncells
        self._F_iminus1jplus1 = [-1] * grid.ncells
        self._func_on_z = [-1] * grid.ncells
        self._equ23_RHS = []
        self._interior_cells = grid.get_interior_cells()
        self._interior_elevs = [-1] * grid.n_interior_cells
        
        for i in self._interior_cells: #What do we do about boundary cells?? Need them to be counted to give enough simultaneous equations.
            #Will need to add exceptions here
            cell_neighbors = grid.get_neighbor_list(i)
            cell_diagonals = grid.get_diagonal_list(i)
            self._z_x[i] = (data.elev[cell_neighbors[0]]-data.elev[cell_neighbors[2]])*0.5*one_over_delta_x
            self._z_y[i] = (data.elev[cell_neighbors[1]]-data.elev[cell_neighbors[3]])*0.5*one_over_delta_y
            self._z_xx[i] = (data.elev[cell_neighbors[0]]-2.*data.elev[i]+data.elev[cell_neighbors[2]])*one_over_delta_x_sqd
            self._z_yy[i] = (data.elev[cell_neighbors[1]]-2.*data.elev[i]+data.elev[cell_neighbors[3]])*one_over_delta_y_sqd
            self._z_xy[i] = (data.elev[cell_diagonals[0]] - data.elev[cell_diagonals[1]] - data.elev[cell_diagonals[3]] + data.elev[cell_diagonals[2]])*0.25*one_over_delta_x*one_over_delta_y
            self._d[i] = 1./(1.-self._b*(self._z_x[i]**2.+self._z_y[i]**2.))
            
            _abd_sqd = self._kappa*self._b*self._d[i]**2.
            self._F_ij[i] = -2.*self._kappa*self._d[i]*(one_over_delta_x_sqd+one_over_delta_y_sqd) - 4.*_abd_sqd*(self._z_x[i]**2.*one_over_delta_x_sqd+self._z_y[i]**2.*one_over_delta_y_sqd)
            self._F_ijminus1[i] = self._kappa*self._d[i]*one_over_delta_x_sqd - _abd_sqd*self._z_x[i]*(self._z_xx[i]+self._z_yy[i])*one_over_delta_x - 4.*_abd_sqd*self._b*self._d[i]*(self._z_x[i]**2.*self._z_xx[i]+self._z_y[i]**2.*self._z_yy[i]+2.*self._z_x[i]*self._z_y[i]*self._z_xy[i])*self._z_x[i]*one_over_delta_x - 2.*_abd_sqd*(self._z_x[i]*self._z_xx[i]*one_over_delta_x-self._z_x[i]**2.*one_over_delta_x_sqd+self._z_y[i]*self._z_xy[i]*one_over_delta_x)
            self._F_ijplus1[i] = self._kappa*self._d[i]*one_over_delta_x_sqd + _abd_sqd*self._z_x[i]*(self._z_xx[i]+self._z_yy[i])*one_over_delta_x + 4.*_abd_sqd*self._b*self._d[i]*(self._z_x[i]**2.*self._z_xx[i]+self._z_y[i]**2.*self._z_yy[i]+2.*self._z_x[i]*self._z_y[i]*self._z_xy[i])*self._z_x[i]*one_over_delta_x + 2.*_abd_sqd*(self._z_x[i]*self._z_xx[i]*one_over_delta_x+self._z_x[i]**2.*one_over_delta_x_sqd+self._z_y[i]*self._z_xy[i]*one_over_delta_x)
            self._F_iminus1j[i] = self._kappa*self._d[i]*one_over_delta_y_sqd - _abd_sqd*self._z_y[i]*(self._z_xx[i]+self._z_yy[i])*one_over_delta_y - 4.*_abd_sqd*self._b*self._d[i]*(self._z_x[i]**2.*self._z_xx[i]+self._z_y[i]**2.*self._z_yy[i]+2.*self._z_x[i]*self._z_y[i]*self._z_xy[i])*self._z_y[i]*one_over_delta_y - 2.*_abd_sqd*(self._z_y[i]*self._z_yy[i]*one_over_delta_y-self._z_y[i]**2.*one_over_delta_y_sqd+self._z_x[i]*self._z_xy[i]*one_over_delta_y)
            self._F_iplus1j[i] = self._kappa*self._d[i]*one_over_delta_y_sqd + _abd_sqd*self._z_y[i]*(self._z_xx[i]+self._z_yy[i])*one_over_delta_y + 4.*_abd_sqd*self._b*self._d[i]*(self._z_x[i]**2.*self._z_xx[i]+self._z_y[i]**2.*self._z_yy[i]+2.*self._z_x[i]*self._z_y[i]*self._z_xy[i])*self._z_y[i]*one_over_delta_y + 2.*_abd_sqd*(self._z_y[i]*self._z_yy[i]*one_over_delta_y+self._z_y[i]**2.*one_over_delta_y_sqd+self._z_x[i]*self._z_xy[i]*one_over_delta_y)
            self._F_iplus1jplus1[i] = _abd_sqd*self._z_x[i]*self._z_y[i]*one_over_delta_x*one_over_delta_y
            self._F_iminus1jminus1[i] = self._F_iplus1jplus1[i]
            self._F_iplus1jminus1[i] = -self._F_iplus1jplus1[i]
            self._F_iminus1jplus1[i] = self._F_iplus1jminus1[i]

            #RHS of equ 6 (see para [20])
            self._func_on_z[i] = self._rock_density/self._sed_density*self._uplift + self._kappa*((self._z_xx[i]+self._z_yy[i])/(1.-(self._z_x[i]**2.+self._z_y[i]**2.)/self._S_crit**2.) + 2.*(self._z_x[i]**2.*self._z_xx[i]+self._z_y[i]**2.*self._z_yy[i]+2.*self._z_x[i]*self._z_y[i]*self._z_xy[i])/(self._S_crit**2.*(1.-(self._z_x[i]**2.+self._z_y[i]**2.)/self._S_crit**2.)**2.))

            self._equ23_RHS.append(data.elev[i] + self._delta_t*(self._func_on_z[i] - (self._F_ij*data.elev[i]+self._F_ijminus1[i]*data.elev[cell_neighbors[2]]+self._F_ijplus1[i]*data.elev[cell_neighbors[0]]+self._F_iminus1j[i]*data.elev[cell_neighbors[3]]+self._F_iplus1j[i]*data.elev[cell_neighbors[1]]+self._F_iminus1jminus1[i]*data.elev[cell_diagonals[2]]+self._F_iplus1jplus1[i]*data.elev[cell_diagonals[0]]+self._F_iplus1jminus1[i]*data.elev[cell_diagonals[1]]+self._F_iminus1jplus1[i]*data.elev[cell_diagonals[3]])))
        
        assert len(self._equ23_RHS) == self._interior_cells

        self._mat_RHS = numpy.matrix(self._equ23_RHS)

        #assemble the operation matrix
        #This is where the BCs get complex...
        for i in range(grid.n_interior_cells): #this is for each cell in the grid_interior_cells[i]
        #Think through whether it's i as index, or self._interior_cells[i]
        #1st row is all boundary
        #core rows are one boundary cell at either end
        #Note we have to add onto the existing value, as we might have already put a value in the cell before the iteration on i gets there
        #Important to adjust the indices STILL
            self._operating_matrix[i,i] = self._operating_matrix[i,i]+1.-self._delta_t*self._F_ij[self._interior_cells[i]]
            if grid.is_interior(self._interior_cells[i]-1):
                self._operating_matrix[i,i-1] = self._operating_matrix[i,i-1]-self._delta_t*self._F_ijminus1[self._interior_cells[i]]
            elif bc.boundary_code[self._interior_cells[i]-1] == bc.TRACKS_CELL_BOUNDARY:
                #Find the cell it's linked to in *this* internal-cells-only matrix:
                j = self._interior_cells.index(bc.tracks_cell[self._interior_cells[i]-1])
                #Put the value in that cell:
                self._operating_matrix[i,j] = self._operating_matrix[i,j]-self._delta_t*self._F_ijminus1[self._interior_cells[i]]
            elif bc.boundary_code[self._interior_cells[i]-1] == bc.FIXED_GRADIENT_BOUNDARY:
                #Which cell is it linked to?
                j = self._interior_cells.index(bc.tracks_cell[self._interior_cells[i]-1])
                #One part goes in this cell...
                self._operating_matrix[i,j] = -self._delta_t*self._F_ijminus1[self._interior_cells[i]]
                #And one part goes on the RHS
                self._mat_RHS[i] = self._mat_RHS[i]+self._delta_t*self._F_ijminus1[self._interior_cells[i]]*grid.dx*bc.gradient[self._interior_cells[i]-1]
            else: #FIXED_VALUE_BOUNDARY
                #Goes to the RHS
                self._mat_RHS[i] = self._mat_RHS[i]+self._delta_t*self._F_ijminus1[self._interior_cells[i]]*data.elev[self._interior_cells[i]-1]
                
            if grid.is_interior(self._interior_cells[i]+1):
                self._operating_matrix[i,i+1] = self._operating_matrix[i,i+1]-self._delta_t*self._F_ijplus1[self._interior_cells[i]]
            elif bc.boundary_code[self._interior_cells[i]+1] == bc.TRACKS_CELL_BOUNDARY:
                j = self._interior_cells.index(bc.tracks_cell[self._interior_cells[i]+1])
                self._operating_matrix[i,j] = self._operating_matrix[i,j]-self._delta_t*self._F_ijplus1[self._interior_cells[i]]
            elif bc.boundary_code[self._interior_cells[i]+1] == bc.FIXED_GRADIENT_BOUNDARY:
                j = self._interior_cells.index(bc.tracks_cell[self._interior_cells[i]+1])
                self._operating_matrix[i,j] = -self._delta_t*self._F_ijplus1[self._interior_cells[i]]
                self._mat_RHS[i] = self._mat_RHS[i]+self._delta_t*self._F_ijplus1[self._interior_cells[i]]*grid.dx*bc.gradient[self._interior_cells[i]+1]
            else:
                self._mat_RHS[i] = self._mat_RHS[i]+self._delta_t*self._F_ijplus1[self._interior_cells[i]]*data.elev[self._interior_cells[i]+1]
                
            if grid.is_interior(self._interior_cells[i]-grid.ncols):
                self._operating_matrix[i,i-grid.ncols] = self._operating_matrix[i,i-grid.ncols]-self._delta_t*self._F_iminus1j[self._interior_cells[i]]
            elif bc.boundary_code[self._interior_cells[i]-grid.ncols] == bc.TRACKS_CELL_BOUNDARY:
                j = self._interior_cells.index(bc.tracks_cell[self._interior_cells[i]-grid.ncols])
                self._operating_matrix[i,j] = self._operating_matrix[i,j]-self._delta_t*self._F_iminus1j[self._interior_cells[i]]
            elif bc.boundary_code[self._interior_cells[i]-grid.ncols] == bc.FIXED_GRADIENT_BOUNDARY:
                j = self._interior_cells.index(bc.tracks_cell[self._interior_cells[i]-grid.ncols])
                self._operating_matrix[i,j] = -self._delta_t*self._F_iminus1j[self._interior_cells[i]]
                self._mat_RHS[i] = self._mat_RHS[i]+self._delta_t*self._F_iminus1j[self._interior_cells[i]]*grid.dx*bc.gradient[self._interior_cells[i]-grid.ncols]
            else:
                self._mat_RHS[i] = self._mat_RHS[i]+self._delta_t*self._F_iminus1j[self._interior_cells[i]]*data.elev[self._interior_cells[i]-grid.ncols]
            
            if grid.is_interior(self._interior_cells[i]+grid.ncols):
                self._operating_matrix[i,i+grid.ncols] = self._operating_matrix[i,i+grid.ncols]-self._delta_t*self._F_iplus1j[self._interior_cells[i]]
            elif bc.boundary_code[self._interior_cells[i]+grid.ncols] == bc.TRACKS_CELL_BOUNDARY:
                j = self._interior_cells.index(bc.tracks_cell[self._interior_cells[i]+grid.ncols])
                self._operating_matrix[i,j] = self._operating_matrix[i,j]-self._delta_t*self._F_iplus1j[self._interior_cells[i]]
            elif bc.boundary_code[self._interior_cells[i]+grid.ncols] == bc.FIXED_GRADIENT_BOUNDARY:
                j = self._interior_cells.index(bc.tracks_cell[self._interior_cells[i]+grid.ncols])
                self._operating_matrix[i,j] = -self._delta_t*self._F_iplus1j[self._interior_cells[i]]
                self._mat_RHS[i] = self._mat_RHS[i]+self._delta_t*self._F_iplus1j[self._interior_cells[i]]*grid.dx*bc.gradient[self._interior_cells[i]+grid.ncols]
            else:
                self._mat_RHS[i] = self._mat_RHS[i]+self._delta_t*self._F_iplus1j[self._interior_cells[i]]*data.elev[self._interior_cells[i]+grid.ncols]
            
            if grid.is_interior(self._interior_cells[i]-grid.ncols-1):
                self._operating_matrix[i,i-grid.ncols-1] = self._operating_matrix[i,i-grid.ncols-1]-self._delta_t*self._F_iminus1jminus1[self._interior_cells[i]]
            elif bc.boundary_code[self._interior_cells[i]-grid.ncols-1] == bc.TRACKS_CELL_BOUNDARY:
                j = self._interior_cells.index(bc.tracks_cell[self._interior_cells[i]-grid.ncols-1])
                self._operating_matrix[i,j] = self._operating_matrix[i,j]-self._delta_t*self._F_iminus1jminus1[self._interior_cells[i]]
            elif bc.boundary_code[self._interior_cells[i]-grid.ncols-1] == bc.FIXED_GRADIENT_BOUNDARY:
                j = self._interior_cells.index(bc.tracks_cell[self._interior_cells[i]-grid.ncols-1])
                self._operating_matrix[i,j] = -self._delta_t*self._F_iminus1jminus1[self._interior_cells[i]]
                self._mat_RHS[i] = self._mat_RHS[i]+self._delta_t*self._F_iminus1jminus1[self._interior_cells[i]]*grid.dx*bc.gradient[self._interior_cells[i]-grid.ncols-1]
            else:
                self._mat_RHS[i] = self._mat_RHS[i]+self._delta_t*self._F_iminus1jminus1[self._interior_cells[i]]*data.elev[self._interior_cells[i]-grid.ncols-1]
            
            if grid.is_interior(self._interior_cells[i]+grid.ncols+1):
                self._operating_matrix[i,i+grid.ncols+1] = self._operating_matrix[i,i+grid.ncols+1]-self._delta_t*self._F_iplus1jplus1[self._interior_cells[i]]
            elif bc.boundary_code[self._interior_cells[i]+grid.ncols+1] == bc.TRACKS_CELL_BOUNDARY:
                j = self._interior_cells.index(bc.tracks_cell[self._interior_cells[i]+grid.ncols+1])
                self._operating_matrix[i,j] = self._operating_matrix[i,j]-self._delta_t*self._F_iplus1jplus1[self._interior_cells[i]]
            elif bc.boundary_code[self._interior_cells[i]+grid.ncols+1] == bc.FIXED_GRADIENT_BOUNDARY:
                j = self._interior_cells.index(bc.tracks_cell[self._interior_cells[i]+grid.ncols+1])
                self._operating_matrix[i,j] = -self._delta_t*self._F_iplus1jplus1[self._interior_cells[i]]
                self._mat_RHS[i] = self._mat_RHS[i]+self._delta_t*self._F_iplus1jplus1[self._interior_cells[i]]*grid.dx*bc.gradient[self._interior_cells[i]+grid.ncols+1]
            else:
                self._mat_RHS[i] = self._mat_RHS[i]+self._delta_t*self._F_iplus1jplus1[self._interior_cells[i]]*data.elev[self._interior_cells[i]+grid.ncols+1]
            
            if grid.is_interior(self._interior_cells[i]+grid.ncols-1):
                self._operating_matrix[i,i+grid.ncols-1] = self._operating_matrix[i,i+grid.ncols-1]-self._delta_t*self._F_iplus1jminus1[self._interior_cells[i]]
            elif bc.boundary_code[self._interior_cells[i]+grid.ncols-1] == bc.TRACKS_CELL_BOUNDARY:
                j = self._interior_cells.index(bc.tracks_cell[self._interior_cells[i]+grid.ncols-1])
                self._operating_matrix[i,j] = self._operating_matrix[i,j]-self._delta_t*self._F_iplus1jminus1[self._interior_cells[i]]
            elif bc.boundary_code[self._interior_cells[i]+grid.ncols-1] == bc.FIXED_GRADIENT_BOUNDARY:
                j = self._interior_cells.index(bc.tracks_cell[self._interior_cells[i]+grid.ncols-1])
                self._operating_matrix[i,j] = -self._delta_t*self._F_iplus1jminus1[self._interior_cells[i]]
                self._mat_RHS[i] = self._mat_RHS[i]+self._delta_t*self._F_iplus1jminus1[self._interior_cells[i]]*grid.dx*bc.gradient[self._interior_cells[i]+grid.ncols-1]
            else:
                self._mat_RHS[i] = self._mat_RHS[i]+self._delta_t*self._F_iplus1jminus1[self._interior_cells[i]]*data.elev[self._interior_cells[i]+grid.ncols-1]
            
            if grid.is_interior(self._interior_cells[i]-grid.ncols+1):
                self._operating_matrix[i,i-grid.ncols+1] = self._operating_matrix[i,i-grid.ncols+1]-self._delta_t*self._F_iminus1jplus1[self._interior_cells[i]]
            elif bc.boundary_code[self._interior_cells[i]-grid.ncols+1] == bc.TRACKS_CELL_BOUNDARY:
                j = self._interior_cells.index(bc.tracks_cell[self._interior_cells[i]-grid.ncols+1])
                self._operating_matrix[i,j] = self._operating_matrix[i,j]-self._delta_t*self._F_iminus1jplus1[self._interior_cells[i]]
            elif bc.boundary_code[self._interior_cells[i]-grid.ncols+1] == bc.FIXED_GRADIENT_BOUNDARY:
                j = self._interior_cells.index(bc.tracks_cell[self._interior_cells[i]-grid.ncols+1])
                self._operating_matrix[i,j] = -self._delta_t*self._F_iminus1jplus1[self._interior_cells[i]]
                self._mat_RHS[i] = self._mat_RHS[i]+self._delta_t*self._F_iminus1jplus1[self._interior_cells[i]]*grid.dx*bc.gradient[self._interior_cells[i]-grid.ncols+1]
            else:
                self._mat_RHS[i] = self._mat_RHS[i]+self._delta_t*self._F_iminus1jplus1[self._interior_cells[i]]*data.elev[self._interior_cells[i]-grid.ncols+1]
            
            
    def update(self, grid, data):
        #Initialize the variables for the step:
        self.initialize(grid, data)
        #Solve interior of grid:
        self._interior_elevs = self._operating_matrix.I * self._mat_RHS
        
        #Load the new elevs into the actual elev data
        j=0 #counter for location in the matrix _interior_elevs
        for i in range(grid.ncells):
            if grid.is_interior(i):
                data.elev[i] = self._interior_elevs[j]
                j = j+1
            elif bc.boundary_code[i] == bc.TRACKS_CELL_BOUNDARY:
                data.elev[i] = self._interior_elevs[self._interior_cells.index(bc.tracks_cell[i])]
            elif bc.boundary_code[i] == bc.FIXED_GRADIENT_BOUNDARY:
                data.elev[i] = self._interior_elevs[self._interior_cells.index(bc.tracks_cell[i])] + bc.gradient[i]*grid.dx
            else:
                continue #The old value of elev is still correct

def main():
    dt = 1.
    nt = 1
    diff = perron_nl_diff()
    cr, mg, vectors = craters.dig_one_crater_smaller()
    for i in range(nt):
        diff.update(mg, vectors)

if __name__ == '__main__':
    main() 
