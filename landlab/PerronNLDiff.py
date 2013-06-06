import numpy

class PerronNLDiff(object):
    '''
    This module uses Taylor Perron's implicit (2011) method to solve the nonlinear hillslope diffusion equation across a rectangular grid for a single timestep. Note it works with the mass flux implicitly, and thus does not actually calculate it. At the moment, this can't handle boundary cells (let's run on a looped grid!)
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
#        self._F_ijless1 = []
#        self._F_ijplus1 = []
#        self._F_iless1j = []
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
        self._operating_matrix = numpy.zeros((grid.ncells, grid.ncells), dtype=float)
        self._delta_x = grid.dx
        self._delta_y = grid.dx
        one_over_delta_x = 1./self._delta_x
        one_over_delta_y = 1./self._delta_y
        one_over_delta_x_sqd = one_over_delta_x**2.
        one_over_delta_y_sqd = one_over_delta_y**2.
        self._delta_t = 100.
        self._b = 1./self._S_crit**2.
        self._z_x = [0] * grid.ncells #ncells, or nactivecells?
        self._z_y = [0] * grid.ncells
        self._z_xx = [0] * grid.ncells
        self._z_yy = [0] * grid.ncells
        self._z_xy = [0] * grid.ncells
        self._d = [0] * grid.ncells
        self._F_ij = [0] * grid.ncells
        self._F_ijless1 = [0] * grid.ncells
        self._F_ijplus1 = [0] * grid.ncells
        self._F_iless1j = [0] * grid.ncells
        self._F_iplus1j = [0] * grid.ncells
        self._F_iplus1jplus1 = [0] * grid.ncells
        self._F_iminus1jminus1 = [0] * grid.ncells
        self._F_iplus1jminus1 = [0] * grid.ncells
        self._F_iminus1jplus1 = [0] * grid.ncells
        self._func_on_z = [0] * grid.ncells
        self._equ23_RHS = [0] * grid.ncells 
        
        for i in grid.get_interior_cells(): #What do we do about boundary cells?? Need them to be counted to give enough simultaneous equations.
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
            self._F_ijless1[i] = self._kappa*self._d[i]*one_over_delta_x_sqd - _abd_sqd*self._z_x[i]*(self._z_xx[i]+self._z_yy[i])*one_over_delta_x - 4.*_abd_sqd*self._b*self._d[i]*(self._z_x[i]**2.*self._z_xx[i]+self._z_y[i]**2.*self._z_yy[i]+2.*self._z_x[i]*self._z_y[i]*self._z_xy[i])*self._z_x[i]*one_over_delta_x - 2.*_abd_sqd*(self._z_x[i]*self._z_xx[i]*one_over_delta_x-self._z_x[i]**2.*one_over_delta_x_sqd+self._z_y[i]*self._z_xy[i]*one_over_delta_x)
            self._F_ijplus1[i] = self._kappa*self._d[i]*one_over_delta_x_sqd + _abd_sqd*self._z_x[i]*(self._z_xx[i]+self._z_yy[i])*one_over_delta_x + 4.*_abd_sqd*self._b*self._d[i]*(self._z_x[i]**2.*self._z_xx[i]+self._z_y[i]**2.*self._z_yy[i]+2.*self._z_x[i]*self._z_y[i]*self._z_xy[i])*self._z_x[i]*one_over_delta_x + 2.*_abd_sqd*(self._z_x[i]*self._z_xx[i]*one_over_delta_x+self._z_x[i]**2.*one_over_delta_x_sqd+self._z_y[i]*self._z_xy[i]*one_over_delta_x)
            self._F_iless1j[i] = self._kappa*self._d[i]*one_over_delta_y_sqd - _abd_sqd*self._z_y[i]*(self._z_xx[i]+self._z_yy[i])*one_over_delta_y - 4.*_abd_sqd*self._b*self._d[i]*(self._z_x[i]**2.*self._z_xx[i]+self._z_y[i]**2.*self._z_yy[i]+2.*self._z_x[i]*self._z_y[i]*self._z_xy[i])*self._z_y[i]*one_over_delta_y - 2.*_abd_sqd*(self._z_y[i]*self._z_yy[i]*one_over_delta_y-self._z_y[i]**2.*one_over_delta_y_sqd+self._z_x[i]*self._z_xy[i]*one_over_delta_y)
            self._F_iplus1j[i] = self._kappa*self._d[i]*one_over_delta_y_sqd + _abd_sqd*self._z_y[i]*(self._z_xx[i]+self._z_yy[i])*one_over_delta_y + 4.*_abd_sqd*self._b*self._d[i]*(self._z_x[i]**2.*self._z_xx[i]+self._z_y[i]**2.*self._z_yy[i]+2.*self._z_x[i]*self._z_y[i]*self._z_xy[i])*self._z_y[i]*one_over_delta_y + 2.*_abd_sqd*(self._z_y[i]*self._z_yy[i]*one_over_delta_y+self._z_y[i]**2.*one_over_delta_y_sqd+self._z_x[i]*self._z_xy[i]*one_over_delta_y)
            self._F_iplus1jplus1[i] = _abd_sqd*self._z_x[i]*self._z_y[i]*one_over_delta_x*one_over_delta_y
            self._F_iminus1jminus1[i] = self._F_iplus1jplus1
            self._F_iplus1jminus1[i] = -self._F_iplus1jplus1
            self._F_iminus1jplus1[i] = self._F_iplus1jminus1

            #RHS of equ 6 (see para [20])
            self._func_on_z[i] = self._rock_density/self._sed_density*self._uplift + self._kappa*((self._z_xx[i]+self._z_yy[i])/(1.-(self._z_x[i]**2.+self._z_y[i]**2.)/self._S_crit**2.) + 2.*(self._z_x[i]**2.*self._z_xx[i]+self._z_y[i]**2.*self._z_yy[i]+2.*self._z_x[i]*self._z_y[i]*self._z_xy[i])/(self._S_crit**2.*(1.-(self._z_x[i]**2.+self._z_y[i]**2.)/self._S_crit**2.)**2.))

            self._equ23_RHS[i] = data.elev[i] + self._delta_t*(self._func_on_z[i] - (self._F_ij*data.elev[i]+self._F_ijless1[i]*data.elev[cell_neighbors[2]]+self._F_ijplus1[i]*data.elev[cell_neighbors[0]]+self._F_iless1j[i]*data.elev[cell_neighbors[3]]+self._F_iplus1j[i]*data.elev[cell_neighbors[1]]+self._F_iminus1jminus1[i]*data.elev[cell_diagonals[2]]+self._F_iplus1jplus1[i]*data.elev[cell_diagonals[0]]+self._F_iplus1jminus1[i]*data.elev[cell_diagonals[1]]+self._F_iminus1jplus1[i]*data.elev[cell_diagonals[3]]))

        self._mat_RHS = numpy.matrix(self._equ23_RHS)

        #assemble the operation matrix
        #...more problems with the boundaries -> make exceptions
        for i in range(grid.ncells): #this is for each cell in the grid
            self._operating_matrix[i,i] = 1.-self._delta_t*self._F_ij[i]
            self._operating_matrix[i,i-1] = -self._delta_t*self._F_ijless1[i]
            self._operating_matrix[i,i+1] = -self._delta_t*self._F_ijplus1[i]
            self._operating_matrix[i,i-grid.ncols] = -self._delta_t*self._F_iless1j[i]
            self._operating_matrix[i,i+grid.ncols] = -self._delta_t*self._F_iplus1j[i]
            self._operating_matrix[i,i-grid.ncols-1] = -self._delta_t*self._F_iminus1jminus1[i]
            self._operating_matrix[i,i+grid.ncols+1] = -self._delta_t*self._F_iplus1jplus1[i]
            self._operating_matrix[i,i+grid.ncols-1] = -self._delta_t*self._F_iplus1jminus1[i]
            self._operating_matrix[i,i-grid.ncols+1] = -self._delta_t*self._F_iminus1jplus1[i]
            
    def update(self):
        data.elev = self._operating_matrix.I * self._mat_RHS

