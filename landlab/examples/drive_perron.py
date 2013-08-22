import numpy
from landlab import RasterModelGrid
from landlab.components.perron_nl_diff import perron_nl_diff_faster as perron_nl_diff
from landlab.components.craters import data
import pylab

class drive_perron(object):
    
    def __init__(self):
        self.nr = 100
        self.nc = 100
        self.dx = 1.
        self.dt = 1.
        self.nt = 10
        self.uplift = 1. #per tstep
        self.build_grid_and_data()
        self.diff = perron_nl_diff(self.mg, self.vectors, self.dt)
       
         
    def build_grid_and_data(self):
        '''
        Builds a nr*nc grid open at the top and bottom edges for testing
        '''
        self.mg = RasterModelGrid()
        self.mg.initialize(self.nr, self.nc, self.dx)
        self.mg.set_noflux_boundaries(False, True, False, True)  # Boundary conditions
        self.vectors = data()
        self.vectors.elev = [0.] * self.mg.ncells
    
    def update(self):
        for i in range(self.mg.ncells):
            self.vectors.elev[i] = self.vectors.elev[i] + self.uplift


def main():
    i = 0
    dp = drive_perron()
    while i<dp.nt:
        dp.update()
        dp.diff.update(dp.mg, dp.vectors)
        i = i+1
        print 'Loop ', i, ' complete'
        pylab.plot(range(dp.nr), dp.vectors.elev[1:(dp.mg.ncells+2-dp.nc):dp.nc])
    pylab.show()
    return dp
    

if __name__=='__main__':
    main()
