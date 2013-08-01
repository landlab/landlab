import landlab.components.craters as craters
from landlab.model_grid import RasterModelGrid
import numpy
from pylab import show, imshow, colorbar, plot

craters=reload(craters)

#craters_data = craters.one_crater_then_degrade()

class data(object):
    '''
        This is where all the whole-grid data lives, as arrays over the various elements of the grid.
        '''
    #Data goes here!!!
    def __init__(self, grid):
        self.elev = grid.create_node_dvector() #some data
        self.flag_already_in_the_list = grid.create_node_dvector()
        self.craters_over_max_radius_not_plotted = grid.create_node_dvector()
        self.impact_sequence = []

#def dig_one_crater_driver(nr, nc, dx, rel_x, rel_y, radius=numpy.nan, angle=numpy.nan):
#    '''
#        This is an ad-hoc script to dig one crater.
#        '''
#    
#    #Setup the grid & impactor object
#    mg = RasterModelGrid()
#    mg.initialize(nr, nc, dx)
#    vectors = data(mg)
#    vectors.elev[:] = 0.
#    cr = craters.impactor()
#
#    #Dig the crater
#    cr.excavate_a_crater_optimized(mg, vectors, forced_radius=radius, forced_angle=angle, forced_pos=[rel_x, rel_y])
#
#    #Finalize
#    elev_raster = mg.node_vector_to_raster(vectors.elev, flip_vertically=True)
#    return elev_raster


#elev_raster = dig_one_crater_driver(120, 120, 0.025, 0.5, 0.5, 0.75, 90.)
#imshow(elev_raster)
#colorbar()
#show()

#This code snippet to produce a landscape w. 500000 impacts, but where the first half of the impacts are 5x bigger on average.
mg = RasterModelGrid()
mg.initialize(1200, 1200, 0.0025)
vectors = data(mg)
vectors.elev[:] = 0.

dict_of_vectors = craters.five_times_reduction(mg, vectors)

