import landlab.components.craters as craters
craters = reload(craters)
from landlab import RasterModelGrid
import numpy
import time
from pylab import show, imshow, colorbar, plot

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

def dig_one_crater_timetest(reps, nr, nc, dx, rel_x, rel_y, radius, switch):
    '''
        This script digs a single crater, either using the crawler method (switch=0), or by sweeping the whole grid (switch=1). It is designed to time the the crater digging scripts in crater.py, with an eye to determining 1. if solving over the whole grid can be faster than the conditional crawler, and 2. if there's a crossover point for small (i.e., typical!) craters.
    '''
    mg = RasterModelGrid()
    mg.initialize(nr, nc, dx)
    vectors = data(mg)
    vectors.elev[:] = 100.
    cr = craters.impactor()

    start_time = time.time()
    for i in range(reps):
        cr._radius = radius
        #print 'Radius: ', cr._radius
        cr.set_depth_from_size()
        #print 'Depth: ', cr._depth
        cr.set_crater_volume()
        cr._xcoord = rel_x*mg.get_grid_xdimension()
        cr._ycoord = rel_y*mg.get_grid_ydimension()
        vertices_array = mg.get_nodes_around_point(cr._xcoord, cr._ycoord)
        distances_to_vertices = []
        #print vertices_array
        #print len(mg.node_x)
        for x in vertices_array:
            distances_to_vertices.append(numpy.sqrt((cr._xcoord-mg.node_x[x])**2. + (cr._ycoord-mg.node_y[x])**2.))
        cr.closest_node_index = vertices_array[numpy.argmin(distances_to_vertices)]
        cr.closest_node_elev = vectors.elev[cr.closest_node_index]
        
        cr._angle_to_horizontal = numpy.pi*0.5*3./3.
        cr._azimuth_of_travel = numpy.pi*1.5
        cr.set_crater_mean_slope_v2(mg, vectors)
        print 'Azimuth of travel: ', cr._azimuth_of_travel
        print 'Angle of ground: ', cr._surface_slope
        print 'Dip direction of ground: ', cr._surface_dip_direction
        if switch:
            cr.set_elev_change_across_grid(mg, vectors)
        else:
            cr.set_elev_change_at_pts(mg, vectors)
        print 'Impact angle to ground normal: ', cr.impactor_angle_to_surface_normal
    
    end_time = time.time()
    print('Elapsed time was %g seconds' % (end_time - start_time))
    
    elev_raster = mg.node_vector_to_raster(vectors.elev, flip_vertically=True)
    vectors.viewing_raster = elev_raster
    return cr, mg, vectors


crater, mg, vectors = dig_one_crater_timetest(1,1200,1200,0.0025, 0.5, 0.5, .75, 1)
imshow(vectors.viewing_raster)
show()

start_time = time.time()
mg, vectors_eroded, profile, xsec = craters.dig_some_craters(mg, vectors, 20000)
end_time = time.time()
print('Elapsed time was %g seconds' % ((end_time - start_time)))

imshow(vectors_eroded.viewing_raster)
show()

