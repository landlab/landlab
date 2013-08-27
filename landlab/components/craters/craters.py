#! /usr/env/python

from random import random
import math
import numpy
from collections import deque
import sys
import time
import scipy.optimize as opt
from landlab import RasterModelGrid #this is the tMesh equivalent module

#these ones only so we can run this module ad-hoc:
from pylab import plot, draw, show, contour, imshow, colorbar
from copy import copy

sys.setrecursionlimit(1500)

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

class impactor(object):
    '''
    This class holds all parameters decribing properties of a single impact structure, and contains methods for recalculating fresh and internally consistent data describing such a impact structure.
    Built DEJH Spring 2013.
    '''
    def __init__(self, min_radius=0.005):
        self._xcoord = -999.
        self._ycoord = -999.
        self._str_regime_cutoff = 0.3 #??? #Holsapple gives 0.1km as the cutoff for R/d scaling, Moore in Housen gives this as the likely Lunar cutoff for true str regime (p. 2495 Housen)
        self._simple_radius_depth_ratio_Pike = 2.55 #Pike thru Holsapple, tho Garvin (2011) gives 2.0 as "typically cited"
        #the ratio for complex craters is called as a method, and not set as such
        self._simple_cr_radius_cutoff = 7.5 #km following Melosh p. 224 - LUNAR
        self._complex_cr_radius_cutoff = 10. #km, folling Melosh, i.e., descending
    #anything in between these values needs a chance of being either - Pike says transition range is R = 6-12.5 km
        self._max_complex_cr_radius = 70. #km, bigger is a "peak ring crater". This is not currently used.
        #self._power_for_complex_crater = 2. #i.e., d = D(r/R)^power. =2 is a parabola. This is no longer used; replaced by get_crater_shape_exp() for all craters simple and complex.
        #we follow Garvin et al 2011 and assume power=1.45 for simple crater
        #NB - if we assume a power law and dictate S=32deg at the rim, we recover n=0.51R/D - gives n~1.3 for Linne (synthetic; if real depth it's ~1) and maxes out ~2 for biggest craters. (Plus infill by lava?) ->This is a nice way to do it. This parameter is presently unused.
        self._minimum_crater = min_radius #km. This is the smallest modelled crater radius. 10m diameter is the strict cutoff for known cr distns
        self.ivanov_a = [-3.0876, -3.557528, 0.781027, 1.021521, -0.156012, -0.444058, 0.019977, 0.08685, -0.005874, -0.006809, 0.000825, 0.0000554] #The coefficients for Ivanov's crater distn function
        self._impactor_angle_to_surface = -999.
        self._angle_to_horizontal = -999.
        self._minimum_ejecta_thickness = 0.000001
        #NB - If this min thickness is changed, the optimization point in excavate_a_crater_optimized() will also need to be changed
        self._beta_factor = 0.5 #this is the arbitrary term that controls how "stretched out" the ejecta field is. <+0.5 prevents "outside the ejecta field" regions forming
        
        self.total_counted_craters = self.ivanov_prod_equ_as_Nequals(self._minimum_crater*2.)
        #Define the Ivanov fn and its derivatives, for use in Newton Raphson optimization. Note this is normalized now.
        self.ivanov_prod_fn = lambda x, N_as_fraction: self.ivanov_prod_equ(x,(N_as_fraction/self.total_counted_craters))
        self.ivanov_prod_fn_1stderiv = lambda x, N_as_fraction: self.ivanov_prod_equ_1stderiv(x, N_as_fraction)
        #self.ivanov_prod_fn_2ndderiv = lambda x: self.ivanov_prod_equ_2ndderiv(x)

    
    def get_complex_radius_depth_ratio(self):
        '''
        This method returns the ratio of radius to max depth for a complex crater, which is itself a function of that crater's radius. This equation is from Pike '77, as quoted in Holsapple '93, p. 358.
        '''
        return self._radius**0.7 #Pike 77, following Holsapple 93, p. 358 - Not clear how I need to call this again to get it to update before use
    
    
    def get_crater_shape_exp(self):
        '''
        This method assumes the max depth and radius of a crater are known. It provides n for a power law of form d = D*(r/R)**n, where D and R are the known values, by assuming the outer edges of the crater sit at angle of repose. This gives very sensible answers; n~2 for big, complex craters (Garvin et al, 2000, p.333: "There is a strong tendency for craters to become more paraboloidal with increasing diameter, independent of location.") and n~1.3 for ~2km simple craters (Garvin following Croft has ~1.18).
        '''
        return 0.51 * self._radius / self._depth
    
    
    def get_complex_peak_radius(self):
        '''
        This method returns the radial distance to the edge of a complex crater central mound, based on the crater radius. This follows an equation presented in Melosh's Planetary Surface Processes.
        '''
        return 0.22 * self._radius #Melosh
    
    
    def get_complex_peak_str_uplift(self):
        '''
        This method returns the radial maximum structural uplift, interpreted as a height, of a complex crater central mound, based on the crater radius. This follows an equation presented in Melosh's Planetary Surface Processes, adjusted for R, not D.
        '''
        return 0.13 * self._radius ** 1.1 #melosh, note correction for R not D
    
    
    #Holsapple & Housen et al 1983 note that in the strength regime, ejecta distribution will NOT scale independently of crater size (it does in the gravity regime) - ejecta will be proportionally FURTHER FROM the rim as craters get smaller (e.g., Housen et al 83, fig 9)
    #...this, of course, makes small craters disproportionately effective diffusers - must find a way of quantifying this!
    #Moore et al 74 cited in Housen et al 83 gives Rce/R as const ~2.35 +0.56 -0.45 std err for Lunar craters 0.3-100km, where Rce is the continuous ejecta limit. Below this, seems to level off abruptly - i.e., Rce ~ R**(<1) (and const will change also to maintain connection at 0.3km)
    #A key point - geomorphic gardening is not totally comparable to regolith gardening. We care about topology - the final crater. Regolith fracture and turnover is more likely to depend on the crater transient depth.
    
    
    def ivanov_prod_equ(self, D, N):
        '''
        This function expresses the crater production function after Ivanov, 1999, in the form 0 = f(D, N). N is the total number density of craters greater than x diameter. Remember, N is in km^-2 Ga^-1. If N = N_counted / N_max, then this is the CPF which will yield the PDF when solved for D given N.
        '''
        sum_terms = numpy.empty(12, dtype=float)
        sum_terms[0] = self.ivanov_a[0] - math.log10(N)
        for i in range(1,12):
            sum_terms[i] = self.ivanov_a[i] * math.log10(D) ** i
        return numpy.sum(sum_terms)
    
    
    #A second version, in the conventional form N = f(D):
    def ivanov_prod_equ_as_Nequals(self, D):
        '''
        This function expresses Ivanov's CPF in the normal form, N = f(D). It returns N (not log10(N)!).
        '''
        sum_terms = numpy.empty(12, dtype=float)
        for i in range (0,12):
            sum_terms[i] = self.ivanov_a[i] * math.log10(D) ** i
        return pow( 10., numpy.sum(sum_terms) )

    
    def ivanov_prod_equ_1stderiv(self, D, N):
        '''
        This is a helper function for the Newton-Raphson optimized solver in solve_ivanov_for_crater_diam. It returns the first derivative for the production function. Note that the N term is a dummy variable used for consistency in the opt.newton(), and does nothing here.
        '''
        sum_terms = numpy.empty(11, dtype=float)
        ln_of_10 = math.log(10)
        for i in range(1,12):
            sum_terms[i-1] = self.ivanov_a[i] * i * math.log(D) ** (i-1.) / ln_of_10**i
        return (numpy.sum(sum_terms) / D)


    def ivanov_prod_equ_2ndderiv(self, D):
        '''
        This is a helper function for the Newton-Raphson optimized solver in solve_ivanov_for_crater_diam. It returns the second derivative for the production function.
        '''
        ln_of_10 = math.log(10)
        sum_terms = numpy.empty(11, dtype=float)
        sum_terms[0] = self.ivanov_a[1] / ln_of_10
        for i in range(2,12):
            sum_terms[i-1] = i * self.ivanov_a[i] * (math.log(D) - i + 1.) * math.log(D)**(i-2.) / ln_of_10**i
        return (-numpy.sum(sum_terms) / (D*D))


    def solve_ivanov_for_crater_diam(self, N_by_fract):
        '''
        Uses a Newton Raphson method to return a random diameter drawn from an Ivanov crater distribution, when provided with a uniformly distributed random number 0->1.
        '''
        args_in = N_by_fract,
        #Estimate zero position by a super ad-hoc Shoemaker distn. Note K is not really right here - we've defined it at D=1; I think it should be equal to Ivanov at D~0.03.
        x_0 = (10.**self.ivanov_a[0] / (self.total_counted_craters * N_by_fract)) ** 0.345
        return opt.newton(func=self.ivanov_prod_fn, x0=x_0, fprime=self.ivanov_prod_fn_1stderiv, args=args_in, tol=0.001)
    
    
    def set_cr_radius_from_shoemaker(self, data):
        '''
        This method is a less accurate (but faster) alternative to the Newton-Raphson-on-Ivanov-distn also available in this object. It takes a random number between 0 and 1, and returns a crater radius based on a py distn N = kD^-2.9, following Shoemaker et al., 1970.
        '''
        self._radius = self._minimum_crater*(random())**-0.345
        if self._radius > self._max_complex_cr_radius:
            data.craters_over_max_radius_not_plotted.append(self._radius)
            print 'Drew a crater above the maximum permitted size. Drawing a new crater...'
            self.set_size(data)


    def set_coords(self, grid, data):
        '''
        This method selects a random location inside the grid onto which to map an impact. It also sets variables for the closest grid node to the impact, and the elevation at that node.
        '''
        #NB - we should be allowing craters OUTSIDE the grid - as long as part of them impinges.
        #This would be relatively easy to implement - allow allocation out to the max crater we expect, then allow runs using these coords on our smaller grid. Can save comp time by checking if there will be impingement before doing the search.
        self._xcoord = random() * (grid.get_grid_xdimension() - grid.dx)
        self._ycoord = random() * (grid.get_grid_ydimension() - grid.dx)
        #Snap impact to grid:
        vertices_array = grid.get_nodes_around_point(self._xcoord, self._ycoord)
        distances_to_vertices = numpy.empty(vertices_array.size)
        counter = 0
        for x in vertices_array:
            #print 'The ID is:', x
            distances_to_vertices[counter] = numpy.sqrt(self._xcoord-grid.node_x[x]**2. + self._ycoord-grid.node_y[x]**2.)
            counter += 1
        self.closest_node_index = vertices_array[numpy.argmin(distances_to_vertices)]
        self.closest_node_elev = data.elev[self.closest_node_index]
        #print 'Closest node elev and index: ', self.closest_node_elev, self.closest_node_index
        # Check we haven't hit a boundary node. Reset if we have:
        #NB - this step can introduce recursion very easily into the program. We need a way out of this that doesn't result in endless recursion. 
        #if not grid.is_interior(self.closest_node_index):
            #self.set_coords(grid, data)
        #NB - snapping to the grid may be quite computationally demanding in a Voronoi. Could dramatically increase speed by picking a node for centre and working out from there (or just leaving it on the node...). Fine for square grid though.
            

    def set_impactor_angles(self):
        '''
        This method sets the angle of impact, assuming the only effect is rotation of the planet under the impactor bombardment (i.e., if the target looks like a circle to the oncoming impactor, there's more limb area there to hit). As long as target is rotating relative to the sun, other (directional) effects should cancel. Angle is given to horizontal. Also sets a random azimuth.
        '''
        self._angle_to_horizontal =  numpy.arccos(random()) #gives sin distn with most values drawn nearer to 0
        self._azimuth_of_travel = random() * 2. * numpy.pi #equal chance of any azimuth
        #Shoemaker 1983 gives a speculative equn for the suppression of D by increasing impact angle (his equ 3), but effect is minor, and it's probably not worth the trouble.
        #Shoemaker 1962 (in book) apparently states low angle impacts are very rare. Need to read this.

    
    def set_size(self, data):
        '''
        This method draws a crater radius at random from the PDF describing crater sizes dictated by the Ivanov distribution, with the probability of obtaining that crater size depending on the relative abundance of that size. Note this method works with the Ivanov distribution according to diameter, but returns a radius.
        '''
        self._radius = 0.5 * self.solve_ivanov_for_crater_diam(random())
        if self._radius > self._max_complex_cr_radius:
            data.craters_over_max_radius_not_plotted.append(self._radius)
            print 'Drew a crater above the maximum permitted size. Drawing a new crater...'
            self.set_size(data)

    
    def set_depth_from_size(self):
        '''
        This method sets a known crater with diameter as either simple or complex, then sets its maximum depth at the center (inferred for complex craters).
        '''
        #Needs to account for simple to complex.
        #Let's ignore the transition at the lowest size scales for now.
        #Need to consider the transition at the top end too!!
        if self._radius < self._simple_cr_radius_cutoff:
            self._crater_type = 0 #simple
        elif self._radius > self._complex_cr_radius_cutoff:
            self._crater_type = 1 #complex
        else:
            py_complex_at_radius = (self._radius - self._simple_cr_radius_cutoff) / (self._complex_cr_radius_cutoff - self._simple_cr_radius_cutoff)
            if random() > py_complex_at_radius:
                self._crater_type = 0 #simple
            else:
                self._crater_type = 1 #complex
                
        if not self._crater_type:
            self._depth = self._radius / self._simple_radius_depth_ratio_Pike
        else:
            self._depth = self._radius / self.get_complex_radius_depth_ratio()
            #Derive the scaling metrics for a complex crater, for subsequent use:
            self._complex_peak_radius = self.get_complex_peak_radius()
            self._complex_peak_str_uplift = self.get_complex_peak_str_uplift()

                
    def set_crater_volume(self):
        '''
        This method uses known crater depth and radius and sets the volume of the excavated cavity.
        '''
        self._cavity_volume = 0.51 * numpy.pi * self._depth * self._radius**3. / (0.51*self._radius + 2.*self._depth)


    def set_crater_mean_slope(self, grid, data):
        '''
        This method takes a crater of known radius, and which has already been "snapped" to the grid through snap_impact_to_grid(mygrid), and returns a spatially averaged value for the local slope of the preexisting topo beneath the cavity footprint. For computational efficiency reasons, the average slope is a mean of local slopes at 8 points 45 degrees apart around the crater rim, plus the center point. This function also  sets the mean surface dip direction.
            NOW DEPRECIATED, USE V2!
        '''
        #This might not be the best way to do this... There are issues with the magnitude of the slope being high if the individual points slope steeply even if they point in different directions.
        half_crater_radius = 0.707 * self._radius
        radial_points = []
        #elevation_array = []
        slope_array = []
        dip_direction_array = []
        #Work round clockwise from 12 o'clock
        slope_pts = [[self._xcoord,self._ycoord+self._radius],
                     [self._xcoord+half_crater_radius,self._ycoord+half_crater_radius],
                     [self._xcoord+self._radius,self._ycoord],
                     [self._xcoord+half_crater_radius,self._ycoord-half_crater_radius],
                     [self._xcoord,self._ycoord-self._radius],
                     [self._xcoord-half_crater_radius,self._ycoord-half_crater_radius],
                     [self._xcoord-self._radius,self._ycoord],
                     [self._xcoord-half_crater_radius,self._ycoord+half_crater_radius]]
        for x,y in slope_pts:
            if x < grid.get_grid_xdimension() and y < grid.get_grid_ydimension():
                radial_points.append(self.snap_coords_to_grid(grid, x, y))
        for a in radial_points:
            if grid.is_interior(a):
                slope, dir = grid.calculate_max_gradient_across_node(data.elev, a)
            #Make sure we stay on the grid!!:
                if not numpy.isnan(slope): #Probably not necessary now; retain for belt and braces
                    slope_array.append(slope)
                    dip_direction_array.append(dir)
                else:
                    print 'Slope is not defined at one of the nodes round the crater rim.'
            else:
                print 'One of the nodes round the crater rim is a boundary cell; no slope calculated here.'
        # Add in the central node:
        slope0, dir0 = grid.calculate_max_gradient_across_node(data.elev, self.closest_node_index)
        if not numpy.isnan(slope0):
            slope_array.append(slope0)
            dip_direction_array.append(dir0)
        self._surface_slope = numpy.mean(slope_array)
        dip_direction_x = []
        dip_direction_y = []
        for a in dip_direction_array:
            dip_direction_x.append(numpy.sin(a))
            dip_direction_y.append(numpy.cos(a))
        dip_direction_x_mean = numpy.mean(dip_direction_x)
        dip_direction_y_mean = numpy.mean(dip_direction_y)
        angle_to_yaxis = numpy.arctan(dip_direction_y_mean/dip_direction_x_mean)
        if dip_direction_x < 0.:
            self._surface_dip_direction = 1.5 * numpy.pi - angle_to_yaxis
        else:
            self._surface_dip_direction = 0.5 * numpy.pi - angle_to_yaxis
        print 'The mean slope is: ', self._surface_slope
    
    
    def set_crater_mean_slope_v2(self, grid, data):
        '''
        This method takes a crater of known radius, and which has already been "snapped" to the grid through snap_impact_to_grid(mygrid), and returns a spatially averaged value for the local slope of the preexisting topo beneath the cavity footprint. This version of the method works by taking four transects across the crater area every 45 degrees around its rim, calculating the slope along each, then setting the slope as the greatest, positive downwards and in the appropriate D8 direction. This function also sets the mean surface dip direction.
        '''
        #It would be trivial to add more divisions, i.e., D16, D32, D64, etc. Code is setup for this; just add to the slope_pts lists, and adjust the set length of the arrays with divisions
        divisions = 4 #len(radial_points1), i.e., half the no. of points == the number of lines across the wheel
        half_crater_radius = 0.707 * self._radius
        radial_points1 = numpy.empty(divisions, dtype=int)
        radial_points2 = numpy.empty(divisions, dtype=int)
        #elevation_array = []
        slope_array = []
        #dip_direction_array = []
        #Work round clockwise from 12 o'clock
        slope_pts1 = [[self._xcoord,self._ycoord+self._radius],
                     [self._xcoord+half_crater_radius,self._ycoord+half_crater_radius],
                     [self._xcoord+self._radius,self._ycoord],
                     [self._xcoord+half_crater_radius,self._ycoord-half_crater_radius]]
        slope_pts2 = [[self._xcoord,self._ycoord-self._radius],
                     [self._xcoord-half_crater_radius,self._ycoord-half_crater_radius],
                     [self._xcoord-self._radius,self._ycoord],
                     [self._xcoord-half_crater_radius,self._ycoord+half_crater_radius]]
        counter=0
        grid_xdimension = grid.get_grid_xdimension()
        grid_ydimension = grid.get_grid_ydimension()
        #this loop tests if the points round the rim are within the grid, then makes a list of nodes closest to that point. If not, it adds the central point to the list.
        for x,y in slope_pts1:
            if 0. < x < grid_xdimension and 0. < y < grid_ydimension:
                radial_points1[counter] = self.snap_coords_to_grid(grid, x, y)
            else:
                radial_points1[counter] = self.closest_node_index
            counter += 1
        counter=0
        for x,y in slope_pts2:
            if 0. < x < grid_xdimension and 0. < y < grid_ydimension:
                radial_points2[counter] = self.snap_coords_to_grid(grid, x, y)
            else:
                radial_points2[counter] = self.closest_node_index
            counter += 1
        for a in range(divisions):
            du = data.elev[radial_points1[a]] - data.elev[radial_points2[a]]
            #If du is negative, it means the surface slopes broadly EAST
            if radial_points1[a] is radial_points2[a]:
                dx = 0.
            elif radial_points1[a] is self.closest_node_index or radial_points2[a] is self.closest_node_index:
                dx = self._radius
            else:
                dx = 2.*self._radius
            if dx:
                slope_array.append(numpy.arctan(du/dx))
            else:
                slope_array.append(numpy.nan)
        slope_array = numpy.array(slope_array)
        #add test to see if crashing is due to problems in slope extraction:
        if numpy.isnan(slope_array).all():
            print 'All slopes beneath area did not evaluate! ABORT ABORT ABORT!!!'
        hi_mag_slope_index = numpy.argmax(numpy.fabs(slope_array))
        hi_mag_slope = slope_array[hi_mag_slope_index]
        #print hi_mag_slope
        if hi_mag_slope > 0.: #i.e., dips WEST (or dead S)
            self._surface_dip_direction = (hi_mag_slope_index/divisions + 1.)*numpy.pi
        elif not hi_mag_slope: #i.e., FLAT, dip dir is arbitrary
            self._surface_dip_direction = 0.
        else: #dips EAST
            self._surface_dip_direction = numpy.pi*hi_mag_slope_index/divisions    
        self._surface_slope = numpy.fabs(hi_mag_slope)
        print 'The slope under the crater cavity footprint is: ', self._surface_slope

        
    def snap_coords_to_grid(self, grid, xcoord, ycoord):
        '''
        This method takes existing coordinates, inside the grid, and returns the ID of the closest grid node.
        '''
            #There should probably be an exception here: return None?
        vertices_array = grid.get_nodes_around_point(self._xcoord, self._ycoord)
        distances_to_vertices = numpy.empty(len(vertices_array),dtype=float)
        counter = 0
        for x in vertices_array:
            distances_to_vertices[counter] = numpy.sqrt((xcoord-grid.node_x[x])**2. + (ycoord-grid.node_y[x])**2.)
            counter += 1
        return vertices_array[numpy.argmin(distances_to_vertices)]


    def set_elev_change_at_pts(self, grid, data):
        '''
        This is a method to take an existing impact properties and a known nearest node to the impact site, and alter the topography to model the impact. It assumes crater radius and depth are known, models cavity shape as a power law where n is a function of R/D, and models ejecta thickness as an exponential decay,sensitive to both ballistic range from tilting and momentum transfer in impact (after Furbish). We DO NOT yet model transition to peak ring craters, or enhanced diffusion by ejecta in the strength regime. Peak ring craters are rejected from the distribution. This version of this method is designed to remove the sheer walls around the edges of craters, and replace them with a true dipping rim.
        '''
        
        #Load in the data, for speed:
        _angle_to_horizontal=self._angle_to_horizontal
        _surface_slope=self._surface_slope
        _surface_dip_direction = self._surface_dip_direction
        _azimuth_of_travel = self._azimuth_of_travel
        pi = numpy.pi
        twopi = 2.*pi
        tan = numpy.tan
        cos = numpy.cos
        sin = numpy.sin
        sqrt = numpy.sqrt
        arctan = numpy.arctan
        _radius = self._radius
        elev = data.elev
        
        #Build a list of nodes to work on. Starts just with the node closest to the center.
        crater_node_list = deque([self.closest_node_index])
        elev_changes = deque()
        #Build an array of flags into the nodelist of the grid to note whether that node has been placed in the list this loop:
        flag_already_in_the_list = numpy.zeros(grid.get_count_of_all_nodes())
        
        #Derive the exponent for the crater shape, shared betw simple & complex:
        crater_bowl_exp = self.get_crater_shape_exp()
        print 'Crater shape exponent: ', crater_bowl_exp
        
        #print 'Radius is: ', _radius
        
        #Derive the effective angle for impact, relative to the surface, beta_eff. The direction is always controlled by the impactor. Draw new impact angles if impact geomtery is impossible.
        while 1:
            #This ugly function applies a correction to the azimuth of travel which will be the actual angle the ejecta is propelled along. This is important since if, e.g., _angle_to_horizontal ~pi/2, then the surface dip direction needs to be the dominant term.
            denominator = cos(_angle_to_horizontal) + sin(_surface_slope)*cos(_surface_dip_direction - _azimuth_of_travel)
            #SOURCE OF THE DIVBYZERO!!!
            #print 'Terms in denominator: ', cos(_angle_to_horizontal),  sin(_surface_slope), cos(_surface_dip_direction - _azimuth_of_travel)
            tan_angle_from_imp_az = sin(_surface_slope)*sin(_surface_dip_direction - _azimuth_of_travel) / denominator
            if denominator > 0.:
                _ejecta_azimuth = (_azimuth_of_travel + arctan(tan_angle_from_imp_az) + twopi)%(twopi)
            else:
                _ejecta_azimuth = (_azimuth_of_travel - arctan(tan_angle_from_imp_az) + pi)%(twopi)
            print 'Angle from impact azimuth: ', _ejecta_azimuth - _azimuth_of_travel
            beta_eff_options = [numpy.nan] * 2
            beta_eff_options[0] = _angle_to_horizontal*cos(_ejecta_azimuth - _azimuth_of_travel) - _surface_slope * cos(_ejecta_azimuth - _surface_dip_direction)
            beta_eff_options[1] = (pi - _angle_to_horizontal*cos(_ejecta_azimuth - _azimuth_of_travel)) + _surface_slope * cos(_azimuth_of_travel - _surface_dip_direction)
            beta_eff = min(beta_eff_options)
            #Need to think thru if this reversal is still required:
            #if beta_eff_options[0] > beta_eff_options[1]:
            #_azimuth_of_travel = (_azimuth_of_travel + pi)%(2.*pi)
            print 'Beta effective: ', beta_eff
            if beta_eff >= 0.:
                break
            else:
                print 'Refreshing the impactor angle'
                self.set_impactor_angles()
                _azimuth_of_travel = self._azimuth_of_travel
                _angle_to_horizontal = self._angle_to_horizontal
        _impactor_angle_to_surface_normal = 0.5*pi - beta_eff
        print 'Impact, ejecta azimuths: ', _azimuth_of_travel, _ejecta_azimuth
        #print 'Beta effective: ', beta_eff
        #rim_height = self._cavity_volume/(3.667*pi*_radius**2.)

        #Derive the new ejecta shape params:
        tan_repose = tan(32.*pi/180.)
        _b = 3.667*_radius*_radius*tan_repose
        _c = self._cavity_volume - 3.667*_radius*_radius*_radius*tan_repose - 0.5*_radius*_radius*tan_repose
        radius_calc = (-_b + sqrt(_b*_b - 2.*tan_repose*_c))/tan_repose
        thickness_at_rim = (_radius - radius_calc)*tan_repose
        #...where thickness(r) = thickness_at_rim*(r/R_true)**-2.75

        #This code crawls out iteratively over the grid under the crater footprint, away from the centerpoint. This could be made much faster in Python if it wasn't a loop. Maybe an alternative would be to define the radius corresponding to the minimum loop, then trim out a box 2*R by 2*R to operate on, and perform all actions on arrays for speed.
        while 1:
            try:
                active_node = crater_node_list.popleft() #i.e., array is not empty
            except:
                break
            else:
                active_node_x_to_center = grid.node_x[active_node] - self._xcoord
                active_node_y_to_center = grid.node_y[active_node] - self._ycoord
                active_node_r_to_center = sqrt(active_node_x_to_center*active_node_x_to_center + active_node_y_to_center*active_node_y_to_center)
                #Special case for if the impact location is right on a gridline:
                if not active_node_x_to_center:
                    if active_node_y_to_center < 0.:
                        active_node_theta = pi
                    else:
                        active_node_theta = 0.
                else:
                    active_node_angle_to_yaxis = arctan(active_node_y_to_center/active_node_x_to_center)
                    if active_node_x_to_center < 0.:
                        active_node_theta = 1.5 * pi - active_node_angle_to_yaxis
                    else:
                        active_node_theta = 0.5 * pi - active_node_angle_to_yaxis

                #In this version, our sweep out from the center will need to account for deposition depth elevating the crater rim, i.e., we need to deposit *before* we cut the cavity. We do this by defining three domains for the node to lie in: 1. r<r_calc, i.e., below the pre-impact surface. No risk of intersecting the surface here. 2. r_calc < r; Th>z_new. this is the domain in the inward sloping rim of the crater ejecta. 3. Th<z_new and beyond. out on the ejecta proper. Note - (1) is not hard & fast rule if the surface dips. Safer is just (Th-lowering)<z_new
                #So, calc the excavation depth for all nodes, just to be on the safe side for strongly tilted geometries:
                if active_node_r_to_center <= _radius:
                    new_z = self.closest_node_elev - self._depth * (1. - (active_node_r_to_center/_radius)**crater_bowl_exp) + thickness_at_rim
                else:
                    new_z = self.closest_node_elev + thickness_at_rim + (active_node_r_to_center-_radius)*tan_repose
                
                #print new_z

                if new_z<elev[active_node]: #Below the original surface
                    #print 'Under surface'
                    elev_changes.append(new_z-elev[active_node])
                    if self._crater_type: #Complex crater, adjust for the peak
                        if active_node_r_to_center <= self._complex_peak_radius: #on the peak
                            new_z  = new_z + self._complex_peak_str_uplift * (1. - active_node_r_to_center/self._complex_peak_radius)
                    #Set the ground elev:
                    elev[active_node] = new_z
                    #Add the neighboring nodes which haven't already been adjusted to the for-processing array:
                    neighbors_active_node = grid.get_neighbor_list(active_node)
                    for x in neighbors_active_node:
                        if not flag_already_in_the_list[x]:
                            if x!=-1: #Not an edge
                                crater_node_list.append(x)
                                flag_already_in_the_list[x] = 1

                else:
                    #Calc the ejecta thickness for a perp. impact. Note it's fine if r<R_true
                    local_flat_thickness = thickness_at_rim*(active_node_r_to_center/_radius)**-2.75
                    #print 'Local flat thickness: ', local_flat_thickness

                    # Modulating functions are f(impactor az, dip dir, impactor angle, dip), the ejecta mass term, and mu(impactor az, dip dir, impactor angle, dip), the distance travelled term. Momentum is the dominant term over geometric effects in mu, by ~an OoM - o just model as a momentum effect.
                    #Thickness = thickness_flat * f_theta / mu
                    #beta_eff (==pi/2 - angle from normal) is already known.
                    theta_eff = _ejecta_azimuth - active_node_theta #This is the angle of the center-to-active-node line to the azimuth along which the ejecta is concentrated
                    tan_beta = tan(_impactor_angle_to_surface_normal*self._beta_factor)
                    sin_theta_sqd = sin(theta_eff) ** 2.
                    cos_theta = cos(theta_eff)
                    #REMEMBER, as tan_beta gets >1, the function describing the ejecta is only valid over ever more restricted ranges of theta!! In other words,
                    if tan_beta and sin_theta_sqd*tan_beta*tan_beta > 1.:
                        thickness_at_active_node = 0.
                        #print 'Outside ejecta cone'
                    else:
                        mu_theta_by_mu0 = tan_beta * cos_theta + sqrt(1. - sin_theta_sqd * tan_beta*tan_beta)
                        f_theta = (tan_beta*tan_beta*(cos_theta*cos_theta-sin_theta_sqd) + 2.*tan_beta*cos_theta*sqrt(1.-tan_beta*tan_beta*sin_theta_sqd) + 1.) / twopi
                        #print f_theta/mu_theta_by_mu0
                        #So, distn_at_angle = distn_vertical_impact*f_theta/mu_theta_by_mu0. Draw the thickness at the active node:
                        #NB-the 2pi is to correct for mismatch in the dimensions of mu and f
                        thickness_at_active_node = f_theta/mu_theta_by_mu0 * twopi * local_flat_thickness
                        if thickness_at_active_node < 0.:
                            thickness_at_active_node = 0.
                    #Now, are we inside or outside the rim?
                    if new_z <= (elev[active_node] + thickness_at_active_node): #inside the rim
                        #print 'Inside the rim, on the ejecta'
                        elev_changes.append(new_z-elev[active_node])
                        elev[active_node] = new_z
                    else: #outside the rim
                        #print 'Outside the rim'
                        elev_changes.append(thickness_at_active_node)
                        elev[active_node] = elev[active_node] + thickness_at_active_node
                    #Add the neighbors to the list, but only if the thickness of the current ejecta layer was nontrivial. Smallest crater is 5m diam, so 1m deep - so its max ejecta thickness at the rim is 10cm! We should probably resolve down to, say, 2.5cm, which is <4 crater radii from the rim of our smallest crater.
                    if thickness_at_active_node > self._minimum_ejecta_thickness:
                        neighbors_active_node = grid.get_neighbor_list(active_node)
                        for x in neighbors_active_node:
                            if not flag_already_in_the_list[x]:
                                if x!=-1:
                                    crater_node_list.append(x)
                                    flag_already_in_the_list[x] = 1
        #print 'f/mu is ', f_theta/mu_theta_by_mu0
        print 'Thickness at rim: ', thickness_at_rim
        
        #Save any data to the higher level:
        self.mass_balance_in_impact = numpy.mean(elev_changes)
        self.ejecta_azimuth = _ejecta_azimuth
        self.impactor_angle_to_surface_normal = _impactor_angle_to_surface_normal

    def set_elev_change_across_grid(self, grid, data):
        '''
        This is a method to take an existing impact properties and a known nearest node to the impact site, and alter the topography to model the impact. It assumes crater radius and depth are known, models cavity shape as a power law where n is a function of R/D, and models ejecta thickness as an exponential decay,sensitive to both ballistic range from tilting and momentum transfer in impact (after Furbish). We DO NOT yet model transition to peak ring craters, or enhanced diffusion by ejecta in the strength regime. Peak ring craters are rejected from the distribution. This version of this method is designed to remove the sheer walls around the edges of craters, and replace them with a true dipping rim.
        This routine differs from set_elev_change_at_pts() as it adjusts elevations at all nodes in the grid, rather than crawling out from the central impact point to a threshold depth. This has the potential to be faster for some craters, as we can operate on the whole grid in C without looping.
        '''

        #Load in the data, for speed:
        _angle_to_horizontal=self._angle_to_horizontal
        _surface_slope=self._surface_slope
        _surface_dip_direction = self._surface_dip_direction
        _azimuth_of_travel = self._azimuth_of_travel
        pi = numpy.pi
        twopi = 2.*pi
        tan = numpy.tan
        cos = numpy.cos
        sin = numpy.sin
        sqrt = numpy.sqrt
        arctan = numpy.arctan
        where = numpy.where
        _radius = self._radius
        elev = data.elev
        
        #Derive the exponent for the crater shape, shared betw simple & complex:
        crater_bowl_exp = self.get_crater_shape_exp()
        
        #Derive the effective angle for impact, relative to the surface, beta_eff. The direction is always controlled by the impactor. Draw new impact angles if impact geomtery is impossible.
        while 1:
            #This ugly function applies a correction to the azimuth of travel which will be the actual angle the ejecta is propelled along. This is important since if, e.g., _angle_to_horizontal ~pi/2, then the surface dip direction needs to be the dominant term.
            denominator = cos(_angle_to_horizontal) + sin(_surface_slope)*cos(_surface_dip_direction - _azimuth_of_travel)
            #print 'Terms in denominator: ', cos(_angle_to_horizontal),  sin(_surface_slope), cos(_surface_dip_direction - _azimuth_of_travel)
            tan_angle_from_imp_az = sin(_surface_slope)*sin(_surface_dip_direction - _azimuth_of_travel) / denominator
            if denominator > 0.:
                _ejecta_azimuth = (_azimuth_of_travel + arctan(tan_angle_from_imp_az) + twopi)%(twopi)
            else:
                _ejecta_azimuth = (_azimuth_of_travel - arctan(tan_angle_from_imp_az) + pi)%(twopi)
            print 'Angle from impact azimuth: ', _ejecta_azimuth - _azimuth_of_travel
            beta_eff_options = [numpy.nan] * 2
            beta_eff_options[0] = _angle_to_horizontal*cos(_ejecta_azimuth - _azimuth_of_travel) - _surface_slope * cos(_ejecta_azimuth - _surface_dip_direction)
            beta_eff_options[1] = (pi - _angle_to_horizontal*cos(_ejecta_azimuth - _azimuth_of_travel)) + _surface_slope * cos(_azimuth_of_travel - _surface_dip_direction)
            beta_eff = min(beta_eff_options)
            #Need to think thru if this reversal is still required:
            #if beta_eff_options[0] > beta_eff_options[1]:
            #_azimuth_of_travel = (_azimuth_of_travel + pi)%(2.*pi)
            print 'Beta effective: ', beta_eff
            if beta_eff >= 0.:
                break
            else:
                print 'Refreshing the impactor angle'
                self.set_impactor_angles()
                _azimuth_of_travel = self._azimuth_of_travel
                _angle_to_horizontal = self._angle_to_horizontal
        _impactor_angle_to_surface_normal = 0.5*pi - beta_eff
        tan_beta = tan(_impactor_angle_to_surface_normal*self._beta_factor)
        tan_beta_sqd = tan_beta*tan_beta
        print 'Impact, ejecta azimuths: ', _azimuth_of_travel, _ejecta_azimuth
        #print 'Beta effective: ', beta_eff
        #rim_height = self._cavity_volume/(3.667*pi*_radius**2.)

        #Derive the new ejecta shape params:
        tan_repose = tan(32.*pi/180.)
        _b = 3.667*_radius*_radius*tan_repose
        _c = self._cavity_volume - 3.667*_radius*_radius*_radius*tan_repose - 0.5*_radius*_radius*tan_repose
        radius_calc = (-_b + sqrt(_b*_b - 2.*tan_repose*_c))/tan_repose
        thickness_at_rim = (_radius - radius_calc)*tan_repose
        #...where thickness(r) = thickness_at_rim*(r/R_true)**-2.75
        
        #Here the code diverges from the crawling version:
        _vec_x_to_center = grid.get_node_x_coords() - self._xcoord
        _vec_y_to_center = grid.get_node_y_coords() - self._ycoord
        _vec_r_to_center = sqrt(_vec_x_to_center*_vec_x_to_center + _vec_y_to_center*_vec_y_to_center)
        try:
            _vec_angle_to_yaxis = arctan(_vec_y_to_center/_vec_x_to_center)
        except: #These cases have the impact right on a gridline.
            _vec_theta = numpy.empty(len(_vec_x_to_center))
            nonzero_nodes = numpy.nonzero(_vec_x_to_center)
            _vec_angle_to_yaxis = arctan(_vec_y_to_center[nonzero_nodes]/_vec_x_to_center[nonzero_nodes])
            _vec_theta[nonzero_nodes] = where(_vec_x_to_center[nonzero_nodes]<0,1.5*pi-_vec_angle_to_yaxis,0.5*pi-_vec_angle_to_yaxis)
            zero_nodes = where(_vec_x_to_center==0.)[0]
            _vec_theta[zero_nodes] = where(_vec_y_to_center[zero_nodes]<0,pi,0.)
        else: #the normal case
            _vec_theta = where(_vec_x_to_center<0,1.5*pi-_vec_angle_to_yaxis,0.5*pi-_vec_angle_to_yaxis)

        #We need to account for deposition depth elevating the crater rim, i.e., we need to deposit *before* we cut the cavity. We do this by defining three domains for the node to lie in: 1. r<r_calc, i.e., below the pre-impact surface. No risk of intersecting the surface here. 2. r_calc < r; Th>z_new. this is the domain in the inward sloping rim of the crater ejecta. 3. Th<z_new and beyond. out on the ejecta proper. Note - (1) is not hard & fast rule if the surface dips. Safer is just (Th-lowering)<z_new
        #So, calc the excavation depth for all nodes, just to be on the safe side for strongly tilted geometries:
        _vec_new_z = where(_vec_r_to_center<=_radius, self.closest_node_elev - self._depth * (1. - (_vec_r_to_center/_radius)**crater_bowl_exp) + thickness_at_rim, self.closest_node_elev + thickness_at_rim + (_vec_r_to_center-_radius)*tan_repose)
    #!!Need to come back to loading data into elev_changes, which keeps track of the mass balance - tho it's less important here.
        nodes_below_surface = where(_vec_new_z<elev)[0]
        nodes_above_surface = where(_vec_new_z>=elev)[0]
        #Check if we need to adjust for a central peak
        if self._crater_type:
            central_peak_pts = where(_vec_r_to_center<=self._complex_peak_radius)
            _vec_new_z[central_peak_pts] = _vec_new_z[central_peak_pts] + self._complex_peak_str_uplift * (1. - _vec_r_to_center[central_peak_pts]/self._complex_peak_radius)
        #Set the ground elev for below ground nodes
        elev[nodes_below_surface] = _vec_new_z[nodes_below_surface]
        #From here on, our new arrays will only be as long as nodes_above_surface
        _vec_flat_thickness_above_surface = thickness_at_rim*(_vec_r_to_center[nodes_above_surface]/_radius)**-2.75
        _vec_theta_eff = _ejecta_azimuth - _vec_theta[nodes_above_surface] #This is the angle of the center-to-active-node line to the azimuth along which the ejecta is concentrated
        _vec_sin_theta_sqd = sin(_vec_theta_eff) ** 2.
        _vec_cos_theta = cos(_vec_theta_eff)
        #REMEMBER, as tan_beta gets >1, the function describing the ejecta is only valid over ever more restricted ranges of theta!! In other words,
        nodes_inside_ejecta = where(_vec_sin_theta_sqd*tan_beta_sqd <= 1.) #these are indices to an array of length nodes_above_surface only
        _vec_thickness = numpy.zeros(nodes_above_surface.size) #i.e., it's zero outside the ejecta
        _vec_mu_theta_by_mu0 = tan_beta * _vec_cos_theta[nodes_inside_ejecta] + sqrt(1. - _vec_sin_theta_sqd[nodes_inside_ejecta] * tan_beta_sqd)
        _vec_f_theta = (tan_beta_sqd*(_vec_cos_theta[nodes_inside_ejecta]**2.-_vec_sin_theta_sqd[nodes_inside_ejecta]) + 2.*tan_beta*_vec_cos_theta[nodes_inside_ejecta]*sqrt(1.-tan_beta_sqd*_vec_sin_theta_sqd[nodes_inside_ejecta]) + 1.) / twopi
        #So, distn_at_angle = distn_vertical_impact*f_theta/mu_theta_by_mu0. Draw the thickness at the active node:
        #NB-the 2pi is to correct for mismatch in the dimensions of mu and f
        thickness_at_nodes_under_ejecta = _vec_f_theta/_vec_mu_theta_by_mu0 * twopi * _vec_flat_thickness_above_surface[nodes_inside_ejecta]
        #Set the thicknesses <0 to 0:
        thickness_at_nodes_under_ejecta = where(thickness_at_nodes_under_ejecta>=0.,thickness_at_nodes_under_ejecta, 0.)
        #make the ejecta thickness map for all nodes above surface:
        _vec_thickness[nodes_inside_ejecta] = thickness_at_nodes_under_ejecta
        #Now, are we inside or outside the rim?
        elev[nodes_above_surface] = where(_vec_new_z[nodes_above_surface]<=elev[nodes_above_surface]+_vec_thickness,_vec_new_z[nodes_above_surface], elev[nodes_above_surface]+_vec_thickness)

        
        #Save any data to the higher level:
        self.mass_balance_in_impact = 0. #this isn't true, but now effectively a flag that we've used this method
        self.ejecta_azimuth = _ejecta_azimuth
        self.impactor_angle_to_surface_normal = _impactor_angle_to_surface_normal


    def excavate_a_crater(self, grid, data):
        '''
        This method executes the most of the other methods of this crater class, and makes the geomorphic changes to a mesh associated with a single bolide impact with randomized properties. It receives parameters of the model grid, and the vector data storage class. It is the primary interface method of this class.
        Unless devtesting, this method has been superceded by excavate_a_crater_optimized().
        '''
        self.set_cr_radius_from_shoemaker(data)
        #self._radius = forced_size
        print 'Radius: ', self._radius
        self.set_depth_from_size()
        self.set_crater_volume()
        self.set_coords(grid, data)
        self.set_impactor_angles()
        self.set_crater_mean_slope_v2(grid, data)
        if numpy.isnan(self._surface_slope):
            print 'Surface slope is not defined for this crater! Is it too big? Crater will not be drawn.'
        else:
            self.set_elev_change_at_pts(grid, data)
            print 'Impactor angle to ground normal: ', self.impactor_angle_to_surface_normal
            print 'Mean mass balance/px: ', self.mass_balance_in_impact
        print '*****'
        #Record the data:
        data.impact_sequence.append({'x': self._xcoord, 'y': self._ycoord, 'r': self._radius, 'volume': self._cavity_volume, 'surface_slope': self._surface_slope, 'normal_angle': self.impactor_angle_to_surface_normal, 'impact_az': self._azimuth_of_travel, 'ejecta_az': self.ejecta_azimuth, 'mass_balance': self.mass_balance_in_impact})

    def excavate_a_crater_whole_grid(self, grid, data):
        '''
            This method executes the most of the other methods of this crater class, and makes the geomorphic changes to a mesh associated with a single bolide impact with randomized properties. It receives parameters of the model grid, and the vector data storage class. It is the primary interface method of this class.
            Unless devtesting, this method has been superceded by excavate_a_crater_optimized().
            '''
        self.set_cr_radius_from_shoemaker(data)
        #self._radius = forced_size
        print 'Radius: ', self._radius
        self.set_depth_from_size()
        self.set_crater_volume()
        self.set_coords(grid, data)
        self.set_impactor_angles()
        self.set_crater_mean_slope_v2(grid, data)
        if numpy.isnan(self._surface_slope):
            print 'Surface slope is not defined for this crater! Is it too big? Crater will not be drawn.'
        else:
            self.set_elev_change_across_grid(grid, data)
            print 'Impactor angle to ground normal: ', self.impactor_angle_to_surface_normal
            print 'Mean mass balance/px: ', self.mass_balance_in_impact
        print '*****'
        #Record the data:
        data.impact_sequence.append({'x': self._xcoord, 'y': self._ycoord, 'r': self._radius, 'volume': self._cavity_volume, 'surface_slope': self._surface_slope, 'normal_angle': self.impactor_angle_to_surface_normal, 'impact_az': self._azimuth_of_travel, 'ejecta_az': self.ejecta_azimuth, 'mass_balance': self.mass_balance_in_impact})

    def excavate_a_crater_optimized(self, grid, data, forced_radius=numpy.nan, forced_angle=numpy.nan, forced_pos=numpy.nan):
        '''
            This method executes the most of the other methods of this crater class, and makes the geomorphic changes to a mesh associated with a single bolide impact with randomized properties. It receives parameters of the model grid, and the vector data storage class. It is the primary interface method of this class.
            This method is optimized to not sweep the whole grid if the crater is small.
            A fixed crater size can be specified with the input variable "forced_radius" (in km), and a fixed impact angle with "forced_angle" (in degrees from vertical - impact azimuth will always be assumed as travel eastwards). Position can be specified with forced_pos, which takes an array-like object with two entries, which are the x and y coordinate in relative position on the grid (e.g., [0.5, 0.5]).
            '''
        if numpy.isnan(forced_radius):
            self.set_cr_radius_from_shoemaker(data)
        else:
            self._radius = forced_radius
        #self._radius = forced_size
        print 'Radius: ', self._radius
        self.set_depth_from_size()
        self.set_crater_volume()
        if numpy.any(forced_pos):
            self.set_coords(grid, data)
        else:
            self._xcoord = forced_pos[0]*(mg.get_grid_xdimension()-mg.dx)
            self._ycoord = forced_pos[1]*(mg.get_grid_ydimension()-mg.dx)
            print self._xcoord, self._ycoord
        if numpy.isnan(forced_angle):
            self.set_impactor_angles()
            #print 'Angle was NaN!'
        else:
            self._angle_to_horizontal = forced_angle/numpy.pi*180.
            self._azimuth_of_travel = 0.5*numpy.pi
        self.set_crater_mean_slope_v2(grid, data)
        if numpy.isnan(self._surface_slope):
            print 'Surface slope is not defined for this crater! Is it too big? Crater will not be drawn.'
        else:
            #NB - this is an empirical optimization. If the minimum crater thickness in set_elev_change_at_pts() changes, so will the optimized value! 0.000001 min thickness == 0.024 km radius here. Normalize by R=0.0025 here => 9.6
            if self._radius/grid.dx < 9.6:
                self.set_elev_change_at_pts(grid, data)
            else:
                self.set_elev_change_across_grid(grid, data)
            print 'Impactor angle to ground normal: ', self.impactor_angle_to_surface_normal
            print 'Mean mass balance/px: ', self.mass_balance_in_impact
        print '*****'
        #Record the data:
        data.impact_sequence.append({'x': self._xcoord, 'y': self._ycoord, 'r': self._radius, 'volume': self._cavity_volume, 'surface_slope': self._surface_slope, 'normal_angle': self.impactor_angle_to_surface_normal, 'impact_az': self._azimuth_of_travel, 'ejecta_az': self.ejecta_azimuth, 'mass_balance': self.mass_balance_in_impact})


#The functions in this segment give control over the execution of this module. Adjusty the params inside the functions to get different effects, and the final line of the file to determine whether you get one crater, or lots of craters.

def dig_some_craters_on_fresh_surface():
    '''
    Ad hoc driver code to make this file run as a standalone:
    '''
    #User-defined params:
    nr = 500
    nc = 500
    dx = 0.0025
    #dt = 1.
    nt = 1000

    #Setup
    mg = RasterModelGrid()
    mg.initialize(nr, nc, dx)
    vectors = data(mg)
    vectors.elev[:] = 100.
    cr = impactor()

    #Update until
    for i in xrange(0,nt):
        print 'Crater number ', i
        cr.excavate_a_crater_optimized(mg, vectors)

    #Finalize
    elev_raster = mg.node_vector_to_raster(vectors.elev, flip_vertically=True)
    #contour(elev_raster)

    imshow(elev_raster)
    colorbar()
    show()
    vectors.viewing_raster = elev_raster
    return cr, mg, vectors

def dig_some_craters(grid, data, nt_in=10000, min_radius=0.005):
    '''
    Takes an existing DTM and peppers it with craters.
    '''
    #dt = 1.

    #Setup
    cr = impactor(min_radius)

    #Update until
    for i in xrange(0,nt_in):
        print 'Crater number ', i
        cr.excavate_a_crater_optimized(grid, data)
    
    #Finalize
    elev_raster = grid.node_vector_to_raster(data.elev, flip_vertically=True)
    #contour(elev_raster)
    
    profile = plot(elev_raster[300,:])
    xsec = plot(elev_raster[:,100])
    #imshow(elev_raster)
    #colorbar()
    #show()
    data.viewing_raster = elev_raster
    return grid, data, profile, xsec
    

def dig_one_crater(nr, nc, dx, rel_x, rel_y, radius):
    '''
    This is an ad-hoc script to dig one crater.
    '''
    #User-defined params:
    #nr = 1200
    #nc = 1200
    #dx = 0.0025
    #dt = 1.
    #nt = 1

    #Setup
    mg = RasterModelGrid()
    mg.initialize(nr, nc, dx)
    vectors = data(mg)
    vectors.elev[:] = 0.
    cr = impactor()

    cr._radius = radius
    print 'Radius: ', cr._radius
    cr.set_depth_from_size()
    print 'Depth: ', cr._depth
    cr.set_crater_volume()
    cr._xcoord = rel_x*mg.get_grid_xdimension()
    cr._ycoord = rel_y*mg.get_grid_ydimension()
    vertices_array = mg.get_nodes_around_point(cr._xcoord, cr._ycoord)
    distances_to_vertices = []
    for x in vertices_array:
        distances_to_vertices.append(numpy.sqrt((cr._xcoord-mg.node_x[x])**2. + (cr._ycoord-mg.node_y[x])**2.))
    cr.closest_node_index = vertices_array[numpy.argmin(distances_to_vertices)]
    cr.closest_node_elev = vectors.elev[cr.closest_node_index]

    cr._angle_to_horizontal = numpy.pi*0.5*2./3.
    cr._azimuth_of_travel = numpy.pi*1.5
    cr.set_crater_mean_slope_v2(mg, vectors)
    print 'Azimuth of travel: ', cr._azimuth_of_travel
    print 'Angle of ground: ', cr._surface_slope
    print 'Dip direction of ground: ', cr._surface_dip_direction
    cr.set_elev_change_across_grid(mg, vectors)
    print 'Impact angle to ground normal: ', cr.impactor_angle_to_surface_normal

    #Finalize
    elev_raster = mg.node_vector_to_raster(vectors.elev, flip_vertically=True)
    #contour(elev_raster)
    #imshow(elev_raster)
    #colorbar()
    #show()
    vectors.viewing_raster = elev_raster
    return cr, mg, vectors


def one_crater_then_degrade():
    #start_time = time.time()
    #cr, mg, vectors = dig_some_craters_on_fresh_surface()
    #cr, mg, vectors = dig_one_crater(120, 120, 0.025, 1., 0.5, 1.)
    #mg_10k, vectors_10k = dig_some_craters(mg, vectors)

    #This code builds a dictionary that contains time slices for each 10k craters hitting a surface:
    #How many times round?
    loops = 50 #500,000 craters
    #Build the dictionary:
    crater_time_sequ = {}
    profile_list = []
    xsec_list = []
    #Initialize the starting condition:
    cr, mg, vectors = dig_one_crater(1200, 1200, 0.0025, 0.5, 0.5, .75)
    #Save the starting conds:
    crater_time_sequ[0] = copy(vectors)
    #Run the loops
    for i in xrange(0,loops):
        mg, vectors, profile, xsec = dig_some_craters(mg, vectors)
        crater_time_sequ[i] = copy(vectors)
        profile_list.append(profile)
        xsec_list.append(xsec)
    show(profile_list)
    #end_time = time.time()
    #print('Elapsed time was %g seconds' % (end_time - start_time))
    return crater_time_sequ

def ten_times_reduction(mg_in, vectors_in, loops=25):
    """
    Depreciated in favour of step_reduce_size.
    """
    crater_time_sequ_50_m_min = {}
    crater_time_sequ_5_m_min = {}
    profile_list = []
    xsec_list = []
    for i in xrange(0,loops):
        mg_in, vectors_in, profile, xsec = dig_some_craters(mg_in, vectors_in, nt_in=10000, min_radius=0.05)
        crater_time_sequ_50_m_min[i] = copy(vectors_in)
    for i in xrange(0, loops):
        mg_in, vectors_in, profile, xsec = dig_some_craters(mg_in, vectors_in, nt_in=10000, min_radius=0.005)
        crater_time_sequ_5_m_min[i] = copy(vectors_in)
    return crater_time_sequ_50_m_min, crater_time_sequ_5_m_min

def step_reduce_size(mg_in, vectors_in, loops=[25,25], interval=10000, min_radius_in=[0.05, 0.005]):
    crater_time_sequ_1st = {}
    crater_time_sequ_2nd = {}
    profile_list = []
    xsec_list = []
    for i in xrange(0,loops[0]):
        mg_in, vectors_in, profile, xsec = dig_some_craters(mg_in, vectors_in, nt_in=interval, min_radius=min_radius_in[0])
        crater_time_sequ_1st[i] = copy(vectors_in)
    for i in xrange(0, loops[1]):
        mg_in, vectors_in, profile, xsec = dig_some_craters(mg_in, vectors_in, nt_in=interval, min_radius=min_radius_in[1])
        crater_time_sequ_2nd[i] = copy(vectors_in)
    return crater_time_sequ_1st, crater_time_sequ_2nd

def plot_hypsometry(plotting_rasters):
    """
    This function plots hypsometry (elev above minimum vs no. of px below it) for each of plotting rasters provided to the function. If plotting_rasters is a dictionary of data, where the plotting rasters are stored as data.plotting_raster, the function will return all entries. If it is a single plotting raster, or a single 'data' object, only that raster will plot. Call 'show()' manually after running this code.
    """
    
    def plot_a_hypsometry_curve(elevs, color=0):
        min_elev=numpy.amin(elevs)
        rel_relief = (elevs.flatten() - min_elev)/(numpy.amax(elevs)-min_elev)
        rel_relief.sort()
        hyps_x_axis = (1.-numpy.array(range(rel_relief.size), dtype=float)/rel_relief.size)
        #print hyps_x_axis
        #print rel_relief
        input_color = color*0.9
        plot(hyps_x_axis, rel_relief, color=str(input_color))
    
    if type(plotting_rasters) == dict:
        for i in range(len(plotting_rasters)):
            fraction = float(i)/len(plotting_rasters)
            plot_a_hypsometry_curve(plotting_rasters[i].viewing_raster, color=fraction)
    elif type(plotting_rasters) == numpy.ndarray:
        plot_a_hypsometry_curve(plotting_rasters)
    else:
        try:
            plot_a_hypsometry_curve(plotting_rasters.viewing_raster)
        except:
            print 'Input type not recognised!'
        
            


#if __name__=='__main__':
#    main()
