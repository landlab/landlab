#! /usr/env/python

import random
import math
import numpy
import sys
import scipy.optimize as opt
from landlab.model_grid import RasterModelGrid #this is the tMesh equivalent module

#these ones only so we can run this module ad-hoc:
from pylab import plot, draw, show, contour, imshow, colorbar
#from copy import copy

sys.setrecursionlimit(1500)

class data(object):
    '''
    This is where all the whole-grid data lives, as arrays over the various elements of the grid.
    '''
    #Data goes here!!!
    def __init__(self):
        self.elev = [] #some data
        self.flag_already_in_the_list = []
        self.craters_over_max_radius_not_plotted = []
        self.impact_sequence = []

class impactor(object):
    '''
    This class holds all parameters decribing properties of a single impact structure, and contains methods for recalculating fresh and internally consistent data describing such a impact structure.
    Built DEJH Spring 2013.
    '''
    def __init__(self):
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
        self._minimum_crater = 0.005 #km. This is the smallest modelled crater radius. 10m diameter is the strict cutoff for known cr distns
        self.ivanov_a = [-3.0876, -3.557528, 0.781027, 1.021521, -0.156012, -0.444058, 0.019977, 0.08685, -0.005874, -0.006809, 0.000825, 0.0000554] #The coefficients for Ivanov's crater distn function
        self._impactor_angle_to_surface = -999.
        self._angle_to_horizontal = -999.
        self._minimum_ejecta_thickness = 0.000001
        self._beta_factor = 0.5
        
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
        sum_terms = [(self.ivanov_a[0] - math.log10(N)),]
        for i in range(1,12):
            sum_terms.append(self.ivanov_a[i] * math.log10(D) ** i)
        return sum(sum_terms)
    
    
    #A second version, in the conventional form N = f(D):
    def ivanov_prod_equ_as_Nequals(self, D):
        '''
        This function expresses Ivanov's CPF in the normal form, N = f(D). It returns N (not log10(N)!).
        '''
        sum_terms = []
        for i in range (0,12):
            sum_terms.append(self.ivanov_a[i] * math.log10(D) ** i)
        return pow( 10., sum(sum_terms) )

    
    def ivanov_prod_equ_1stderiv(self, D, N):
        '''
        This is a helper function for the Newton-Raphson optimized solver in solve_ivanov_for_crater_diam. It returns the first derivative for the production function. Note that the N term is a dummy variable used for consistency in the opt.newton(), and does nothing here.
        '''
        sum_terms = []
        ln_of_10 = math.log(10)
        for i in range(1,12):
            sum_terms.append(self.ivanov_a[i] * i * math.log(D) ** (i-1.) / ln_of_10**i)
        return (sum(sum_terms) / D)


    def ivanov_prod_equ_2ndderiv(self, D):
        '''
        This is a helper function for the Newton-Raphson optimized solver in solve_ivanov_for_crater_diam. It returns the second derivative for the production function.
        '''
        ln_of_10 = math.log(10)
        sum_terms = [self.ivanov_a[1] / ln_of_10,]
        for i in range(2,12):
            sum_terms.append(i * self.ivanov_a[i] * (math.log(D) - i + 1.) * math.log(D)**(i-2.) / ln_of_10**i )
        return (-sum(sum_terms) / D**2.)


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
        self._radius = self._minimum_crater*(random.random())**-0.345
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
        self._xcoord = random.random() * (grid.get_grid_xdimension()-grid.dx)
        self._ycoord = random.random() * (grid.get_grid_ydimension()-grid.dx)
        #Snap impact to grid:
        vertices_array = grid.get_nodes_around_point(self._xcoord, self._ycoord)
        distances_to_vertices = []
        for x in vertices_array:
            #print 'The ID is:', x
            distances_to_vertices.append(numpy.sqrt((self._xcoord-grid.x(x))**2. + (self._ycoord-grid.y(x))**2.))
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
        self._angle_to_horizontal =  numpy.arccos(random.random()) #gives sin distn with most values drawn nearer to 0
        self._azimuth_of_travel = random.random() * 2. * numpy.pi #equal chance of any azimuth
        #Shoemaker 1983 gives a speculative equn for the suppression of D by increasing impact angle (his equ 3), but effect is minor, and it's probably not worth the trouble.
        #Shoemaker 1962 (in book) apparently states low angle impacts are very rare. Need to read this.

    
    def set_size(self, data):
        '''
        This method draws a crater radius at random from the PDF describing crater sizes dictated by the Ivanov distribution, with the probability of obtaining that crater size depending on the relative abundance of that size. Note this method works with the Ivanov distribution according to diameter, but returns a radius.
        '''
        self._radius = 0.5 * self.solve_ivanov_for_crater_diam(random.random())
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
            if random.random() > py_complex_at_radius:
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
            if x < grid.get_grid_xdimension and y < grid.get_grid_ydimension:
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
        #It would be trivial to add more divisions, i.e., D16, D32, D64, etc. Code is setup for this; just add to the slope_pts lists.
        half_crater_radius = 0.707 * self._radius
        radial_points1 = []
        radial_points2 = []
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
        for x,y in slope_pts1:
            if 0. < x < grid.get_grid_xdimension and 0. < y < grid.get_grid_ydimension:
                radial_points1.append(self.snap_coords_to_grid(grid, x, y))
            else:
                radial_points1.append(self.closest_node_index)
        for x,y in slope_pts2:
            if 0. < x < grid.get_grid_xdimension and 0. < y < grid.get_grid_ydimension:
                radial_points2.append(self.snap_coords_to_grid(grid, x, y))
            else:
                radial_points2.append(self.closest_node_index)
        divisions = len(radial_points1)
        for a in range(0, divisions):
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
        hi_mag_slope_index = numpy.argmax(numpy.fabs(slope_array))
        hi_mag_slope = slope_array[hi_mag_slope_index]
        #print hi_mag_slope
        if hi_mag_slope > 0: #i.e., dips WEST (or dead S)
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
        distances_to_vertices = []
        for x in vertices_array:
            distances_to_vertices.append(numpy.sqrt((xcoord-grid.x(x))**2. + (ycoord-grid.y(x))**2.))
        return vertices_array[numpy.argmin(distances_to_vertices)]


    def set_elev_change_at_pts(self, grid, data):
        '''
        This method takes existing impact properties and a known nearest node to the impact site, and alters the topography to model the impact. It assumes crater radius and depth are known, models cavity shape as a power law where n is a function of R/D, and models ejecta thickness as an exponential decay,sensitive to both ballistic range from tilting and momentum transfer in impact (after Furbish). We DO NOT yet model transition to peak ring craters, or enhanced diffusion by ejecta in the strength regime. Peak ring craters are rejected from the distribution.
        '''
        #Build a list of nodes to work on. Starts just with the node closest to the center.
        crater_node_list = [self.closest_node_index]
        #Build an array of flags into the nodelist of the grid to note whether that node has been placed in the list this loop:
        data.flag_already_in_the_list = [0] * grid.ncells
        #Could we have the node properties as dictionaries of arrays?:
        #grid.nodes[adjusted] = zeros(grid.ncells)
        
        #Derive the exponent for the crater shape, shared betw simple & complex:
        crater_bowl_exp = self.get_crater_shape_exp()
        print 'Crater shape exponent: ', crater_bowl_exp

        #print 'Radius is: ', self._radius
        
        #Derive the effective angle for impact, relative to the surface, beta_eff. The direction is always controlled by the impactor. Draw new impact angles if impact geomtery is impossible.
        while 1:
            #This ugly function applies a correction to the azimuth of travel which will be the actual angle the ejecta is propelled along. This is important since if, e.g., _angle_to_horizontal ~pi/2, then the surface dip direction needs to be the dominant term.
            denominator = numpy.cos(self._angle_to_horizontal) + numpy.sin(self._surface_slope)*numpy.cos(self._surface_dip_direction - self._azimuth_of_travel)
            tan_angle_from_imp_az = numpy.sin(self._surface_slope)*numpy.sin(self._surface_dip_direction - self._azimuth_of_travel) / denominator
            if denominator > 0.:
                self._ejecta_azimuth = (self._azimuth_of_travel + numpy.arctan(tan_angle_from_imp_az) + 2.*numpy.pi)%(2.*numpy.pi)
            else:
                self._ejecta_azimuth = (self._azimuth_of_travel - numpy.arctan(tan_angle_from_imp_az) + numpy.pi)%(2.*numpy.pi)
            print 'Angle from impact azimuth: ', self._ejecta_azimuth - self._azimuth_of_travel
            beta_eff_options = [numpy.nan] * 2
            beta_eff_options[0] = self._angle_to_horizontal*numpy.cos(self._ejecta_azimuth - self._azimuth_of_travel) - self._surface_slope * numpy.cos(self._ejecta_azimuth - self._surface_dip_direction)
            beta_eff_options[1] = (numpy.pi - self._angle_to_horizontal*numpy.cos(self._ejecta_azimuth - self._azimuth_of_travel)) + self._surface_slope * numpy.cos(self._azimuth_of_travel - self._surface_dip_direction)
            beta_eff = min(beta_eff_options)
            #Need to think thru if this reversal is still required:
            #if beta_eff_options[0] > beta_eff_options[1]:
                #self._azimuth_of_travel = (self._azimuth_of_travel + numpy.pi)%(2.*numpy.pi)
            print 'Beta effective: ', beta_eff
            if beta_eff >= 0.:
                break
            else:
                print 'Refreshing the impactor angle'
                self.set_impactor_angles()
        self._impactor_angle_to_surface_normal = 0.5*numpy.pi - beta_eff
        print 'Impact, ejecta azimuths: ', self._azimuth_of_travel, self._ejecta_azimuth
        #print 'Beta effective: ', beta_eff
        rim_height = self._cavity_volume/(3.667*numpy.pi*self._radius**2.)
        
        while len(crater_node_list): #i.e., array is not empty
            low_cell_flag = False
            active_node = crater_node_list.pop(0)
            active_node_x_to_center = grid.x(active_node) - self._xcoord
            active_node_y_to_center = grid.y(active_node) - self._ycoord
            active_node_r_to_center = numpy.sqrt(active_node_x_to_center**2. + active_node_y_to_center**2.)
            #Special case for if the impact location is right on a gridline:
            if not active_node_x_to_center:
                if active_node_y_to_center < 0.:
                    active_node_theta = numpy.pi
                else:
                    active_node_theta = 0.
            else:
                active_node_angle_to_yaxis = numpy.arctan(active_node_y_to_center/active_node_x_to_center)
                if active_node_x_to_center < 0.:
                    active_node_theta = 1.5 * numpy.pi - active_node_angle_to_yaxis
                else:
                    active_node_theta = 0.5 * numpy.pi - active_node_angle_to_yaxis
        
            #flag the node as adjusted this timestep:
            #data.flag_adjusted_this_tstep[active_node] = 1
            
            if active_node_r_to_center <= self._radius: #In the cavity
                #print 'In the crater, r to center: ', active_node_r_to_center
                #For inside crater, needs to be an absolute set relative to the deepest point, not rel to surface.
                #We assume the elevation of the nearest node to the center is an adequate proxy for the elevation of the actual impact spot. We need to add the rim height too!
                new_z = self.closest_node_elev - self._depth * (1. - (active_node_r_to_center/self._radius)**crater_bowl_exp) + rim_height
                if numpy.isnan(new_z):
                    print 'New_z is nan!!!'
                    raw_input()
                if new_z > data.elev[active_node]:
                    new_z = data.elev[active_node]
                    low_cell_flag = True
                if self._crater_type: #Complex crater, adjust for the peak
                    if active_node_r_to_center <= self._complex_peak_radius: #on the peak
                        new_z  = new_z + self._complex_peak_str_uplift * (1. - active_node_r_to_center/self._complex_peak_radius)
                #Make sure I check these geometries!!
                #Set the ground elev:
                data.elev[active_node] = new_z
                #Add the neighboring nodes which haven't already been adjusted to the for-processing array:
                neighbors_active_node = grid.get_neighbor_list(active_node)
                for x in neighbors_active_node:
                    if not data.flag_already_in_the_list[x]:
                        if x!=-1: #Not an edge
                            crater_node_list.append(x)
                            data.flag_already_in_the_list[x] = 1
                
                #print crater_node_list
                #raw_input()
        
            if active_node_r_to_center > self._radius or low_cell_flag: #this is less efficient than else, but whatever
            #Outside the cavity, but less than the ejecta limit, or within the circle of the cavity, but surface is weirdly low (e.g., stongly sloping surface)
                #This is the tough part! For pure ejection,
                # Thickness = V / (3.667*3.141*R^2.) * (r/R)^-2.75
                #Note the exponent is dodgy, though Housen 83 puts theoretical bounds, 2.5<n<3.0 (his Table 2)
                # Modulating functions are f(impactor az, dip dir, impactor angle, dip), the ejecta mass term, and mu(impactor az, dip dir, impactor angle, dip), the distance travelled term. Momentum is the dominant term over geometric effects in mu, by ~an OoM - o just model as a momentum effect.
                #Thickness = thickness_flat * f_theta / mu
                #beta_eff is already known.
                theta_eff = self._ejecta_azimuth - active_node_theta #This is the angle of the center-to-active-node line to the azimuth along which the ejecta is concentrated
                tan_beta = numpy.tan(self._impactor_angle_to_surface_normal)
                sin_theta_sqd = numpy.sin(theta_eff) ** 2.
                cos_theta = numpy.cos(theta_eff)
                #REMEMBER, as tan_beta gets >1, the function describing the ejecta is only valid over ever more restricted ranges of theta!! In other words,
                if tan_beta and sin_theta_sqd > 1./tan_beta**2:
                    thickness_at_active_node = 0.
                else:
                    mu_theta_by_mu0 = tan_beta * cos_theta + numpy.sqrt(1. - sin_theta_sqd * tan_beta**2.)
                    f_theta = (tan_beta**2.*(cos_theta**2.-sin_theta_sqd) + 2.*tan_beta*cos_theta*numpy.sqrt(1.-tan_beta**2.*sin_theta_sqd) + 1.) / (2.*numpy.pi)
                    #print f_theta/mu_theta_by_mu0
                    #So, distn_at_angle = distn_vertical_impact*f_theta/mu_theta_by_mu0. Draw the thickness at the active node:
                    thickness_at_active_node = f_theta/mu_theta_by_mu0 * 2.*numpy.pi*rim_height*(active_node_r_to_center/self._radius)**-2.75
                if numpy.isnan(thickness_at_active_node):
                    print 'thickness is nan!!!'
                    raw_input()
                #Adjust the elev of the node:
                data.elev[active_node] = data.elev[active_node] + thickness_at_active_node
                #Add the neighbors to the list, but only if the thickness of the current ejecta layer was nontrivial. Smallest crater is 5m diam, so 1m deep - so its max ejecta thickness at the rim is 10cm! We should probably resolve down to, say, 2.5cm, which is <4 crater radii from the rim of our smallest crater.
                if thickness_at_active_node > 0.001:
                    neighbors_active_node = grid.get_neighbor_list(active_node)
                    for x in neighbors_active_node:
                        if not data.flag_already_in_the_list[x]:
                            if x!=-1:
                                crater_node_list.append(x)
                                data.flag_already_in_the_list[x] = 1
        #print 'f/mu is ', f_theta/mu_theta_by_mu0


    def set_elev_change_at_pts_v2(self, grid, data):
        '''
        This is an alternative method to take an existing impact properties and a known nearest node to the impact site, and alter the topography to model the impact. It assumes crater radius and depth are known, models cavity shape as a power law where n is a function of R/D, and models ejecta thickness as an exponential decay,sensitive to both ballistic range from tilting and momentum transfer in impact (after Furbish). We DO NOT yet model transition to peak ring craters, or enhanced diffusion by ejecta in the strength regime. Peak ring craters are rejected from the distribution. This version of this method is designed to remove the sheer walls around the edges of craters, and replace them with a true dipping rim.
        '''
        #Build a list of nodes to work on. Starts just with the node closest to the center.
        crater_node_list = [self.closest_node_index]
        self.elev_changes = []
        #Build an array of flags into the nodelist of the grid to note whether that node has been placed in the list this loop:
        data.flag_already_in_the_list = [0] * grid.ncells
        #Could we have the node properties as dictionaries of arrays?:
        #grid.nodes[adjusted] = zeros(grid.ncells)
        
        #Derive the exponent for the crater shape, shared betw simple & complex:
        crater_bowl_exp = self.get_crater_shape_exp()
        print 'Crater shape exponent: ', crater_bowl_exp
        
        #print 'Radius is: ', self._radius
        
        #Derive the effective angle for impact, relative to the surface, beta_eff. The direction is always controlled by the impactor. Draw new impact angles if impact geomtery is impossible.
        while 1:
            #This ugly function applies a correction to the azimuth of travel which will be the actual angle the ejecta is propelled along. This is important since if, e.g., _angle_to_horizontal ~pi/2, then the surface dip direction needs to be the dominant term.
            denominator = numpy.cos(self._angle_to_horizontal) + numpy.sin(self._surface_slope)*numpy.cos(self._surface_dip_direction - self._azimuth_of_travel)
            print 'Terms in denominator: ', numpy.cos(self._angle_to_horizontal),  numpy.sin(self._surface_slope), numpy.cos(self._surface_dip_direction - self._azimuth_of_travel)
            tan_angle_from_imp_az = numpy.sin(self._surface_slope)*numpy.sin(self._surface_dip_direction - self._azimuth_of_travel) / denominator
            if denominator > 0.:
                self._ejecta_azimuth = (self._azimuth_of_travel + numpy.arctan(tan_angle_from_imp_az) + 2.*numpy.pi)%(2.*numpy.pi)
            else:
                self._ejecta_azimuth = (self._azimuth_of_travel - numpy.arctan(tan_angle_from_imp_az) + numpy.pi)%(2.*numpy.pi)
            print 'Angle from impact azimuth: ', self._ejecta_azimuth - self._azimuth_of_travel
            beta_eff_options = [numpy.nan] * 2
            beta_eff_options[0] = self._angle_to_horizontal*numpy.cos(self._ejecta_azimuth - self._azimuth_of_travel) - self._surface_slope * numpy.cos(self._ejecta_azimuth - self._surface_dip_direction)
            beta_eff_options[1] = (numpy.pi - self._angle_to_horizontal*numpy.cos(self._ejecta_azimuth - self._azimuth_of_travel)) + self._surface_slope * numpy.cos(self._azimuth_of_travel - self._surface_dip_direction)
            beta_eff = min(beta_eff_options)
            #Need to think thru if this reversal is still required:
            #if beta_eff_options[0] > beta_eff_options[1]:
            #self._azimuth_of_travel = (self._azimuth_of_travel + numpy.pi)%(2.*numpy.pi)
            print 'Beta effective: ', beta_eff
            if beta_eff >= 0.:
                break
            else:
                print 'Refreshing the impactor angle'
                self.set_impactor_angles()
        self._impactor_angle_to_surface_normal = 0.5*numpy.pi - beta_eff
        print 'Impact, ejecta azimuths: ', self._azimuth_of_travel, self._ejecta_azimuth
        #print 'Beta effective: ', beta_eff
        #rim_height = self._cavity_volume/(3.667*numpy.pi*self._radius**2.)

        #Derive the new ejecta shape params:
        tan_repose = numpy.tan(32.*numpy.pi/180.)
        _b = 3.667*self._radius**2.*tan_repose
        _c = self._cavity_volume - 3.667*self._radius**3.*tan_repose - 0.5*self._radius**2.*tan_repose
        radius_calc = (-_b + numpy.sqrt(_b**2. - 2.*tan_repose*_c))/tan_repose
        thickness_at_rim = (self._radius - radius_calc)*tan_repose
        #...where thickness(r) = thickness_at_rim*(r/R_true)**-2.75

        while len(crater_node_list): #i.e., array is not empty
            active_node = crater_node_list.pop(0)
            active_node_x_to_center = grid.x(active_node) - self._xcoord
            active_node_y_to_center = grid.y(active_node) - self._ycoord
            active_node_r_to_center = numpy.sqrt(active_node_x_to_center**2. + active_node_y_to_center**2.)
            #Special case for if the impact location is right on a gridline:
            if not active_node_x_to_center:
                if active_node_y_to_center < 0.:
                    active_node_theta = numpy.pi
                else:
                    active_node_theta = 0.
            else:
                active_node_angle_to_yaxis = numpy.arctan(active_node_y_to_center/active_node_x_to_center)
                if active_node_x_to_center < 0.:
                    active_node_theta = 1.5 * numpy.pi - active_node_angle_to_yaxis
                else:
                    active_node_theta = 0.5 * numpy.pi - active_node_angle_to_yaxis

            #In this version, our sweep out from the center will need to account for deposition depth elevating the crater rim, i.e., we need to deposit *before* we cut the cavity. We do this by defining three domains for the node to lie in: 1. r<r_calc, i.e., below the pre-impact surface. No risk of intersecting the surface here. 2. r_calc < r; Th>z_new. this is the domain in the inward sloping rim of the crater ejecta. 3. Th<z_new and beyond. out on the ejecta proper. Note - (1) is not hard & fast rule if the surface dips. Safer is just (Th-lowering)<z_new
            #So, calc the excavation depth for all nodes, just to be on the safe side for strongly tilted geometries:
            if active_node_r_to_center <= self._radius:
                new_z = self.closest_node_elev - self._depth * (1. - (active_node_r_to_center/self._radius)**crater_bowl_exp) + thickness_at_rim
            else:
                new_z = self.closest_node_elev + thickness_at_rim + (active_node_r_to_center-self._radius)*tan_repose
            
            #print new_z

            if new_z<data.elev[active_node]: #Below the original surface
                #print 'Under surface'
                self.elev_changes.append(new_z-data.elev[active_node])
                if self._crater_type: #Complex crater, adjust for the peak
                    if active_node_r_to_center <= self._complex_peak_radius: #on the peak
                        new_z  = new_z + self._complex_peak_str_uplift * (1. - active_node_r_to_center/self._complex_peak_radius)
                #Set the ground elev:
                data.elev[active_node] = new_z
                #Add the neighboring nodes which haven't already been adjusted to the for-processing array:
                neighbors_active_node = grid.get_neighbor_list(active_node)
                for x in neighbors_active_node:
                    if not data.flag_already_in_the_list[x]:
                        if x!=-1: #Not an edge
                            crater_node_list.append(x)
                            data.flag_already_in_the_list[x] = 1

            else:
                #Calc the ejecta thickness for a perp. impact. Note it's fine if r<R_true
                local_flat_thickness = thickness_at_rim*(active_node_r_to_center/self._radius)**-2.75
                #print 'Local flat thickness: ', local_flat_thickness

                # Modulating functions are f(impactor az, dip dir, impactor angle, dip), the ejecta mass term, and mu(impactor az, dip dir, impactor angle, dip), the distance travelled term. Momentum is the dominant term over geometric effects in mu, by ~an OoM - o just model as a momentum effect.
                #Thickness = thickness_flat * f_theta / mu
                #beta_eff (==pi/2 - angle from normal) is already known.
                theta_eff = self._ejecta_azimuth - active_node_theta #This is the angle of the center-to-active-node line to the azimuth along which the ejecta is concentrated
                tan_beta = numpy.tan(self._impactor_angle_to_surface_normal*self._beta_factor)
                sin_theta_sqd = numpy.sin(theta_eff) ** 2.
                cos_theta = numpy.cos(theta_eff)
                #REMEMBER, as tan_beta gets >1, the function describing the ejecta is only valid over ever more restricted ranges of theta!! In other words,
                if tan_beta and sin_theta_sqd*tan_beta**2. > 1.:
                    thickness_at_active_node = 0.
                    #print 'Outside ejecta cone'
                else:
                    mu_theta_by_mu0 = tan_beta * cos_theta + numpy.sqrt(1. - sin_theta_sqd * tan_beta**2.)
                    f_theta = (tan_beta**2.*(cos_theta**2.-sin_theta_sqd) + 2.*tan_beta*cos_theta*numpy.sqrt(1.-tan_beta**2.*sin_theta_sqd) + 1.) / (2.*numpy.pi)
                    #print f_theta/mu_theta_by_mu0
                    #So, distn_at_angle = distn_vertical_impact*f_theta/mu_theta_by_mu0. Draw the thickness at the active node:
                    #NB-the 2pi is to correct for mismatch in the dimensions of mu and f
                    thickness_at_active_node = f_theta/mu_theta_by_mu0 * 2.*numpy.pi*local_flat_thickness
                    if thickness_at_active_node < 0.:
                        thickness_at_active_node = 0.
                #Now, are we inside or outside the rim?
                if new_z <= (data.elev[active_node] + thickness_at_active_node): #inside the rim
                    #print 'Inside the rim, on the ejecta'
                    self.elev_changes.append(new_z-data.elev[active_node])
                    data.elev[active_node] = new_z
                else: #outside the rim
                    #print 'Outside the rim'
                    self.elev_changes.append(thickness_at_active_node)
                    data.elev[active_node] = data.elev[active_node] + thickness_at_active_node
                #Add the neighbors to the list, but only if the thickness of the current ejecta layer was nontrivial. Smallest crater is 5m diam, so 1m deep - so its max ejecta thickness at the rim is 10cm! We should probably resolve down to, say, 2.5cm, which is <4 crater radii from the rim of our smallest crater.
                if thickness_at_active_node > self._minimum_ejecta_thickness:
                    neighbors_active_node = grid.get_neighbor_list(active_node)
                    for x in neighbors_active_node:
                        if not data.flag_already_in_the_list[x]:
                            if x!=-1:
                                crater_node_list.append(x)
                                data.flag_already_in_the_list[x] = 1
        #print 'f/mu is ', f_theta/mu_theta_by_mu0
        print 'Thickness at rim: ', thickness_at_rim



    def excavate_a_crater(self, grid, data):
        '''
        This method executes the most of the other methods of this crater class, and makes the geomorphic changes to a mesh associated with a single bolide impact with randomized properties. It receives parameters of the model grid, and the vector data storage class. It is the primary interface method of this class.
        '''
        self.set_cr_radius_from_shoemaker(data)
        print 'Radius: ', self._radius
        self.set_depth_from_size()
        self.set_crater_volume()
        self.set_coords(grid, data)
        self.set_impactor_angles()
        self.set_crater_mean_slope_v2(grid, data)
        if numpy.isnan(self._surface_slope):
            print 'Surface slope is not defined for this crater! Is it too big? Crater will not be drawn.'
        else:
            self.set_elev_change_at_pts_v2(grid, data)
            print 'Impactor angle to ground normal: ', self._impactor_angle_to_surface_normal
            print 'Mean mass balance/px: ', numpy.mean(self.elev_changes)
        print '*****'
        #Record the data:
        data.impact_sequence.append({'x': self._xcoord, 'y': self._ycoord, 'r': self._radius, 'volume': self._cavity_volume, 'surface_slope': self._surface_slope, 'normal_angle': self._impactor_angle_to_surface_normal, 'impact_az': self._azimuth_of_travel, 'ejecta_az': self._ejecta_azimuth, 'mass_balance': numpy.mean(self.elev_changes)})


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
    vectors = data()
    vectors.elev = [100.] * mg.ncells
    cr = impactor()

    #Update until
    for i in range(0,nt):
        print 'Crater number ', i
        cr.excavate_a_crater(mg, vectors)

    #Finalize
    elev_raster = mg.cell_vector_to_raster(vectors.elev)
    #contour(elev_raster)
    flipped_elev_raster = numpy.empty_like(elev_raster)
    for i in range(0,nr):
        flipped_elev_raster[i,:] = elev_raster[(nr-i-1),:]

    imshow(flipped_elev_raster)
    colorbar()
    show()
    vectors.viewing_raster = flipped_elev_raster
    return cr, mg, vectors

def dig_some_craters(grid, data):
    '''
    Takes an existing DTM and peppers it with craters.
    '''
    #dt = 1.
    nt = 10000

    #Setup
    cr = impactor()

    #Update until
    for i in range(0,nt):
        print 'Crater number ', i
        cr.excavate_a_crater(grid, data)
    
    #Finalize
    elev_raster = grid.cell_vector_to_raster(data.elev)
    #contour(elev_raster)
    flipped_elev_raster = numpy.empty_like(elev_raster)
    for i in range(0,grid.nrows):
        flipped_elev_raster[i,:] = elev_raster[(grid.nrows-i-1),:]
    
    profile = plot(elev_raster[600,:])
    xsec = plot(elev_raster[:,1100])
    #imshow(flipped_elev_raster)
    #colorbar()
    #show()
    data.viewing_raster = flipped_elev_raster
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
    vectors = data()
    vectors.elev = []
    for i in range(0, nr):
        vectors.elev = vectors.elev + [100.]*nc
        #vectors.elev = vectors.elev + [100.-i*0.003]*nc
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
        distances_to_vertices.append(numpy.sqrt((cr._xcoord-mg.x(x))**2. + (cr._ycoord-mg.y(x))**2.))
    cr.closest_node_index = vertices_array[numpy.argmin(distances_to_vertices)]
    cr.closest_node_elev = vectors.elev[cr.closest_node_index]

    cr._angle_to_horizontal = numpy.pi*0.5*3./3.
    cr._azimuth_of_travel = numpy.pi*1.5
    cr.set_crater_mean_slope_v2(mg, vectors)
    print 'Azimuth of travel: ', cr._azimuth_of_travel
    print 'Angle of ground: ', cr._surface_slope
    print 'Dip direction of ground: ', cr._surface_dip_direction
    cr.set_elev_change_at_pts_v2(mg, vectors)
    print 'Impact angle to ground normal: ', cr._impactor_angle_to_surface_normal

    #Finalize
    elev_raster = mg.cell_vector_to_raster(vectors.elev)
    #contour(elev_raster)
    flipped_elev_raster = numpy.empty_like(elev_raster)
    for i in range(0,nr):
        flipped_elev_raster[i,:] = elev_raster[(nr-i-1),:]
    imshow(flipped_elev_raster)
    colorbar()
    show()
    vectors.viewing_raster = flipped_elev_raster
    return cr, mg, vectors


def main():
#   cr, mg, vectors = dig_some_craters_on_fresh_surface()
    cr, mg, vectors = dig_one_crater(120, 120, 0.025, 0.5, 0.5, 1.)
#   mg_10k, vectors_10k = dig_some_craters(mg, vectors)

#    #This code builds a dictionary that contains time slices for each 10k craters hitting a surface:
#    #How many times round?
#    loops = 5 #500,000 craters!!
#    #Build the dictionary:
#    crater_time_sequ = {}
#    profile_list = []
#    xsec_list = []
#    #Initialize the starting condition:
#    cr, mg, vectors = dig_one_crater(1200, 1200, 0.0025, 0.5, 0.5, 1.)
#    #Save the starting conds:
#    crater_time_sequ[0] = copy(vectors)
#    #Run the loops
#    for i in range(0,loops):
#        mg, vectors, profile, xsec = dig_some_craters(mg, vectors)
#        crater_time_sequ[i] = copy(vectors)
#        profile_list.append(profile)
#        xsec_list.append(xsec)
#    show(profile_list)

if __name__=='__main__':
    main()