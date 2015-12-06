#! /usr/env/python
# -*- coding: utf-8 -*-

"""
.. deprecated:: 0.4
This code is no longer supported.
Use 'dig_craters.py" instead.
"""

from random import random
import math
import numpy
from collections import deque
import sys
import time
import scipy.optimize as opt
from sympy import Symbol
from sympy.solvers import solve
from sympy.utilities.lambdify import lambdify
import six

from landlab import RasterModelGrid #this is the tMesh equivalent module

#these ones only so we can run this module ad-hoc:
from pylab import plot, draw, show, contour, imshow, colorbar
from copy import deepcopy as copy

#sys.setrecursionlimit(1500)

class data(object):
    '''
    This is where all the whole-grid data lives, as arrays over the various elements of the grid.
    '''
    #Data goes here!!!
    def __init__(self, grid):
        self.elev = grid.zeros(centering='node') #some data
        self.flag_already_in_the_list = grid.zeros(centering='node')
        self.craters_over_max_radius_not_plotted = grid.zeros(centering='node')
        self.impact_sequence = []

class impactor(object):
    '''
    This class holds all parameters decribing properties of a single impact structure, and contains methods for recalculating fresh and internally consistent data describing such a impact structure.
    Built DEJH Spring 2013.
    '''
    def __init__(self, min_radius=0.005):
        self._xcoord = -999.
        self._ycoord = -999.
        self.tan_repose = numpy.tan(32.*numpy.pi/180.)
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
        self._angle_to_vertical = -999.
        self._minimum_ejecta_thickness = 0.00000001
        #NB - If this min thickness is changed, the optimization point in excavate_a_crater_optimized() will also need to be changed
        self._beta_factor = 0.5 #this is the arbitrary term that controls how "stretched out" the ejecta field is. <+0.5 prevents "outside the ejecta field" regions forming

        self.total_counted_craters = self.ivanov_prod_equ_as_Nequals(self._minimum_crater*2.)
        #Define the Ivanov fn and its derivatives, for use in Newton Raphson optimization. Note this is normalized now.
        self.ivanov_prod_fn = lambda x, N_as_fraction: self.ivanov_prod_equ(x,(N_as_fraction/self.total_counted_craters))
        self.ivanov_prod_fn_1stderiv = lambda x, N_as_fraction: self.ivanov_prod_equ_1stderiv(x, N_as_fraction)
        #self.ivanov_prod_fn_2ndderiv = lambda x: self.ivanov_prod_equ_2ndderiv(x)

        self.V = Symbol('V') #Crater cavity vol
        self.r0 = Symbol('r0') #Crater rim radius
        self.T = Symbol('T') #Crater rim ejecta thickness
        self.r = Symbol('r') #actual dist from crater center of a given pt
        self.solution_for_rim_thickness = solve(8./3.*self.T*numpy.pi*self.r0**2 + 0.33333*numpy.pi*self.T*(self.r0**2+(self.r0-self.T/self.tan_repose)**2+self.r0*(self.r0-self.T/self.tan_repose)) - self.V, self.T)
        #...gives a list of 3 sympy expressions, f(V,r), for T
        self.expression_for_local_thickness = self.T*(self.r/self.r0)**-2.75

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
            six.print_('Drew a crater above the maximum permitted size. Drawing a new crater...')
            self.set_size(data)


    def set_coords(self, grid, data):
        '''
        This method selects a random location inside the grid onto which to map an impact. It also sets variables for the closest grid node to the impact, and the elevation at that node.
        '''
        #NB - we should be allowing craters OUTSIDE the grid - as long as part of them impinges.
        #This would be relatively easy to implement - allow allocation out to the max crater we expect, then allow runs using these coords on our smaller grid. Can save comp time by checking if there will be impingement before doing the search.
        self._xcoord = random() * grid.get_grid_xdimension()
        self._ycoord = random() * grid.get_grid_ydimension()
        #print (self._xcoord, self._ycoord)
        #print grid.dx
        #print grid.number_of_node_columns
        #Snap impact to grid:
        self.closest_node_index = grid.snap_coords_to_grid(self._xcoord, self._ycoord)
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
        self._angle_to_vertical =  numpy.arcsin(random()) #gives sin distn with most values drawn nearer to 0
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
            six.print_('Drew a crater above the maximum permitted size. Drawing a new crater...')
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

        six.print_('Depth: ', self._depth)


    def set_crater_volume(self):
        '''
        This method uses known crater depth and radius and sets the volume of the excavated cavity.
        '''
        self._cavity_volume = 0.51 * numpy.pi * self._depth * self._radius**3. / (0.51*self._radius + 2.*self._depth)


    def create_lambda_fn_for_ejecta_thickness(self):
        """
        This method takes the complicated equation that relates "flat" ejecta thickness (symmetrical, with impact angle=0) to radius and cavity volume which is set in __init__(), and solves it for a given pair of impact specific parameters, V_cavity & crater radius.
        Both the cavity volume and crater radius need to have been set before this method is called.
        Method returns a lambda function for the radially symmetrical ejecta thickness distribution as a function of distance from crater center, r. i.e., call unique_expression_for_local_thickness(r) to calculate a thickness.
        Added DEJH Sept 2013.
        """
        local_solution_for_rim_thickness = self.solution_for_rim_thickness[0].subs(self.V, self._cavity_volume)
        unique_expression_for_local_thickness = self.expression_for_local_thickness.subs({self.r0:self._radius, self.T:local_solution_for_rim_thickness})
        unique_expression_for_local_thickness = lambdify(self.r, unique_expression_for_local_thickness)
        return unique_expression_for_local_thickness


    def set_crater_mean_slope_v2(self, grid, data):
        '''
        This method takes a crater of known radius, and which has already been "snapped" to the grid through snap_impact_to_grid(mygrid), and returns a spatially averaged value for the local slope of the preexisting topo beneath the cavity footprint. This version of the method works by taking four transects across the crater area every 45 degrees around its rim, calculating the slope along each, then setting the slope as the greatest, positive downwards and in the appropriate D8 direction. This function also sets the mean surface dip direction.
        In here, we start to assume a convex and structured grid, such that if pts N and W on the rim are in the grid, so is the point NW.
        This version is vectorized, and so hopefully faster.
        DEJH, Sept 2013.
        '''
        #It would be trivial to add more divisions, i.e., D16, D32, D64, etc. Code is setup for this; just add to the slope_pts lists, and adjust the set length of the arrays with divisions
        divisions = 4 #len(radial_points1), i.e., half the no. of points == the number of lines across the wheel
        half_crater_radius = 0.707 * self._radius
        #Work round clockwise from 12 o'clock
        slope_pts1 = numpy.array([[self._xcoord,self._ycoord+self._radius],
                     [self._xcoord+half_crater_radius,self._ycoord+half_crater_radius],
                     [self._xcoord+self._radius,self._ycoord],
                     [self._xcoord+half_crater_radius,self._ycoord-half_crater_radius]]) #(4,2)
        slope_pts2 = numpy.array([[self._xcoord,self._ycoord-self._radius],
                     [self._xcoord-half_crater_radius,self._ycoord-half_crater_radius],
                     [self._xcoord-self._radius,self._ycoord],
                     [self._xcoord-half_crater_radius,self._ycoord+half_crater_radius]]) #(4,2)
        inbounds_test1 = grid.is_point_on_grid(slope_pts1[:,0], slope_pts1[:,1])
        inbounds_test2 = grid.is_point_on_grid(slope_pts2[:,0], slope_pts2[:,1])
        distance_array = (inbounds_test1.astype(float) + inbounds_test2.astype(float))*self._radius
        radial_points1 = numpy.where(inbounds_test1, (slope_pts1[:,0],slope_pts1[:,1]), [[self._xcoord], [self._ycoord]])
        radial_points2 = numpy.where(inbounds_test2, (slope_pts2[:,0],slope_pts2[:,1]), [[self._xcoord], [self._ycoord]])
        radial_points1 = grid.snap_coords_to_grid(radial_points1[0,:], radial_points1[1,:])
        radial_points2 = grid.snap_coords_to_grid(radial_points2[0,:], radial_points2[1,:])
        #print 'radial pts arrays: ', radial_points1, radial_points2
        #print 'On grid?: ', inbounds_test1, inbounds_test2
        #print 'Dist array: ', distance_array
        slope_array = numpy.where(distance_array, (data.elev[radial_points1]-data.elev[radial_points2])/distance_array, numpy.nan)
        slope_array = numpy.arctan(slope_array)
        #if slope  is negative, it means the surface slopes broadly EAST
        try:
            hi_mag_slope_index = numpy.nanargmax(numpy.fabs(slope_array))
            hi_mag_slope = slope_array[hi_mag_slope_index]
        except:
            self.surface_slope = 1.e-10
            self.surface_dip_direction = self._azimuth_of_travel
            six.print_('Unable to assign crater slope by this method. Is crater of size comparable with grid?')
            six.print_('Setting slope to zero')
        else:
            #print 'Slope array: ', slope_array
            if hi_mag_slope > 0.: #i.e., dips WEST (or dead S)
                self._surface_dip_direction = (hi_mag_slope_index/float(divisions) + 1.)*numpy.pi
            elif not hi_mag_slope: #i.e., FLAT, dip dir is arbitrary; set to the travel direction of the impactor
                self._surface_dip_direction = self._azimuth_of_travel
            else: #dips EAST
                self._surface_dip_direction = numpy.pi*hi_mag_slope_index/float(divisions)
            self._surface_slope = numpy.fabs(hi_mag_slope)
            six.print_('The slope under the crater cavity footprint is: ', self._surface_slope)

    def set_elev_change_crawler(self, grid, data):
        """
        """
        #Load in the data, for speed and conciseness:
        _angle_to_vertical=self._angle_to_vertical
        _surface_slope=self._surface_slope
        _surface_dip_direction = self._surface_dip_direction
        _azimuth_of_travel = self._azimuth_of_travel
        pi = numpy.pi
        twopi = 2.*pi
        tan = numpy.tan
        cos = numpy.cos
        sin = numpy.sin
        sqrt = numpy.sqrt
        #arctan = numpy.arctan
        arccos = numpy.arccos
        #where = numpy.where
        _radius = self._radius
        tan_repose = self.tan_repose
        self.mass_balance_in_impact = 0.
        crater_vol_below_ground = 0.

        #Build a list of nodes to work on. Starts just with the node closest to the center.
        crater_node_list = deque([self.closest_node_index])
        #elev_changes = deque()
        #Build an array of flags into the nodelist of the grid to note whether that node has been placed in the list this loop:
        flag_already_in_the_list = numpy.zeros(grid.number_of_nodes)

        #Derive the exponent for the crater shape, shared betw simple & complex:
        crater_bowl_exp = self.get_crater_shape_exp()
        six.print_('Crater shape exponent: ', crater_bowl_exp)

        #Derive the effective angle for impact, relative to the surface normal, beta_eff. The direction is always controlled by the impactor. Draw new impact angles if impact geomtery is impossible.
        while 1:
            #epsilon is the angle between the surface normal and the impactor angle to vertical, projected along the line of travel.
            rake_in_surface_plane = _surface_dip_direction - _azimuth_of_travel
            six.print_('_surface_dip_direction: ', _surface_dip_direction)
            six.print_('_azimuth_of_travel: ', _azimuth_of_travel)
            six.print_('_angle_to_vertical: ', _angle_to_vertical)
            absolute_rake = numpy.fabs(rake_in_surface_plane)
            if not _surface_slope:
                epsilon = 0.
                beta_eff = _angle_to_vertical
            else:
                if absolute_rake in (0.5*pi, 1.5*pi):
                    epsilon = 0.
                else:
                    epsilon = arccos(cos(rake_in_surface_plane)*cos(_surface_slope)*sqrt(1.+(tan(rake_in_surface_plane)/cos(_surface_slope))**2.))
                #Make the necessary adjustment to the angle to vertical to reflect dipping surface:
                if absolute_rake <= 0.5*pi or absolute_rake >= 1.5*pi:
                    beta_eff = _angle_to_vertical + epsilon
                else:
                    beta_eff = _angle_to_vertical + epsilon - pi
            six.print_('Beta effective: ', beta_eff)
            #Make correction to ejecta direction needed if angle to normal is small and slope is large in the opposite direction:
            if 0. <= beta_eff <= 90.:
                _ejecta_azimuth = _azimuth_of_travel
                break
            elif beta_eff < 0.:
                #reverse the azimuth, make beta positive again
                beta_eff = -beta_eff
                _ejecta_azimuth = (_azimuth_of_travel+pi)%twopi
                break
            else:
                six.print_('Impact geometry was not possible! Refreshing the impactor angle...')
                self.set_impactor_angles()
                _azimuth_of_travel = self._azimuth_of_travel
                _angle_to_vertical = self._angle_to_vertical

        #apply correction to beta to suppress nonlinear BC problem and make ejecta patterns "look like" they actually do:
        tan_beta = tan(beta_eff*self._beta_factor)
        tan_beta_sqd = tan_beta*tan_beta
        #print 'Impact, ejecta azimuths: ', _azimuth_of_travel, _ejecta_azimuth

        unique_expression_for_local_thickness = self.create_lambda_fn_for_ejecta_thickness()
        thickness_at_rim = unique_expression_for_local_thickness(_radius)
        six.print_('thickness_at_rim: ', thickness_at_rim)

        #This code crawls out iteratively over the grid under the crater footprint, away from the centerpoint.
        while 1:
            try:
                active_node = crater_node_list.popleft() #i.e., array is not empty
            except:
                break
            else:
                pre_elev = data.elev[active_node]
                _r_to_center, _theta = grid.get_distances_of_nodes_to_point((self._xcoord,self._ycoord), get_az='angles', node_subset=active_node)

                ##We need to account for deposition depth elevating the crater rim, i.e., we need to deposit *before* we cut the cavity. We do this by defining three domains for the node to lie in: 1. r<r_calc, i.e., below the pre-impact surface. No risk of intersecting the surface here. 2. r_calc < r; Th>z_new. this is the domain in the inward sloping rim of the crater ejecta. 3. Th<z_new and beyond. out on the ejecta proper. Note - (1) is not hard & fast rule if the surface dips. Safer is just (Th-lowering)<z_new
                ##So, calc the excavation depth for all nodes, just to be on the safe side for strongly tilted geometries:
                #(_vec_new_z is the depth that would be excavated inside the cavity, including projected depths ouside the cavity.
                _new_z = self.closest_node_elev+thickness_at_rim-self._depth
                if _r_to_center<=_radius:
                    _new_z += self._depth * (_r_to_center/_radius)**crater_bowl_exp
                else:
                    _new_z += self._depth + (_r_to_center-_radius)*tan_repose
                #_nodes_below_surface = _vec_new_z<elev[footprint_nodes]
                #_nodes_above_surface = numpy.logical_not(_nodes_below_surface)
                #Check if we need to adjust for a central peak
                #Set the ground elev for below ground nodes
                if _new_z<pre_elev:
                    if self._crater_type:
                        if _r_to_center<=self._complex_peak_radius:
                            _new_z = _new_z + self._complex_peak_str_uplift * (1. - _r_to_center/self._complex_peak_radius)
                    elev = _new_z
                    depth_excavated = (pre_elev-_new_z)
                    crater_vol_below_ground += depth_excavated
                    self.mass_balance_in_impact -= depth_excavated
                    neighbors_active_node = grid.get_active_neighbors_at_node(active_node)
                    for x in neighbors_active_node:
                        if not flag_already_in_the_list[x]:
                            if x!=-1: #Not an edge
                                crater_node_list.append(x)
                                flag_already_in_the_list[x] = 1

                else: #Above ground nodes
                    _flat_thickness_above_surface = unique_expression_for_local_thickness(_r_to_center)
                    _theta_eff = _ejecta_azimuth - _theta #This is the angle of the center-to-active-node line to the azimuth along which the ejecta is concentrated
                    _sin_theta_sqd = sin(_theta_eff) ** 2.
                    _cos_theta = cos(_theta_eff)

                    #This material is not necessary as part of this footprint-based method, as we already forbid beta_factor>0.5
                    ##REMEMBER, as tan_beta gets >1, the function describing the ejecta is only valid over ever more restricted ranges of theta!! In other words,
                    #nodes_inside_ejecta = where(_vec_sin_theta_sqd*tan_beta_sqd <= 1.) #these are indices to an array of length nodes_above_surface only
                    #_vec_thickness = numpy.zeros(nodes_above_surface.size) #i.e., it's zero outside the ejecta
                    _vec_mu_theta_by_mu0 = tan_beta * _cos_theta + sqrt(1. - _sin_theta_sqd * tan_beta_sqd)
                    _vec_f_theta = (tan_beta_sqd*(_cos_theta**2.-_sin_theta_sqd) + 2.*tan_beta*_cos_theta*sqrt(1.-tan_beta_sqd*_sin_theta_sqd) + 1.) / twopi
                    #So, distn_at_angle = distn_vertical_impact*f_theta/mu_theta_by_mu0. Draw the thickness at the active node:
                    #NB-the 2pi is to correct for mismatch in the dimensions of mu and f
                    _thickness = _vec_f_theta/_vec_mu_theta_by_mu0 * twopi * _flat_thickness_above_surface
                    #Set the thicknesses <0 to 0:
                    if _thickness<0.:
                        _thickness = 0.
                    #Now, are we inside or outside the rim?
                    potential_ejecta_thickness = pre_elev + _thickness
                    if _new_z<=potential_ejecta_thickness:
                        elev = _new_z
                        self.mass_balance_in_impact += elev - pre_elev
                    else:
                        elev = potential_ejecta_thickness
                        self.mass_balance_in_impact += _thickness

                    if _thickness > self._minimum_ejecta_thickness:
                        neighbors_active_node = grid.get_active_neighbors_at_node(active_node)
                        for x in neighbors_active_node:
                            if not flag_already_in_the_list[x]:
                                if x!=-1:
                                    crater_node_list.append(x)
                                    flag_already_in_the_list[x] = 1

                #Save any data to the higher level:
                data.elev[active_node] = elev #the refs get broken somewhere...

        self.mass_balance_in_impact = self.mass_balance_in_impact/crater_vol_below_ground
        self.ejecta_azimuth = _ejecta_azimuth
        self.impactor_angle_to_surface_normal = beta_eff #note in this case this is the *effective* angle (in the direction of travel), not the actual angle to the surface.
        #Uncomment this line to see which nodes are under the footprint:
        #data.elev[footprint_nodes] = 10.

        return data.elev


    def set_elev_change_only_beneath_footprint(self, grid, data):
        '''
        This is a method to take an existing impact properties and a known nearest node to the impact site, and alter the topography to model the impact. It assumes crater radius and depth are known, models cavity shape as a power law where n is a function of R/D, and models ejecta thickness as an exponential decay,sensitive to both ballistic range from tilting and momentum transfer in impact (after Furbish). We DO NOT yet model transition to peak ring craters, or enhanced diffusion by ejecta in the strength regime. Peak ring craters are rejected from the distribution. This version of this method is designed to remove the sheer walls around the edges of craters, and replace them with a true dipping rim.
        This routine differs from other set_elev_change...() as it adjusts elevations only on nodes which fall w/i a certain footprint, rather than crawling out from the central impact point to a threshold depth or solving the whole grid.
        This version of the code does NOT correct for slope dip direction - because Furbish showed momentum almost always wins, and these impactors have a lot of momentum!
        NB - this function ASSUMES that the "beta factor" in the model is <=0.5, i.e., nonlinearities can't develop in the ejecta field, and the impact point is always within the (circular) ejecta footprint.
        Created DEJH Sept 2013.
        '''
        #Load in the data, for speed and conciseness:
        _angle_to_vertical=self._angle_to_vertical
        _surface_slope=self._surface_slope
        _surface_dip_direction = self._surface_dip_direction
        _azimuth_of_travel = self._azimuth_of_travel
        pi = numpy.pi
        twopi = 2.*pi
        tan = numpy.tan
        cos = numpy.cos
        sin = numpy.sin
        sqrt = numpy.sqrt
        #arctan = numpy.arctan
        arccos = numpy.arccos
        where = numpy.where
        _radius = self._radius
        elev = data.elev
        tan_repose = self.tan_repose

        #Derive the exponent for the crater shape, shared betw simple & complex:
        crater_bowl_exp = self.get_crater_shape_exp()

        #Derive the effective angle for impact, relative to the surface normal, beta_eff. The direction is always controlled by the impactor. Draw new impact angles if impact geomtery is impossible.
        while 1:
            #epsilon is the angle between the surface normal and the impactor angle to vertical, projected along the line of travel.
            rake_in_surface_plane = _surface_dip_direction - _azimuth_of_travel
            six.print_('_surface_dip_direction: ', _surface_dip_direction)
            six.print_('_azimuth_of_travel: ', _azimuth_of_travel)
            six.print_('_angle_to_vertical: ', _angle_to_vertical)
            absolute_rake = numpy.fabs(rake_in_surface_plane)
            if not _surface_slope:
                epsilon = 0.
                beta_eff = _angle_to_vertical
            else:
                if absolute_rake in (0.5*pi, 1.5*pi):
                    epsilon = 0.
                else:
                    epsilon = arccos(cos(rake_in_surface_plane)*cos(_surface_slope)*sqrt(1.+(tan(rake_in_surface_plane)/cos(_surface_slope))**2.))
                #Make the necessary adjustment to the angle to vertical to reflect dipping surface:
                if absolute_rake <= 0.5*pi or absolute_rake >= 1.5*pi:
                    beta_eff = _angle_to_vertical + epsilon
                else:
                    beta_eff = _angle_to_vertical + epsilon - pi
            six.print_('Beta effective: ', beta_eff)
            #Make correction to ejecta direction needed if angle to normal is small and slope is large in the opposite direction:
            if 0. <= beta_eff <= 90.:
                _ejecta_azimuth = _azimuth_of_travel
                break
            elif beta_eff < 0.:
                #reverse the azimuth, make beta positive again
                beta_eff = -beta_eff
                _ejecta_azimuth = (_azimuth_of_travel+pi)%twopi
                break
            else:
                six.print_('Impact geometry was not possible! Refreshing the impactor angle...')
                self.set_impactor_angles()
                _azimuth_of_travel = self._azimuth_of_travel
                _angle_to_vertical = self._angle_to_vertical

        #apply correction to beta to suppress nonlinear BC problem and make ejecta patterns "look like" they actually do:
        tan_beta = tan(beta_eff*self._beta_factor)
        tan_beta_sqd = tan_beta*tan_beta
        #print 'Impact, ejecta azimuths: ', _azimuth_of_travel, _ejecta_azimuth

        unique_expression_for_local_thickness = self.create_lambda_fn_for_ejecta_thickness()
        thickness_at_rim = unique_expression_for_local_thickness(_radius)

        #The center is transposed by rmax*tan(beta) in the direction of impact:
        #Might be an issue here - we use the rmax, so doesn't the *position* of the max thickness then also depend on rmax??
        #Should it actually just be r? -> no, use the MEAN: see Furbish para 27? For now, just use the max radius - this will be conservative, at least.
        #solve the ejecta radial thickness equ w. the min thickness to get the max radius:
        max_radius_ejecta_on_flat = _radius * (thickness_at_rim/self._minimum_ejecta_thickness)**0.3636
        six.print_('thickness_at_rim: ', thickness_at_rim)
        six.print_('max_radius_ejecta: ', max_radius_ejecta_on_flat)
        footprint_center_x = self._xcoord+sin(_azimuth_of_travel)*max_radius_ejecta_on_flat*tan_beta
        footprint_center_y = self._ycoord+cos(_azimuth_of_travel)*max_radius_ejecta_on_flat*tan_beta
        distances_to_footprint_center = grid.get_distances_of_nodes_to_point((footprint_center_x,footprint_center_y))

        #There's currently issues doing this "properly" - and would also extend anyway as by moving the ejecta field, it's also possible to chop off part of the actual crater.
        #Resolve by also doing a full distance map for the actual impact site, and making a union of the two distance thresholds as the footprint:
        distances_to_crater_center, azimuths_to_crater_center = grid.get_distances_of_nodes_to_point((self._xcoord,self._ycoord), get_az=1)
        footprint_nodes = numpy.logical_or(distances_to_footprint_center<=max_radius_ejecta_on_flat, distances_to_crater_center<=max_radius_ejecta_on_flat)

        #Define a shortcut identity for the elev[footprint_nodes] patch, so we don't have to keep looking it up:
        elevs_under_footprint = elev[footprint_nodes]
        #make a local copy of the elevations under the footprint, so we can do a mass balance at the end:
        old_elevs_under_footprint = numpy.copy(elevs_under_footprint)

        _vec_r_to_center = distances_to_crater_center[footprint_nodes]
        _vec_theta = azimuths_to_crater_center[footprint_nodes]
        #print distances_to_crater_center.shape, _vec_r_to_center.shape

        ##We need to account for deposition depth elevating the crater rim, i.e., we need to deposit *before* we cut the cavity. We do this by defining three domains for the node to lie in: 1. r<r_calc, i.e., below the pre-impact surface. No risk of intersecting the surface here. 2. r_calc < r; Th>z_new. this is the domain in the inward sloping rim of the crater ejecta. 3. Th<z_new and beyond. out on the ejecta proper. Note - (1) is not hard & fast rule if the surface dips. Safer is just (Th-lowering)<z_new
        ##So, calc the excavation depth for all nodes, just to be on the safe side for strongly tilted geometries:
        #(_vec_new_z is the depth that would be excavated inside the cavity, including projected depths ouside the cavity.
        _vec_new_z = numpy.empty_like(_vec_r_to_center)
        _nodes_within_crater = _vec_r_to_center<=_radius
        _nodes_outside_crater = numpy.logical_not(_nodes_within_crater)
        _vec_new_z[:] = self.closest_node_elev + thickness_at_rim - self._depth
        _vec_new_z[_nodes_within_crater] += self._depth * (_vec_r_to_center[_nodes_within_crater]/_radius)**crater_bowl_exp
        _vec_new_z[_nodes_outside_crater] += self._depth + (_vec_r_to_center[_nodes_outside_crater]-_radius)*tan_repose
        _nodes_below_surface = _vec_new_z<elev[footprint_nodes]
        _nodes_above_surface = numpy.logical_not(_nodes_below_surface)
        #Check if we need to adjust for a central peak
        if self._crater_type:
            central_peak_pts = _vec_r_to_center<=self._complex_peak_radius
            _vec_new_z[central_peak_pts] = _vec_new_z[central_peak_pts] + self._complex_peak_str_uplift * (1. - _vec_r_to_center[central_peak_pts]/self._complex_peak_radius)
        #Set the ground elev for below ground nodes
        elevs_under_footprint[_nodes_below_surface] = _vec_new_z[_nodes_below_surface]
        #From here on, our new arrays will only be as long as nodes_above_surface
        _vec_flat_thickness_above_surface = unique_expression_for_local_thickness(_vec_r_to_center[_nodes_above_surface])
        _vec_theta_eff = _ejecta_azimuth - _vec_theta[_nodes_above_surface] #This is the angle of the center-to-active-node line to the azimuth along which the ejecta is concentrated
        _vec_sin_theta_sqd = sin(_vec_theta_eff) ** 2.
        _vec_cos_theta = cos(_vec_theta_eff)

        #This material is not necessary as part of this footprint-based method, as we already forbid beta_factor>0.5
        ##REMEMBER, as tan_beta gets >1, the function describing the ejecta is only valid over ever more restricted ranges of theta!! In other words,
        #nodes_inside_ejecta = where(_vec_sin_theta_sqd*tan_beta_sqd <= 1.) #these are indices to an array of length nodes_above_surface only
        #_vec_thickness = numpy.zeros(nodes_above_surface.size) #i.e., it's zero outside the ejecta
        _vec_mu_theta_by_mu0 = tan_beta * _vec_cos_theta + sqrt(1. - _vec_sin_theta_sqd * tan_beta_sqd)
        _vec_f_theta = (tan_beta_sqd*(_vec_cos_theta**2.-_vec_sin_theta_sqd) + 2.*tan_beta*_vec_cos_theta*sqrt(1.-tan_beta_sqd*_vec_sin_theta_sqd) + 1.) / twopi
        #So, distn_at_angle = distn_vertical_impact*f_theta/mu_theta_by_mu0. Draw the thickness at the active node:
        #NB-the 2pi is to correct for mismatch in the dimensions of mu and f
        _vec_thickness = _vec_f_theta/_vec_mu_theta_by_mu0 * twopi * _vec_flat_thickness_above_surface
        #Set the thicknesses <0 to 0:
        _vec_thickness = where(_vec_thickness>=0.,_vec_thickness, 0.)
        #Now, are we inside or outside the rim?
        elevs_under_footprint[_nodes_above_surface] = where(_vec_new_z[_nodes_above_surface]<=(elevs_under_footprint[_nodes_above_surface]+_vec_thickness),_vec_new_z[_nodes_above_surface], elevs_under_footprint[_nodes_above_surface]+_vec_thickness)

        #Save any data to the higher level:
        data.elev[footprint_nodes] = elevs_under_footprint #the refs get broken somewhere...
        elev_diff = elevs_under_footprint-old_elevs_under_footprint
        self.mass_balance_in_impact = numpy.sum(elev_diff)/-numpy.sum(elev_diff[elev_diff<0.]) #positive is mass gain, negative is mass loss. This is currently a mass fraction, given relative to volume (h*px#) excavated from below the original surface.
        #whole_grid = numpy.zeros_like(data.elev)
        #whole_grid[footprint_nodes]=elev_diff
        #whole_grid = grid.node_vector_to_raster(whole_grid, flip_vertically=True)
        #imshow(whole_grid), colorbar()
        self.ejecta_azimuth = _ejecta_azimuth
        self.impactor_angle_to_surface_normal = beta_eff #note in this case this is the *effective* angle (in the direction of travel), not the actual angle to the surface.
        six.print_('Vol of crater cavity: ', self._cavity_volume)
        #print 'Vol below ground: ', -numpy.sum(elev_diff[elev_diff<0.])
        #print 'Vol above ground: ', numpy.sum(elev_diff[elev_diff>0.])
        #print 'Total mass balance: ', numpy.sum(elev_diff)

        #Uncomment this line to see which nodes are under the footprint:
        #data.elev[footprint_nodes] = 10.

        return data.elev


    def excavate_a_crater(self, grid, data, **kwds):
        '''
            This method executes the most of the other methods of this crater class, and makes the geomorphic changes to a mesh associated with a single bolide impact with randomized properties. It receives parameters of the model grid, and the vector data storage class. It is the primary interface method of this class.
            This method is optimized to only calculate the elevation changes for an impact within its ejecta footprint.
            A fixed crater size can be specified with the input variable "forced_radius" (in km), and a fixed impact angle with "forced_angle" (in degrees from vertical - impact azimuth will always be assumed as travel eastwards). Position can be specified with forced_pos, which takes an array-like object with two entries, which are the x and y coordinate in relative position on the grid (e.g., [0.5, 0.5]).
        '''
        try:
            self._radius = kwds['forced_radius']
        except:
            six.print_('Randomly generating impact radius...')
            self.set_cr_radius_from_shoemaker(data)
        six.print_('Radius: ', self._radius)
        self.set_depth_from_size()
        self.set_crater_volume()
        try:
            self._xcoord = kwds['forced_pos'][0]*grid.get_grid_xdimension()
        except:
            six.print_('Randomly generating impact site...')
            self.set_coords(grid, data)
        else:
            try:
                self._ycoord = kwds['forced_pos'][1]*grid.get_grid_ydimension()
                six.print_(self._xcoord, self._ycoord)
                self.closest_node_index = grid.snap_coords_to_grid(self._xcoord, self._ycoord)
                self.closest_node_elev = data.elev[self.closest_node_index]
            except:
                six.print_('Could not set specified position. Was a 2 item iterable provided?')
        try:
            self._angle_to_vertical = kwds['forced_angle']*numpy.pi/180.
        except:
            six.print_('Randomly generating impactor angle...')
            self.set_impactor_angles()
        else:
            #assert 0. <= self._angle_to_vertical <= 90.
            self._azimuth_of_travel = 0.5*numpy.pi
        try:
            self._minimum_crater = kwds['minimum_radius']
        except:
            pass

        self.set_crater_mean_slope_v2(grid, data)
        if numpy.isnan(self._surface_slope):
            six.print_('Surface slope is not defined for this crater! Is it too big? Crater will not be drawn.')
        else:
            self.set_elev_change_only_beneath_footprint(grid, data)
            #self.set_elev_change_crawler(grid, data)
        six.print_('Impactor angle to ground normal: ', self.impactor_angle_to_surface_normal)
        six.print_('Mass balance in impact: ', self.mass_balance_in_impact)
        six.print_('*****')
        #Record the data:
        #Is this making copies, or just by reference? Check output.
        data.impact_sequence.append({'x': self._xcoord, 'y': self._ycoord, 'r': self._radius, 'volume': self._cavity_volume, 'surface_slope': self._surface_slope, 'normal_angle': self.impactor_angle_to_surface_normal, 'impact_az': self._azimuth_of_travel, 'ejecta_az': self.ejecta_azimuth, 'mass_balance': self.mass_balance_in_impact})


#The functions in this segment give control over the execution of this module. Adjust the params inside the functions to get different effects, and the final line of the file to determine whether you get one crater, or lots of craters.
#These should really be in a separate driver file.

def dig_some_craters(use_existing_grid=0, grid_dimension_in=1000, dx_in=0.0025, dy_in=0.0025, n_craters=1, surface_slope=0., **kwds):
    '''
    Ad hoc driver code to make this file run as a standalone.
    If a surface_slope is specified, it should be in degrees, and the resulting surface will dip west.
    If use_existing_grid is set, it should (for now) be a tuple containing (grid, data).
    If force_crater_properties is specified, it should be keywords for excavate_a_crater(), comprising as many as desired of: forced_radius=value, forced_angle=value, forced_pos=(rel_x,rel_y).
    '''
    #User-defined params:
    nr = grid_dimension_in
    nc = grid_dimension_in
    dx = dx_in
    dy = dy_in
    #dt = 1.
    nt = n_craters

    #Setup
    if not use_existing_grid:
        mg = RasterModelGrid((nr, nc), spacing=(dy, dx))
        vectors = data(mg)
        if not surface_slope:
            vectors.elev[:] = 1.
        else:
            vertical_rise_in_one_node_step = dx*numpy.tan(surface_slope*numpy.pi/180.)
            for i in range(nr):
                vectors.elev = numpy.tile(numpy.array(range(nr))*vertical_rise_in_one_node_step, nr)
        #add some noise:
        vectors.elev += numpy.random.uniform(0.,0.00001,vectors.elev.shape)
    else:
        try:
            mg = use_existing_grid[0]
            vectors = use_existing_grid[1]
        except:
            six.print_('Could not set variables for existing grid!')

    if not 'cr' in locals():
        cr = impactor()

    #Update until
    for i in xrange(0,nt):
        six.print_('Crater number ', i)
        cr.excavate_a_crater(mg, vectors, **kwds)

    #Finalize
    elev_raster = mg.node_vector_to_raster(vectors.elev, flip_vertically=True)
    #contour(elev_raster)

    #imshow(elev_raster)
    #colorbar()
    #show()
    vectors.viewing_raster = copy(elev_raster)
    return cr, mg, vectors

def dig_one_crater_then_degrade(loops=1, step=500):
    #Build the dictionary:
    crater_time_sequ = {}
    #Initialize the starting condition:
    cr, mg, vectors = dig_some_craters(grid_dimension_in=1000, dx_in=0.002, dy_in=0.002, n_craters=1, forced_radius = 0.5, forced_angle=0., forced_pos=(0.5,0.5))
    #Save the starting conds:
    crater_time_sequ[0] = copy(vectors.impact_sequence)
    numpy.savetxt('saved_elevs0', vectors.viewing_raster)
    #Run the loops
    for i in xrange(0,loops):
        cr, mg, vectors = dig_some_craters(use_existing_grid=(mg,vectors), n_craters=step)
        crater_time_sequ[i+1] = copy(vectors.impact_sequence)
        numpy.savetxt('saved_elevs'+str(i+1), vectors.viewing_raster)
    #end_time = time.time()
    #six.print_('Elapsed time was %g seconds' % (end_time - start_time))
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

def step_reduce_size(mg_in, vectors_in, loops=(25, 25), interval=10000,
                     min_radius_in=(0.05, 0.005)):
    crater_time_sequ_1st = {}
    crater_time_sequ_2nd = {}
    profile_list = []
    xsec_list = []
    for i in xrange(0,loops[0]):
        mg_in, vectors_in, profile, xsec = dig_some_craters(
            mg_in, vectors_in, nt_in=interval, min_radius=min_radius_in[0])
        crater_time_sequ_1st[i] = copy(vectors_in)
    for i in xrange(0, loops[1]):
        mg_in, vectors_in, profile, xsec = dig_some_craters(
            mg_in, vectors_in, nt_in=interval, min_radius=min_radius_in[1])
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
            six.print_('Input type not recognised!')


def mass_balance_tests():
    """
    This script applies tests to try to pin down the cause of the weird mass balance issues existing in this module.
    """



if __name__=='__main__':
    dig_one_crater_then_degrade(loops=4, step=20)
