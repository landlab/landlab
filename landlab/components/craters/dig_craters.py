#! /usr/env/python
"""
This component supercedes an older version, "craters.py".
This component excavates impact craters across a surface. The properties of the
craters broadly follow those observed for the Moon. Craters obey realistic size-
frequency distributions; at the moment this is forced to follow the Shoemaker
scaling, N = k*D^-2.9.
Importantly, this module incorporates a dependence of  ejecta distribution on
impact momentum.
At the moment, this component does not produce complex craters, largely because
it is optimized to run for craters in the meter-km scale. Craters are dug
perpendicular to the local surface, unless this would steepen local slopes so
much as to severely violate mass balance for the crater,or if slopes would
exceed 90 degrees. In other respects, strength effects are not well
incorporated.
As yet, this module does not have a true time interface - a single call of
excavate_a_crater() will dig a single crater, and adjust its size appropriately
for the Shoemaker distribution. It can't yet adjust itself to imposed timesteps
other than this natural timestep.
At the moment, this component requires you to run on a square regular grid.
"""

from random import random
import numpy
from sympy import Symbol
from sympy.solvers import solve
from sympy.utilities.lambdify import lambdify
#import pylab
from time import sleep
import pandas as pd
import six
from six.moves import zip

from landlab import ModelParameterDictionary

class impactor(object):
    '''
    This class holds all parameters decribing properties of a single impact
    structure, and contains methods for recalculating fresh and internally
    consistent data describing such a impact structure.
    Built DEJH Winter 2013, after an earlier version, Spring 2013.
    '''

    def __init__(self, grid, input_stream):
        self.grid = grid
        inputs = ModelParameterDictionary(input_stream)

        #test the necessary fields are all already present:
        try:
            self.elev = grid.at_node['topographic__elevation']
        except:
            six.print_('elevations not found in grid!')

        #User sets:
        self._minimum_crater = inputs.read_float('min_radius') #km. This is the smallest modelled crater radius. 10m diameter is the strict cutoff for known cr distns
        self._minimum_ejecta_thickness = inputs.read_float('min_ejecta_thickness') #0.00000001
        self._record_impacts_flag = inputs.read_int('record_impacts') #Not actually used; they're recorded by default
        #These parameters are optional
        try:
            self._radius = inputs.read_float('forced_radius')
        except:
            six.print_('Impact radii will be randomly generated.')
            self.radius_auto_flag = 1
            self.set_cr_radius_from_shoemaker()
        else:
            self.radius_auto_flag = 0
        try:
            self._xcoord = 0.5*grid.dx + inputs.read_float('x_position')*(grid.get_grid_xdimension()-grid.dx)
            self._ycoord = 0.5*grid.dy + inputs.read_float('y_position')*(grid.get_grid_ydimension()-grid.dy)
        except:
            six.print_('Impact sites will be randomly generated.')
            self.position_auto_flag = 1
            self.set_coords()
        else:
            self.position_auto_flag = 0
            self.closest_node_index = grid.find_nearest_node((self._xcoord, self._ycoord))
            self.closest_node_elev = self.elev[self.closest_node_index]
        try:
            self._angle_to_vertical = inputs.read_float('forced_angle')*numpy.pi/180.
            assert self._angle_to_vertical <= 0.5*numpy.pi
        except:
            six.print_('Impactor angles will be randomly generated.')
            self.angle_auto_flag = 1
            self.set_impactor_angles()
        else:
            self.angle_auto_flag = 0
            self._azimuth_of_travel = 0.5*numpy.pi

        self.cheater_flag = 0

        #The user has no control over the following inbuilt parameters describing crater scaling:
        self.tan_repose = numpy.tan(32.*numpy.pi/180.)
        self._beta_factor = 0.5 #this is the arbitrary term that controls how "stretched out" the ejecta field is. <+0.5 prevents "outside the ejecta field" regions forming
        self._simple_radius_depth_ratio_Pike = 2.55 #Pike thru Holsapple, tho Garvin (2011) gives 2.0 as "typically cited"

        self.V = Symbol('V') #Crater cavity vol
        self.r0 = Symbol('r0') #Crater rim radius
        self.T = Symbol('T') #Crater rim ejecta thickness
        self.r = Symbol('r') #actual dist from crater center of a given pt
        self.solution_for_rim_thickness = solve(8./3.*self.T*numpy.pi*self.r0**2 + 0.33333*numpy.pi*self.T*(self.r0**2+(self.r0-self.T/self.tan_repose)**2+self.r0*(self.r0-self.T/self.tan_repose)) - self.V, self.T)
        #...gives a list of 3 sympy expressions, f(V,r), for T
        self.expression_for_local_thickness = self.T*(self.r/self.r0)**-2.75
        self.loop_dict = {} #this holds information on how looped BCs execute - see set_slope()

        #perform a BC condition check:
        if not numpy.all(numpy.equal(grid.status_at_node[numpy.nonzero(grid.status_at_node)], 3)):
            self.looped_BCs = False
            six.print_('*****-----*****-----*****')
            six.print_('This module is designed to run with looped boundary conditions.')
            six.print_('Proceed at your own risk!')
            six.print_('Significant mass leaks are likely to develop.')
            six.print_('*****-----*****-----*****')
            sleep(3.)
        else:
            self.looped_BCs = True

        self.grid = grid
        self.twod_node_store = numpy.empty((self.grid.number_of_node_rows,self.grid.number_of_node_columns), dtype=int)
        self.dummy_1 = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype=float)
        self.dummy_2 = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype=float)
        self.dummy_3 = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype=float)
        self.dummy_4 = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype=float)
        self.dummy_5 = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype=float)
        self.dummy_6 = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype=float)
        self.dummy_int = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype=int)
        self.dummy_bool = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype=bool)
        self.double_dummy = numpy.empty((2,self.grid.number_of_node_rows*self.grid.number_of_node_columns), dtype=float)
        #this is a dummy parameter to store intermediate node id refs in, notably in create_square_footsix.print_(). It is designed to stop memory management blowing up
        self.crater_footprint_max = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype = int)
        self.footprint_max = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype = int)
        self._vec_r_to_center = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype = float)
        self._vec_theta = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype = float)
        self.slope_offsets_rel_to_center = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype = float)
        self._vec_new_z = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype = float)
        self.z_difference = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype = float)
        self.elevs_ground_less_new_z = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype = float)
        self._vec_flat_thickness_above_surface = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype = float)
        self._vec_mu_theta_by_mu0 = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype = float)
        self._vec_f_theta = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype = float)
        self._vec_thickness = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype = float)
        self.final_elev_diffs = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype = float)
        self.elev_diff = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype = float)
        self.pre_impact_elev = self.elev.copy()

        #full length arrays for memory management

        #Build the permanent maps for the distances and azimuths between nodes:
        ###The resulting matrix is too big in practical cases (GB of memory use). Need to take another apporach - see set_elevation_change_only_beneath_footsix.print_()
        #print 'Building the distances map for the grid... may take some time...'
        #self.all_node_distances_map, self.all_node_azimuths_map = grid.build_all_node_distances_azimuths_maps()
        #print '...Done.'

        self.impact_property_dict = {}

        six.print_('Craters component setup complete!')

    def draw_new_parameters(self):
        '''
        This method updates the core properties of radius, position and angle to
        vertical for a new impact crater.
        '''
        #Need to ensure the elevs and grid have been updated before this call...
        if self.radius_auto_flag == 1:
            self.set_cr_radius_from_shoemaker()
        if self.position_auto_flag == 1:
            self.set_coords()
        if self.angle_auto_flag == 1:
            self.set_impactor_angles()

    def get_crater_shape_exp(self):
        '''
        This method assumes the max depth and radius of a crater are known.
        It provides n for a power law of form d = D*(r/R)**n, where D and R are
        the known values, by assuming the outer edges of the crater sit at angle
        of repose. This gives very sensible answers; n~2 for big, complex
        craters (Garvin et al, 2000, p.333: "There is a strong tendency for
        craters to become more paraboloidal with increasing diameter,
        independent of location.") and n~1.3 for ~2km simple craters (Garvin,
        following Croft has ~1.18).
        '''
        return 0.51 * self._radius / self._depth

    #Holsapple & Housen et al 1983 note that in the strength regime, ejecta
    #distribution will NOT scale independently of crater size (it does in the
    #gravity regime) - ejecta will be proportionally FURTHER FROM the rim as craters
    #get smaller (e.g., Housen et al 83, fig 9).
    #...this, of course, makes small craters disproportionately effective diffusers
    #- must find a way of quantifying this!
    #Moore et al 74 cited in Housen et al 83 gives Rce/R as const ~2.35 +0.56 -0.45
    #std err for Lunar craters 0.3-100km, where Rce is the continuous ejecta limit.
    #Below this, seems to level off abruptly - i.e., Rce ~ R**(<1) (and const will
    #change also to maintain connection at 0.3km).
    #A key point - geomorphic gardening is not totally comparable to regolith
    #gardening. We care about topology - the final crater. Regolith fracture and
    #turnover is more likely to depend on the crater transient depth.

    def set_cr_radius_from_shoemaker(self):
        '''
        This method takes a random number between 0 and 1, and returns a crater
        radius based on a py distn N = kD^-2.9, following Shoemaker et al., 1970.
        '''
        self._radius = self._minimum_crater*(random())**-0.345

    def set_coords(self):
        '''
        This method selects a random location inside the grid onto which to map
        an impact. It also sets variables for the closest grid node to the
        impact, and the elevation at that node.
        '''
        #NB - we should be allowing craters OUTSIDE the grid - as long as part of them impinges.
        #This would be relatively easy to implement - allow allocation out to the max crater we expect, then allow runs using these coords on our smaller grid. Can save comp time by checking if there will be impingement before doing the search.
        grid = self.grid
        self._xcoord = 0.5*grid.dx + random() * (grid.get_grid_xdimension() - grid.dx)
        self._ycoord = 0.5*grid.dy + random() * (grid.get_grid_ydimension() - grid.dy)
        #Snap impact to grid:
        self.closest_node_index = grid.find_nearest_node((self._xcoord, self._ycoord))
        self.closest_node_elev = self.elev[self.closest_node_index]
        #NB - snapping to the grid may be quite computationally demanding in a Voronoi.

    def check_coords_and_angles_for_grazing(self):
        '''
        This method migrates the coords of a given impactor if its normal angle
        means it would clip other topo before striking home.
        It assumes both set_impactor_angles() and set_coords() have already both
        been called.
        '''
        ###At the moment, this method is introducing subtle bias: low angle impactors
        ###are more likely to get rejected, as they disappear off the grid edges
        ###more readily.
        ###i.e., need to introduce looping.
        sin_az = numpy.sin(self._azimuth_of_travel)
        cos_az = numpy.cos(self._azimuth_of_travel)
        alpha = self._azimuth_of_travel - numpy.pi
        cos_alpha = numpy.cos(alpha)
        sin_alpha = numpy.sin(alpha)
        grid = self.grid
        gridx = grid.get_grid_xdimension()
        gridy = grid.get_grid_ydimension()
        dx = grid.dx
        dy - grid.dy
        x = self._xcoord
        y = self._ycoord
        height = self.closest_node_elev
        not_done = True
        impact_at_edge = False
        counter = 0
        while not_done and counter<30:
            six.print_('...')
            if sin_alpha<0.: #travelling east
                line_horiz = -(x-0.5*dx) #...so line extends TO THE WEST
            else: #travelling W
                line_horiz = gridx - x - 0.5*dx
            if cos_alpha<0.: #travelling N
                line_vert = -(y-0.5*dy)
            else: #travelling S
                line_vert = gridy - y - 0.5*dy
            #How many divisions?
            hyp_line_vert = numpy.fabs(line_vert/cos_alpha)
            hyp_line_horiz = numpy.fabs(line_horiz/sin_alpha)
            #print hyp_line_vert, hyp_line_horiz
            num_divisions = int(min(hyp_line_vert,hyp_line_horiz)//dx) #Spacing set to dx, ALONG LINE (so spacing always <dx in x,y)
            if num_divisions > 0:
                self.dummy_1[:(num_divisions)] = xrange(num_divisions)
                numpy.add(self.dummy_1[:(num_divisions)], 1., out=self.dummy_2[:(num_divisions)])
                numpy.multiply(self.dummy_2[:(num_divisions)], dx, out=self.dummy_1[:(num_divisions)])
                #line_points = (numpy.arange(num_divisions)+1.)*dx #[dx,2dx,3dx...]; arbitrary reduction in length at end to keep clear of dodgy grid edge
            else:
                self.dummy_1[:1] = numpy.array((0.,))
                num_divisions = 1
                impact_at_edge = True
                six.print_('AD HOC FIX')
            numpy.multiply(self.dummy_1[:(num_divisions)], sin_alpha, out=self.dummy_2[:(num_divisions)])
            numpy.add(self.dummy_2[:(num_divisions)], x, out=self.dummy_5[:(num_divisions)])
            #line_xcoords = x + sin_alpha*line_points
            numpy.multiply(self.dummy_1[:(num_divisions)], cos_alpha, out=self.dummy_2[:(num_divisions)])
            numpy.add(self.dummy_2[:(num_divisions)], y, out=self.dummy_6[:(num_divisions)])
            #line_ycoords = y + cos_alpha*line_points #negative dimensions should sort themselves out
            assert numpy.all(self.dummy_5[:(num_divisions)]>=0) or impact_at_edge
            self.dummy_2[:(num_divisions)] = grid.find_nearest_node((self.dummy_5[:(num_divisions)], self.dummy_6[:(num_divisions)]))
            #snapped_pts_along_line = grid.snap_coords_to_grid(line_xcoords,line_ycoords) #closest to to furthest from impact
            #print snapped_pts_along_line
            try:
                numpy.divide(self.dummy_1[:(num_divisions)],numpy.tan(self._angle_to_vertical), out=self.dummy_3[:(num_divisions)])
                numpy.add(self.dummy_3[:(num_divisions)], height, out=self.dummy_4[:(num_divisions)])
                #impactor_elevs_along_line = line_points/numpy.tan(self._angle_to_vertical) + height
            except ZeroDivisionError: #vertical impactor, no need to migrate anything
                not_done = False
            else:
                try: #start at edge, work in
                    numpy.less_equal(self.dummy_4[(num_divisions-1)::-1], self.elev[self.dummy_2[(num_divisions-1)::-1].astype(int, copy=False)], out=self.dummy_int[:(num_divisions)])
                    #points_under_surface_reversed = impactor_elevs_along_line[::-1]<=self.elev[snapped_pts_along_line[::-1]]
                    reversed_index = numpy.argmax(self.dummy_int[:(num_divisions)])
                    #reversed_index = numpy.argmax(points_under_surface_reversed)
                    intersect_pt_node_status = self.grid.status_at_node[self.dummy_2[(num_divisions-1)::-1][reversed_index]]
                    #intersect_pt_node_status = self.grid.status_at_node[snapped_pts_along_line[::-1][reversed_index]]
                except IndexError: #only one item in array
                    six.print_(len(self.dummy_4[:(num_divisions)]))
                    assert len(self.dummy_4[:(num_divisions)]) == 1
                    numpy.less_equal(self.dummy_4[:(num_divisions)], self.elev[self.dummy_2[:(num_divisions)].astype(int, copy=False)], out=self.dummy_int[:(num_divisions)])
                    #points_under_surface_reversed = impactor_elevs_along_line<=self.elev[snapped_pts_along_line]
                    reversed_index = numpy.argmax(self.dummy_int[:(num_divisions)])
                    #reversed_index = numpy.argmax(points_under_surface_reversed)
                    intersect_pt_node_status = self.grid.status_at_node[self.dummy_2[:(num_divisions)]]
                    #intersect_pt_node_status = self.grid.status_at_node[snapped_pts_along_line]
                if numpy.any(self.dummy_int[:(num_divisions)]):
                #if numpy.any(points_under_surface_reversed):
                    #print 'points under surface...', numpy.sum(points_under_surface_reversed)
                    #reverse the array order as impactor comes from far and approaches the impact point
                    #If reversed_index is 0, the impact NEVER makes it above ground on the grid, and needs to be discarded.
                    #BUT, if boundaries are looped, we need to keep going!!
                    if (reversed_index and intersect_pt_node_status != 3) and not impact_at_edge: #not the final point (i.e., closest to edge), and not a looped boundary node
                        six.print_('migrating...')
                        index_of_impact = num_divisions - reversed_index
                        self._xcoord = self.dummy_5[:(num_divisions)][index_of_impact]
                        self._ycoord = self.dummy_6[:(num_divisions)][index_of_impact]
                        #self._xcoord = line_xcoords[index_of_impact]
                        #self._ycoord = line_ycoords[index_of_impact]
                        self.closest_node_index = self.dummy_2[:(num_divisions)][index_of_impact]
                        #self.closest_node_index = snapped_pts_along_line[index_of_impact]
                        self.closest_node_elev = self.elev[self.closest_node_index]
                        not_done = False
                    else:
                        if self.looped_BCs:
                            six.print_('Need to loop the impact trajectory...')
                            shorter_edge = numpy.argmin(numpy.array([hyp_line_vert,hyp_line_horiz]))
                            if shorter_edge: #the horizontal line is shorter, line makes contact with a side
                                if sin_alpha>0.: #the RHS
                                    horiz_coord = 0.50001*dx #note the switch to the other side - here's the loop! Make sure we're not right on the grid edge, to avoid snapping problems
                                    vert_coord = y + (gridx-x-0.5*dx)*cos_alpha/sin_alpha
                                    height += numpy.sqrt((gridx-x-0.5*dx)**2. + (vert_coord-y)**2.)/numpy.tan(self._angle_to_vertical)
                                else: #the LHS
                                    horiz_coord = gridx - 0.50001*dx
                                    vert_coord = y - (x-0.5*dx)*cos_alpha/sin_alpha
                                    height += numpy.sqrt((x-0.5*dx)**2. + (vert_coord-y)**2.)/numpy.tan(self._angle_to_vertical)
                            else: #vertical is shorter
                                if cos_alpha>0.: #travelling S, line extends N
                                    vert_coord = 0.50001*dy #ditto on switch
                                    horiz_coord = x + (gridy-y-0.5*dy)*sin_alpha/cos_alpha
                                    height += numpy.sqrt((gridy-y-0.5*dy)**2. + (horiz_coord-x)**2.)/numpy.tan(self._angle_to_vertical)
                                else:
                                    vert_coord = gridy - 0.50001*dy
                                    horiz_coord = x - (y-0.5*dy)*sin_alpha/cos_alpha
                                    height += numpy.sqrt((y-0.5*dy)**2. + (horiz_coord-x)**2.)/numpy.tan(self._angle_to_vertical)
                            #...& no break!
                            x = horiz_coord
                            y = vert_coord
                            six.print_(x/gridx,y/gridy)
                            counter += 1
                            if counter >= 30:
                                #This is a duff crater; kill it by setting r=0
                                self._radius = 0.000001
                                six.print_('Aborted this crater. Its trajectory was not physically plausible!')
                        elif self.position_auto_flag == 1:
                            self.draw_new_parameters()
                            self.check_coords_and_angles_for_grazing()
                            not_done = False
                        else:
                            #This is a duff crater; kill it by setting r=0
                            self._radius = 0.000001
                            six.print_('Aborted this crater. Its trajectory was not physically plausible!')
                            not_done = False
                else:
                    not_done = False



    def set_impactor_angles(self):
        '''
        This method sets the angle of impact, assuming the only effect is
        rotation of the planet under the impactor bombardment (i.e., if the
        target looks like a circle to the oncoming impactor, there's more limb
        area there to hit). As long as target is rotating relative to the sun,
        other (directional) effects should cancel. i.e., it draws from a sine
        distribution.
        Angle is given to vertical.
        Also sets a random azimuth.
        '''
        self._angle_to_vertical =  abs(numpy.arcsin(random())) #gives sin distn with most values drawn nearer to 0
        self._azimuth_of_travel = random() * 2. * numpy.pi #equal chance of any azimuth
        #Shoemaker 1983 gives a speculative equn for the suppression of D by increasing impact angle (his equ 3), but effect is minor, and it's probably not worth the trouble.
        #Shoemaker 1962 (in book) apparently states low angle impacts are very rare.

    def set_depth_from_size(self):
        '''
        This method sets the maximum depth at the center for a crater of known
        (i.e., already set) radius.
        '''
        #Let's ignore the strength transition at the lowest size scales for now.
        self._depth = self._radius / self._simple_radius_depth_ratio_Pike

    def set_crater_volume(self):
        '''
        This method uses known crater depth and radius and sets the volume of
        the excavated cavity.
        Note this is is the cavity volume, not the subsurface excavation vol.
        The true excavated volume is used to set the ejecta volumes in the
        ..._BAND_AID methods.
        '''
        radius = self._radius
        depth = self._depth
        self._cavity_volume = 0.51 * numpy.pi * depth * radius*radius*radius / (0.51*radius + 2.*depth)

    def create_lambda_fn_for_ejecta_thickness(self):
        """
        This method takes the complicated equation that relates "flat" ejecta
        thickness (symmetrical, with impact angle=0) to radius and cavity volume
        which is set in __init__(), and solves it for a given pair of impact
        specific parameters, V_cavity & crater radius.
        Both the cavity volume and crater radius need to have been set before
        this method is called.
        Method returns a lambda function for the radially symmetrical ejecta
        thickness distribution as a function of distance from crater center, r.
        i.e., call unique_expression_for_local_thickness(r) to calculate a
        thickness.
        Added DEJH Sept 2013.
        """
        local_solution_for_rim_thickness = self.solution_for_rim_thickness[0].subs(self.V, self._cavity_volume)
        unique_expression_for_local_thickness = self.expression_for_local_thickness.subs({self.r0:self._radius, self.T:local_solution_for_rim_thickness})
        unique_expression_for_local_thickness = lambdify(self.r, unique_expression_for_local_thickness)
        return unique_expression_for_local_thickness

    def create_lambda_fn_for_ejecta_thickness_BAND_AID(self, excavated_vol):
        """
        This method takes the complicated equation that relates "flat" ejecta
        thickness (symmetrical, with impact angle=0) to radius and cavity volume
        which is set in __init__(), and solves it for a given pair of impact
        specific parameters, V_cavity & crater radius.
        Both the cavity volume and crater radius need to have been set before
        this method is called.
        Method returns a lambda function for the radially symmetrical ejecta
        thickness distribution as a function of distance from crater center, r.
        i.e., call unique_expression_for_local_thickness(r) to calculate a
        thickness.
        This method is exactly equivalent to the non _BAND_AID equivalent,
        except it takes the volume to use as a parameter rather than reading it
        from the impactor object.
        Added DEJH Sept 2013.
        """
        local_solution_for_rim_thickness = self.solution_for_rim_thickness[0].subs(self.V, excavated_vol)
        unique_expression_for_local_thickness = self.expression_for_local_thickness.subs({self.r0:self._radius, self.T:local_solution_for_rim_thickness})
        unique_expression_for_local_thickness = lambdify(self.r, unique_expression_for_local_thickness)
        return unique_expression_for_local_thickness


    def set_crater_mean_slope_v2(self):
        '''
        This method takes a crater of known radius, and which has already been
        "snapped" to the grid through snap_impact_to_grid(mygrid), and returns a
        spatially averaged value for the local slope of the preexisting topo
        beneath the cavity footprint. This version of the method works by taking
        four transects across the crater area every 45 degrees around its rim,
        calculating the slope along each, then setting the slope as the greatest,
        positive downwards and in the appropriate D8 direction. This function
        also sets the mean surface dip direction.
        In here, we start to assume a convex and structured grid, such that if
        pts N and W on the rim are in the grid, so is the point NW.
        This version is vectorized, and so hopefully faster.
        This method is largely superceded by set_crater_mean_slope_v3(), which
        uses an GIS-style routine to set slopes. However, it may still be
        preferable for grid-marginal craters.
        DEJH, Sept 2013.
        '''
        grid = self.grid
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
        radial_points1 = grid.find_nearest_node((radial_points1[0,:], radial_points1[1,:]))
        radial_points2 = grid.find_nearest_node((radial_points2[0,:], radial_points2[1,:]))
        #print 'radial pts arrays: ', radial_points1, radial_points2
        #print 'On grid?: ', inbounds_test1, inbounds_test2
        #print 'Dist array: ', distance_array
        slope_array = numpy.where(distance_array, (self.elev[radial_points1]-self.elev[radial_points2])/distance_array, numpy.nan)
        slope_array = numpy.arctan(slope_array)
        #if slope  is negative, it means the surface slopes broadly EAST
        try:
            hi_mag_slope_index = numpy.nanargmax(numpy.fabs(slope_array))
            hi_mag_slope = slope_array[hi_mag_slope_index]
        except:
            self._surface_slope = 1.e-10
            #self_direction = self._azimuth_of_travel
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
            #print 'The slope under the crater cavity footprint is: ', self._surface_slope


    def set_crater_mean_slope_v3(self):
        '''
        Runs on a square which encapsulates the crater.
        If some of the nodes are off the grid AND the boundaries aren't looped,
        it falls back on v2.
        If they are, it will freely loo the nodes back onto the grid, and return
        a real slope.
        This version uses the Horn, 1981 algorithm, the same one used by many
        GIS packages.
        '''
        grid = self.grid
        elev = self.elev
        r = 0.7071*self._radius
        x = self._xcoord
        y = self._ycoord
        #dx = grid.dx
        slope_pts =numpy.array([[x-r,y-r],[x,y-r],[x+r,y-r],[x-r,y],[x,y],[x+r,y],[x-r,y+r],[x,y+r],[x+r,y+r]])
        pts_on_grid = grid.is_point_on_grid(slope_pts[:,0],slope_pts[:,1]) #needs to be on **interior** grid
        if not numpy.all(pts_on_grid) and not self.looped_BCs:
            slope_coords_ongrid = slope_pts[pts_on_grid]
            slope_pts_ongrid = grid.find_nearest_node((slope_coords_ongrid[:,0],slope_coords_ongrid[:,1]))
            cardinal_elevs = elev[slope_pts_ongrid]
            self.closest_node_index = grid.find_nearest_node((self._xcoord, self._ycoord))
            self.set_crater_mean_slope_v2()
        else:
            slope_pts %= numpy.array([self.grid.get_grid_xdimension(), self.grid.get_grid_ydimension()]) #added new, to remove boundaries
            slope_pts_ongrid = grid.find_nearest_node((slope_pts[:,0],slope_pts[:,1]))
            self.closest_node_index = slope_pts_ongrid[4]
            cardinal_elevs = elev[slope_pts_ongrid]
            #Now the Horn '81 algorithm for weighted max slope: (careful w signs! altered to give down as +ve)
            S_we = ((cardinal_elevs[6]+2*cardinal_elevs[3]+cardinal_elevs[0])-(cardinal_elevs[8]+2*cardinal_elevs[5]+cardinal_elevs[2]))/(8.*r)
            S_sn = ((cardinal_elevs[0]+2*cardinal_elevs[1]+cardinal_elevs[2])-(cardinal_elevs[6]+2*cardinal_elevs[7]+cardinal_elevs[8]))/(8.*r)
            self._surface_slope = numpy.sqrt(S_we*S_we + S_sn*S_sn)
            if not S_we:
                if S_sn<0.:
                    self._surface_dip_direction = numpy.pi
                else:
                    self._surface_dip_direction = 0.
            else: #general case
                angle_to_xaxis = numpy.arctan(S_sn/S_we) #+ve is CCW rotation from x axis
                self._surface_dip_direction = ((1.-numpy.sign(S_we))*0.5)*numpy.pi + (0.5*numpy.pi-angle_to_xaxis)
        self.closest_node_elev = numpy.mean(cardinal_elevs)
        #self.closest_node_elev = cardinal_elevs[4]




#    #@profile
#    def set_elev_change_only_beneath_footprint_BAND_AID(self):
#        '''
#        This is a method to take an existing impact properties and a known
#        nearest node to the impact site, and alter the topography to model the
#        impact. It assumes crater radius and depth are known, models cavity
#        shape as a power law where n is a function of R/D, and models ejecta
#        thickness as an exponential decay,sensitive to both ballistic range from
#        tilting and momentum transfer in impact (after Furbish). We DO NOT yet
#        model transition to peak ring craters, or enhanced diffusion by ejecta
#        in the strength regime. All craters are dug perpendicular to the geoid,
#        not the surface.
#        This version of the code does NOT correct for slope dip direction -
#        because Furbish showed momentum almost always wins, and these impactors
#        have a lot of momentum!
#        NB - this function ASSUMES that the "beta factor" in the model is <=0.5,
#        i.e., nonlinearities can't develop in the ejecta field, and the impact
#        point is always within the (circular) ejecta footprint.
#        This version of this method ("_band_aid"!) uses a quick and dirty fix
#        which substitutes the actual excavated volume into the equn to derive
#        ejecta thicknesses. It also digs craters perpendicular to the local
#        surface, mimicking some aspects of a "strength dominated" impact -
#        unless doing so would create extremely strongly tilted craters, excavate
#        large sheets of material, or otherwise destabilize the mass balance of
#        the component.
#        It pays no regard for the inefficiency of doing that!
#        Created DEJH Dec 2013
#        '''
#        #Load in the data, for speed and conciseness:
#        _angle_to_vertical=self._angle_to_vertical
#        _surface_slope=self._surface_slope
#        _surface_dip_direction = self._surface_dip_direction
#        _azimuth_of_travel = self._azimuth_of_travel
#        pi = numpy.pi
#        twopi = 2.*pi
#        tan = numpy.tan
#        cos = numpy.cos
#        sin = numpy.sin
#        sqrt = numpy.sqrt
#        #arctan = numpy.arctan
#        arccos = numpy.arccos
#        where = numpy.where
#        _radius = self._radius
#        elev = self.elev
#        tan_repose = self.tan_repose
#        grid = self.grid
#        self.cheater_flag = 0
#        pre_impact_elev = elev.copy()
#
#        #Derive the exponent for the crater shape, shared betw simple & complex:
#        crater_bowl_exp = self.get_crater_shape_exp()
#
#        #Derive the effective angle for impact, relative to the surface normal, beta_eff. The direction is always controlled by the impactor. Draw new impact angles if impact geomtery is impossible.
#        while 1:
#            #epsilon is the angle between the surface normal and the impactor angle to vertical, projected along the line of travel.
#            rake_in_surface_plane = _surface_dip_direction - _azimuth_of_travel
#            absolute_rake = numpy.fabs(rake_in_surface_plane)
#            if not _surface_slope:
#                epsilon = 0.
#                beta_eff = _angle_to_vertical
#            else:
#                if absolute_rake in (0.5*pi, 1.5*pi):
#                    epsilon = 0.
#                else:
#                    epsilon = arccos(cos(rake_in_surface_plane)*cos(_surface_slope)*sqrt(1.+(tan(rake_in_surface_plane)/cos(_surface_slope))**2.))
#                #Make the necessary adjustment to the angle to vertical to reflect dipping surface:
#                if absolute_rake <= 0.5*pi or absolute_rake >= 1.5*pi:
#                    beta_eff = _angle_to_vertical + epsilon
#                else:
#                    beta_eff = _angle_to_vertical + epsilon - pi
#            #print 'Beta effective: ', beta_eff
#            #Make correction to ejecta direction needed if angle to normal is small and slope is large in the opposite direction:
#            if 0. <= beta_eff <= 0.5*pi:
#                _ejecta_azimuth = _azimuth_of_travel
#                break
#            elif beta_eff < 0.:
#                #reverse the azimuth, make beta positive again
#                beta_eff = -beta_eff
#                _ejecta_azimuth = (_azimuth_of_travel+pi)%twopi
#                break
#            else:
#                #Note this situation should now be forbidden as we correct for grazing impacts
#                print 'Impact geometry was not possible! Refreshing the impactor angle...'
#                self.set_impactor_angles()
#                _azimuth_of_travel = self._azimuth_of_travel
#                _angle_to_vertical = self._angle_to_vertical
#
#        mass_bal_corrector_for_slope = self.correct_for_slope()
#        mass_bal_corrector_for_angle_to_vertical = self.correct_for_angle_to_vertical()
#
#        #apply correction to beta to suppress nonlinear BC problem and make ejecta patterns "look like" they actually do:
#        tan_beta = tan(beta_eff*self._beta_factor)
#        tan_beta_sqd = tan_beta*tan_beta
#        #print 'Impact, ejecta azimuths: ', _azimuth_of_travel, _ejecta_azimuth
#
#        ###Here's the oh-so-ugly ad-hoc fix to address the mass balance problems
#        #We're calculating the approx excavated volume directly to feed the ejecta thickness calculator...
#        crater_edge_type_flag, crater_num_repeats = self.footprint_edge_type((self._xcoord,self._ycoord),4.*_radius)
#        crater_footprint_iterator = self.create_square_footsix.print_((self._xcoord,self._ycoord),4.*_radius, crater_edge_type_flag) #Arbitrary increase in radius to try to catch the additional excavated material
#        #sharp lips w/i the grid are CATASTROPHIC (=>get super-steep slopes) - this 4* really has to be sufficient.
#        crater_center_offset_iterator = self.center_offset_list_for_looped_BCs((self._xcoord,self._ycoord), crater_edge_type_flag, crater_num_repeats)
#        excavated_volume = 0.
#        for i,j in izip(crater_footprint_iterator,crater_center_offset_iterator):
#            _vec_r_to_center_excav,_vec_theta_excav = grid.get_distances_of_nodes_to_point(j, get_az='angles', node_subset=i)
#            slope_offsets_rel_to_center_excav = -_vec_r_to_center_excav*numpy.tan(self._surface_slope)*numpy.cos(_vec_theta_excav-self._surface_dip_direction) #node leading negative - downslopes are +ve!!
#            _vec_new_z_excav = slope_offsets_rel_to_center_excav + self.closest_node_elev - 1.*self._depth #Arbitrary assumed scaling: rim is 0.1 of total depth
#            _vec_new_z_excav += self._depth * (_vec_r_to_center_excav/_radius)**crater_bowl_exp #note there's another fix here - we've relaxed the constraint that slopes outside one radius are at repose. It will keep curling up, so poke out quicker
#            z_difference = self.elev[i] - _vec_new_z_excav #+ve when excavating
#            excavated_volume += numpy.sum(numpy.where(z_difference>0.,z_difference,0.))*grid.dx*grid.dx
#        #This all assumes we're representing the surface slope accurately now. If we are, then the real radius will remain roughly as the calculated value as we shear the crater into the surface. If not, could have problems.
#
#        unique_expression_for_local_thickness = self.create_lambda_fn_for_ejecta_thickness_BAND_AID(excavated_volume*mass_bal_corrector_for_slope*mass_bal_corrector_for_angle_to_vertical)
#        thickness_at_rim = unique_expression_for_local_thickness(_radius)
#        if thickness_at_rim<0.:
#            thickness_at_rim = 0. #this shouldn't be a problem anymore
#
#        self.excavated_volume = excavated_volume
#        self.rim_thickness = thickness_at_rim
#
#        #The center is transposed by rmax*tan(beta) in the direction of impact:
#        #Might be an issue here - we use the rmax, so doesn't the *position* of the max thickness then also depend on rmax??
#        #Should it actually just be r? -> no, use the MEAN: see Furbish para 27? For now, just use the max radius - this will be conservative, at least.
#        #solve the ejecta radial thickness equ w. the min thickness to get the max radius:
#        max_radius_ejecta_on_flat = _radius * (thickness_at_rim/self._minimum_ejecta_thickness)**0.3636 #+ 0.5*grid.dx #0.5dx added as a conservative buffer, as we are going to use the closest node as the center... Now removed as we don't use this method
#        #print 'thickness_at_rim: ', thickness_at_rim
#        #print 'max_radius_ejecta: ', max_radius_ejecta_on_flat
#
#        displacement_distance = max_radius_ejecta_on_flat*tan_beta
#        #need a variable for if the exacavated radius is ever outside the locus defined by the MREOF around the footprint center:
#        x_impact_offset = sin(_azimuth_of_travel)*displacement_distance
#        y_impact_offset = cos(_azimuth_of_travel)*displacement_distance
#        #enlarge the area all round if the crater itself peeks out of the ejecta area at low angles:
#        if x_impact_offset+4*self._radius > max_radius_ejecta_on_flat:
#            print 'A low-angle crater!'
#            max_radius_ejecta_on_flat = x_impact_offset+4.*self._radius
#        if y_impact_offset+4*self._radius > max_radius_ejecta_on_flat:
#            print 'A low-angle crater!'
#            max_radius_ejecta_on_flat = y_impact_offset+4.*self._radius
#
#        footprint_center_x = self._xcoord+x_impact_offset
#        footprint_center_y = self._ycoord+y_impact_offset
#        edge_type_flag, num_repeats = self.footprint_edge_type((footprint_center_x,footprint_center_y),max_radius_ejecta_on_flat)
#        footprint_iterator = self.create_square_footsix.print_((footprint_center_x,footprint_center_y),max_radius_ejecta_on_flat, edge_type_flag)
#        center_offset_iterator = self.center_offset_list_for_looped_BCs((self._xcoord,self._ycoord), edge_type_flag, num_repeats)
#
#        for footprint_nodes,center_tuple in izip(footprint_iterator, center_offset_iterator):
#            elev = self.elev
#            _vec_r_to_center, _vec_theta = grid.get_distances_of_nodes_to_point(center_tuple, get_az='angles') #, node_subset=footprint_nodes)
#            _vec_r_to_center = _vec_r_to_center[footprint_nodes]
#            _vec_theta = _vec_theta[footprint_nodes]
#
#            ##We need to account for deposition depth elevating the crater rim, i.e., we need to deposit *before* we cut the cavity. We do this by defining three domains for the node to lie in: 1. r<r_calc, i.e., below the pre-impact surface. No risk of intersecting the surface here. 2. r_calc < r; Th>z_new. this is the domain in the inward sloping rim of the crater ejecta. 3. Th<z_new and beyond. out on the ejecta proper. Note - (1) is not hard & fast rule if the surface dips. Safer is just (Th-lowering)<z_new
#            _vec_new_z = numpy.empty_like(_vec_theta)
#            _vec_new_z.fill(self.closest_node_elev + thickness_at_rim - self._depth)
#
#            _vec_new_z += self._depth * (_vec_r_to_center/_radius)**crater_bowl_exp
#
#            _vec_theta_eff = _ejecta_azimuth - _vec_theta #This is the angle of the center-to-active-node line to the azimuth along which the ejecta is concentrated
#            _vec_sin_theta_sqd = sin(_vec_theta_eff) ** 2.
#            _vec_cos_theta = cos(_vec_theta_eff)
#
#            #need to be way more careful in the way we handle filling within the cavity. Try this:
#            inner_radius = _radius #- thickness_at_rim/tan_repose
#            _nodes_below_theoretical_ground_level = _vec_r_to_center<=inner_radius
#            elevs_cavity_less_ground = _vec_new_z[_nodes_below_theoretical_ground_level]-pre_impact_elev[footprint_nodes][_nodes_below_theoretical_ground_level]
#            #Assume everything in the inner radius gets filled to the crater level:
#            elevs_ground_less_new_z = pre_impact_elev[footprint_nodes] - _vec_new_z
#            volume_to_fill_inner_crater_divots = 0.#numpy.sum(numpy.where(elevs_cavity_less_ground>0., elevs_cavity_less_ground, 0.))*grid.dx*grid.dx
#            volume_to_remove_highs = numpy.sum(numpy.where(elevs_ground_less_new_z>0., elevs_ground_less_new_z, 0.))*grid.dx*grid.dx
#            #print 'lows ', volume_to_fill_inner_crater_divots
#            #print 'highs', volume_to_remove_highs
#            #...then one more iteration on the volumes:
#            unique_expression_for_local_thickness = self.create_lambda_fn_for_ejecta_thickness_BAND_AID((volume_to_remove_highs-volume_to_fill_inner_crater_divots)*mass_bal_corrector_for_slope*mass_bal_corrector_for_angle_to_vertical)
#            #thickness_at_rim = unique_expression_for_local_thickness(_radius)
#            #if thickness_at_rim<0.:
#            #    thickness_at_rim = 0. #this shouldn't be a problem anymore
#            self.excavated_volume = excavated_volume
#            #self.rim_thickness = thickness_at_rim
#
#            _vec_flat_thickness_above_surface = unique_expression_for_local_thickness(_vec_r_to_center)
#
#            #This material is not necessary as part of this footprint-based method, as we already forbid beta_factor>0.5
#            ##REMEMBER, as tan_beta gets >1, the function describing the ejecta is only valid over ever more restricted ranges of theta!! In other words,
#            _vec_mu_theta_by_mu0 = tan_beta * _vec_cos_theta + sqrt(1. - _vec_sin_theta_sqd * tan_beta_sqd)
#            _vec_f_theta = (tan_beta_sqd*(_vec_cos_theta**2.-_vec_sin_theta_sqd) + 2.*tan_beta*_vec_cos_theta*sqrt(1.-tan_beta_sqd*_vec_sin_theta_sqd) + 1.) / twopi
#            #So, distn_at_angle = distn_vertical_impact*f_theta/mu_theta_by_mu0. Draw the thickness at the active node:
#            #NB-the 2pi is to correct for mismatch in the dimensions of mu and f
#            _vec_thickness = _vec_f_theta/_vec_mu_theta_by_mu0 * twopi * _vec_flat_thickness_above_surface
#            #Set the thicknesses <0 to 0:
#            _vec_thickness_positive = where(_vec_thickness>=0.,_vec_thickness, 0.)
#            #Now, are we inside or outside the rim?
#            absolute_ejecta_elevations = pre_impact_elev[footprint_nodes]+_vec_thickness_positive
#            final_elevs = numpy.amin(numpy.vstack((absolute_ejecta_elevations, _vec_new_z)), axis=0)
#            final_elevs -= pre_impact_elev[footprint_nodes] #need to make it a difference, not absolute, if we're stacking tiles
#            #can't use the init_mass_balance to assess if this is necessary if we're looping now...
#            self.elev[footprint_nodes] += final_elevs
#            #Save any data to the higher level:
#            if edge_type_flag == 'C' or edge_type_flag == 'X':
#                elev_diff = self.elev[footprint_nodes] - pre_impact_elev[footprint_nodes]
#        if edge_type_flag != 'C' and edge_type_flag != 'X':
#            elev_diff = self.elev - pre_impact_elev
#        self.mass_balance_in_impact = numpy.sum(elev_diff)/-numpy.sum(elev_diff[elev_diff<0.]) #positive is mass gain, negative is mass loss. This is currently a mass fraction, given relative to volume (h*px#) excavated from below the original surface.
#
#        print 'sum elev_diff', numpy.sum(elev_diff)
#        print 'sum elev_diff<0', -numpy.sum(elev_diff[elev_diff<0.])
#        #now perform the check on mass balance to try to iron out the mass balance issues.
#
#        self.ejecta_azimuth = _ejecta_azimuth
#        self.impactor_angle_to_surface_normal = beta_eff #note in this case this is the *effective* angle (in the direction of travel), not the actual angle to the surface.
#        #print 'Vol of crater cavity: ', self._cavity_volume
#
#        #set this to see where the footprint is, and to see values:
#        #self.elev[footprint_nodes] = _vec_r_to_center
#
#        return self.elev


    def set_elev_change_only_beneath_footprint_BAND_AID_memory_save(self):
        '''
        This is a method to take an existing impact properties and a known
        nearest node to the impact site, and alter the topography to model the
        impact. It assumes crater radius and depth are known, models cavity
        shape as a power law where n is a function of R/D, and models ejecta
        thickness as an exponential decay,sensitive to both ballistic range from
        tilting and momentum transfer in impact (after Furbish). We DO NOT yet
        model transition to peak ring craters, or enhanced diffusion by ejecta
        in the strength regime. All craters are dug perpendicular to the geoid,
        not the surface.
        This version of the code does NOT correct for slope dip direction -
        because Furbish showed momentum almost always wins, and these impactors
        have a lot of momentum!
        NB - this function ASSUMES that the "beta factor" in the model is <=0.5,
        i.e., nonlinearities can't develop in the ejecta field, and the impact
        point is always within the (circular) ejecta footprint.
        This version of this method ("_band_aid"!) uses a quick and dirty fix
        which substitutes the actual excavated volume into the equn to derive
        ejecta thicknesses. It also digs craters perpendicular to the local
        surface, mimicking some aspects of a "strength dominated" impact -
        unless doing so would create extremely strongly tilted craters, excavate
        large sheets of material, or otherwise destabilize the mass balance of
        the component.
        It pays no regard for the inefficiency of doing that!
        Created DEJH Dec 2013
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
        elev = self.elev
        tan_repose = self.tan_repose
        grid = self.grid
        self.cheater_flag = 0
        elev_diff_below_ground = 0.
        total_elev_diff = 0.
        #pre_impact_elev = elev.copy() #####DON'T DO THIS
        all_the_footprint_nodes = numpy.empty(self.grid.number_of_node_rows*self.grid.number_of_node_columns, dtype = int) #this will store the union of footprint across all loops over the grid
        number_of_total_footprint_nodes = 0.

        #Derive the exponent for the crater shape, shared betw simple & complex:
        crater_bowl_exp = self.get_crater_shape_exp()

        #Derive the effective angle for impact, relative to the surface normal, beta_eff. The direction is always controlled by the impactor. Draw new impact angles if impact geomtery is impossible.
        counter = 0
        while 1:
            #epsilon is the angle between the surface normal and the impactor angle to vertical, projected along the line of travel.
            rake_in_surface_plane = _surface_dip_direction - _azimuth_of_travel
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
            if 0. <= beta_eff <= 0.5*numpy.pi:
                _ejecta_azimuth = _azimuth_of_travel
                break
            elif -0.5*numpy.pi < beta_eff < 0.:
                #reverse the azimuth, make beta positive again
                beta_eff = -beta_eff
                _ejecta_azimuth = (_azimuth_of_travel+pi)%twopi
                break
            else:
                #Note this situation should now be forbidden as we correct for grazing impacts
                six.print_('Impact geometry was not possible! Refreshing the impactor angle...')
                self.set_impactor_angles()
                _azimuth_of_travel = self._azimuth_of_travel
                _angle_to_vertical = self._angle_to_vertical
                if counter > 20:
                    beta_eff = 0.
                    _ejecta_azimuth = self._azimuth_of_travel
                    break
                else:
                    counter+=1

        mass_bal_corrector_for_slope = self.correct_for_slope()
        mass_bal_corrector_for_angle_to_vertical = self.correct_for_angle_to_vertical()
        mass_bal_corrector_for_size = self.correct_for_crater_size()

        #apply correction to beta to suppress nonlinear BC problem and make ejecta patterns "look like" they actually do:
        tan_beta = tan(beta_eff*self._beta_factor)
        #if abs(tan_beta) > 1.:
        #    print 'tan_beta is ', tan_beta
        #    #Q&D fix for tan_beta>1... but why does this arise?
        #    if tan_beta > 1.:
        #        tan_beta = 1.
        #    elif tan_beta < -1.:
        #        tan_beta = -1.
        tan_beta_sqd = tan_beta*tan_beta
        #print 'Impact, ejecta azimuths: ', _azimuth_of_travel, _ejecta_azimuth

        ###Here's the oh-so-ugly ad-hoc fix to address the mass balance problems
        #We're calculating the approx excavated volume directly to feed the ejecta thickness calculator...
        crater_edge_type_flag, crater_num_repeats = self.footprint_edge_type((self._xcoord,self._ycoord),4.*_radius)
        crater_footprint_iterator = self.create_square_footsix.print_((self._xcoord,self._ycoord),4.*_radius, crater_edge_type_flag, self.crater_footprint_max) #Arbitrary increase in radius to try to catch the additional excavated material
        #sharp lips w/i the grid are CATASTROPHIC (=>get super-steep slopes) - this 4* really has to be sufficient.
        crater_center_offset_iterator = self.center_offset_list_for_looped_BCs((self._xcoord,self._ycoord), crater_edge_type_flag, crater_num_repeats)
        excavated_volume = 0.
        for footprint_tuple,j in six.moves.zip(crater_footprint_iterator,crater_center_offset_iterator):
            i = footprint_tuple[0]
            self.crater_footprint_max = footprint_tuple[1]
            grid.get_distances_of_nodes_to_point(j, get_az='angles', node_subset=self.crater_footprint_max[:i], out_distance=self._vec_r_to_center[:i], out_azimuth=self._vec_theta[:i])
            numpy.subtract(self._vec_theta[:i], self._surface_dip_direction, out=self.dummy_1[:i])
            numpy.cos(self.dummy_1[:i], out=self.dummy_2[:i])
            numpy.multiply(self.dummy_2[:i], self._vec_r_to_center[:i], out=self.dummy_1[:i])
            numpy.multiply(self.dummy_1[:i], -numpy.tan(self._surface_slope), out=self.slope_offsets_rel_to_center[:i])
            #self.slope_offsets_rel_to_center[:i] = -self._vec_r_to_center[:i]*numpy.tan(self._surface_slope)*numpy.cos(self._vec_theta[:i]-self._surface_dip_direction) #node leading negative - downslopes are +ve!!
            numpy.add(self.slope_offsets_rel_to_center[:i], self.closest_node_elev - 0.9*self._depth, out=self.dummy_1[:i])
            #self._vec_new_z[:i] = self.slope_offsets_rel_to_center[:i] + self.closest_node_elev - 0.9*self._depth #Arbitrary assumed scaling: rim is 0.1 of total depth... this was 0, so our correction fns will need updating 7/3
            numpy.divide(self._vec_r_to_center[:i], _radius, out=self.dummy_2[:i])
            numpy.power(self.dummy_2[:i], crater_bowl_exp, out=self.dummy_3[:i])
            numpy.multiply(self.dummy_3[:i], self._depth, out=self.dummy_2[:i])
            numpy.add(self.dummy_2[:i], self.dummy_1[:i], out=self._vec_new_z[:i])
            #self._vec_new_z[:i] += self._depth * (self._vec_r_to_center[:i]/_radius)**crater_bowl_exp #note there's another fix here - we've relaxed the constraint that slopes outside one radius are at repose. It will keep curling up, so poke out quicker
            numpy.subtract(self.elev[self.crater_footprint_max[:i]], self._vec_new_z[:i], out=self.z_difference[:i])
            #self.z_difference[:i] = self.elev[self.crater_footprint_max[:i]] - self._vec_new_z[:i] #+ve when excavating
            numpy.less_equal(self.z_difference[:i], 0., self.dummy_bool[:i])
            self.z_difference[:i][self.dummy_bool[:i]] = 0.
            excavated_volume += numpy.sum(self.z_difference[:i])*grid.dx*grid.dy
            #excavated_volume += numpy.sum(numpy.where(self.z_difference[:i]>0.,self.z_difference[:i],0.))*grid.dx*grid.dx
        #This all assumes we're representing the surface slope accurately now. If we are, then the real radius will remain roughly as the calculated value as we shear the crater into the surface. If not, could have problems.

        unique_expression_for_local_thickness = self.create_lambda_fn_for_ejecta_thickness_BAND_AID(excavated_volume*mass_bal_corrector_for_slope*mass_bal_corrector_for_angle_to_vertical*mass_bal_corrector_for_size)
        thickness_at_rim = unique_expression_for_local_thickness(_radius)
        if thickness_at_rim<0.:
            #raise ValueError
            thickness_at_rim = 0. #this shouldn't be a problem anymore
        six.print_('First pass thickness: ', thickness_at_rim)
        six.print_('dist above ground: ', thickness_at_rim - self._depth)

        #now repeat the loop to improve precision on mass balance:
        if crater_edge_type_flag=='C' or crater_edge_type_flag=='X': #cases with only 1 loop; can reuse previous variables for speed
            self._vec_new_z[:i].fill(self.closest_node_elev + thickness_at_rim - self._depth)
            six.print_('dist above ground: ', thickness_at_rim - self._depth)
            numpy.add(self.dummy_2[:i], self._vec_new_z[:i], out=self._vec_new_z[:i])
            numpy.subtract(self.pre_impact_elev[self.crater_footprint_max[:i]], self._vec_new_z[:i], out=self.elevs_ground_less_new_z[:i])
            volume_to_fill_inner_crater_divots = 0.
            numpy.less_equal(self.elevs_ground_less_new_z[:i], 0., self.dummy_bool[:i])
            self.elevs_ground_less_new_z[:i][self.dummy_bool[:i]] = 0.
            volume_to_remove_highs = numpy.sum(self.elevs_ground_less_new_z[:i])*grid.dx*grid.dy
            unique_expression_for_local_thickness = self.create_lambda_fn_for_ejecta_thickness_BAND_AID((volume_to_remove_highs-volume_to_fill_inner_crater_divots)*mass_bal_corrector_for_slope*mass_bal_corrector_for_angle_to_vertical*mass_bal_corrector_for_size)
            thickness_at_rim = unique_expression_for_local_thickness(_radius)
        else:
            crater_edge_type_flag, crater_num_repeats = self.footprint_edge_type((self._xcoord,self._ycoord),4.*_radius)
            crater_footprint_iterator = self.create_square_footsix.print_((self._xcoord,self._ycoord),4.*_radius, crater_edge_type_flag, self.crater_footprint_max) #Arbitrary increase in radius to try to catch the additional excavated material
            crater_center_offset_iterator = self.center_offset_list_for_looped_BCs((self._xcoord,self._ycoord), crater_edge_type_flag, crater_num_repeats)
            volume_to_remove_highs = 0.
            for footprint_tuple,center_tuple in zip(crater_footprint_iterator,crater_center_offset_iterator):
#####is this right to do every loop?
                six.print_('recalc values...')
                i = footprint_tuple[0]
                self.crater_footprint_max = footprint_tuple[1]
                footprint_nodes = self.crater_footprint_max[:i]
                grid.get_distances_of_nodes_to_point(center_tuple, get_az='angles', node_subset=footprint_nodes, out_distance=self._vec_r_to_center[:i], out_azimuth=self._vec_theta[:i])
                self._vec_new_z[:i].fill(self.closest_node_elev + thickness_at_rim - self._depth)
                numpy.divide(self._vec_r_to_center[:i], _radius, out=self.dummy_2[:i])
                numpy.power(self.dummy_2[:i], crater_bowl_exp, out=self.dummy_3[:i])
                numpy.multiply(self.dummy_3[:i], self._depth, out=self.dummy_2[:i])
                numpy.add(self.dummy_2[:i], self._vec_new_z[:i], out=self._vec_new_z[:i])
                numpy.subtract(self.pre_impact_elev[footprint_nodes], self._vec_new_z[:i], out=self.elevs_ground_less_new_z[:i])
                volume_to_fill_inner_crater_divots = 0.#numpy.sum(numpy.where(elevs_cavity_less_ground>0., elevs_cavity_less_ground, 0.))*grid.dx*grid.dx
                numpy.less_equal(self.elevs_ground_less_new_z[:i], 0., self.dummy_bool[:i])
                self.elevs_ground_less_new_z[:i][self.dummy_bool[:i]] = 0.
                volume_to_remove_highs += numpy.sum(self.elevs_ground_less_new_z[:i])*grid.dx*grid.dy
            unique_expression_for_local_thickness = self.create_lambda_fn_for_ejecta_thickness_BAND_AID((volume_to_remove_highs-volume_to_fill_inner_crater_divots)*mass_bal_corrector_for_slope*mass_bal_corrector_for_angle_to_vertical*mass_bal_corrector_for_size)
            thickness_at_rim = unique_expression_for_local_thickness(_radius)

        six.print_(volume_to_remove_highs)
        six.print_(volume_to_fill_inner_crater_divots)
        if thickness_at_rim<0.:
            thickness_at_rim = 0. #this shouldn't be a problem anymore
        self.rim_thickness = thickness_at_rim

        #The center is transposed by rmax*tan(beta) in the direction of impact:
        #Might be an issue here - we use the rmax, so doesn't the *position* of the max thickness then also depend on rmax??
        #Should it actually just be r? -> no, use the MEAN: see Furbish para 27? For now, just use the max radius - this will be conservative, at least.
        #solve the ejecta radial thickness equ w. the min thickness to get the max radius:
        max_radius_ejecta_on_flat = _radius * (thickness_at_rim/self._minimum_ejecta_thickness)**0.3636 #+ 0.5*grid.dx #0.5dx added as a conservative buffer, as we are going to use the closest node as the center... Now removed as we don't use this method
        #print 'thickness_at_rim: ', thickness_at_rim
        #print 'max_radius_ejecta: ', max_radius_ejecta_on_flat

        displacement_distance = max_radius_ejecta_on_flat*tan_beta
        #need a variable for if the exacavated radius is ever outside the locus defined by the MREOF around the footprint center:
        x_impact_offset = sin(_azimuth_of_travel)*displacement_distance
        y_impact_offset = cos(_azimuth_of_travel)*displacement_distance
        #enlarge the area all round if the crater itself peeks out of the ejecta area at low angles:
        if x_impact_offset+4*self._radius > max_radius_ejecta_on_flat:
            six.print_('A low-angle crater!')
            max_radius_ejecta_on_flat = x_impact_offset+4.*self._radius
        if y_impact_offset+4*self._radius > max_radius_ejecta_on_flat:
            six.print_('A low-angle crater!')
            max_radius_ejecta_on_flat = y_impact_offset+4.*self._radius

        footprint_center_x = self._xcoord+x_impact_offset
        footprint_center_y = self._ycoord+y_impact_offset
        edge_type_flag, num_repeats = self.footprint_edge_type((footprint_center_x,footprint_center_y),max_radius_ejecta_on_flat)
        six.print_(edge_type_flag)
        footprint_iterator = self.create_square_footsix.print_((footprint_center_x,footprint_center_y),max_radius_ejecta_on_flat, edge_type_flag, self.crater_footprint_max)
        center_offset_iterator = self.center_offset_list_for_looped_BCs((self._xcoord,self._ycoord), edge_type_flag, num_repeats)
        total_elev_diffs = 0.
        elev_diffs_below_ground = 0.

        for footprint_tuple,center_tuple in zip(footprint_iterator, center_offset_iterator):
            i = footprint_tuple[0]
            self.crater_footprint_max = footprint_tuple[1]
            six.print_('looping... effective center is at ', center_tuple)
            footprint_nodes = self.crater_footprint_max[:i]
            #print footprint_nodes
            #self._vec_r_to_center[:i], self._vec_theta[:i] = grid.get_distances_of_nodes_to_point(center_tuple, get_az='angles', node_subset=footprint_nodes)
            grid.get_distances_of_nodes_to_point(center_tuple, get_az='angles', node_subset=footprint_nodes, out_distance=self._vec_r_to_center[:i], out_azimuth=self._vec_theta[:i])
            #_vec_r_to_center = _vec_r_to_center[footprint_nodes]
            #_vec_theta = _vec_theta[footprint_nodes]

            ##We need to account for deposition depth elevating the crater rim, i.e., we need to deposit *before* we cut the cavity. We do this by defining three domains for the node to lie in: 1. r<r_calc, i.e., below the pre-impact surface. No risk of intersecting the surface here. 2. r_calc < r; Th>z_new. this is the domain in the inward sloping rim of the crater ejecta. 3. Th<z_new and beyond. out on the ejecta proper. Note - (1) is not hard & fast rule if the surface dips. Safer is just (Th-lowering)<z_new
            #_vec_new_z = numpy.empty_like(_vec_theta)
            #self._vec_new_z[:i].fill(self.closest_node_elev - 0.9*self._depth)
            self._vec_new_z[:i].fill(self.closest_node_elev + thickness_at_rim - self._depth)

            numpy.divide(self._vec_r_to_center[:i], _radius, out=self.dummy_2[:i])
            numpy.power(self.dummy_2[:i], crater_bowl_exp, out=self.dummy_3[:i])
            numpy.multiply(self.dummy_3[:i], self._depth, out=self.dummy_2[:i])
            numpy.add(self.dummy_2[:i], self._vec_new_z[:i], out=self._vec_new_z[:i])
            #self._vec_new_z[:i] += self._depth * (self._vec_r_to_center[:i]/_radius)**crater_bowl_exp

            #self._vec_theta_eff[:i] = _ejecta_azimuth - self._vec_theta[:i] #This is the angle of the center-to-active-node line to the azimuth along which the ejecta is concentrated
            #_vec_sin_theta_sqd = sin(_ejecta_azimuth - self._vec_theta[:i]) ** 2.
            #_vec_cos_theta = cos(_ejecta_azimuth - self._vec_theta[:i])

            #need to be way more careful in the way we handle filling within the cavity. Try this:
            #inner_radius = _radius #- thickness_at_rim/tan_repose
            #_nodes_below_theoretical_ground_level = self._vec_r_to_center[:i]<=inner_radius
            #elevs_cavity_less_ground = self._vec_new_z[:i][_nodes_below_theoretical_ground_level]-pre_impact_elev[footprint_nodes][_nodes_below_theoretical_ground_level]
            #Assume everything in the inner radius gets filled to the crater level:
            ##numpy.subtract(self.pre_impact_elev[footprint_nodes], self._vec_new_z[:i], out=self.elevs_ground_less_new_z[:i])
            #self.elevs_ground_less_new_z[:i] = pre_impact_elev[footprint_nodes] - self._vec_new_z[:i]
            ##volume_to_fill_inner_crater_divots = 0.#numpy.sum(numpy.where(elevs_cavity_less_ground>0., elevs_cavity_less_ground, 0.))*grid.dx*grid.dx
            ##self.dummy_1[:i] = numpy.where(self.elevs_ground_less_new_z[:i]>0., self.elevs_ground_less_new_z[:i], 0.)
            ##volume_to_remove_highs += numpy.sum(self.dummy_1[:i])*grid.dx*grid.dx
            #volume_to_remove_highs = numpy.sum(numpy.where(self.elevs_ground_less_new_z[:i]>0., self.elevs_ground_less_new_z[:i], 0.))*grid.dx*grid.dx
            #print 'lows ', volume_to_fill_inner_crater_divots
            #print 'highs', volume_to_remove_highs
            #...then one more iteration on the volumes:
            ##unique_expression_for_local_thickness = self.create_lambda_fn_for_ejecta_thickness_BAND_AID((volume_to_remove_highs-volume_to_fill_inner_crater_divots)*mass_bal_corrector_for_slope*mass_bal_corrector_for_angle_to_vertical*mass_bal_corrector_for_size)
            #thickness_at_rim = unique_expression_for_local_thickness(_radius)
            #if thickness_at_rim<0.:
            #    thickness_at_rim = 0. #this shouldn't be a problem anymore
            ##self.excavated_volume = excavated_volume
            #self.rim_thickness = thickness_at_rim

            self._vec_flat_thickness_above_surface[:i] = unique_expression_for_local_thickness(self._vec_r_to_center[:i])

            #This material is not necessary as part of this footprint-based method, as we already forbid beta_factor>0.5
            ##REMEMBER, as tan_beta gets >1, the function describing the ejecta is only valid over ever more restricted ranges of theta!! In other words,
            numpy.subtract(_ejecta_azimuth, self._vec_theta[:i], out=self.dummy_1[:i])
            numpy.cos(self.dummy_1[:i], out=self.dummy_2[:i])
            numpy.sin(self.dummy_1[:i], out=self.dummy_3[:i])
            numpy.square(self.dummy_3[:i], out=self.dummy_1[:i])
            numpy.multiply(self.dummy_1[:i], tan_beta_sqd, out=self.dummy_4[:i])
            if tan_beta_sqd > 1.:
                six.print_('tan_beta_sqd error!', beta_eff, tan_beta_sqd, epsilon)
                raise ValueError
            numpy.subtract(1., self.dummy_4[:i], out=self.dummy_3[:i])
            if numpy.any(numpy.less(self.dummy_3[:i], 0.)):
                six.print_('about to crash...', numpy.sum(numpy.less(self.dummy_3[:i], 0.)), self.dummy_3[:i][numpy.where(numpy.less(self.dummy_3[:i], 0.))])
                raise ValueError
            numpy.sqrt(self.dummy_3[:i], out=self.dummy_4[:i])
            numpy.multiply(self.dummy_2[:i], tan_beta, out=self.dummy_3[:i])
            numpy.add(self.dummy_3[:i], self.dummy_4[:i], out=self._vec_mu_theta_by_mu0[:i])
            #self._vec_mu_theta_by_mu0[:i] = tan_beta * cos(_ejecta_azimuth - self._vec_theta[:i]) + sqrt(1. - sin(_ejecta_azimuth - self._vec_theta[:i]) ** 2. * tan_beta_sqd)
            numpy.square(self.dummy_2[:i], out=self.dummy_3[:i])
            numpy.subtract(self.dummy_3[:i], self.dummy_1[:i], out=self.dummy_4[:i])
            numpy.multiply(self.dummy_4[:i], tan_beta_sqd, out=self.dummy_3[:i]) #the LH fragment
            numpy.multiply(self.dummy_1[:i], tan_beta_sqd, out=self.dummy_4[:i])
            numpy.subtract(1., self.dummy_4[:i], out=self.dummy_1[:i])
            numpy.sqrt(self.dummy_1[:i], out=self.dummy_4[:i])
            numpy.multiply(self.dummy_2[:i], self.dummy_4[:i], out=self.dummy_1[:i])
            numpy.multiply(self.dummy_1[:i], 2.*tan_beta, out=self.dummy_4[:i])
            numpy.add(self.dummy_3[:i], self.dummy_4[:i], out=self.dummy_1[:i]) #the RH fragment
            numpy.add(self.dummy_1[:i], 1., out=self.dummy_3[:i])
            numpy.divide(self.dummy_3[:i], twopi, out=self._vec_f_theta[:i])
            #self._vec_f_theta[:i] = (tan_beta_sqd*(cos(_ejecta_azimuth - self._vec_theta[:i])**2.-sin(_ejecta_azimuth - self._vec_theta[:i]) ** 2.) + 2.*tan_beta*cos(_ejecta_azimuth - self._vec_theta[:i])*sqrt(1.-tan_beta_sqd*sin(_ejecta_azimuth - self._vec_theta[:i]) ** 2.) + 1.) / twopi
            #So, distn_at_angle = distn_vertical_impact*f_theta/mu_theta_by_mu0. Draw the thickness at the active node:
            #NB-the 2pi is to correct for mismatch in the dimensions of mu and f
            numpy.divide(self._vec_f_theta[:i], self._vec_mu_theta_by_mu0[:i], out=self.dummy_1[:i])
            numpy.multiply(self.dummy_1[:i], self._vec_flat_thickness_above_surface[:i], out=self.dummy_2[:i])
            numpy.multiply(self.dummy_2[:i], twopi, out=self._vec_thickness[:i])
            #self._vec_thickness[:i] = self._vec_f_theta[:i]/self._vec_mu_theta_by_mu0[:i] * twopi * self._vec_flat_thickness_above_surface[:i]
            #Set the thicknesses <0 to 0:
            numpy.less(self._vec_thickness[:i], 0., self.dummy_bool[:i])
            self._vec_thickness[:i][self.dummy_bool[:i]] = 0.
            #Now, are we inside or outside the rim?

            numpy.add(self.pre_impact_elev[footprint_nodes], self._vec_thickness[:i], out=self.dummy_1[:i]) #the ejecta thicknesses
            self.double_dummy[0,:i] = self.dummy_1[:i]
            self.double_dummy[1,:i] = self._vec_new_z[:i] #the excavation depths
            numpy.amin(self.double_dummy[:,:i], axis=0, out=self.dummy_2[:i])
            numpy.subtract(self.dummy_2[:i], self.pre_impact_elev[footprint_nodes], out=self.final_elev_diffs[:i])
            #can't use the init_mass_balance to assess if this is necessary if we're looping now...
            numpy.add(self.final_elev_diffs[:i],self.elev[footprint_nodes], out=self.dummy_1[:i])
            self.elev[footprint_nodes] = self.dummy_1[:i]
            #Update the mass balance:
            total_elev_diffs += numpy.sum(self.final_elev_diffs[:i])
            numpy.less(self.final_elev_diffs[:i], 0., out=self.dummy_bool[:i])
            elev_diffs_below_ground += -numpy.sum(self.final_elev_diffs[:i][self.dummy_bool[:i]])

        six.print_(total_elev_diffs)
        six.print_(elev_diffs_below_ground)
        self.mass_balance_in_impact = total_elev_diffs / elev_diffs_below_ground #positive is mass gain, negative is mass loss. This is currently a mass fraction, given relative to volume (h*px#) excavated from below the original surface.

        self.ejecta_azimuth = _ejecta_azimuth
        self.impactor_angle_to_surface_normal = beta_eff #note in this case this is the *effective* angle (in the direction of travel), not the actual angle to the surface.
        #print 'Vol of crater cavity: ', self._cavity_volume

        return self.elev


    def footprint_edge_type(self, center, eff_radius):
        '''
        Returns the edge type of a given node footprint to be build with create_
        square_footsix.print_()  around the given center.
        Returns 2 values:
        * One of N,S,E,W,NW,NE,SW,SE, if the footprint overlaps a single edge or
          corner; 'C' for an entirely enclosed footprint; 'I' if it intersects
          opposite sides of the grid (i.e., bigger than the grid); 'X' if looped
          boundaries aren't set.
        * An integer. If the first return is 'I', this is the number of "whole"
          grids the footprint could enclose (minimum 1). If it's something else,
          returns 0.
        '''
        assert type(center) == tuple
        assert len(center) == 2
        grid_x = self.grid.get_grid_xdimension()-self.grid.dx #as the edge nodes are looped!
        grid_y = self.grid.get_grid_ydimension()-self.grid.dy
        six.print_('center, r, x', center[0], eff_radius, grid_x)
        left_repeats = -int((center[0]-eff_radius)//grid_x)
        right_repeats = int((center[0]+eff_radius)//grid_x)
        top_repeats = int((center[1]+eff_radius)//grid_y)
        bottom_repeats = -int((center[1]-eff_radius)//grid_y)
        if (left_repeats and right_repeats) or (top_repeats and bottom_repeats):
            big_foot = True
        else:
            big_foot = False
        if big_foot and self.looped_BCs:
            return 'I', int(max([left_repeats,right_repeats,top_repeats,bottom_repeats]))
        elif self.looped_BCs:
            if left_repeats:
                if top_repeats:
                    flag = 'NW'
                elif bottom_repeats:
                    flag = 'SW'
                else:
                    flag = 'W'
            elif right_repeats:
                if top_repeats:
                    flag = 'NE'
                elif bottom_repeats:
                    flag = 'SE'
                else:
                    flag = 'E'
            elif top_repeats:
                flag = 'N'
            elif bottom_repeats:
                flag = 'S'
            else:
                flag = 'C'
            return flag, 0
        else:
            return 'X', 0

    def create_square_footprint(self, center, eff_radius, footprint_edge_type, array_in=None):
        '''
        This method creates a square footprint of nodes around a given center
        point, with a specified halfwidth.
        It is designed to avoid the need to actually search the whole grid
        in order to establish the footprint, to accelerate the craters
        module.
        "Center" is a tuple, (x,y), the footprint's center.
        eff_radius is the footprint halfwidth.
        footprint_edge_type and whole_grid_repeats are the outputs from
        footprint_edge_type() for this footprint.
        The function is a generator, designed to interface with a loop on the
        grids it outputs. The method generates intelligently:
        -If loops AREN'T active, it returns just a single array and terminates.
        Note we can never return the IDs of the boundary nodes - we shouldn't
        ever be operating directly on them in the meat of this module.
        -If loops ARE active, it differentiates between big ('I') footprints
        and the others-
            Edge type 'C' - the entire footprint is within the grid, and one
              array is output.
            Edge type 'N', 'E', 'S', 'W' - the footprint overlaps with one
              edge only. Two arrays follow; the center grid nodes, then the
              edge nodes.
            Edge type 'NE','SE','SW','NW' - the footprint overlaps with a
              corner and two edges. three arrays follow; the center tile, then
              the edges (clockwise), then the corner.

          If it's 'I', it means the footprint is "large" compared to the
          grid. The method then uses the integer stored in whole_grid_repeats.
          this works as 1: one whole grid plus edges needed.
          2: a 3x3 grid surrounded by edges. 3: a 5x5 grid surrounded by edges.
          The intention is that these will be rare enough just using whole grids
          is not too inefficient.


          ______________________________________
          |                                     | 'I'
          |      * * * * * * * * * * *          |
          |      *  ___              *          |
          |      * |   | 'C'        _*___       |
          |      * |___|           | *   | 'E'  |
          |      *                 | *   |      |
          |    __*__               |_*___|      |
          |   |  *  | 'SW'           *          |
          |   |  * * * * * * * * * * *          |
          |   |_____|                           |
          |                                     |
          |_____________________________________|


        '''
        if array_in is not None:
            x = numpy.empty(self.grid.number_of_node_columns-2)
            y_column = numpy.empty(self.grid.number_of_node_rows-2).reshape((self.grid.number_of_node_rows-2,1))
        assert type(center) == tuple
        assert len(center) == 2
        if footprint_edge_type == 'I':
            if array_in is not None:
                number_of_footprint_nodes = self.grid.number_of_interior_nodes
                array_in = self.grid.core_nodes
                while 1:
                    yield (number_of_footprint_nodes, array_in)
            else:
                while 1:
                    yield self.grid.core_nodes
        else:
            center_array = numpy.array(center)
            dx = self.grid.dx
            dy = self.grid.dy
            max_cols = self.grid.number_of_node_columns - 2
            max_rows = self.grid.number_of_node_rows - 2
            max_dims_array = numpy.array([max_cols,max_rows])
            left_bottom = ((center_array - eff_radius)//dx).astype(int) + 1 #the (leftmost, bottommost) node row/col included within the footprint
            right_top = ((center_array + eff_radius)//dx).astype(int)
            right_top_nonzero = right_top.copy()
            left_bottom_nonzero = numpy.where(left_bottom<1,1,left_bottom) #we can NEVER return IDs of the boundary nodes
            if right_top_nonzero[0]>max_cols:
                right_top_nonzero[0] = max_cols
            if right_top_nonzero[1]>max_rows:
                right_top_nonzero[1] = max_rows
            if array_in is not None:
                len_x = right_top_nonzero[0]-left_bottom_nonzero[0]+1
                len_y = right_top_nonzero[1]-left_bottom_nonzero[1]+1
                x[:len_x] = numpy.arange(right_top_nonzero[0]-left_bottom_nonzero[0]+1) + left_bottom_nonzero[0]
                y_column[:len_y] = numpy.arange(right_top_nonzero[1]-left_bottom_nonzero[1]+1).reshape((len_y,1)) + left_bottom_nonzero[1]
                self.twod_node_store[:len_y,:len_x] = x[:len_x] + y_column[:len_y]*self.grid.number_of_node_columns
                array_in[:((len_x)*(len_y))] = self.twod_node_store[:len_y,:len_x].flatten()
            else:
                x = numpy.arange(right_top_nonzero[0]-left_bottom_nonzero[0]+1) + left_bottom_nonzero[0]
                y = numpy.arange(right_top_nonzero[1]-left_bottom_nonzero[1]+1) + left_bottom_nonzero[1]
                y_column = y.reshape((y.shape[0],1))
                footprint_nodes_2dim = x + y_column*self.grid.number_of_node_columns
                center_nodes = footprint_nodes_2dim.flatten()
            flag = footprint_edge_type

            if flag != 'C' and self.looped_BCs:
                if array_in is not None:
                    yield (len_x*len_y, array_in)
                else:
                    yield center_nodes
                left_bottom_nonzero_edge = (max_dims_array + left_bottom - 1)%max_dims_array + 1
                right_top_nonzero_edge = (right_top - max_dims_array - 1)%max_dims_array + 1 #this fiddly adding and subtracting unity serves to force boundary interior nodes to appear in their original positions, not their "ghost" boundary positions in the edge grids
                if flag == 'N':
                    if array_in is not None:
                        len_x = right_top_nonzero_edge[0]-left_bottom_nonzero_edge[0]+1
                        len_y = right_top_nonzero_edge[1]
                        x[:len_x] = numpy.arange(right_top_nonzero_edge[0]-left_bottom_nonzero_edge[0]+1) + left_bottom_nonzero_edge[0]
                        y_column[:len_y] = numpy.arange(right_top_nonzero_edge[1]).reshape((len_y,1)) + 1
                        self.twod_node_store[:len_y,:len_x] = x[:len_x] + y_column[:len_y]*self.grid.number_of_node_columns
                        array_in[:((len_x)*(len_y))] = self.twod_node_store[:len_y,:len_x].flatten()
                        yield (len_x*len_y, array_in)
                        return
                    else:
                        x = numpy.arange(right_top_nonzero_edge[0]-left_bottom_nonzero_edge[0]+1) + left_bottom_nonzero_edge[0]
                        y = numpy.arange(right_top_nonzero_edge[1]) + 1
                        y_column = y.reshape((y.shape[0],1))
                        footprint_nodes_2dim = x + y_column*self.grid.number_of_node_columns
                        yield footprint_nodes_2dim.flatten()
                        return
                elif flag == 'S':
                    if array_in is not None:
                        len_x = right_top_nonzero_edge[0]-left_bottom_nonzero_edge[0]+1
                        len_y = max_rows-left_bottom_nonzero_edge[1]+1
                        x[:len_x] = numpy.arange(right_top_nonzero_edge[0]-left_bottom_nonzero_edge[0]+1) + left_bottom_nonzero_edge[0]
                        y_column[:len_y] = numpy.arange(max_rows-left_bottom_nonzero_edge[1]+1).reshape((len_y,1)) + left_bottom_nonzero_edge[1]
                        self.twod_node_store[:len_y,:len_x] = x[:len_x] + y_column[:len_y]*self.grid.number_of_node_columns
                        array_in[:((len_x)*(len_y))] = self.twod_node_store[:len_y,:len_x].flatten()
                        yield (len_x*len_y, array_in)
                        return
                    else:
                        x = numpy.arange(right_top_nonzero_edge[0]-left_bottom_nonzero_edge[0]+1) + left_bottom_nonzero_edge[0]
                        y = numpy.arange(max_rows-left_bottom_nonzero_edge[1]+1) + left_bottom_nonzero_edge[1]
                        y_column = y.reshape((y.shape[0],1))
                        yield footprint_nodes_2dim.flatten()
                        return
                elif flag == 'E':
                    if array_in is not None:
                        len_x = right_top_nonzero_edge[0]
                        len_y = right_top_nonzero_edge[1]-left_bottom_nonzero_edge[1]+1
                        x[:len_x] = numpy.arange(right_top_nonzero_edge[0]) + 1
                        y_column[:len_y] = numpy.arange(right_top_nonzero_edge[1]-left_bottom_nonzero_edge[1]+1).reshape((len_y,1)) + left_bottom_nonzero_edge[1]
                        self.twod_node_store[:len_y,:len_x] = x[:len_x] + y_column[:len_y]*self.grid.number_of_node_columns
                        array_in[:((len_x)*(len_y))] = self.twod_node_store[:len_y,:len_x].flatten()
                        yield (len_x*len_y, array_in)
                        return
                    else:
                        x = numpy.arange(right_top_nonzero_edge[0]) + 1
                        y = numpy.arange(right_top_nonzero_edge[1]-left_bottom_nonzero_edge[1]+1) + left_bottom_nonzero_edge[1]
                        y_column = y.reshape((y.shape[0],1))
                        footprint_nodes_2dim = x + y_column*self.grid.number_of_node_columns
                        yield footprint_nodes_2dim.flatten()
                        return
                elif flag == 'W':
                    if array_in is not None:
                        len_x = max_cols-left_bottom_nonzero_edge[0]+1
                        len_y = right_top_nonzero_edge[1]-left_bottom_nonzero_edge[1]+1
                        x[:len_x] = numpy.arange(max_cols-left_bottom_nonzero_edge[0]+1) + left_bottom_nonzero_edge[0]
                        y_column[:len_y] = numpy.arange(right_top_nonzero_edge[1]-left_bottom_nonzero_edge[1]+1).reshape((len_y,1)) + left_bottom_nonzero_edge[1]
                        self.twod_node_store[:len_y,:len_x] = x[:len_x] + y_column[:len_y]*self.grid.number_of_node_columns
                        array_in[:((len_x)*(len_y))] = self.twod_node_store[:len_y,:len_x].flatten()
                        yield (len_x*len_y, array_in)
                        return
                    else:
                        x = numpy.arange(max_cols-left_bottom_nonzero_edge[0]+1) + left_bottom_nonzero_edge[0]
                        y = numpy.arange(right_top_nonzero_edge[1]-left_bottom_nonzero_edge[1]+1) + left_bottom_nonzero_edge[1]
                        y_column = y.reshape((y.shape[0],1))
                        footprint_nodes_2dim = x + y_column*self.grid.number_of_node_columns
                        yield footprint_nodes_2dim.flatten()
                        return
                elif flag == 'NW':
                    if array_in is not None:
                        len_x = max_cols-left_bottom_nonzero_edge[0]+1
                        len_y = max_rows-left_bottom_nonzero_edge[1]+1
                        x[:len_x] = numpy.arange(max_cols-left_bottom_nonzero_edge[0]+1) + left_bottom_nonzero_edge[0]
                        y_column[:len_y] = numpy.arange(max_rows-left_bottom_nonzero_edge[1]+1).reshape((len_y,1)) + left_bottom_nonzero_edge[1]
                        self.twod_node_store[:len_y,:len_x] = x[:len_x] + y_column[:len_y]*self.grid.number_of_node_columns
                        array_in[:((len_x)*(len_y))] = self.twod_node_store[:len_y,:len_x].flatten()
                        yield (len_x*len_y, array_in)
                        len_x = right_top_nonzero_edge[0]
                        len_y = right_top_nonzero_edge[1]
                        x[:len_x] = numpy.arange(right_top_nonzero_edge[0]) + 1
                        y_column[:len_y] = numpy.arange(right_top_nonzero_edge[1]).reshape((len_y,1)) + 1
                        self.twod_node_store[:len_y,:len_x] = x[:len_x] + y_column[:len_y]*self.grid.number_of_node_columns
                        array_in[:((len_x)*(len_y))] = self.twod_node_store[:len_y,:len_x].flatten()
                        yield (len_x*len_y, array_in)
                        len_x = max_cols-left_bottom_nonzero_edge[0]+1
                        len_y = right_top_nonzero_edge[1]
                        x[:len_x] = numpy.arange(max_cols-left_bottom_nonzero_edge[0]+1) + left_bottom_nonzero_edge[0]
                        y_column[:len_y] = numpy.arange(right_top_nonzero_edge[1]).reshape((len_y,1)) + 1
                        self.twod_node_store[:len_y,:len_x] = x[:len_x] + y_column[:len_y]*self.grid.number_of_node_columns
                        array_in[:((len_x)*(len_y))] = self.twod_node_store[:len_y,:len_x].flatten()
                        yield (len_x*len_y, array_in)
                        return
                    else:
                        x = numpy.arange(max_cols-left_bottom_nonzero_edge[0]+1) + left_bottom_nonzero_edge[0]
                        y = numpy.arange(max_rows-left_bottom_nonzero_edge[1]+1) + left_bottom_nonzero_edge[1]
                        y_column = y.reshape((y.shape[0],1))
                        footprint_nodes_2dim = x + y_column*self.grid.number_of_node_columns
                        west = footprint_nodes_2dim.flatten()
                        yield west
                        x = numpy.arange(right_top_nonzero_edge[0]) + 1
                        y = numpy.arange(right_top_nonzero_edge[1]) + 1
                        y_column = y.reshape((y.shape[0],1))
                        footprint_nodes_2dim = x + y_column*self.grid.number_of_node_columns
                        north = footprint_nodes_2dim.flatten()
                        yield north
                        x = numpy.arange(max_cols-left_bottom_nonzero_edge[0]+1) + left_bottom_nonzero_edge[0]
                        y = numpy.arange(right_top_nonzero_edge[1]) + 1
                        y_column = y.reshape((y.shape[0],1))
                        footprint_nodes_2dim = x + y_column*self.grid.number_of_node_columns
                        corner = footprint_nodes_2dim.flatten()
                        yield corner
                        return
                elif flag == 'NE':
                    if array_in is not None:
                        len_x = max_cols-left_bottom_nonzero_edge[0]+1
                        len_y = right_top_nonzero_edge[1]
                        x[:len_x] = numpy.arange(max_cols-left_bottom_nonzero_edge[0]+1) + left_bottom_nonzero_edge[0]
                        y_column[:len_y] = numpy.arange(right_top_nonzero_edge[1]).reshape((len_y,1)) + 1
                        self.twod_node_store[:len_y,:len_x] = x[:len_x] + y_column[:len_y]*self.grid.number_of_node_columns
                        array_in[:((len_x)*(len_y))] = self.twod_node_store[:len_y,:len_x].flatten()
                        yield (len_x*len_y, array_in)
                        len_x = right_top_nonzero_edge[0]
                        len_y = max_rows-left_bottom_nonzero_edge[1]+1
                        x[:len_x] = numpy.arange(right_top_nonzero_edge[0]) + 1
                        y_column[:len_y] = numpy.arange(max_rows-left_bottom_nonzero_edge[1]+1).reshape((len_y,1)) + left_bottom_nonzero_edge[1]
                        self.twod_node_store[:len_y,:len_x] = x[:len_x] + y_column[:len_y]*self.grid.number_of_node_columns
                        array_in[:((len_x)*(len_y))] = self.twod_node_store[:len_y,:len_x].flatten()
                        yield (len_x*len_y, array_in)
                        len_x = right_top_nonzero_edge[0]
                        len_y = right_top_nonzero_edge[1]
                        x[:len_x] = numpy.arange(right_top_nonzero_edge[0]) + 1
                        y_column[:len_y] = numpy.arange(right_top_nonzero_edge[1]).reshape((len_y,1)) + 1
                        self.twod_node_store[:len_y,:len_x] = x[:len_x] + y_column[:len_y]*self.grid.number_of_node_columns
                        array_in[:((len_x)*(len_y))] = self.twod_node_store[:len_y,:len_x].flatten()
                        yield (len_x*len_y, array_in)
                        return
                    else:
                        x = numpy.arange(max_cols-left_bottom_nonzero_edge[0]+1) + left_bottom_nonzero_edge[0]
                        y = numpy.arange(right_top_nonzero_edge[1]) + 1
                        y_column = y.reshape((y.shape[0],1))
                        footprint_nodes_2dim = x + y_column*self.grid.number_of_node_columns
                        north = footprint_nodes_2dim.flatten()
                        yield north
                        x = numpy.arange(right_top_nonzero_edge[0]) + 1
                        y = numpy.arange(max_rows-left_bottom_nonzero_edge[1]+1) + left_bottom_nonzero_edge[1]
                        y_column = y.reshape((y.shape[0],1))
                        footprint_nodes_2dim = x + y_column*self.grid.number_of_node_columns
                        east = footprint_nodes_2dim.flatten()
                        yield east
                        x = numpy.arange(right_top_nonzero_edge[0]) + 1
                        y = numpy.arange(right_top_nonzero_edge[1]) + 1
                        y_column = y.reshape((y.shape[0],1))
                        footprint_nodes_2dim = x + y_column*self.grid.number_of_node_columns
                        corner = footprint_nodes_2dim.flatten()
                        yield corner
                        return
                elif flag == 'SE':
                    if array_in is not None:
                        len_x = right_top_nonzero_edge[0]
                        len_y = right_top_nonzero_edge[1]
                        x[:len_x] = numpy.arange(right_top_nonzero_edge[0]) + 1
                        y_column[:len_y] = numpy.arange(right_top_nonzero_edge[1]).reshape((len_y,1)) + 1
                        self.twod_node_store[:len_y,:len_x] = x[:len_x] + y_column[:len_y]*self.grid.number_of_node_columns
                        array_in[:((len_x)*(len_y))] = self.twod_node_store[:len_y,:len_x].flatten()
                        yield (len_x*len_y, array_in)
                        len_x = max_cols-left_bottom_nonzero_edge[0]+1
                        len_y = max_rows-left_bottom_nonzero_edge[1]+1
                        x[:len_x] = numpy.arange(max_cols-left_bottom_nonzero_edge[0]+1) + left_bottom_nonzero_edge[0]
                        y_column[:len_y] = numpy.arange(max_rows-left_bottom_nonzero_edge[1]+1).reshape((len_y,1)) + left_bottom_nonzero_edge[1]
                        self.twod_node_store[:len_y,:len_x] = x[:len_x] + y_column[:len_y]*self.grid.number_of_node_columns
                        array_in[:((len_x)*(len_y))] = self.twod_node_store[:len_y,:len_x].flatten()
                        yield (len_x*len_y, array_in)
                        len_x = right_top_nonzero_edge[0]
                        len_y = max_rows-left_bottom_nonzero_edge[1]+1
                        x[:len_x] = numpy.arange(right_top_nonzero_edge[0]) + 1
                        y_column[:len_y] = numpy.arange(max_rows-left_bottom_nonzero_edge[1]+1).reshape((len_y,1)) + left_bottom_nonzero_edge[1]
                        self.twod_node_store[:len_y,:len_x] = x[:len_x] + y_column[:len_y]*self.grid.number_of_node_columns
                        array_in[:((len_x)*(len_y))] = self.twod_node_store[:len_y,:len_x].flatten()
                        yield (len_x*len_y, array_in)
                        return
                    else:
                        x = numpy.arange(right_top_nonzero_edge[0]) + 1
                        y = numpy.arange(right_top_nonzero_edge[1]) + 1
                        y_column = y.reshape((y.shape[0],1))
                        footprint_nodes_2dim = x + y_column*self.grid.number_of_node_columns
                        east = footprint_nodes_2dim.flatten()
                        yield east
                        x = numpy.arange(max_cols-left_bottom_nonzero_edge[0]+1) + left_bottom_nonzero_edge[0]
                        y = numpy.arange(max_rows-left_bottom_nonzero_edge[1]+1) + left_bottom_nonzero_edge[1]
                        y_column = y.reshape((y.shape[0],1))
                        footprint_nodes_2dim = x + y_column*self.grid.number_of_node_columns
                        south = footprint_nodes_2dim.flatten()
                        yield south
                        x = numpy.arange(right_top_nonzero_edge[0]) + 1
                        y = numpy.arange(max_rows-left_bottom_nonzero_edge[1]+1) + left_bottom_nonzero_edge[1]
                        y_column = y.reshape((y.shape[0],1))
                        footprint_nodes_2dim = x + y_column*self.grid.number_of_node_columns
                        corner = footprint_nodes_2dim.flatten()
                        yield corner
                        return
                elif flag == 'SW':
                    if array_in is not None:
                        len_x = right_top_nonzero_edge[0]
                        len_y = max_rows-left_bottom_nonzero_edge[1]+1
                        x[:len_x] = numpy.arange(right_top_nonzero_edge[0]) + 1
                        y_column[:len_y] = numpy.arange(max_rows-left_bottom_nonzero_edge[1]+1).reshape((len_y,1)) + left_bottom_nonzero_edge[1]
                        self.twod_node_store[:len_y,:len_x] = x[:len_x] + y_column[:len_y]*self.grid.number_of_node_columns
                        array_in[:((len_x)*(len_y))] = self.twod_node_store[:len_y,:len_x].flatten()
                        yield (len_x*len_y, array_in)
                        len_x = max_cols-left_bottom_nonzero_edge[0]+1
                        len_y = right_top_nonzero_edge[1]
                        x[:len_x] = numpy.arange(max_cols-left_bottom_nonzero_edge[0]+1) + left_bottom_nonzero_edge[0]
                        y_column[:len_y] = numpy.arange(right_top_nonzero_edge[1]).reshape((len_y,1)) + 1
                        self.twod_node_store[:len_y,:len_x] = x[:len_x] + y_column[:len_y]*self.grid.number_of_node_columns
                        array_in[:((len_x)*(len_y))] = self.twod_node_store[:len_y,:len_x].flatten()
                        yield (len_x*len_y, array_in)
                        len_x = max_cols-left_bottom_nonzero_edge[0]+1
                        len_y = max_rows-left_bottom_nonzero_edge[1]+1
                        x[:len_x] = numpy.arange(max_cols-left_bottom_nonzero_edge[0]+1) + left_bottom_nonzero_edge[0]
                        y_column[:len_y] = numpy.arange(max_rows-left_bottom_nonzero_edge[1]+1).reshape((len_y,1)) + left_bottom_nonzero_edge[1]
                        self.twod_node_store[:len_y,:len_x] = x[:len_x] + y_column[:len_y]*self.grid.number_of_node_columns
                        array_in[:((len_x)*(len_y))] = self.twod_node_store[:len_y,:len_x].flatten()
                        yield (len_x*len_y, array_in)
                        return
                    else:
                        x = numpy.arange(right_top_nonzero_edge[0]) + 1
                        y = numpy.arange(max_rows-left_bottom_nonzero_edge[1]+1) + left_bottom_nonzero_edge[1]
                        y_column = y.reshape((y.shape[0],1))
                        footprint_nodes_2dim = x + y_column*self.grid.number_of_node_columns
                        south = footprint_nodes_2dim.flatten()
                        yield south
                        x = numpy.arange(max_cols-left_bottom_nonzero_edge[0]+1) + left_bottom_nonzero_edge[0]
                        y = numpy.arange(right_top_nonzero_edge[1]) + 1
                        y_column = y.reshape((y.shape[0],1))
                        footprint_nodes_2dim = x + y_column*self.grid.number_of_node_columns
                        west = footprint_nodes_2dim.flatten()
                        yield west
                        x = numpy.arange(max_cols-left_bottom_nonzero_edge[0]+1) + left_bottom_nonzero_edge[0]
                        y = numpy.arange(max_rows-left_bottom_nonzero_edge[1]+1) + left_bottom_nonzero_edge[1]
                        y_column = y.reshape((y.shape[0],1))
                        footprint_nodes_2dim = x + y_column*self.grid.number_of_node_columns
                        corner = footprint_nodes_2dim.flatten()
                        yield corner
                        return
                else:
                    raise IndexError('Corner flag not set correctly!')
            elif flag == 'C' and self.looped_BCs:
                if array_in is not None:
                    yield (len_x*len_y, array_in)
                    return
                else:
                    yield center_nodes
                    return
            else:
                if array_in is not None:
                    yield (len_x*len_y, array_in)
                    return
                else:
                    yield center_nodes
                    return


    def center_offset_list_for_looped_BCs(self, center_tuple, flag_from_footprint_edge_type, whole_grid_repeats_from_fet):
        '''
        This helper method takes the output from the function self.footprint_
        edge_type() and returns a list of center tuples paired with the arrays
        in that output so that iterating on those arrays is easier.
        If the output flag is 'I' (a big footprint), it instead returns a
        number_of_tiles long list (e.g., 9, 25, 49...) of the relevant
        offsets, in effective ID order (i.e., from bottom left, working across
        line-by-line), as (x, y) tuples.

        It needs to be supplied with the second output from self.create_square
        _footprint() also.

        This method is a generator.
        '''
        grid_x = self.grid.get_grid_xdimension()-self.grid.dx #as the edge nodes are looped!
        grid_y = self.grid.get_grid_ydimension()-self.grid.dy
        assert type(center_tuple) == tuple

        if flag_from_footprint_edge_type == 'I':
            assert type(whole_grid_repeats_from_fet) == int
            for i in xrange(2*whole_grid_repeats_from_fet+1):
                x_offset = center_tuple[0] + (i-whole_grid_repeats_from_fet) * grid_x
                for j in xrange(2*whole_grid_repeats_from_fet+1):
                    y_offset = center_tuple[1] + (j-whole_grid_repeats_from_fet) * grid_y
                    yield (x_offset, y_offset)
        else:
            yield center_tuple
            if flag_from_footprint_edge_type == 'C' or flag_from_footprint_edge_type == 'X':
                return
            elif flag_from_footprint_edge_type == 'N':
                yield (center_tuple[0],center_tuple[1]-grid_y)
                return
            elif flag_from_footprint_edge_type == 'S':
                yield (center_tuple[0],center_tuple[1]+grid_y)
                return
            elif flag_from_footprint_edge_type == 'E':
                yield (center_tuple[0]-grid_x,center_tuple[1])
                return
            elif flag_from_footprint_edge_type == 'W':
                yield (center_tuple[0]+grid_x,center_tuple[1])
                return
            elif flag_from_footprint_edge_type == 'NE':
                yield (center_tuple[0],center_tuple[1]-grid_y)
                yield (center_tuple[0]-grid_x,center_tuple[1])
                yield (center_tuple[0]-grid_x,center_tuple[1]-grid_y)
                return
            elif flag_from_footprint_edge_type == 'SE':
                yield (center_tuple[0]-grid_x,center_tuple[1])
                yield (center_tuple[0],center_tuple[1]+grid_y)
                yield (center_tuple[0]-grid_x,center_tuple[1]+grid_y)
                return
            elif flag_from_footprint_edge_type == 'SW':
                yield (center_tuple[0],center_tuple[1]+grid_y)
                yield (center_tuple[0]+grid_x,center_tuple[1])
                yield (center_tuple[0]+grid_x,center_tuple[1]+grid_y)
                return
            elif flag_from_footprint_edge_type == 'NW':
                yield (center_tuple[0]+grid_x,center_tuple[1])
                yield (center_tuple[0],center_tuple[1]-grid_y)
                yield (center_tuple[0]+grid_x,center_tuple[1]-grid_y)
                return
            else:
                raise IndexError('boundary type not recognised!')


    def correct_for_slope(self):
        '''
        This method uses an empirically observed correlation between surface
        slope and mass_balance to correct the mass balance in the impact back
        towards zero.
        The calibration is performed only for vertical impacts. A second,
        subsequent calibration is needed to correct for impact angle losses (a
        different geometrical effect).
        It returns 1/(mass_bal+1), which is the multiplier to use on
        '''
        coeffs = numpy.array([  1.90725634e+03,  -2.22353553e+03,   9.94701549e+02,
                               -2.04517399e+02,   2.04931537e+01,  -5.84367364e-01,
                               -3.53244682e-01])
        powers = numpy.array([6.,5.,4.,3.,2.,1.,0.])
        synthetic_mass_balance = coeffs*self._surface_slope**powers
        return 1./(numpy.sum(synthetic_mass_balance) + 1.)
        #return 1.


    def correct_for_angle_to_vertical(self):

        coeffs = numpy.array([ -1.52913024e-11,   3.86656244e-09,  -3.81596690e-07,
                                1.76792349e-05,  -4.04787246e-04,   3.04029462e-03,
                                4.28318928e-02])
        powers = numpy.array([6.,5.,4.,3.,2.,1.,0.])
        synthetic_mass_balance = coeffs*(self._angle_to_vertical*180./numpy.pi)**powers
        return 1./(numpy.sum(synthetic_mass_balance) + 1.)
        #return 1.

    def correct_for_crater_size(self):
        """
        NB: the largest crater in the size calibration was 320m in radius.
        Larger craters may have poor mass balances. Though, the mass balance
        problems related to size become much less severe as the craters get
        larger.
        """
        coeffs = numpy.array([ -2.41642810e+04,   2.25573864e+04,  -8.13732097e+03,
                                1.44581478e+03,  -1.32719840e+02,   5.97063876e+00,
                               -4.17789847e-02])
        powers = numpy.array([6.,5.,4.,3.,2.,1.,0.])
        synthetic_mass_balance = coeffs*self._radius**powers
        return 1./(numpy.sum(synthetic_mass_balance) + 1.)
        #return 1.


    #@profile
    def excavate_a_crater_noangle(self, grid):
        '''
        This method executes the most of the other methods of this crater
        class, and makes the geomorphic changes to a mesh associated with a
        single bolide impact with randomized properties. It receives and works
        on the data fields attached to the model grid.
        This version calls the _no_angle() method above, and thus produces
        always radially symmetric distributions.
        ***This is one of the primary interface method of this class.***
        '''
        self.grid = grid
        self.elev = grid.at_node['topographic__elevation']
        self.draw_new_parameters()
        #These get updated in set_crater_mean_slope_v3()
        self.closest_node_index = grid.find_nearest_node((self._xcoord, self._ycoord))
        self.closest_node_elev = self.elev[self.closest_node_index]
        self.check_coords_and_angles_for_grazing()
        self.set_crater_mean_slope_v3()
        self.set_depth_from_size()
        self.set_crater_volume()

        if self._radius !=1e-6:
            if numpy.isnan(self._surface_slope):
                six.print_('Surface slope is not defined for this crater! Is it too big? Crater will not be drawn.')
            else:
                self.set_elev_change_only_beneath_footprint_no_angular_BAND_AID()
            #print 'Impactor angle to ground normal: ', self.impactor_angle_to_surface_normal
            six.print_('Mass balance in impact: ', self.mass_balance_in_impact)
            #print '*****'
            #Record the data:
            #Is this making copies, or just by reference? Check output. Should be copies, as these are floats, not more complex objects.
            self.impact_property_dict = {'x': self._xcoord, 'y': self._ycoord, 'r': self._radius, 'volume': self._cavity_volume, 'surface_slope': self._surface_slope, 'normal_angle': self.impactor_angle_to_surface_normal, 'impact_az': self._azimuth_of_travel, 'ejecta_az': self.ejecta_azimuth, 'mass_balance': self.mass_balance_in_impact, 'redug_crater': self.cheater_flag}
        else:
            self.impact_property_dict = {'x': -1., 'y': -1., 'r': -1., 'volume': -1., 'surface_slope': -1., 'normal_angle': -1., 'impact_az': -1., 'ejecta_az': -1., 'mass_balance': -1., 'redug_crater': -1}
        return self.grid


    def excavate_a_crater_furbish(self, grid, single_process=False):
        '''
        This method executes the most of the other methods of this crater
        class, and makes the geomorphic changes to a mesh associated with a
        single bolide impact with randomized properties. It receives and works
        on the data fields attached to the model grid.
        This version implements the full, angle dependent versions of the
        methods, following Furbish et al., 2007.
        ***This is one of the primary interface method of this class.***
        '''
        self.grid = grid
        self.elev = grid.at_node['topographic__elevation']
        self.draw_new_parameters()
        #These get updated in set_crater_mean_slope_v3()
        self.closest_node_index = grid.find_nearest_node((self._xcoord, self._ycoord))
        self.closest_node_elev = self.elev[self.closest_node_index]
        self.check_coords_and_angles_for_grazing()
        self.set_crater_mean_slope_v3()
        self.set_depth_from_size()
        self.set_crater_volume()
        if not single_process:
            self.pre_impact_elev = self.elev.copy()
        #print 'Surface slope: ', self._surface_slope
        #print 'Dip dir: ', self._surface_dip_direction

        if self._radius !=1e-6:
            if numpy.isnan(self._surface_slope):
                six.print_('Surface slope is not defined for this crater! Is it too big? Crater will not be drawn.')
            else:
                self.set_elev_change_only_beneath_footprint_BAND_AID_memory_save()
                six.print_('Rim thickness: ', self.rim_thickness)
            #print 'Impactor angle to ground normal: ', self.impactor_angle_to_surface_normal
            six.print_('Mass balance in impact: ', self.mass_balance_in_impact)
            if self.mass_balance_in_impact<-0.9 or numpy.isnan(self.mass_balance_in_impact):
                six.print_('radius: ', self._radius)
                six.print_('location: ', self._xcoord, self._ycoord)
                six.print_('surface dip dir: ', self._surface_dip_direction)
                six.print_('surface slope: ', self._surface_slope)
                six.print_('travel azimuth: ', self._azimuth_of_travel)
                six.print_('angle to normal: ', self.impactor_angle_to_surface_normal)
                six.print_('Cavity volume: ', self._cavity_volume)
                #print 'Excavated volume: ', self.excavated_volume
                six.print_('Rim thickness: ', self.rim_thickness)
                #pylab.imshow(self.grid.node_vector_to_raster(self.elev))
                #pylab.colorbar()
                #pylab.show()
            six.print_('*****')
            #Record the data:
            #Is this making copies, or just by reference? Check output. Should be copies, as these are floats, not more complex objects.
        #    self.impact_property_dict = {'x': self._xcoord, 'y': self._ycoord, 'r': self._radius, 'volume': self._cavity_volume, 'surface_slope': self._surface_slope, 'normal_angle': self.impactor_angle_to_surface_normal, 'impact_az': self._azimuth_of_travel, 'ejecta_az': self.ejecta_azimuth, 'mass_balance': self.mass_balance_in_impact, 'redug_crater': self.cheater_flag}
        #else:
        #    self.impact_property_dict = {'x': -1., 'y': -1., 'r': -1., 'volume': -1., 'surface_slope': -1., 'normal_angle': -1., 'impact_az': -1., 'ejecta_az': -1., 'mass_balance': -1., 'redug_crater': -1}
        return self.grid

    def crawl_roughness(self, window_dimension, accelerate=True):
        '''
        Takes a dimension, D, then passes a DxD window over the elevation grid.
        Returns the mean and standard error (SD/sqrt(n)) of the local relief in
        all windows of that size possible on the grid.
        self.elev_r, the rasterized form of the grid, must already exist.
        '''
        #for clarity, a shortcut:
        D = window_dimension
        #make the storing list:
        reliefs = []
        #iterate in x, then y dimensions
        if accelerate and (D/min([self.grid.shape[0], self.grid.shape[1]])<0.1):
            for i in xrange(0,self.grid.shape[1]-D, D):
                for j in xrange(0,self.grid.shape[0]-D, D):
                    min_elev = numpy.amin(self.elev_r[j:(j+D),i:(i+D)])
                    max_elev = numpy.amax(self.elev_r[j:(j+D),i:(i+D)])
                    reliefs.append(max_elev-min_elev)
        else:
            #this will take significantly longer, but does not subsample the grid
            for i in xrange(self.grid.shape[1]-D):
                for j in xrange(self.grid.shape[0]-D):
                    min_elev = numpy.amin(self.elev_r[j:(j+D),i:(i+D)])
                    max_elev = numpy.amax(self.elev_r[j:(j+D),i:(i+D)])
                    reliefs.append(max_elev-min_elev)
        #calc the stats
        reliefs = numpy.array(reliefs)
        mean_relief = numpy.mean(reliefs)
        stderr_relief = numpy.std(reliefs)/numpy.sqrt(len(reliefs))

        return mean_relief, stderr_relief


    def calculate_scale_roughness_dependence(self, elevations, min_window_size=2, max_window_size=100, step=5):
        self.elev = elevations
        #create the 2D grid
        self.elev_r = self.grid.node_vector_to_raster(elevations)
        #max_window_size = min([self.grid.shape[0], self.grid.shape[1]])
        assert max_window_size < min([self.grid.shape[0], self.grid.shape[1]])
        #make the storage items
        window_sizes = numpy.arange(min_window_size,(max_window_size+1), step, dtype=int)
        means = []
        stderrs = []
        for i in window_sizes:
            six.print_('calculating window size ', i)
            mean_relief, stderr_relief = self.crawl_roughness(i)
            means.append(mean_relief)
            stderrs.append(stderr_relief)

        return window_sizes, numpy.array(means), numpy.array(stderrs)


    def calculate_smoothed_ffts(self, elevations, smoothing_window=10000):
        '''
        Takes the elevations on the grid, and returns the real part of the fast
        Fourier transform of those elevations, smoothed at the supplied scale.
        Smoothing is applied on a log scale, but the returned values are both in
        their natural forms (i.e., not logged).

        Returns the frequency domain values, and the fft real part at each of
        those frequencies.
        '''
        self.elev = elevations
        fft = numpy.log(numpy.fft.fft(elevations))
        fftfreq = numpy.fft.fftfreq(len(elevations))
        fftreal = fft.real[numpy.logical_not(numpy.isnan(fftfreq))]
        fftfreq = fftfreq[numpy.logical_not(numpy.isnan(fftfreq))]
        fftfreq_nonan = fftfreq[numpy.logical_not(numpy.isnan(fftreal))]
        fftreal_nonan = fftreal[numpy.logical_not(numpy.isnan(fftreal))]
        fftreal_series = pd.Series(fftreal_nonan)
        fftreal_smoothed = pd.rolling_mean(fftreal_series, smoothing_window)

        return fftfreq_nonan, numpy.exp(fftreal_smoothed)

#differences = []
#for i in xrange(len(scale)):
#    difference = []
#    for j in xrange(len(means)-1):
#        difference.append(means[j+1][i]-means[j][i])
#    differences.append(difference)


    @property
    def crater_radius(self):
        return self._radius

    @property
    def impact_xy_location(self):
        return (self._xcoord, self._ycoord)

    @property
    def surface_slope_beneath_crater(self):
        return self._surface_slope

    @property
    def surface_dip_direction_beneath_crater(self):
        self._surface_dip_direction

    @property
    def impactor_travel_azimuth(self):
        return self._azimuth_of_travel

    @property
    def impact_angle_to_normal(self):
        return self.impactor_angle_to_surface_normal

    @property
    def cavity_volume(self):
        return self._cavity_volume

    @property
    def ejecta_direction_azimuth(self):
        return self.ejecta_azimuth

    @property
    def mass_balance(self):
        '''
        Mass balance in the impact. Negative means mass loss in the impact.
        '''
        return self.mass_balance_in_impact

    @property
    def slope_insensitive(self):
        return self.cheater_flag

