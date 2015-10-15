""" Landlab component that simulates overland flow.

 This component simulates overland flow using
 the 2-D numerical model of shallow-water flow
 over topography using the de Almeida et al., 2012
 algorithm for storage-cell inundation modeling.

Written by Jordan Adams, based on code written by Greg Tucker.
Last updated: July 17, 2015

"""
from landlab import Component, ModelParameterDictionary
import numpy as np
import os
import warnings
from landlab.grid.structured_quad import links

class OverlandFlow(Component):
    """  Landlab component that simulates overland flow using the de Almeida
    et al., 2012 approximations of the 1D shallow water equations to be used
    for 2D flood inundation modeling.

    This component calculates discharge, depth and shear stress after some
    precipitation event across any raster grid. Default input file is named
    "overland_flow_input.txt' and is contained in the
    landlab.components.overland_flow folder.

        Inputs
        ------
        grid : Requires a RasterGridModel instance

        input_file : Contains necessary and optional inputs. If not given,
            default input file and values is used.
            - Manning's n is needed, default value of 0.01.
            - Storm duration is needed IF rainfall_duration is not passed in
                the initialization
            - Rainfall intensity is needed IF rainfall_intensity is not passed
                in the initialization
            - Model run time can be provided in initialization. If not it is
                set to the storm duration

        Constants
        ---------
        h_init : float
            Some initial depth in the channels. Default = 0.001 m
        g : float
            Gravitational acceleration, \x0crac{m}{s^2}
        alpha : float
            Non-dimensional time step factor from Bates et al., (2010)
        rho : integer
            Density of water, \x0crac{kg}{m^3}
        ten_thirds : float
            Precalculated value of \x0crac{10}{3} which is used in the implicit
            shallow water equation.


        >>> DEM_name = 'DEM_name.asc'
        >>> (rg, z) = read_esri_ascii(DEM_name) # doctest: +SKIP
        >>> of = OverlandFlow(rg) # doctest: +SKIP

    """
    _name = 'OverlandFlow'

    _input_var_names = set(['water_depth', 'topographic__elevation'])

    _output_var_names = set(['water_depth',
     'water_discharge',
     'shear_stress',
     'water_discharge_at_nodes',
     'water_surface_slope_at_nodes'])

    _var_units = {'water_depth': 'm',
     'water_discharge': 'm3/s',
     'shear_stress': 'Pa',
     'water_discharge_at_nodes': 'm3/s',
     'water_surface_slope_at_nodes': 'm/m',
     'topographic__elevation': 'm'}

    _var_mapping = {'water_depth': 'node',
     'topographic__elevtation': 'node',
     'water_discharge': 'active_link',
     'shear_stress': 'node',
     'water_discharge_at_nodes': 'node',
     'water_surface_slope_at_nodes': 'node'}

    _var_mapping = {'water_depth': 'The depth of water at each node.',
     'topographic__elevtation': 'The land surface elevation.',
     'water_discharge': 'The discharge of water on active links.',
     'shear_stress': 'The calculated shear stress at each node.',
     'water_discharge_at_nodes':
         'The water discharge from surrounding links mapped onto nodes.',
     'water_surface_slope_at_nodes':
         'The slope of the water surface at each node.'}

    def __init__(self, grid, input_file = None, use_fixed_links=False, **kwds):

        super(OverlandFlow, self).__init__(grid, **kwds)

        # First we copy our grid
        self._grid = grid

        # Then, we look for a input file...
        if input_file is not None:
            inputs = ModelParameterDictionary(input_file)
        else:
            warnings.warn(
            "No input file provided! Default file and default values "
            "will be used")
            _DEFAULT_INPUT_FILE = os.path.join(os.path.dirname(__file__),
                                               'overland_flow_input.txt')
            input_file = _DEFAULT_INPUT_FILE
            inputs = ModelParameterDictionary(input_file)

        # And here we look to see what parameters are set in the input file.
        # If a parameter is not found, default values are set below.

        #This is an initial thin layer of water to prevent divide by zero errors
        try:
            self.h_init = inputs.read_float('h_init')
        except:
            self.h_init = 0.001

        # This is the time step coeffcient, described in Bates et al., 2010 and
        # de Almeida et al., 2012
        try:
            self.alpha = inputs.read_float('alpha')
        except:
            self.alpha = 0.7

        # Manning's roughness coefficient or Manning's n
        try:
            self.mannings_n = inputs.read_float('Mannings_n')
        except:
            self.mannings_n = 0.01

        # Gravitational acceleration in L/T^2
        try:
            self.g = inputs.read_float('g')
        except:
            self.g = 9.8

        # Weighting factor from de Almeida et al., 2012
        try:
            self.theta = inputs.read_float('theta')
        except:
            self.theta = 0.8

        # Rainfall intensity
        try:
            self.rainfall_intensity = inputs.read_float('rainfall_intensity')
        except:
            self.rainfall_intensity = 0.0

        # Setting up all fields found at nodes.
        for name in self._input_var_names:
            if name not in self._grid.at_node:
                self._grid.add_zeros('node', name, units=self._var_units[name])

        for name in self._output_var_names:
            if name not in self._grid.at_node:
                self._grid.add_zeros('node', name, units=self._var_units[name])

        # Now setting up fields at the links...
        # For water discharge
        self.water_discharge = grid.add_zeros(
            'link', 'water_discharge',
            units=self._var_units['water_discharge'])

        # For water depths calculated at links
        self.h_links = grid.add_zeros('link', 'water_depth',
                                      units=self._var_units['water_depth'])
        self.h_links += self.h_init

        # For water surface slopes at links
        self.slope = grid.add_zeros('link', 'water_surface_slope')

        # Pre-calculated values included for speed.
        self.ten_thirds = 10.0 / 3.0
        self.seven_over_three = 7.0 / 3.0
        self.mannings_n_squared = self.mannings_n * self.mannings_n

        # Start time of simulation is at 1.0 s
        self.elapsed_time = 1.0

        self.dt = None
        self.dhdt = grid.create_node_array_zeros()

        # When we instantiate the class we recognize that neighbors have not
        # been found. After the user either calls self.set_up_neighbor_array
        # or self.overland_flow this will be set to True. This is done so
        # that every iteration of self.overland_flow does NOT need to
        # reinitalize the neighbors and saves computation time.
        self.neighbor_flag = False

        # When looking for neighbors, we automatically ignore inactive links
        # by default. However, what about when we want to look at fixed links
        # too? By default, we ignore these, but if they are important to your
        # model and will be updated in your driver loop, they can be used by
        # setting the flag in the initialization of  the class to 'True'
        self.use_fixed_links = use_fixed_links

        # Assigning a class variable to the water depth field and adding
        # the initial thin water depth
        self.h = self._grid['node']['water_depth'] = (
            self._grid['node']['water_depth'] + self.h_init)

        # Assigning a class variable to the water discharge field.
        self.q = self._grid['link']['water_discharge']

        # Assiging a class variable to the elevation field.
        self.z = self._grid.at_node['topographic__elevation']

    def gear_time_step(self, grid):

        # Adaptive time stepper from Bates et al., 2010 and
        # de Almeida et al., 2012
        self.dt = (self.alpha * self._grid.dx /
            np.sqrt(self.g * np.amax(self._grid.at_node['water_depth'])))

        return self.dt

    def set_up_neighbor_arrays(self, grid):

        # This function gets arrays of neighboring horizontal and vertical
        # links which are needed for the de Almeida solution

        # First we identify all active links
        self.active_ids = links.active_link_ids(grid.shape, grid.status_at_node)

        # And then find all horizontal link IDs (Active and Inactive)
        self.horizontal_ids = links.horizontal_link_ids(grid.shape)

        # And make the array 1-D
        self.horizontal_ids = self.horizontal_ids.flatten()

        # Find all horizontal active link ids
        self.horizontal_active_link_ids = links.horizontal_active_link_ids(
            grid.shape, self.active_ids)

        # Now we repeat this process for the vertical links.
        # First find the vertical link ids and reshape it into a 1-D array
        self.vertical_ids = links.vertical_link_ids(grid.shape).flatten()

        # Find the *active* verical link ids
        self.vertical_active_link_ids = links.vertical_active_link_ids(
            grid.shape, self.active_ids)

        if self.use_fixed_links==True:
            fixed_link_ids = links.fixed_link_ids(
                grid.shape, grid.status_at_node)
            fixed_horizontal_links = links.horizontal_fixed_link_ids(
                grid.shape, fixed_link_ids)
            fixed_vertical_links = links.vertical_fixed_link_ids(
                grid.shape, fixed_link_ids)
            self.horizontal_active_link_ids= np.maximum(
                self.horizontal_active_link_ids, fixed_horizontal_links)
            self.vertical_active_link_ids= np.maximum(
                self.vertical_active_link_ids, fixed_vertical_links)
            self.active_neighbors = find_active_neighbors_for_fixed_links(grid)


        # Using the active vertical link ids we can find the north
        # and south vertical neighbors
        self.north_neighbors = links.vertical_north_link_neighbor(
            grid.shape, self.vertical_active_link_ids)
        self.south_neighbors = links.vertical_south_link_neighbor(
            grid.shape, self.vertical_active_link_ids)

        # Using the horizontal active link ids, we can find the west and
        # east neighbors
        self.west_neighbors = links.horizontal_west_link_neighbor(
            grid.shape, self.horizontal_active_link_ids)
        self.east_neighbors = links.horizontal_east_link_neighbor(
            grid.shape, self.horizontal_active_link_ids)

        # Set up arrays for discharge in the horizontal & vertical directions.
        self.q_horizontal = np.zeros(links.number_of_horizontal_links(
            grid.shape))
        self.q_vertical = np.zeros(links.number_of_vertical_links(
            grid.shape))
        # Once the neighbor arrays are set up, we change the flag to True!
        self.neighbor_flag = True


    def overland_flow(self, grid, dt = None, **kwds):
        """Generate overland flow across a grid.

        For one time step, this generates 'overland flow' across a given grid
        by calculating discharge at each node.

        Using the depth slope product, shear stress is calculated at every
        node.

        Outputs water depth, discharge and shear stress values through time at
        every point in the input grid.


        Parameters
        ------
        grid : Requires a RasterGridModel instance

        dt : either set when called or the fxn will do it for you.

        """

        # First, we check and see if the neighbor arrays have been initialized
        if self.neighbor_flag is False:
            self.set_up_neighbor_arrays(grid)

        # In case another component has added data to the fields, we just
        # reset our water depths, topographic elevations and water discharge
        # variables to the fields.
        self.h = self._grid['node']['water_depth']
        self.z = self._grid['node']['topographic__elevation']
        self.q = self._grid['link']['water_discharge']
        self.h_links = self._grid['link']['water_depth']

        # Here we identify the core nodes and active link ids for later use.
        self.core_nodes = self._grid.core_nodes
        self.active_links = self._grid.active_links

        # Per Bates et al., 2010, this solution needs to find the difference
        # between the highest water surface in the two cells and the
        # highest bed elevation
        zmax = self._grid.max_of_link_end_node_values(self.z)
        w = self.h + self.z
        wmax = self._grid.max_of_link_end_node_values(w)
        hflow = wmax - zmax

        # Insert this water depth into an array of water depths at the links.
        self.h_links[self.active_links] = hflow

        # Now we calculate the slope of the water surface elevation at
        # active links
        water_surface_slope = self._grid.calculate_gradients_at_active_links(w)

        # And insert these values into an array of all links
        self.slope[self.active_links] = water_surface_slope

        # If the user chooses to set boundary links to the neighbor value, we
        # set the discharge array to have the boundary links set to their
        # neighbor value
        if self.use_fixed_links == True:
            self.q[grid.fixed_links] = self.q[self.active_neighbors]

        # Now we can calculate discharge. To handle links with neighbors that
        # do not exist, we will do a fancy indexing trick. Non-existent links
        # or inactive links have an index of '-1', which in Python, looks to
        # the end of a list or array. To accommodate these '-1' indices, we
        # will simply insert an value of 0.0 discharge (in units of L^2/T)
        # to the end of the discharge array.
        self.q = np.append(self.q, [0])

        # Now we calculate discharge in the horizontal direction
        self.q_horizontal = ((self.theta * self.q_horizontal + (1 - self.theta)
            / 2 * (self.q[self.west_neighbors] + self.q[self.east_neighbors]) -
            self.g * self.h_links[self.horizontal_ids] * self.dt *
            self.slope[self.horizontal_ids]) / (1 + self.g * self.dt *
            self.mannings_n_squared * abs(self.q_horizontal) /
            self.h_links[self.horizontal_ids] ** self.seven_over_three))

        # ... and in the vertical direction
        self.q_vertical = ((self.theta * self.q_vertical + (1 - self.theta) /
            2 * (self.q[self.north_neighbors] + self.q[self.south_neighbors]) -
            self.g * self.h_links[self.vertical_ids] * self.dt *
            self.slope[self.vertical_ids]) / (1 + self.g * self.dt *
            self.mannings_n_squared * abs(self.q_vertical) /
            self.h_links[self.vertical_ids] ** self.seven_over_three))

        # Now to return the array to its original length (length of number of
        # all links), we delete the extra 0.0 value from the end of the array.
        self.q = np.delete(self.q, len(self.q) - 1)

        # And put the horizontal and vertical arrays back together, to create
        # the discharge array.
        self.q = np.concatenate((self.q_vertical, self.q_horizontal), axis=0)

        # Updating the discharge array to have the boundary links set to
        # their neighbor
        if self.use_fixed_links == True:
            self.q[grid.fixed_links] = self.q[self.active_neighbors]

        # To prevent water from draining too fast for our time steps...
        # Our Froude number.
        Fr = 0.8
        # Our two limiting factors, the froude number and courant number.
        # Looking a calculated q to be compared to our Fr number.
        calculated_q = (self.q / self.h_links) / np.sqrt(self.g * self.h_links)

        # Looking at our calculated q and comparing it to our Courant number,
        q_courant = self.q*self.dt/grid.dx

        # Water depth split equally between four links..
        water_div_4 = self.h_links/4.

        # IDs where water discharge is positive...
        (positive_q, ) = np.where(self.q>0)

        # ... and negative.
        (negative_q, ) = np.where(self.q<0)

        # Where does our calculated q exceed the Froude number? If q does
        # exceed the Froude number, we are getting supercritical flow and
        # discharge needs to be reduced to maintain stability.
        (Froude_logical, ) = np.where((calculated_q) > Fr)
        (Froude_abs_logical, ) = np.where(abs(calculated_q) > Fr)

        # Where does our calculated q exceed the Courant number and water
        # depth divided amongst 4 links? If the calculated q exceeds the
        # Courant number and is greater than the water depth divided by 4
        # links, we reduce discharge to maintain stability.
        (water_logical, ) = np.where(q_courant>water_div_4)
        (water_abs_logical, ) = np.where(abs(q_courant)>water_div_4)

        # Where are these conditions met? For positive and negative q, there
        # are specific rules to reduce q. This step finds where the discharge
        # values are positive or negative and where discharge exceeds the
        # Froude or Courant number.
        self.if_statement_1 = np.intersect1d(positive_q, Froude_logical)
        self.if_statement_2 = np.intersect1d(negative_q, Froude_abs_logical)
        self.if_statement_3 = np.intersect1d(positive_q, water_logical)
        self.if_statement_4 = np.intersect1d(negative_q, water_abs_logical)

        # Rules 1 and 2 reduce discharge by the Froude number.
        self.q[self.if_statement_1] = (self.h_links[self.if_statement_1] *
            (np.sqrt(self.g * self.h_links[self.if_statement_1]) * Fr))

        self.q[self.if_statement_2] = (0 - (self.h_links[self.if_statement_2] *
            np.sqrt(self.g * self.h_links[self.if_statement_2]) * Fr))

        # Rules 3 and 4 reduce discharge by the Courant number.
        self.q[self.if_statement_3] = (((self.h_links[self.if_statement_3] *
            grid.dx) / 5.) / self.dt)

        self.q[self.if_statement_4] = (0 - (self.h_links[self.if_statement_4] *
            grid.dx / 5.) / self.dt)

        # Once stability has been restored, we calculate the change in water
        # depths on all core nodes by finding the difference between the inputs
        # (rainfall) and the inputs/outputs (flux divergence of discharge)
        self.dhdt = (self.rainfall_intensity -
        self._grid.calculate_flux_divergence_at_nodes(self.q[self.active_links]))

        # Updating our water depths...
        self.h[self.core_nodes] = (self.h[self.core_nodes] +
            self.dhdt[self.core_nodes] * self.dt)

        # To prevent divide by zero errors, a minimum threshold water depth
        # must be maintained. To reduce mass imbalances, this is set to
        # find locations where water depth is smaller than h_init (default is
        # is 0.001) and the new value is self.h_init * 10^-3. This was set as
        # it showed the smallest amount of mass creation in the grid during
        # testing.
        self.h[np.where(self.h < self.h_init)] = (self.h_init * (10.0**-3))

        # And reset our field values with the newest water depth and discharge.
        self._grid.at_node['water_depth'] = self.h
        self._grid.at_link['water_discharge'] = self.q


    @property
    def input_var_names(self):
        return self._input_var_names

    @property
    def output_var_names(self):
        return self._output_var_names

    @property
    def var_units(self):
        return self._var_units

    @property
    def var_mapping(self):
        return self._var_mapping

def find_active_neighbors_for_fixed_links(grid):
    '''
    Specialized link ID function used to ID the active links that neighbor
    fixed links in the vertical and horizontal directions.

    If the user wants to assign fixed gradients or values to the fixed
    links dynamically, this function identifies the nearest active_link
    neighbor.

    Each fixed link can either have 0 or 1 active neighbor. This function
    finds if and where that active neighbor is and stores those IDs in
    an array.
    '''

    shape = grid.shape
    status_at_node = grid.status_at_node

    # First, we identify fixed links using node status
    fixed_links = links.fixed_link_ids(shape, status_at_node)

    # Identifying *just* fixed links IDs.
    fixed_ids_only = fixed_links[np.where(fixed_links > -1)]

    # Identifying active link IDs.
    active_links = links.active_link_ids(shape, status_at_node)

    # Identifying vertical active link IDs.
    vertical_active_links = (links.vertical_active_link_ids(shape,
                                                            active_links))

    # Identifying horizontal active link IDs.
    horizontal_active_links = (links.horizontal_active_link_ids(shape,
                                                            active_links))

    # Identifying north vertical active link IDs.
    north_vert = (links.vertical_north_link_neighbor(shape,
                                                    vertical_active_links))

    # Identifying south verical active link IDs.
    south_vert = (links.vertical_south_link_neighbor(shape,
                                                    vertical_active_links))

    # Identifying horizontal east active link IDs.
    east_hori = (links.horizontal_east_link_neighbor(shape,
                                                horizontal_active_links))

    # Identifying horizontal west active link IDs.
    west_hori = (links.horizontal_west_link_neighbor(shape,
                                                horizontal_active_links))

    # Because each fixed link can have at most 1 active neighbor, there
    # is at least one "BAD_INDEX_VALUE" link neighbor (-1). The maximum
    # ID value will be the active neighbor. This finds the N/S vertical
    # active neighbor and the E/W horizontal neighbor.
    max_vertical_neighbor = np.maximum(north_vert, south_vert)
    max_horizontal_neighbor = np.maximum(east_hori, west_hori)

    # Concatenating the vertical and horizontal arrays to get one
    # neighbor array of len(all links)
    all_active_neighbors = (np.concatenate((max_vertical_neighbor,
                                        max_horizontal_neighbor), axis=0))

    # Getting JUST the active neighbor IDs for fixed links. This
    # sets the array to a new length - that of len(fixed_links)
    all_active_neighbors = all_active_neighbors[fixed_ids_only]

    return all_active_neighbors

