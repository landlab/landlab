"""Landlab component that simulates overland flow.

This component simulates overland flow using the 2-D numerical model of
shallow-water flow over topography using the de Almeida et al., 2012
algorithm for storage-cell inundation modeling.

.. codeauthor:: Jordan Adams
"""

from landlab import Component, ModelParameterDictionary, MissingKeyError
import numpy as np
import warnings
import os
from landlab.grid.structured_quad import links


class OverlandFlow(Component):

    """Simulate overland flow using the de Almeida et al. approximations.

    Landlab component that simulates overland flow using the de Almeida
    et al., 2012 approximations of the 1D shallow water equations to be used
    for 2D flood inundation modeling.

    This component calculates discharge, depth and shear stress after some
    precipitation event across any raster grid. Default input file is named
    "overland_flow_input.txt' and is contained in the
    landlab.components.overland_flow folder.

    The :class:`OverlandFlow` component contains necessary and optional
    inputs (read from a :class:`~landlab.ModelParameterDictionary` file). If
    not given, default input file and values is used.

    - Manning's n is needed, default value of 0.01 (`MANNINGS_N`).
    - Storm duration is needed IF rainfall_duration is not passed in the
      initialization (`STORM_DURATION`).
    - Rainfall intensity is needed *if rainfall_intensity is not passed in the
      initialization* (`RAINFALL_INTENSITY`).
    - Model run time can be provided in initialization. If not it is set to
      the storm duration.

    Parameters
    ----------
    grid : RasterModelGrid
        A RasterGridModel.
    input_file : str, optional
        Name of component input file (as a
        :class:`~landlab.ModelParameterDictionary`).  If not provided,
        default values will be used.
    use_fixed_links : boolean
        Use fixed links.

    Notes
    -----
    Some constants used in the :class:`OverlandFlow` component:

    - ``h_init`` (float): Initial depth in the channels. Default = 0.001 m
    - ``g`` (float): Gravitational acceleration [m s^2]
    - ``alpha`` (float): Non-dimensional time step factor from Bates et al.,
      (2010)
    - ``rho`` (integer): Density of water, [kg m^3]
    - ``ten_thirds`` (float): Precalculated value of [10 / 3] which is
      used in the implicit shallow water equation.

    Examples
    --------
    >>> DEM_name = 'DEM_name.asc'
    >>> (rg, z) = read_esri_ascii(DEM_name) # doctest: +SKIP
    >>> of = OverlandFlow(rg) # doctest: +SKIP
    """
    _name = 'OverlandFlow'

    _input_var_names = set(['water_depth', 'topographic__elevation'])

    _output_var_names = set([
        'water_depth',
        'water_discharge',
        'shear_stress',
        'water_discharge_at_nodes',
        'water_surface_slope_at_nodes'])

    _var_units = {
        'water_depth': 'm',
        'water_discharge': 'm3/s',
        'shear_stress': 'Pa',
        'water_discharge_at_nodes': 'm3/s',
        'water_surface_slope_at_nodes': 'm/m',
        'topographic__elevation': 'm'}

    _var_mapping = {
        'water_depth': 'node',
        'topographic__elevtation': 'node',
        'water_discharge': 'active_link',
        'shear_stress': 'node',
        'water_discharge_at_nodes': 'node',
        'water_surface_slope_at_nodes': 'node'}

    _var_mapping = {
        'water_depth': 'The depth of water at each node.',
        'topographic__elevtation': 'The land surface elevation.',
        'water_discharge': 'The discharge of water on active links.',
        'shear_stress': 'The calculated shear stress at each node.',
        'water_discharge_at_nodes':
            'The water discharge from surrounding links mapped onto nodes.',
        'water_surface_slope_at_nodes':
            'The slope of the water surface at each node.'}

    def __init__(self, grid, input_file=None, use_fixed_links=False, **kwds):
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

        # This is an initial thin layer of water to prevent divide by zero
        # errors
        try:
            self.h_init = inputs.read_float('h_init')
        except MissingKeyError:
            self.h_init = 0.001

        # This is the time step coeffcient, described in Bates et al., 2010 and
        # de Almeida et al., 2012
        try:
            self.alpha = inputs.read_float('alpha')
        except MissingKeyError:
            self.alpha = 0.7

        # Manning's roughness coefficient or Manning's n
        try:
            self.mannings_n = inputs.read_float('Mannings_n')
        except MissingKeyError:
            self.mannings_n = 0.01

        # Gravitational acceleration in L/T^2
        try:
            self.g = inputs.read_float('g')
        except MissingKeyError:
            self.g = 9.8

        # Weighting factor from de Almeida et al., 2012
        try:
            self.theta = inputs.read_float('theta')
        except MissingKeyError:
            self.theta = 0.8

        # Rainfall intensity
        try:
            self.rainfall_intensity = inputs.read_float('rainfall_intensity')
        except MissingKeyError:
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

        # When we instantiate the class we recognize that neighbors have not
        # been found. After the user either calls self.set_up_neighbor_array or
        # self.overland_flow this will be set to True. This is done so that
        # every iteration of self.overland_flow does NOT need to reinitalize
        # the neighbors and saves computation time.
        self.neighbor_flag = False

        # When looking for neighbors, we automatically ignore inactive links
        # by default.  However, what about when we want to look at fixed links
        # too? By default, we ignore these, but if they are important to your
        # model and will be updated in your driver loop, they can be used by
        # setting the flag in the initialization of the class to 'True'

        self.use_fixed_links = use_fixed_links

        # Assigning a class variable to the water depth field and adding the
        # initial thin water depth
        self.h = self._grid['node']['water_depth'] = (
            self._grid['node']['water_depth'] + self.h_init)

        # Assigning a class variable to the water discharge field.
        self.q = self._grid['link']['water_discharge']

        # Assiging a class variable to the elevation field.
        self.z = self._grid.at_node['topographic__elevation']

    def gear_time_step(self, grid):
        # Adaptive time stepper from Bates et al., 2010 and de Almeida et al.,
        # 2012
        dt = (self.alpha * self._grid.dx /
              np.sqrt(self.g * np.amax(self._grid.at_node['water_depth'])))

        return dt

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

        if self.use_fixed_links:
            fixed_link_ids = links.fixed_link_ids(grid.shape, grid.status_at_node)
            fixed_horizontal_links = links.horizontal_fixed_link_ids(
                grid.shape, fixed_link_ids)
            fixed_vertical_links = links.vertical_fixed_link_ids(
                grid.shape, fixed_link_ids)
            self.horizontal_active_link_ids = np.maximum(
                self.horizontal_active_link_ids, fixed_horizontal_links)
            self.vertical_active_link_ids = np.maximum(
                self.vertical_active_link_ids, fixed_vertical_links)

        # Using the active vertical link ids we can find the north and south
        # vertical neighbors
        self.north_neighbors = links.vertical_north_link_neighbor(
            grid.shape, self.vertical_active_link_ids)
        self.south_neighbors = links.vertical_south_link_neighbor(
            grid.shape, self.vertical_active_link_ids)

        # Using the horizontal active link ids, we can find the west and east
        # neighbors
        self.west_neighbors = links.horizontal_west_link_neighbor(
            grid.shape, self.horizontal_active_link_ids)
        self.east_neighbors = links.horizontal_east_link_neighbor(
            grid.shape, self.horizontal_active_link_ids)

        # Set up arrays for discharge in the horizontal and vertical
        # directions.
        self.q_horizontal = np.zeros(
            links.number_of_horizontal_links(grid.shape))
        self.q_vertical = np.zeros(
            links.number_of_vertical_links(grid.shape))

        # Once the neighbor arrays are set up, we change the flag to True!
        self.neighbor_flag = True

    def overland_flow(self, grid, dt=None, **kwds):
        """Generate overland flow across a grid.

        For one time step, this generates 'overland flow' across a given grid
        by calculating discharge at each node.

        Using the depth slope product, shear stress is calculated at every
        node.

        Outputs water depth, discharge and shear stress values through time at
        every point in the input grid.

        Parameters
        ----------
        grid : RasterModelGrid
            A RasterModelGrid
        dt : float, optional
            Time step for overland flow.
        """
        # If no dt is provided, one will be calculated using
        # self.gear_time_step()
        if dt is None:
            dt = self.gear_time_step(grid)

        # Next, we check and see if the neighbor arrays have been initialized
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
        # between the highest water surface in the two cells and the highest
        # bed elevation
        zmax = self._grid.max_of_link_end_node_values(self.z)
        w = self.h + self.z
        wmax = self._grid.max_of_link_end_node_values(w)
        hflow = wmax - zmax

        # Insert this water depth into an array of water depths at the links.
        self.h_links[self.active_links] = hflow

        # Now we calculate the slope of the water surface elevation at active
        # links
        water_surface_slope = self._grid.calculate_gradients_at_active_links(w)

        # And insert these values into an array of all links
        self.slope[self.active_links] = water_surface_slope

        # Now we can calculate discharge. To handle links with neighbors that
        # do not exist, we will do a fancy indexing trick. Non-existent links
        # or inactive links have an index of '-1', which in Python, looks to
        # the end of a list or array. To accommodate these '-1' indices, we
        # will simply insert an value of 0.0 discharge (in units of L^2/T) to
        # the end of the discharge array.
        self.q = np.append(self.q, [0])

        # Now we calculate discharge in the horizontal direction
        self.q_horizontal = (self.theta * self.q_horizontal + (1 - self.theta) / 2 * (self.q[self.west_neighbors] + self.q[self.east_neighbors]) - self.g * self.h_links[self.horizontal_ids] * dt * self.slope[self.horizontal_ids]) / (1 + self.g * dt * self.mannings_n_squared * abs(self.q_horizontal) / self.h_links[self.horizontal_ids] ** self.seven_over_three)

        # ... and in the vertical direction
        self.q_vertical = (self.theta * self.q_vertical + (1 - self.theta) / 2 * (self.q[self.north_neighbors] + self.q[self.south_neighbors]) - self.g * self.h_links[self.vertical_ids] * dt * self.slope[self.vertical_ids]) / (1 + self.g * dt * self.mannings_n_squared * abs(self.q_vertical) / self.h_links[self.vertical_ids] ** self.seven_over_three)

        # Now to return the array to its original length (length of number of
        # all links), we delete the extra 0.0 value from the end of the array.
        self.q = np.delete(self.q, len(self.q) - 1)

        # And put the horizontal and vertical arrays back together, to create
        # the discharge array.
        self.q = np.concatenate((self.q_vertical, self.q_horizontal), axis=0)

        # Update our water depths
        dhdt = (self.rainfall_intensity -
                self._grid.calculate_flux_divergence_at_nodes(
                    self.q[self.active_links]))

        self.h[self.core_nodes] = (self.h[self.core_nodes] +
                                   dhdt[self.core_nodes] * dt)

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
