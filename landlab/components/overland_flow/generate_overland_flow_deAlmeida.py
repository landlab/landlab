"""Landlab component that simulates overland flow.

This component simulates overland flow using the 2-D numerical model of
shallow-water flow over topography using the de Almeida et al., 2012
algorithm for storage-cell inundation modeling.

.. codeauthor:: Jordan Adams

Examples
--------
>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> from landlab.components.overland_flow import OverlandFlow

Create a grid on which to calculate overland flow.

>>> grid = RasterModelGrid((4, 5))

The grid will need some data to provide the overland flow component. To
check the names of the fields that provide input to the overland flow
component use the *input_var_names* class property.

>>> OverlandFlow.input_var_names
('water__depth', 'topographic__elevation')

Create fields of data for each of these input variables.

>>> grid.at_node['topographic__elevation'] = np.array([
...     0., 0., 0., 0., 0.,
...     1., 1., 1., 1., 1.,
...     2., 2., 2., 2., 2.,
...     3., 3., 3., 3., 3.])
>>> grid.at_node['water__depth'] = np.array([
...     0. , 0. , 0. , 0. , 0. ,
...     0. , 0. , 0. , 0. , 0. ,
...     0. , 0. , 0. , 0. , 0. ,
...     0.1, 0.1, 0.1, 0.1, 0.1])

Instantiate the `OverlandFlow` component to work on this grid, and run it.

>>> of = OverlandFlow(grid, steep_slopes=True)
>>> of.overland_flow()

After calculating the overland flow, new fields have been added to the
grid. Use the *output_var_names* property to see the names of the fields that
have been changed.

>>> of.output_var_names
('water__depth', 'water__discharge', 'water_surface__gradient')

The `water__depth` field is defined at nodes.

>>> of.var_loc('water__depth')
'node'
>>> grid.at_node['water__depth'] # doctest: +NORMALIZE_WHITESPACE
array([  1.00000000e-05,   1.00000000e-05,   1.00000000e-05,
         1.00000000e-05,   1.00000000e-05,   1.00000000e-05,
         1.00000000e-05,   1.00000000e-05,   1.00000000e-05,
         1.00000000e-05,   1.00000000e-05,   2.00100000e-02,
         2.00100000e-02,   2.00100000e-02,   1.00000000e-05,
         1.00010000e-01,   1.00010000e-01,   1.00010000e-01,
         1.00010000e-01,   1.00010000e-01])

The `water__discharge` field is defined at links. Because our initial
topography was a dipping plane, there is no water discharge in the horizontal
direction, only toward the bottom of the grid.

>>> of.var_loc('water__discharge')
'link'
>>> q = grid.at_link['water__discharge'] # doctest: +NORMALIZE_WHITESPACE
>>> np.all(q[grid.horizontal_links] == 0.)
True
>>> np.all(q[grid.vertical_links] <= 0.)
True

The *water_surface__gradient* is also defined at links.

>>> of.var_loc('water_surface__gradient')
'link'
>>> grid.at_link['water_surface__gradient'] # doctest: +NORMALIZE_WHITESPACE
array([ 0. ,  0. ,  0. ,  0. ,
        0. ,  1. ,  1. ,  1. ,  0. ,
        0. ,  0. ,  0. ,  0. ,
        0. ,  1. ,  1. ,  1. ,  0. ,
        0. ,  0. ,  0. ,  0. ,
        0. ,  1.1,  1.1,  1.1,  0. ,
        0. ,  0. ,  0. ,  0. ])
"""
from landlab import Component, FieldError
import numpy as np
from landlab.grid.structured_quad import links
from landlab.utils.decorators import use_file_name_or_kwds


_SEVEN_OVER_THREE = 7.0 / 3.0


class OverlandFlow(Component):

    """Simulate overland flow using de Almeida approximations.

    Landlab component that simulates overland flow using the de Almeida
    et al., 2012 approximations of the 1D shallow water equations to be used
    for 2D flood inundation modeling.

    This component calculates discharge, depth and shear stress after some
    precipitation event across any raster grid. Default input file is named
    "overland_flow_input.txt' and is contained in the
    landlab.components.overland_flow folder.

    Parameters
    ----------
    grid : RasterModelGrid
        A landlab grid.
    h_init : float, optional
        Thicknes of initial thin layer of water to prevent divide by zero
        errors (m).
    alpha : float, optional
        Time step coeffcient, described in Bates et al., 2010 and
        de Almeida et al., 2012.
    mannings_n : float, optional
        Manning's roughness coefficient.
    g : float, optional
        Acceleration due to gravity (m/s^2).
    theta : float, optional
        Weighting factor from de Almeida et al., 2012.
    rainfall_intensity : float, optional
        Rainfall intensity.
    """

    _name = 'OverlandFlow'

    _input_var_names = (
        'water__depth',
        'topographic__elevation',
    )

    _output_var_names = (
        'water__depth',
        'water__discharge',
        'water_surface__gradient',
    )

    _var_units = {
        'water__depth': 'm',
        'water__discharge': 'm3/s',
        'topographic__elevation': 'm',
        'water_surface__gradient': '-',
    }

    _var_mapping = {
        'water__depth': 'node',
        'topographic__elevtation': 'node',
        'water__discharge': 'link',
        'water_surface__gradient': 'link',
    }

    _var_doc = {
        'water__depth': 'The depth of water at each node.',
        'topographic__elevtation': 'The land surface elevation.',
        'water__discharge': 'The discharge of water on active links.',
        'water_surface__gradient': 'Downstream gradient of the water surface.',
    }

    @use_file_name_or_kwds
    def __init__(self, grid, use_fixed_links=False, h_init=0.00001, alpha=0.7,
                 mannings_n=0.03, g=9.81, theta=0.8, rainfall_intensity=0.0,
                 steep_slopes = False, **kwds):
        """Create a overland flow component.

        Parameters
        ----------
        grid : RasterModelGrid
            A landlab grid.
        h_init : float, optional
            Thicknes of initial thin layer of water to prevent divide by zero
            errors (m).
        alpha : float, optional
            Time step coeffcient, described in Bates et al., 2010 and
            de Almeida et al., 2012.
        mannings_n : float, optional
            Manning's roughness coefficient.
        g : float, optional
            Acceleration due to gravity (m/s^2).
        theta : float, optional
            Weighting factor from de Almeida et al., 2012.
        rainfall_intensity : float, optional
            Rainfall intensity.
        """
        super(OverlandFlow, self).__init__(grid, **kwds)

        # First we copy our grid
        self._grid = grid

        self.h_init = h_init
        self.alpha = alpha
        self.mannings_n = mannings_n
        self.g = g
        self.theta = theta
        self.rainfall_intensity = rainfall_intensity
        self.steep_slopes = steep_slopes

        # Now setting up fields at the links...
        # For water discharge
        try:
            self.water__discharge = grid.add_zeros(
                'water__discharge', at='link',
                units=self._var_units['water__discharge'])
        except FieldError:
            # Field was already set; still, fill it with zeros
            self.water__discharge = grid.at_link['water__discharge']
            self.water__discharge.fill(0.)

        # For water depths calculated at links
        try:
            self.h_links = grid.add_zeros('water__depth', at='link',
                                          units=self._var_units[
                                              'water__depth'])
        except FieldError:
            self.h_links = grid.at_link['water__depth']
            self.h_links.fill(0.)
        self.h_links += self.h_init

        # For water surface slopes at links
        try:
            self.slope = grid.add_zeros('water_surface__gradient', at='link')
        except FieldError:
            self.slope = grid.at_link['water_surface__gradient']
            self.slope.fill(0.)

        # Start time of simulation is at 1.0 s
        self.elapsed_time = 1.0

        self.dt = None
        self.dhdt = grid.zeros()

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
        self.h = self._grid['node']['water__depth'] = (
            self._grid['node']['water__depth'] + self.h_init)

        # Assigning a class variable to the water discharge field.
        self.q = self._grid['link']['water__discharge']

        # Assiging a class variable to the elevation field.
        self.z = self._grid.at_node['topographic__elevation']

    def calc_time_step(self):
        """Calculate time step.

        Adaptive time stepper from Bates et al., 2010 and de Almeida
        et al., 2012
        """
        self.dt = (self.alpha * self._grid.dx / np.sqrt(self.g * np.amax(
            self._grid.at_node['water__depth'])))

        return self.dt

    def set_up_neighbor_arrays(self):
        """Create and initialize link neighbor arrays.

        Set up arrays of neighboring horizontal and vertical links that are
        needed for the de Almeida solution.
        """
        # First we identify all active links
        self.active_ids = links.active_link_ids(self.grid.shape,
                                                self.grid.status_at_node)

        # And then find all horizontal link IDs (Active and Inactive)
        self.horizontal_ids = links.horizontal_link_ids(self.grid.shape)

        # And make the array 1-D
        self.horizontal_ids = self.horizontal_ids.flatten()

        # Find all horizontal active link ids
        self.horizontal_active_link_ids = links.horizontal_active_link_ids(
            self.grid.shape, self.active_ids)

        # Now we repeat this process for the vertical links.
        # First find the vertical link ids and reshape it into a 1-D array
        self.vertical_ids = links.vertical_link_ids(self.grid.shape).flatten()

        # Find the *active* verical link ids
        self.vertical_active_link_ids = links.vertical_active_link_ids(
            self.grid.shape, self.active_ids)

        if self.use_fixed_links is True:
            fixed_link_ids = links.fixed_link_ids(
                self.grid.shape, self.grid.status_at_node)
            fixed_horizontal_links = links.horizontal_fixed_link_ids(
                self.grid.shape, fixed_link_ids)
            fixed_vertical_links = links.vertical_fixed_link_ids(
                self.grid.shape, fixed_link_ids)
            self.horizontal_active_link_ids = np.maximum(
                self.horizontal_active_link_ids, fixed_horizontal_links)
            self.vertical_active_link_ids = np.maximum(
                self.vertical_active_link_ids, fixed_vertical_links)
            self.active_neighbors = find_active_neighbors_for_fixed_links(
                self.grid)

        # Using the active vertical link ids we can find the north
        # and south vertical neighbors
        self.north_neighbors = links.vertical_north_link_neighbor(
            self.grid.shape, self.vertical_active_link_ids)
        self.south_neighbors = links.vertical_south_link_neighbor(
            self.grid.shape, self.vertical_active_link_ids)

        # Using the horizontal active link ids, we can find the west and
        # east neighbors
        self.west_neighbors = links.horizontal_west_link_neighbor(
            self.grid.shape, self.horizontal_active_link_ids)
        self.east_neighbors = links.horizontal_east_link_neighbor(
            self.grid.shape, self.horizontal_active_link_ids)

        # Set up arrays for discharge in the horizontal & vertical directions.
        self.q_horizontal = np.zeros(links.number_of_horizontal_links(
            self.grid.shape))
        self.q_vertical = np.zeros(links.number_of_vertical_links(
            self.grid.shape))

        # Once the neighbor arrays are set up, we change the flag to True!
        self.neighbor_flag = True

    def overland_flow(self, dt=None):
        """Generate overland flow across a grid.

        For one time step, this generates 'overland flow' across a given grid
        by calculating discharge at each node.

        Using the depth slope product, shear stress is calculated at every
        node.

        Outputs water depth, discharge and shear stress values through time at
        every point in the input grid.
        """
        if dt is None:
            self.calc_time_step()

        # First, we check and see if the neighbor arrays have been initialized
        if self.neighbor_flag is False:
            self.set_up_neighbor_arrays()

        # In case another component has added data to the fields, we just
        # reset our water depths, topographic elevations and water discharge
        # variables to the fields.
        self.h = self.grid['node']['water__depth']
        self.z = self.grid['node']['topographic__elevation']
        self.q = self.grid['link']['water__discharge']
        self.h_links = self.grid['link']['water__depth']

        # Here we identify the core nodes and active link ids for later use.
        self.core_nodes = self.grid.core_nodes
        self.active_links = self.grid.active_links

        # Per Bates et al., 2010, this solution needs to find the difference
        # between the highest water surface in the two cells and the
        # highest bed elevation
        zmax = self._grid.map_max_of_link_nodes_to_link(self.z)
        w = self.h + self.z
        wmax = self._grid.map_max_of_link_nodes_to_link(w)
        hflow = wmax[self._grid.active_links] - zmax[self._grid.active_links]

        # Insert this water depth into an array of water depths at the links.
        self.h_links[self.active_links] = hflow

        # Now we calculate the slope of the water surface elevation at
        # active links
        water_surface_gradient = (
            self.grid.calc_grad_at_link(w)[self.grid.active_links])

        # And insert these values into an array of all links
        self.slope[self.active_links] = water_surface_gradient

        # If the user chooses to set boundary links to the neighbor value, we
        # set the discharge array to have the boundary links set to their
        # neighbor value
        if self.use_fixed_links is True:
            self.q[self.grid.fixed_links] = self.q[self.active_neighbors]

        # Now we can calculate discharge. To handle links with neighbors that
        # do not exist, we will do a fancy indexing trick. Non-existent links
        # or inactive links have an index of '-1', which in Python, looks to
        # the end of a list or array. To accommodate these '-1' indices, we
        # will simply insert an value of 0.0 discharge (in units of L^2/T)
        # to the end of the discharge array.
        self.q = np.append(self.q, [0])

        horiz = self.horizontal_ids
        vert = self.vertical_ids
        # Now we calculate discharge in the horizontal direction
        self.q[horiz] = ((
            self.theta * self.q[horiz] + (1 - self.theta) /
            2 * (self.q[self.west_neighbors] + self.q[self.east_neighbors]) -
            self.g * self.h_links[self.horizontal_ids] * self.dt *
            self.slope[self.horizontal_ids]) / (
                1 + self.g * self.dt * self.mannings_n ** 2. *
                abs(self.q[horiz]) / self.h_links[self.horizontal_ids] **
                _SEVEN_OVER_THREE))

        # ... and in the vertical direction
        self.q[vert] = ((
            self.theta * self.q[vert] + (1 - self.theta) /
            2 * (self.q[self.north_neighbors] + self.q[self.south_neighbors]) -
            self.g * self.h_links[self.vertical_ids] * self.dt *
            self.slope[self.vertical_ids]) / (
                1 + self.g * self.dt * self.mannings_n ** 2. *
                abs(self.q[vert]) / self.h_links[self.vertical_ids] **
                _SEVEN_OVER_THREE))

        # Now to return the array to its original length (length of number of
        # all links), we delete the extra 0.0 value from the end of the array.
        self.q = np.delete(self.q, len(self.q) - 1)

        # And put the horizontal and vertical arrays back together, to create
        # the discharge array.
        # self.q = np.concatenate((self.q_vertical, self.q_horizontal), axis=0)

        # Updating the discharge array to have the boundary links set to
        # their neighbor
        if self.use_fixed_links is True:
            self.q[self.grid.fixed_links] = self.q[self.active_neighbors]

        if self.steep_slopes is True:
            # To prevent water from draining too fast for our time steps...
            # Our Froude number.
            Fr = 0.8
            # Our two limiting factors, the froude number and courant number.
            # Looking a calculated q to be compared to our Fr number.
            calculated_q = (self.q / self.h_links) / np.sqrt(self.g *
                                                             self.h_links)

            # Looking at our calculated q and comparing it to our Courant no.,
            q_courant = self.q * self.dt / self.grid.dx

            # Water depth split equally between four links..
            water_div_4 = self.h_links / 4.

            # IDs where water discharge is positive...
            (positive_q, ) = np.where(self.q > 0)

            # ... and negative.
            (negative_q, ) = np.where(self.q < 0)

            # Where does our calculated q exceed the Froude number? If q does
            # exceed the Froude number, we are getting supercritical flow and
            # discharge needs to be reduced to maintain stability.
            (Froude_logical, ) = np.where((calculated_q) > Fr)
            (Froude_abs_logical, ) = np.where(abs(calculated_q) > Fr)

            # Where does our calculated q exceed the Courant number and water
            # depth divided amongst 4 links? If the calculated q exceeds the
            # Courant number and is greater than the water depth divided by 4
            # links, we reduce discharge to maintain stability.
            (water_logical, ) = np.where(q_courant > water_div_4)
            (water_abs_logical, ) = np.where(abs(q_courant) > water_div_4)

            # Where are these conditions met? For positive and negative q,
            # there are specific rules to reduce q. This step finds where the
            # discharge values are positive or negative and where discharge
            # exceeds the Froude or Courant number.
            self.if_statement_1 = np.intersect1d(positive_q, Froude_logical)
            self.if_statement_2 = np.intersect1d(negative_q,
                                                 Froude_abs_logical)
            self.if_statement_3 = np.intersect1d(positive_q, water_logical)
            self.if_statement_4 = np.intersect1d(negative_q, water_abs_logical)

            # Rules 1 and 2 reduce discharge by the Froude number.
            self.q[self.if_statement_1] = (
                self.h_links[self.if_statement_1] *
                (np.sqrt(self.g * self.h_links[self.if_statement_1]) * Fr))

            self.q[self.if_statement_2] = (
                0. - (self.h_links[self.if_statement_2] *
                      np.sqrt(self.g * self.h_links[self.if_statement_2]) *
                      Fr))

            # Rules 3 and 4 reduce discharge by the Courant number.
            self.q[self.if_statement_3] = (((
                self.h_links[self.if_statement_3] * self.grid.dx) / 5.) /
                self.dt)

            self.q[self.if_statement_4] = (
                0. - (self.h_links[self.if_statement_4] * self.grid.dx / 5.) /
                self.dt)

#        # Once stability has been restored, we calculate the change in water
#        # depths on all core nodes by finding the difference between inputs
#        # (rainfall) and the inputs/outputs (flux divergence of discharge)
        self.dhdt = (self.rainfall_intensity - self.grid.calc_flux_div_at_node(
                                                                    self.q))

        # Updating our water depths...
        self.h[self.core_nodes] = (self.h[self.core_nodes] +
                                   self.dhdt[self.core_nodes] * self.dt)

#        # To prevent divide by zero errors, a minimum threshold water depth
#        # must be maintained. To reduce mass imbalances, this is set to
#        # find locations where water depth is smaller than h_init (default is
#        # is 0.001) and the new value is self.h_init * 10^-3. This was set as
#        # it showed the smallest amount of mass creation in the grid during
#        # testing.
        if self.steep_slopes is True:
            self.h[np.where(self.h < self.h_init)] = self.h_init * 10.0 ** -3

        # And reset our field values with the newest water depth and discharge.
        self.grid.at_node['water__depth'] = self.h
        self.grid.at_link['water__discharge'] = self.q


def find_active_neighbors_for_fixed_links(grid):
    """Find active link neighbors for every fixed link.

    Specialized link ID function used to ID the active links that neighbor
    fixed links in the vertical and horizontal directions.

    If the user wants to assign fixed gradients or values to the fixed
    links dynamically, this function identifies the nearest active_link
    neighbor.

    Each fixed link can either have 0 or 1 active neighbor. This function
    finds if and where that active neighbor is and stores those IDs in
    an array.

    Parameters
    ----------
    grid : RasterModelGrid
        A landlab grid.

    Returns
    -------
    ndarray of int, shape `(*, )`
        Flat array of links.


    Examples
    --------
    >>> from landlab.grid.structured_quad.links import neighbors_at_link
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.overland_flow.generate_overland_flow_deAlmeida import find_active_neighbors_for_fixed_links

    >>> from landlab import RasterModelGrid, FIXED_GRADIENT_BOUNDARY

    >>> grid = RasterModelGrid((4, 5))
    >>> grid.status_at_node[:5] = FIXED_GRADIENT_BOUNDARY
    >>> grid.status_at_node[::5] = FIXED_GRADIENT_BOUNDARY
    >>> grid.status_at_node # doctest: +NORMALIZE_WHITESPACE
    array([2, 2, 2, 2, 2,
           2, 0, 0, 0, 1,
           2, 0, 0, 0, 1,
           2, 1, 1, 1, 1], dtype=int8)

    >>> grid.fixed_links
    array([ 5,  6,  7,  9, 18])
    >>> grid.active_links
    array([10, 11, 12, 14, 15, 16, 19, 20, 21, 23, 24, 25])

    >>> find_active_neighbors_for_fixed_links(grid)
    array([14, 15, 16, 10, 19])

    >>> rmg = RasterModelGrid((4, 7))

    >>> rmg.at_node['topographic__elevation'] = rmg.zeros(at='node')
    >>> rmg.at_link['topographic__slope'] = rmg.zeros(at='link')

    >>> rmg.set_fixed_link_boundaries_at_grid_edges(True, True, True, True)
    >>> find_active_neighbors_for_fixed_links(rmg)
    array([20, 21, 22, 23, 24, 14, 17, 27, 30, 20, 21, 22, 23, 24])
    """
    neighbors = links.neighbors_at_link(grid.shape, grid.fixed_links).flat
    return neighbors[np.in1d(neighbors, grid.active_links)]
