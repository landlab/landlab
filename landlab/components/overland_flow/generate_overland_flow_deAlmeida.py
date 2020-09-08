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
('surface_water__depth', 'topographic__elevation')

Create fields of data for each of these input variables.

>>> grid.at_node['topographic__elevation'] = np.array([
...     0., 0., 0., 0., 0.,
...     1., 1., 1., 1., 1.,
...     2., 2., 2., 2., 2.,
...     3., 3., 3., 3., 3.])
>>> grid.at_node['surface_water__depth'] = np.array([
...     0. , 0. , 0. , 0. , 0. ,
...     0. , 0. , 0. , 0. , 0. ,
...     0. , 0. , 0. , 0. , 0. ,
...     0.1, 0.1, 0.1, 0.1, 0.1])

Instantiate the `OverlandFlow` component to work on this grid, and run it.

>>> of = OverlandFlow(grid, steep_slopes=True)
>>> of.run_one_step()

After calculating the overland flow, new fields have been added to the
grid. Use the *output_var_names* property to see the names of the fields that
have been changed.

>>> of.output_var_names
('surface_water__depth', 'surface_water__discharge', 'water_surface__gradient')

The `surface_water__depth` field is defined at nodes.

>>> of.var_loc('surface_water__depth')
'node'
>>> grid.at_node['surface_water__depth'] # doctest: +NORMALIZE_WHITESPACE
array([  1.00000000e-05,   1.00000000e-05,   1.00000000e-05,
         1.00000000e-05,   1.00000000e-05,   1.00000000e-05,
         1.00000000e-05,   1.00000000e-05,   1.00000000e-05,
         1.00000000e-05,   1.00000000e-05,   2.00100000e-02,
         2.00100000e-02,   2.00100000e-02,   1.00000000e-05,
         1.00010000e-01,   1.00010000e-01,   1.00010000e-01,
         1.00010000e-01,   1.00010000e-01])

The `surface_water__discharge` field is defined at links. Because our initial
topography was a dipping plane, there is no water discharge in the horizontal
direction, only toward the bottom of the grid.

>>> of.var_loc('surface_water__discharge')
'link'
>>> q = grid.at_link['surface_water__discharge'] # doctest: +NORMALIZE_WHITESPACE
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
import numpy as np
import scipy.constants

from landlab import Component, FieldError

from . import _links as links

_SEVEN_OVER_THREE = 7.0 / 3.0


def _active_links_at_node(grid, *args):
    """_active_links_at_node([node_ids]) Active links of a node.

    .. note::

        This function returns links that are in *clockwise* order,
        rather than the standard *counterclockwise* ordering that
        landlab uses everywhere else.

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.
    node_ids : int or list of ints
        ID(s) of node(s) for which to find connected active links

    Returns
    -------
    (4, N) ndarray
        The ids of active links attached to grid nodes with
        *node_ids*. If *node_ids* is not given, return links for all of the
        nodes in the grid. Link ids are listed in clockwise order starting
        with the south link. Diagonal links are never returned.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.overland_flow.generate_overland_flow_deAlmeida import _active_links_at_node
    >>> grid = RasterModelGrid((3, 4))
    >>> grid.links_at_node[5]
    array([ 8, 11,  7,  4])
    >>> _active_links_at_node(grid, (5, 6))
    array([[ 4,  5],
           [ 7,  8],
           [11, 12],
           [ 8,  9]])
    >>> _active_links_at_node(grid)
    array([[-1, -1, -1, -1, -1,  4,  5, -1, -1, 11, 12, -1],
           [-1, -1, -1, -1, -1,  7,  8,  9, -1, -1, -1, -1],
           [-1,  4,  5, -1, -1, 11, 12, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1,  7,  8,  9, -1, -1, -1, -1, -1]])

    LLCATS: DEPR LINF NINF
    """
    active_links_at_node = grid.links_at_node.copy()
    active_links_at_node[grid.active_link_dirs_at_node == 0] = -1
    active_links_at_node = active_links_at_node[:, (3, 2, 1, 0)]

    if len(args) == 0:
        return active_links_at_node.T
    elif len(args) == 1:
        node_ids = np.broadcast_arrays(args[0])[0]
        return active_links_at_node[node_ids, :].T
    else:
        raise ValueError("only zero or one arguments accepted")


class OverlandFlow(Component):

    """Simulate overland flow using de Almeida approximations.

    Landlab component that simulates overland flow using the de Almeida
    et al., 2012 approximations of the 1D shallow water equations to be used
    for 2D flood inundation modeling.

    This component calculates discharge, depth and shear stress after some
    precipitation event across any raster grid. Default input file is named
    "overland_flow_input.txt' and is contained in the
    landlab.components.overland_flow folder.

    The primary method of this class is :func:`run_one_step`.

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Adams, J., Gasparini, N., Hobley, D., Tucker, G., Hutton, E., Nudurupati,
    S., Istanbulluoglu, E. (2017). The Landlab v1. 0 OverlandFlow component:
    a Python tool for computing shallow-water flow across watersheds.
    Geoscientific Model Development  10(4), 1645.
    https://dx.doi.org/10.5194/gmd-10-1645-2017

    **Additional References**

    de Almeida, G., Bates, P., Freer, J., Souvignet, M. (2012). Improving the
    stability of a simple formulation of the shallow water equations for 2-D
    flood modeling. Water Resources Research 48(5)
    https://dx.doi.org/10.1029/2011wr011570

    """

    _name = "OverlandFlow"

    _unit_agnostic = False

    _cite_as = """@article{adams2017landlab,
        title={The Landlab v1. 0 OverlandFlow component: a Python
            tool for computing shallow-water flow across watersheds},
        author={Adams, Jordan M and Gasparini, Nicole M and
            Hobley, Daniel EJ and Tucker, Gregory E and
            Hutton, Eric WH and Nudurupati, Sai S and
            Istanbulluoglu, Erkan},
        journal={Geoscientific Model Development},
        volume={10},
        number={4},
        pages={1645},
        year={2017},
        publisher={Copernicus GmbH}
        }
    """

    _info = {
        "surface_water__depth": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of water on the surface",
        },
        "surface_water__discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m3/s",
            "mapping": "link",
            "doc": "Volumetric discharge of surface water",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "water_surface__gradient": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "link",
            "doc": "Downstream gradient of the water surface.",
        },
    }

    def __init__(
        self,
        grid,
        default_fixed_links=False,
        h_init=0.00001,
        alpha=0.7,
        mannings_n=0.03,
        g=scipy.constants.g,
        theta=0.8,
        rainfall_intensity=0.0,
        steep_slopes=False,
    ):
        """Create an overland flow component.

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
            Rainfall intensity. Default is zero.
        steep_slopes : bool, optional
            Modify the algorithm to handle steeper slopes at the expense of
            speed. If model runs become unstable, consider setting to True.
        """
        super().__init__(grid)

        # First we copy our grid

        self._h_init = h_init
        self._alpha = alpha

        if isinstance(mannings_n, str):
            self._mannings_n = self._grid.at_link[mannings_n]
        else:
            self._mannings_n = mannings_n

        self._g = g
        self._theta = theta
        self.rainfall_intensity = rainfall_intensity
        self._steep_slopes = steep_slopes

        # Now setting up fields at the links...
        # For water discharge
        try:
            self._q = grid.add_zeros(
                "surface_water__discharge",
                at="link",
                units=self._info["surface_water__discharge"]["units"],
            )

        except FieldError:
            # Field was already set; still, fill it with zeros
            self._q = grid.at_link["surface_water__discharge"]
            self._q.fill(0.0)

        # For water depths calculated at links
        try:
            self._h_links = grid.add_zeros(
                "surface_water__depth",
                at="link",
                units=self._info["surface_water__depth"]["units"],
            )
        except FieldError:
            self._h_links = grid.at_link["surface_water__depth"]
            self._h_links.fill(0.0)
        self._h_links += self._h_init

        self._h = grid.at_node["surface_water__depth"]
        self._h += self._h_init

        # For water surface slopes at links
        try:
            self._water_surface_slope = grid.add_zeros(
                "water_surface__gradient", at="link"
            )
        except FieldError:
            self._water_surface_slope = grid.at_link["water_surface__gradient"]
            self._water_surface_slope.fill(0.0)

        # Start time of simulation is at 1.0 s
        self._elapsed_time = 1.0

        self._dt = None
        self._dhdt = grid.zeros()

        # When we instantiate the class we recognize that neighbors have not
        # been found. After the user either calls self.set_up_neighbor_array
        # or self.overland_flow this will be set to True. This is done so
        # that every iteration of self.overland_flow does NOT need to
        # reinitalize the neighbors and saves computation time.
        self._neighbor_flag = False

        # When looking for neighbors, we automatically ignore inactive links
        # by default. However, what about when we want to look at fixed links
        # too? By default, we ignore these, but if they are important to your
        # model and will be updated in your driver loop, they can be used by
        # setting the flag in the initialization of  the class to 'True'
        self._default_fixed_links = default_fixed_links

        # Assiging a class variable to the elevation field.
        self._z = self._grid.at_node["topographic__elevation"]

    @property
    def h(self):
        """The depth of water at each node."""
        return self._h

    @property
    def dt(self):
        """dt: Component timestep."""
        return self._dt

    @dt.setter
    def dt(self, dt):
        assert dt > 0
        self._dt = dt

    @property
    def rainfall_intensity(self):
        """rainfall_intensity: the rainfall rate [m/s]

        Must be positive.
        """
        return self._rainfall_intensity

    @rainfall_intensity.setter
    def rainfall_intensity(self, rainfall_intensity):
        if rainfall_intensity >= 0:
            self._rainfall_intensity = rainfall_intensity
        else:
            raise ValueError("Rainfall intensity must be positive")

    def calc_time_step(self):
        """Calculate time step.

        Adaptive time stepper from Bates et al., 2010 and de Almeida et
        al., 2012
        """
        self._dt = (
            self._alpha
            * self._grid.dx
            / np.sqrt(self._g * np.amax(self._grid.at_node["surface_water__depth"]))
        )

        return self._dt

    def set_up_neighbor_arrays(self):
        """Create and initialize link neighbor arrays.

        Set up arrays of neighboring horizontal and vertical links that
        are needed for the de Almeida solution.
        """
        # First we identify all active links

        self._active_ids = links.active_link_ids(
            self._grid.shape, self._grid.status_at_node
        )

        self._active_links_at_open_bdy = _active_links_at_node(
            self.grid, self.grid.open_boundary_nodes
        ).transpose()

        self._active_links_at_open_bdy = self._active_links_at_open_bdy[
            np.where(self._active_links_at_open_bdy > -1)
        ]

        # And then find all horizontal link IDs (Active and Inactive)
        self._horizontal_ids = links.horizontal_link_ids(self._grid.shape)

        # And make the array 1-D
        self._horizontal_ids = self._horizontal_ids.flatten()

        # Find all horizontal active link ids
        self._horizontal_active_link_ids = links.horizontal_active_link_ids(
            self._grid.shape, self._active_ids
        )

        # Now we repeat this process for the vertical links.
        # First find the vertical link ids and reshape it into a 1-D array
        self._vertical_ids = links.vertical_link_ids(self._grid.shape).flatten()

        # Find the *active* verical link ids
        self._vertical_active_link_ids = links.vertical_active_link_ids(
            self._grid.shape, self._active_ids
        )

        if self._default_fixed_links is True:
            fixed_link_ids = links.fixed_link_ids(
                self._grid.shape, self._grid.status_at_node
            )
            fixed_horizontal_links = links.horizontal_fixed_link_ids(
                self._grid.shape, fixed_link_ids
            )
            fixed_vertical_links = links.vertical_fixed_link_ids(
                self._grid.shape, fixed_link_ids
            )
            self._horizontal_active_link_ids = np.maximum(
                self._horizontal_active_link_ids, fixed_horizontal_links
            )
            self._vertical_active_link_ids = np.maximum(
                self._vertical_active_link_ids, fixed_vertical_links
            )
            self._active_neighbors = find_active_neighbors_for_fixed_links(self._grid)

        self._vert_bdy_ids = self._active_links_at_open_bdy[
            links.is_vertical_link(self._grid.shape, self._active_links_at_open_bdy)
        ]

        self._vert_bdy_ids = links.nth_vertical_link(
            self._grid.shape, self._vert_bdy_ids
        )

        self._horiz_bdy_ids = self._active_links_at_open_bdy[
            links.is_horizontal_link(self._grid.shape, self._active_links_at_open_bdy)
        ]

        self._horiz_bdy_ids = links.nth_horizontal_link(
            self._grid.shape, self._horiz_bdy_ids
        )

        # Using the active vertical link ids we can find the north
        # and south vertical neighbors
        self._north_neighbors = links.vertical_north_link_neighbor(
            self._grid.shape, self._vertical_active_link_ids
        )
        self._south_neighbors = links.vertical_south_link_neighbor(
            self._grid.shape, self._vertical_active_link_ids
        )

        # Using the horizontal active link ids, we can find the west and
        # east neighbors
        self._west_neighbors = links.horizontal_west_link_neighbor(
            self._grid.shape, self._horizontal_active_link_ids
        )
        self._east_neighbors = links.horizontal_east_link_neighbor(
            self._grid.shape, self._horizontal_active_link_ids
        )

        # replace bdy condition links
        (ids,) = np.where(self._west_neighbors[self._horiz_bdy_ids] == -1)
        ids = self._horiz_bdy_ids[ids]
        self._west_neighbors[ids] = self._horizontal_active_link_ids[ids]

        (ids,) = np.where(self._east_neighbors[self._horiz_bdy_ids] == -1)
        ids = self._horiz_bdy_ids[ids]
        self._east_neighbors[ids] = self._horizontal_active_link_ids[ids]

        (ids,) = np.where(self._north_neighbors[self._vert_bdy_ids] == -1)
        ids = self._vert_bdy_ids[ids]
        self._north_neighbors[ids] = self._vertical_active_link_ids[ids]

        (ids,) = np.where(self._south_neighbors[self._vert_bdy_ids] == -1)
        ids = self._vert_bdy_ids[ids]
        self._south_neighbors[ids] = self._vertical_active_link_ids[ids]

        # Set up arrays for discharge in the horizontal & vertical directions.
        self._q_horizontal = np.zeros(
            links.number_of_horizontal_links(self._grid.shape)
        )
        self._q_vertical = np.zeros(links.number_of_vertical_links(self._grid.shape))

        # Once the neighbor arrays are set up, we change the flag to True!
        self._neighbor_flag = True

    def overland_flow(self, dt=None):
        """Generate overland flow across a grid.

        For one time step, this generates 'overland flow' across a given grid
        by calculating discharge at each node.

        Using the depth slope product, shear stress is calculated at every
        node.

        Outputs water depth, discharge and shear stress values through time at
        every point in the input grid.
        """
        # DH adds a loop to enable an imposed tstep while maintaining stability
        local_elapsed_time = 0.0
        if dt is None:
            dt = np.inf  # to allow the loop to begin
        while local_elapsed_time < dt:
            dt_local = self.calc_time_step()
            # Can really get into trouble if nothing happens but we still run:
            if not dt_local < np.inf:
                break
            if local_elapsed_time + dt_local > dt:
                dt_local = dt - local_elapsed_time
            self._dt = dt_local

            # First, we check and see if the neighbor arrays have been
            # initialized
            if self._neighbor_flag is False:
                self.set_up_neighbor_arrays()

            # In case another component has added data to the fields, we just
            # reset our water depths, topographic elevations and water
            # discharge variables to the fields.
            self._h = self._grid["node"]["surface_water__depth"]
            self._z = self._grid["node"]["topographic__elevation"]
            self._q = self._grid["link"]["surface_water__discharge"]
            self._h_links = self._grid["link"]["surface_water__depth"]

            # Here we identify the core nodes and active links for later use.
            self._core_nodes = self._grid.core_nodes
            self._active_links = self._grid.active_links

            # Per Bates et al., 2010, this solution needs to find difference
            # between the highest water surface in the two cells and the
            # highest bed elevation
            zmax = self._grid.map_max_of_link_nodes_to_link(self._z)
            w = self._h + self._z
            wmax = self._grid.map_max_of_link_nodes_to_link(w)
            hflow = wmax[self._grid.active_links] - zmax[self._grid.active_links]

            # Insert this water depth into an array of water depths at the
            # links.
            self._h_links[self._active_links] = hflow

            # Now we calculate the slope of the water surface elevation at
            # active links
            self._water_surface__gradient = self._grid.calc_grad_at_link(w)[
                self._grid.active_links
            ]

            # And insert these values into an array of all links
            self._water_surface_slope[
                self._active_links
            ] = self._water_surface__gradient
            # If the user chooses to set boundary links to the neighbor value,
            # we set the discharge array to have the boundary links set to
            # their neighbor value
            if self._default_fixed_links is True:
                self._q[self._grid.fixed_links] = self._q[self._active_neighbors]

            # Now we can calculate discharge. To handle links with neighbors
            # that do not exist, we will do a fancy indexing trick. Non-
            # existent links or inactive links have an index of '-1', which in
            # Python, looks to the end of a list or array. To accommodate these
            # '-1' indices, we will simply insert an value of 0.0 discharge (in
            # units of L^2/T) to the end of the discharge array.
            self._q = np.append(self._q, [0])

            horiz = self._horizontal_ids
            vert = self._vertical_ids
            # Now we calculate discharge in the horizontal direction
            try:
                self._q[horiz] = (
                    self._theta * self._q[horiz]
                    + (1.0 - self._theta)
                    / 2.0
                    * (self._q[self._west_neighbors] + self._q[self._east_neighbors])
                    - self._g
                    * self._h_links[horiz]
                    * self._dt
                    * self._water_surface_slope[horiz]
                ) / (
                    1
                    + self._g
                    * self._dt
                    * self._mannings_n ** 2.0
                    * abs(self._q[horiz])
                    / self._h_links[horiz] ** _SEVEN_OVER_THREE
                )

                # ... and in the vertical direction
                self._q[vert] = (
                    self._theta * self._q[vert]
                    + (1 - self._theta)
                    / 2.0
                    * (self._q[self._north_neighbors] + self._q[self._south_neighbors])
                    - self._g
                    * self._h_links[vert]
                    * self._dt
                    * self._water_surface_slope[vert]
                ) / (
                    1
                    + self._g
                    * self._dt
                    * self._mannings_n ** 2.0
                    * abs(self._q[vert])
                    / self._h_links[vert] ** _SEVEN_OVER_THREE
                )

            except ValueError:
                self._mannings_n = self._grid["link"]["mannings_n"]
                # if manning's n in a field
                # calc discharge in horizontal
                self._q[horiz] = (
                    self._theta * self._q[horiz]
                    + (1.0 - self._theta)
                    / 2.0
                    * (self._q[self._west_neighbors] + self._q[self._east_neighbors])
                    - self._g
                    * self._h_links[horiz]
                    * self._dt
                    * self._water_surface_slope[horiz]
                ) / (
                    1
                    + self._g
                    * self._dt
                    * self._mannings_n[horiz] ** 2.0
                    * abs(self._q[horiz])
                    / self._h_links[horiz] ** _SEVEN_OVER_THREE
                )

                # ... and in the vertical direction
                self._q[vert] = (
                    self._theta * self._q[vert]
                    + (1 - self._theta)
                    / 2.0
                    * (self._q[self._north_neighbors] + self._q[self._south_neighbors])
                    - self._g
                    * self._h_links[vert]
                    * self._dt
                    * self._water_surface_slope[self._vertical_ids]
                ) / (
                    1
                    + self._g
                    * self._dt
                    * self._mannings_n[vert] ** 2.0
                    * abs(self._q[vert])
                    / self._h_links[vert] ** _SEVEN_OVER_THREE
                )

            # Now to return the array to its original length (length of number
            # of all links), we delete the extra 0.0 value from the end of the
            # array.
            self._q = np.delete(self._q, len(self._q) - 1)

            # Updating the discharge array to have the boundary links set to
            # their neighbor
            if self._default_fixed_links is True:
                self._q[self._grid.fixed_links] = self._q[self._active_neighbors]

            if self._steep_slopes is True:
                # To prevent water from draining too fast for our time steps...
                # Our Froude number.
                Fr = 1.0
                # Our two limiting factors, the froude number and courant
                # number.
                # Looking a calculated q to be compared to our Fr number.
                calculated_q = (self._q / self._h_links) / np.sqrt(
                    self._g * self._h_links
                )

                # Looking at our calculated q and comparing it to Courant no.,
                q_courant = self._q * self._dt / self._grid.dx

                # Water depth split equally between four links..
                water_div_4 = self._h_links / 4.0

                # IDs where water discharge is positive...
                (positive_q,) = np.where(self._q > 0)

                # ... and negative.
                (negative_q,) = np.where(self._q < 0)

                # Where does our calculated q exceed the Froude number? If q
                # does exceed the Froude number, we are getting supercritical
                # flow and discharge needs to be reduced to maintain stability.
                (Froude_logical,) = np.where((calculated_q) > Fr)
                (Froude_abs_logical,) = np.where(abs(calculated_q) > Fr)

                # Where does our calculated q exceed the Courant number and
                # water depth divided amongst 4 links? If the calculated q
                # exceeds the Courant number and is greater than the water
                # depth divided by 4 links, we reduce discharge to maintain
                # stability.
                (water_logical,) = np.where(q_courant > water_div_4)
                (water_abs_logical,) = np.where(abs(q_courant) > water_div_4)

                # Where are these conditions met? For positive and negative q,
                # there are specific rules to reduce q. This step finds where
                # the discharge values are positive or negative and where
                # discharge exceeds the Froude or Courant number.
                self._if_statement_1 = np.intersect1d(positive_q, Froude_logical)
                self._if_statement_2 = np.intersect1d(negative_q, Froude_abs_logical)
                self._if_statement_3 = np.intersect1d(positive_q, water_logical)
                self._if_statement_4 = np.intersect1d(negative_q, water_abs_logical)

                # Rules 1 and 2 reduce discharge by the Froude number.
                self._q[self._if_statement_1] = self._h_links[self._if_statement_1] * (
                    np.sqrt(self._g * self._h_links[self._if_statement_1]) * Fr
                )

                self._q[self._if_statement_2] = 0.0 - (
                    self._h_links[self._if_statement_2]
                    * np.sqrt(self._g * self._h_links[self._if_statement_2])
                    * Fr
                )

                # Rules 3 and 4 reduce discharge by the Courant number.
                self._q[self._if_statement_3] = (
                    (self._h_links[self._if_statement_3] * self._grid.dx) / 5.0
                ) / self._dt

                self._q[self._if_statement_4] = (
                    0.0
                    - (self._h_links[self._if_statement_4] * self._grid.dx / 5.0)
                    / self._dt
                )

            # Once stability has been restored, we calculate the change in
            # water depths on all core nodes by finding the difference between
            # inputs (rainfall) and the inputs/outputs (flux divergence of
            # discharge)
            self._dhdt = self._rainfall_intensity - self._grid.calc_flux_div_at_node(
                self._q
            )

            # Updating our water depths...
            self._h[self._core_nodes] = (
                self._h[self._core_nodes] + self._dhdt[self._core_nodes] * self._dt
            )

            # To prevent divide by zero errors, a minimum threshold water depth
            # must be maintained. To reduce mass imbalances, this is set to
            # find locations where water depth is smaller than h_init (default
            # is 0.001) and the new value is self._h_init * 10^-3. This was set
            # as it showed the smallest amount of mass creation in the grid
            # during testing.
            if self._steep_slopes is True:
                self._h[self._h < self._h_init] = self._h_init * 10.0 ** -3

            # And reset our field values with the newest water depth and
            # discharge.
            self._grid.at_node["surface_water__depth"] = self._h
            self._grid.at_link["surface_water__discharge"] = self._q
            #
            #
            #            self._helper_q = self._grid.map_upwind_node_link_max_to_node(self._q)
            #            self._helper_s = self._grid.map_upwind_node_link_max_to_node(
            #                                                    self._water_surface_slope)
            #
            #            self._helper_q = self._grid.map_max_of_link_nodes_to_link(self._helper_q)
            #            self._helper_s = self._grid.map_max_of_link_nodes_to_link(self._helper_s)
            #
            #            self._grid['link']['surface_water__discharge'][
            #                     self._active_links_at_open_bdy] = self._helper_q[
            #                     self._active_links_at_open_bdy]
            #
            #            self._grid['link']['water_surface__gradient'][
            #                self._active_links_at_open_bdy] = self._helper_s[
            #                self._active_links_at_open_bdy]
            # Update nodes near boundary locations - nodes adjacent to
            # boundaries may have discharge and water surface slopes
            # artifically reduced due to boundary effects. This step removes
            # those errors.

            if dt is np.inf:
                break
            local_elapsed_time += self._dt

    def run_one_step(self, dt=None):
        """Generate overland flow across a grid.

        For one time step, this generates 'overland flow' across a given grid
        by calculating discharge at each node.

        Using the depth slope product, shear stress is calculated at every
        node.

        Outputs water depth, discharge and shear stress values through time at
        every point in the input grid.
        """
        self.overland_flow(dt=dt)

    def discharge_mapper(self, input_discharge, convert_to_volume=False):
        """Maps discharge value from links onto nodes.

        This method takes the discharge values on links and determines the
        links that are flowing INTO a given node. The fluxes moving INTO a
        given node are summed.

        This method ignores all flow moving OUT of a given node.

        This takes values from the OverlandFlow component (by default) in
        units of [L^2/T]. If the convert_to_cms flag is raised as True, this
        method converts discharge to units [L^3/T] - as of Aug 2016, only
        operates for square RasterModelGrid instances.

        The output array is of length grid.number_of_nodes and can be used
        with the Landlab imshow_grid plotter.

        Returns a numpy array (discharge_vals)
        """

        discharge_vals = np.zeros(self._grid.number_of_links)
        discharge_vals[:] = input_discharge[:]

        if convert_to_volume:
            discharge_vals *= self._grid.dx

        discharge_vals = (
            discharge_vals[self._grid.links_at_node] * self._grid.link_dirs_at_node
        )

        discharge_vals = discharge_vals.flatten()

        discharge_vals[np.where(discharge_vals < 0)] = 0.0

        discharge_vals = discharge_vals.reshape(self._grid.number_of_nodes, 4)

        discharge_vals = discharge_vals.sum(axis=1)

        return discharge_vals


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
    >>> from landlab import NodeStatus, RasterModelGrid
    >>> from landlab.components.overland_flow.generate_overland_flow_deAlmeida import find_active_neighbors_for_fixed_links

    >>> grid = RasterModelGrid((4, 5))
    >>> grid.status_at_node[:5] = NodeStatus.FIXED_GRADIENT
    >>> grid.status_at_node[::5] = NodeStatus.FIXED_GRADIENT
    >>> grid.status_at_node # doctest: +NORMALIZE_WHITESPACE
    array([2, 2, 2, 2, 2,
           2, 0, 0, 0, 1,
           2, 0, 0, 0, 1,
           2, 1, 1, 1, 1], dtype=uint8)

    >>> grid.fixed_links
    array([ 5,  6,  7,  9, 18])
    >>> grid.active_links
    array([10, 11, 12, 14, 15, 16, 19, 20, 21, 23, 24, 25])

    >>> find_active_neighbors_for_fixed_links(grid)
    array([14, 15, 16, 10, 19])

    >>> rmg = RasterModelGrid((4, 7))

    >>> rmg.at_node['topographic__elevation'] = rmg.zeros(at='node')
    >>> rmg.at_link['topographic__slope'] = rmg.zeros(at='link')
    >>> rmg.status_at_node[rmg.perimeter_nodes] = rmg.BC_NODE_IS_FIXED_GRADIENT
    >>> find_active_neighbors_for_fixed_links(rmg)
    array([20, 21, 22, 23, 24, 14, 17, 27, 30, 20, 21, 22, 23, 24])
    """
    neighbors = links.neighbors_at_link(grid.shape, grid.fixed_links).flat
    return neighbors[np.in1d(neighbors, grid.active_links)]
