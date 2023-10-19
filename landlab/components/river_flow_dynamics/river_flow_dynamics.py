"""river_flow_dynamics.py

This component implements a semi-implicit, semi-Lagrangian finite-volume approximation of
the depth-averaged shallow water equations originally proposed by Casulli and Cheng in 1992,
and subsequent related work.

Written by Sebastian Bernal and Angel Monsalve.

Last updated: October 18, 2023

Examples
--------

River Flow Dynamics Simulation Example

This example demonstrates how to simulate river flow dynamics using the Landlab library.

First, import necessary libraries such as NumPy, Matplotlib, and Landlab components.

>>> import numpy as np
>>> import matplotlib.pyplot as plt
>>> from landlab import RasterModelGrid
>>> from landlab.components import river_flow_dynamics
>>> from landlab.io import read_esri_ascii
>>> from landlab.plot.imshow import imshow_grid

Create a rectangular grid for flow dynamics calculations with specified dimensions
and cell size.

>>> nRows = 20
>>> nCols = 60
>>> cellSize = 0.1
>>> grid = RasterModelGrid((nRows,nCols), xy_spacing=(cellSize,cellSize))

Defining the Topography: Set up the initial topographic elevation for the grid,
creating a basic rectangular channel with a slope of 0.01.

>>> te = grid.add_zeros("topographic__elevation", at="node")
>>> te += 0.059 - 0.01*grid.x_of_node
>>> te[grid.y_of_node > 1.5] = 1.0
>>> te[grid.y_of_node < 0.5] = 1.0

We could visualize the elevation profile using 'imshow_grid'.
imshow_grid(grid, "topographic__elevation")
Explore the grid's middle longitudinal section.

>>> middleBedProfile = np.reshape(te,(nRows,nCols))[10,:]
>>> np.round(middleBedProfile,3)
array([ 0.059,  0.058,  0.057,  0.056,  0.055,  0.054,  0.053,  0.052,
        0.051,  0.05 ,  0.049,  0.048,  0.047,  0.046,  0.045,  0.044,
        0.043,  0.042,  0.041,  0.04 ,  0.039,  0.038,  0.037,  0.036,
        0.035,  0.034,  0.033,  0.032,  0.031,  0.03 ,  0.029,  0.028,
        0.027,  0.026,  0.025,  0.024,  0.023,  0.022,  0.021,  0.02 ,
        0.019,  0.018,  0.017,  0.016,  0.015,  0.014,  0.013,  0.012,
        0.011,  0.01 ,  0.009,  0.008,  0.007,  0.006,  0.005,  0.004,
        0.003,  0.002,  0.001, -0.   ])

Instantiating the Component: Initialize the river_flow_dynamics component with
specified parameters, including the time step and Manning's roughness coefficient.
The grid will need some data to run the river_flow_dynamics component.
To check the names of the required inputs, use the 'input_var_names' class property.

>>> river_flow_dynamics.input_var_names
('surface_water__depth',
 'surface_water__elevation',
 'surface_water__velocity',
 'topographic__elevation')

To determine where these fields are mapped, use 'var_mapping'.

>>> river_flow_dynamics.var_mapping
(('surface_water__depth', 'node'),
 ('surface_water__elevation', 'node'),
 ('surface_water__velocity', 'link'),
 ('topographic__elevation', 'node'))

Create fields of data for each of these input variables. The channel is initially empty,
so we create the water depth field.

>>> h = grid.add_zeros("surface_water__depth", at="node")

Water velocity is zero everywhere since there is no water yet.

>>> vel = grid.add_zeros("surface_water__velocity", at="link")

Calculate the initial water surface elevation from water depth and topographic elevation.

>>> wse = grid.add_zeros("surface_water__elevation", at="node")
>>> wse += h + te

Specify the nodes at which water enters the domain, and also the associated links.
These will serve as the inlet boundary conditions for water depth and velocity.
In this case, water flows from left to right at a depth of 0.5 meters with a velocity
of 0.45 m/s.

>>> fixed_entry_nodes = np.arange(300,910,60)
>>> fixed_entry_links = grid.links_at_node[fixed_entry_nodes][:,0]

Set the fixed values for these entry nodes/links.

>>> entry_nodes_h_values   = np.full(11, 0.5)
>>> entry_links_vel_values = np.full(11, 0.45)

Instantiate 'river_flow_dynamics' with the previously defined arguments.

>>> rfd = river_flow_dynamics(grid, dt=0.1, mannings_n=0.012,
...     fixed_entry_nodes=fixed_entry_nodes,
...     fixed_entry_links=fixed_entry_links,entry_nodes_h_values=entry_nodes_h_values,
...     entry_links_vel_values=entry_links_vel_values)

Run the simulation for 100 timesteps (equivalent to 10 seconds).

>>> n_timesteps = 100
>>> for timestep in range(n_timesteps):
...     rfd.run_one_step()

Examine the flow depth at the center of the channel after 10 seconds.

>>> np.reshape(grid['node']["surface_water__depth"],(nRows,nCols))[10,:]
array([ 0.5       ,  0.49119496,  0.48070117,  0.47309464,  0.46764464,
        0.46370528,  0.46076177,  0.45840315,  0.45619896,  0.45419369,
        0.45250159,  0.45106348,  0.44972981,  0.4484113 ,  0.44702824,
        0.445566  ,  0.44414898,  0.44279566,  0.44144646,  0.44007219,
        0.43868299,  0.43726816,  0.43582451,  0.4344022 ,  0.43300387,
        0.43159008,  0.43012984,  0.42867334,  0.42723278,  0.42577057,
        0.42426757,  0.42277662,  0.42127323,  0.41977445,  0.41825745,
        0.41674194,  0.41522539,  0.41368181,  0.41214349,  0.41060916,
        0.40905859,  0.40748906,  0.40591859,  0.40435914,  0.40275767,
        0.40113684,  0.39954249,  0.39795401,  0.39635525,  0.39474439,
        0.39313631,  0.39155279,  0.38993334,  0.38828883,  0.38664499,
        0.3849828 ,  0.38333126,  0.381655  ,  0.37977436,  0.37881807])

And the velocity at links along the center of the channel

>>> linksAtCenter = grid.links_at_node[np.array(np.arange(600,660))][:-1,0]
>>> grid['link']["surface_water__velocity"][linksAtCenter]
array([ 0.45      ,  0.59396861,  0.69574417,  0.75731793,  0.79261482,
        0.8123394 ,  0.82700801,  0.84314289,  0.85567049,  0.85979031,
        0.85838515,  0.85658338,  0.85716582,  0.85944389,  0.86308335,
        0.86497462,  0.86463743,  0.86367451,  0.86314567,  0.86565779,
        0.86981418,  0.8726041 ,  0.87338014,  0.87403681,  0.87494614,
        0.87549414,  0.87571505,  0.87786193,  0.87953196,  0.8794075 ,
        0.88030366,  0.88296806,  0.88376619,  0.88519116,  0.88610576,
        0.88744898,  0.88909536,  0.8898242 ,  0.89135825,  0.89318975,
        0.89358086,  0.89405077,  0.89630327,  0.89764043,  0.8976985 ,
        0.89895516,  0.90065785,  0.90219923,  0.90227029,  0.90168375,
        0.90460332,  0.90538973,  0.9036968 ,  0.90546049,  0.90666316,
        0.90682276,  0.90754984,  0.90761995,  0.91305765])

"""

import numpy as np
import scipy as sp

from landlab import Component, FieldError


class river_flow_dynamics(Component):
    """Simulate surface fluid flow based on Casulli and Cheng (1992).

    This Landlab component simulates surface fluid flow using the approximations of the
    2D shallow water equations developed by Casulli and Cheng in 1992. It calculates water
    depth and velocity across the raster grid, given a specific input discharge.

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    Casulli, V., Cheng, R.T. (1992). “Semi-implicit finite difference methods for
    three-dimensional shallow water flow”. International Journal for Numerical Methods
    in Fluids. 15: 629-648.
    https://doi.org/10.1002/fld.1650150602

    """

    _name = "river_flow_dynamics"

    _unit_agnostic = False

    _info = {
        "surface_water__depth": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of water on the surface",
        },
        "surface_water__velocity": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m/s",
            "mapping": "link",
            "doc": "Speed of water flow above the surface",
        },
        "surface_water__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Water surface elevation at time N",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
    }

    def __init__(
        self,
        grid,
        dt=0.01,  # Sets the time step (s)
        eddy_viscosity=1e-4,  # ddy viscosity coefficient
        mannings_n=0.012,  # Manning's n
        threshold_depth=0.01,  # Sets the wet/dry threshold
        theta=0.5,  # Degree of 'implicitness' of the solution
        fixed_entry_nodes=None,  # Node IDs where flow enters the domain
        fixed_entry_links=None,  # Link IDs where flow enters the domain
        entry_nodes_h_values=None,  # Water depth at nodes where flow enters the domain
        entry_links_vel_values=None,  # Water velocity at links where flow enters the domain
        pcg_tolerance=1e-05,  # Preconditioned Conjugate Gradient convergence tolerance
        pcg_max_iterations=None,  # Preconditioned Conjugate Gradient max iterations
        surface_water__elevation_at_N_1=None,  # Surf water elev at prev. time
        surface_water__elevation_at_N_2=None,  # Surf water elev at prev prev time
        surface_water__velocity_at_N_1=None,  # Speed of water at prev time
    ):
        """Simulate the vertical-averaged surface fluid flow

        Simulate vertical-averaged surface fluid flow using the Casulli and Cheng (1992)
        approximations of the 2D shallow water equations. This Landlab component calculates
        water depth and velocity across the raster grid based on a given input discharge.

        Parameters
        ----------
        grid : RasterModelGrid
            A grid.
        dt : float, optional
            Time step in seconds. If not provided, it is calculated from the CFL condition.
        eddy_viscosity : float, optional
            Eddy viscosity coefficient. Default = 1e-4 :math:`m^2 / s`
        mannings_n : float or array_like, optional
            Manning's roughness coefficient. Default = 0.012 :math:`s / m^1/3`
        threshold_depth : float, optional
            Threshold at which a cell is considered wet. Default = 0.01 m
        theta : float, optional
            Degree of 'implicitness' of the solution, ranging between 0.5 and 1.0.
            Default: 0.5. When set to 0.5, the approximation is centered in time;
            when set to 1.0, it is fully implicit.
        fixed_entry_nodes : array_like or None, optional
            Node IDs where flow enters the domain (Dirichlet boundary condition).
            If not provided, existing water in the domain is not renewed.
        fixed_entry_links : array_like or None, optional
            Link IDs where flow enters the domain (Dirichlet boundary condition).
            If not provided, existing water in the domain is not renewed.
        entry_nodes_h_values : array_like, optional
            Water depth values at nodes where flow enters the domain
            (Dirichlet boundary condition).
            If not provided, existing water in the domain is not renewed.
        entry_links_vel_values : array_like, optional
            Water velocity values at links where flow enters the domain
            (Dirichlet boundary condition).
            If not provided, existing water in the domain is not renewed.
        pcg_tolerance : float, optional
            Tolerance for convergence in the Preconditioned Conjugate Gradient
            method. Default: 1e-05.
        pcg_max_iterations : integer, optional
            Maximum number of iterations for the Preconditioned Conjugate Gradient
            method. Iteration stops after maxiter steps, even if the specified
            tolerance is not achieved. Default: None.
        surface_water__elevation_at_N_1: float, optional
            Water surface elevation at time N-1.
            units: m, mapping: node
        surface_water__elevation_at_N_2: float, optional
            Water surface elevation at time N-2.
            units: m, mapping: node
        surface_water__velocity_at_N_1: float, optional
            Speed of water flow above the surface at time N-1",
            units: m/2, mapping: link
        """
        super().__init__(grid)

        # User inputs
        self._dt = dt
        self._eddy_viscosity = eddy_viscosity
        self._g = sp.constants.g
        self._mannings_n = mannings_n
        self._threshold_depth = threshold_depth
        self._theta = theta
        self._pcg_tolerance = pcg_tolerance
        self._pcg_max_iterations = pcg_max_iterations

        # Getting topography for further calculations
        self._additional_z = 10  # To set the virtual reference elevation (z=0)
        self._max_elevation = self._grid["node"]["topographic__elevation"].max()
        self._z = self._grid["node"]["topographic__elevation"]
        self._z = (self._z.max() + self._additional_z) - self._z

        if fixed_entry_nodes is None:
            fixed_entry_nodes = []
        self._fixed_entry_nodes = fixed_entry_nodes

        if fixed_entry_links is None:
            fixed_entry_links = []
        self._fixed_entry_links = fixed_entry_links

        if entry_nodes_h_values is None:
            entry_nodes_h_values = []
        self._entry_nodes_h_values = entry_nodes_h_values

        if entry_links_vel_values is None:
            entry_links_vel_values = []
        self._entry_links_vel_values = entry_links_vel_values

        # Creating fields if they don't exist
        try:
            self._grid["node"]["surface_water__depth"] = grid.add_zeros(
                "surface_water__depth",
                at="node",
                units=self._info["surface_water__depth"]["units"],
            )
        except FieldError:
            self._grid["node"]["surface_water__depth"] = grid.at_node[
                "surface_water__depth"
            ]

        try:
            self._grid["link"]["surface_water__velocity"] = grid.add_zeros(
                "surface_water__velocity",
                at="link",
                units=self._info["surface_water__velocity"]["units"],
            )
        except FieldError:
            self._grid["link"]["surface_water__velocity"] = grid.at_link[
                "surface_water__velocity"
            ]

        try:
            self._grid["node"]["surface_water__elevation"] = grid.add_zeros(
                "surface_water__elevation",
                at="node",
                units=self._info["surface_water__elevation"]["units"],
            )
            self._grid["node"]["surface_water__elevation"] = (
                self._grid["node"]["surface_water__depth"] - self._z
            )
        except FieldError:
            self._grid["node"]["surface_water__elevation"] = grid.at_node[
                "surface_water__elevation"
            ]

        if surface_water__elevation_at_N_1 is None:
            self._surface_water__elevation_at_N_1 = np.zeros(self._grid.number_of_nodes)
        else:
            if surface_water__elevation_at_N_1.size > 0:
                if (
                    surface_water__elevation_at_N_1.shape[0]
                    == self._grid.number_of_nodes
                ):
                    self._surface_water__elevation_at_N_1 = (
                        surface_water__elevation_at_N_1
                    )
                else:
                    raise ValueError(
                        "surface_water__elevation_at_N_1 \
                        does not have the same dimensions of the grid's nodes"
                    )

        if surface_water__elevation_at_N_2 is None:
            self._surface_water__elevation_at_N_2 = np.zeros(self._grid.number_of_nodes)
        else:
            if surface_water__elevation_at_N_2.size > 0:
                if (
                    surface_water__elevation_at_N_2.shape[0]
                    == self._grid.number_of_nodes
                ):
                    self._surface_water__elevation_at_N_2 = (
                        surface_water__elevation_at_N_2
                    )
                else:
                    raise ValueError(
                        "surface_water__elevation_at_N_2 \
                        does not have the same dimensions of the grid's nodes"
                    )

        if surface_water__velocity_at_N_1 is None:
            self._surface_water__velocity_at_N_1 = np.zeros(self._grid.number_of_links)
        else:
            if surface_water__velocity_at_N_1.size > 0:
                if (
                    surface_water__velocity_at_N_1.shape[0]
                    == self._grid.number_of_links
                ):
                    self._surface_water__velocity_at_N_1 = (
                        surface_water__velocity_at_N_1
                    )
                else:
                    raise ValueError(
                        "surface_water__velocity_at_N_1 \
                        does not have the same dimensions of the grid's links"
                    )

        # Assigning a class variable to the fields
        self._h = self._grid["node"]["surface_water__depth"]
        self._vel = self._grid["link"]["surface_water__velocity"]
        self._vel_at_N_1 = self._surface_water__velocity_at_N_1
        self._eta = self._grid["node"]["surface_water__elevation"] - (
            self._max_elevation + self._additional_z
        )
        self._eta_at_N_1 = self._surface_water__elevation_at_N_1 - (
            self._max_elevation + self._additional_z
        )
        self._eta_at_N_2 = self._surface_water__elevation_at_N_2 - (
            self._max_elevation + self._additional_z
        )

        # Getting nodes and links at each edge
        self._nodes_at_right_edge = self._grid.nodes_at_right_edge
        self._nodes_at_top_edge = self._grid.nodes_at_top_edge
        self._nodes_at_left_edge = self._grid.nodes_at_left_edge
        self._nodes_at_bottom_edge = self._grid.nodes_at_bottom_edge

        self._nx = len(self._nodes_at_top_edge)
        self._ny = len(self._nodes_at_left_edge)

        self._dx = self._grid.dx
        self._dy = self._grid.dy

        # Open boundary conditions
        # water can leave the domain at everywhere, only limited by topography
        self._grid.status_at_node[
            self._nodes_at_left_edge
        ] = self._grid.BC_NODE_IS_FIXED_VALUE
        self._grid.status_at_node[
            self._nodes_at_right_edge
        ] = self._grid.BC_NODE_IS_FIXED_VALUE
        self._grid.status_at_node[
            self._nodes_at_bottom_edge
        ] = self._grid.BC_NODE_IS_FIXED_VALUE
        self._grid.status_at_node[
            self._nodes_at_top_edge
        ] = self._grid.BC_NODE_IS_FIXED_VALUE

        # Identifying node and link ids for later use.
        self._core_nodes = self._grid.core_nodes
        self._corner_nodes = self._grid.corner_nodes
        self._active_links = self._grid.active_links
        self._horizontal_links = self._grid.horizontal_links
        self._vertical_links = self._grid.vertical_links
        self._number_of_nodes = self._grid.number_of_nodes

        self._adjacent_nodes_at_corner_nodes = np.array(
            [
                [
                    self._nodes_at_top_edge[-2],
                    self._nodes_at_right_edge[-2],
                ],  # Top right
                [self._nodes_at_top_edge[1], self._nodes_at_left_edge[-2]],  # Top left
                [
                    self._nodes_at_left_edge[1],
                    self._nodes_at_bottom_edge[1],
                ],  # Bottom left
                [
                    self._nodes_at_right_edge[1],
                    self._nodes_at_bottom_edge[-2],
                ],  # Bottom right
            ]
        )

        # Updating open boundary nodes/links
        self._open_boundary_nodes = self._grid.boundary_nodes
        self._open_boundary_links = np.unique(
            self._grid.links_at_node[self._open_boundary_nodes]
        )

        self._open_boundary_nodes = np.setdiff1d(
            self._open_boundary_nodes, self._fixed_entry_nodes
        )
        self._open_boundary_links = np.setdiff1d(
            self._open_boundary_links, self._fixed_entry_links
        )

        self._fixed_corner_nodes = np.setdiff1d(
            self._corner_nodes, self._open_boundary_nodes
        )
        self._open_corner_nodes = np.setdiff1d(
            self._corner_nodes, self._fixed_corner_nodes
        )

        self._open_boundary_nodes = np.setdiff1d(
            self._open_boundary_nodes, self._open_corner_nodes
        )

        # Using fixed entry nodes/links only when they exist
        if len(self._fixed_entry_nodes) > 0:
            self._fixed_nodes_exist = True
        else:
            self._fixed_nodes_exist = False

        if len(self._fixed_entry_links) > 0:
            self._fixed_links_exist = True
        else:
            self._fixed_links_exist = False

        # Updating grid fixed values according to the user input
        if self._fixed_nodes_exist is True:
            self._h[self._fixed_entry_nodes] = entry_nodes_h_values
            self._eta[self._fixed_entry_nodes] = (
                entry_nodes_h_values - self._z[self._fixed_entry_nodes]
            )
        if self._fixed_links_exist is True:
            self._vel[self._fixed_entry_links] = entry_links_vel_values

        # Mapping node values at links
        self._z_at_links = self._grid.map_mean_of_link_nodes_to_link(self._z)
        self._h_at_links = self._grid.map_mean_of_link_nodes_to_link(self._h)
        self._eta_at_links = self._h_at_links - self._z_at_links

        # Passing values to the time step N
        self._h_at_N = self._h.copy()
        self._h_at_N_at_links = self._grid.map_mean_of_link_nodes_to_link(self._h_at_N)

        self._vel_at_N = self._vel.copy()
        self._eta_at_N = self._eta.copy()

        # Boolean for wet nodes/links
        self._wet_nodes = np.where(self._h_at_N >= self._threshold_depth, True, False)
        self._wet_links = np.where(
            self._h_at_N_at_links >= self._threshold_depth, True, False
        )

    # Defining some functions
    def find_nearest_link(self, x_coordinates, y_coordinates, objective_links="all"):
        """Link nearest a point.

        Find the index to the link nearest the given x, y coordinates.
        Returns the indices of the links nearest the given coordinates.

        """
        # Defining the set of links that are going to be used
        if objective_links == "all":
            objective_links = np.arange(self._grid.number_of_links)
        elif objective_links == "horizontal":
            objective_links = self._horizontal_links
        elif objective_links == "vertical":
            objective_links = self._vertical_links
        # if (objective_links == "all") END

        # Coordinates of all the RasterModelGrid links
        x_of_objective_links = np.unique(self._grid.xy_of_link[objective_links][:, 0])
        y_of_objective_links = np.unique(self._grid.xy_of_link[objective_links][:, 1])

        # Getting the closest link-coordinate to the exit point
        tempCalc1 = np.repeat(x_coordinates, len(x_of_objective_links)).reshape(
            len(x_coordinates), len(x_of_objective_links)
        )
        tempCalc2 = np.tile(x_of_objective_links, len(x_coordinates)).reshape(
            len(x_coordinates), len(x_of_objective_links)
        )
        indices = abs(tempCalc2 - tempCalc1).argmin(axis=1)
        nearest_x = x_of_objective_links[indices]

        tempCalc1 = np.repeat(y_coordinates, len(y_of_objective_links)).reshape(
            len(y_coordinates), len(y_of_objective_links)
        )
        tempCalc2 = np.tile(y_of_objective_links, len(y_coordinates)).reshape(
            len(y_coordinates), len(y_of_objective_links)
        )
        indices = abs(tempCalc2 - tempCalc1).argmin(axis=1)
        nearest_y = y_of_objective_links[indices]

        # Getting the closest link to link
        tempCalc1 = np.repeat(
            nearest_x, len(self._grid.xy_of_link[objective_links][:, 0])
        ).reshape(len(nearest_x), len(self._grid.xy_of_link[objective_links][:, 0]))
        tempCalc2 = np.tile(
            self._grid.xy_of_link[objective_links][:, 0], len(x_coordinates)
        ).reshape(len(x_coordinates), len(self._grid.xy_of_link[objective_links][:, 0]))
        tempB1 = tempCalc1 == tempCalc2

        tempCalc1 = np.repeat(
            nearest_y, len(self._grid.xy_of_link[objective_links][:, 1])
        ).reshape(len(nearest_y), len(self._grid.xy_of_link[objective_links][:, 1]))
        tempCalc2 = np.tile(
            self._grid.xy_of_link[objective_links][:, 1], len(y_coordinates)
        ).reshape(len(y_coordinates), len(self._grid.xy_of_link[objective_links][:, 1]))
        tempB2 = tempCalc1 == tempCalc2

        tempCalc3 = (
            np.repeat(objective_links, len(x_coordinates))
            .reshape((len(objective_links), len(y_coordinates)))
            .T
        )
        nearest_link = tempCalc3[tempB1 * tempB2]

        return nearest_link.astype(int)

    def find_adjacent_links_at_link(self, current_link, objective_links="horizontal"):
        """Get adjacent links to the link.

        This function finds the links at right, above, left and below the given link.
        Similar purpose to the "adjacent_nodes_at_node" function.
        Return the adjacent links in as many rows as given links.
        Link IDs are returned as columns in clock-wise order starting from East (E, N, W, S).

        """

        # Defining the set of links that are going to be used
        if objective_links == "horizontal":
            objective_links = self._horizontal_links
            reshape_pair = (self._ny, self._nx - 1)
        elif objective_links == "vertical":
            objective_links = self._vertical_links
            reshape_pair = (self._ny - 1, self._nx)
        # if (objective_links == "horizontal") END

        # Coordinates of the current link
        x_of_current_links = self._grid.xy_of_link[current_link][:, 0]
        y_of_current_links = self._grid.xy_of_link[current_link][:, 1]

        # Coordinates of all the RasterModelGrid links
        x_of_objective_links = np.unique(self._grid.xy_of_link[objective_links][:, 0])
        y_of_objective_links = np.unique(self._grid.xy_of_link[objective_links][:, 1])

        # Getting links that share the same y-coordinates
        # The following matrices are built to be compared to each other.
        # tempCalc1 repeats "y_of_current_links" for every x-coordinate in
        # "objective_links": cols = "y_of_current_links"
        # tempCalc2 repeats "y_of_objective_links" for every x-coordinate in
        # the "current_link": rows = "y_of_objective_links"
        # tempCalc3 give us the index to extract all the objective links that
        # are located in the same row than the current link: rows = [0, 1, 2, ...]
        tempCalc1 = np.repeat(y_of_current_links, len(y_of_objective_links)).reshape(
            len(y_of_current_links), len(y_of_objective_links)
        )
        tempCalc2 = np.tile(y_of_objective_links, len(y_of_current_links)).reshape(
            len(y_of_current_links), len(y_of_objective_links)
        )
        tempCalc3 = (
            np.repeat(np.arange(len(y_of_objective_links)), len(y_of_current_links))
            .reshape(len(y_of_objective_links), len(y_of_current_links))
            .T
        )

        indices = tempCalc3[tempCalc1 == tempCalc2]
        links_at_same_rows = objective_links.reshape(reshape_pair)[indices, :]
        links_at_same_rows = np.append(
            np.array([-np.ones_like(current_link)]).T, links_at_same_rows, axis=1
        )
        links_at_same_rows = np.append(
            links_at_same_rows, np.array([-np.ones_like(current_link)]).T, axis=1
        )

        # Getting links that share the same x-coordinates
        # The following matrices are built to be compared to each other.
        # tempCalc1 repeats "x_of_current_links" for every x-coordinate in
        # "objective_links": cols = "x_of_current_links"
        # tempCalc2 repeats "x_of_objective_links" for every x-coordinate in
        # the "current_link": rows = "x_of_objective_links"
        # tempCalc3 give us the index to extract all the objective links that
        # are located in the same row than the current link: rows = [0, 1, 2, ...]
        tempCalc1 = np.repeat(x_of_current_links, len(x_of_objective_links)).reshape(
            len(x_of_current_links), len(x_of_objective_links)
        )
        tempCalc2 = np.tile(x_of_objective_links, len(x_of_current_links)).reshape(
            len(x_of_current_links), len(x_of_objective_links)
        )
        tempCalc3 = (
            np.repeat(np.arange(len(x_of_objective_links)), len(x_of_current_links))
            .reshape(len(x_of_objective_links), len(x_of_current_links))
            .T
        )

        indices = tempCalc3[tempCalc1 == tempCalc2]
        links_at_same_cols = objective_links.reshape(reshape_pair)[:, indices].T
        links_at_same_cols = np.append(
            np.array([-np.ones_like(current_link)]).T, links_at_same_cols, axis=1
        )
        links_at_same_cols = np.append(
            links_at_same_cols, np.array([-np.ones_like(current_link)]).T, axis=1
        )

        # Extracing the adjacent links to current link (E,N,W,S)
        adjacent_links_at_link = np.zeros((current_link.shape[0], 4))

        # Rows (E,W)
        tempCalc1 = np.repeat(current_link, links_at_same_rows.shape[1]).reshape(
            current_link.shape[0], links_at_same_rows.shape[1]
        )
        tempCalc2 = (
            np.repeat(np.arange(links_at_same_rows.shape[1]), current_link.shape[0])
            .reshape(links_at_same_rows.shape[1], current_link.shape[0])
            .T
        )
        tempCalc3 = tempCalc2[tempCalc1 == links_at_same_rows]

        adjacent_links_at_link[:, 0] = links_at_same_rows[
            (range(links_at_same_rows.shape[0])), (tempCalc3 + 1)
        ]
        adjacent_links_at_link[:, 2] = links_at_same_rows[
            (range(links_at_same_rows.shape[0])), (tempCalc3 - 1)
        ]

        # Cols (N,S)
        tempCalc1 = np.repeat(current_link, links_at_same_cols.shape[1]).reshape(
            current_link.shape[0], links_at_same_cols.shape[1]
        )
        tempCalc2 = (
            np.repeat(np.arange(links_at_same_cols.shape[1]), current_link.shape[0])
            .reshape(links_at_same_cols.shape[1], current_link.shape[0])
            .T
        )
        tempCalc3 = tempCalc2[tempCalc1 == links_at_same_cols]

        adjacent_links_at_link[:, 1] = links_at_same_cols[
            (range(links_at_same_cols.shape[0])), (tempCalc3 + 1)
        ]
        adjacent_links_at_link[:, 3] = links_at_same_cols[
            (range(links_at_same_cols.shape[0])), (tempCalc3 - 1)
        ]

        return adjacent_links_at_link.astype(int)

    def path_line_tracing(self):
        """ " Path line tracing algorithm.

        This function implements the semi-analytical path line tracing method
        of Pollock (1988).

        The semi-analytical path line tracing method was developed for particle
        tracking in ground water flow models. The assumption that each directional
        velocity component varies linearly in its coordinate directions within
        each computational volume or cell underlies the method.
        Linear variation allows the derivation of an analytical expression for
        the path line of a particle across a volume.

        Given an initial point located at each volume faces of the domain, particle
        trayectories are traced backwards on time. Then, this function returns
        the departure point of the particle at the beginning of the time step.
        """

        # Calculating the partial time-step TAUx, TAUy, dt - sum_partial_times
        sum_partial_times = np.zeros_like(self._u_vel_of_particle)
        remaining_time = self._dt - sum_partial_times
        keep_tracing = np.where(remaining_time > 0, True, False)

        while np.any(remaining_time > 0):
            # Using the previous exit point as the new entry point
            self._x_of_particle = np.where(
                keep_tracing, self._x_at_exit_point, self._x_of_particle
            )
            self._y_of_particle = np.where(
                keep_tracing, self._y_at_exit_point, self._y_of_particle
            )

            # Checking if the particles departs (backwards) from a link position (True)
            tempBx = np.isin(
                self._x_of_particle, self._grid.xy_of_link[self._active_links][:, 0]
            )  # Particles located on horizontal-links/vertical-faces
            tempBy = np.isin(
                self._y_of_particle, self._grid.xy_of_link[self._active_links][:, 1]
            )  # Particles located on vertical-links/horizontal-faces

            # True, particles depart from link positions.
            # False, particles depart from random locations inside a cell
            tempBxy = tempBx + tempBy

            # Getting surrounding links for particles located inside a cell
            tempCalc1 = np.append(
                np.array([self._x_of_particle]), np.array([self._y_of_particle]), axis=0
            )
            tempCalc2 = self._grid.find_nearest_node(tempCalc1, mode="raise")
            temp_links_from_node = self._grid.links_at_node[tempCalc2]
            nodes_from_particle = tempCalc2

            # Getting surrounding links for particles located at link positions
            tempCalc1 = np.where(
                self._u_vel_of_particle >= 0,
                np.array([self._x_of_particle]) - self._dx / 10,
                np.array([self._x_of_particle]) + self._dx / 10,
            )
            tempCalc2 = np.where(
                self._v_vel_of_particle >= 0,
                np.array([self._y_of_particle]) - self._dy / 10,
                np.array([self._y_of_particle]) + self._dy / 10,
            )
            tempCalc3 = np.append(tempCalc1, tempCalc2, axis=0)
            tempCalc4 = self._grid.find_nearest_node(tempCalc3, mode="raise")
            temp_links_from_link = self._grid.links_at_node[tempCalc4]
            nodes_from_particle = np.where(tempBxy, tempCalc4, nodes_from_particle)

            # Getting links around particle
            tempBxy = np.tile(tempBxy, 4).reshape(4, len(tempBxy)).T
            links_at_particle = np.where(
                tempBxy, temp_links_from_link, temp_links_from_node
            )

            # Defining links based on velocity direction
            link_at_x2 = np.where(
                self._u_vel_of_particle >= 0,
                links_at_particle[:, 0],
                links_at_particle[:, 2],
            )
            link_at_x1 = np.where(
                self._u_vel_of_particle >= 0,
                links_at_particle[:, 2],
                links_at_particle[:, 0],
            )
            link_at_y2 = np.where(
                self._v_vel_of_particle >= 0,
                links_at_particle[:, 1],
                links_at_particle[:, 3],
            )
            link_at_y1 = np.where(
                self._v_vel_of_particle >= 0,
                links_at_particle[:, 3],
                links_at_particle[:, 1],
            )

            x_at_x2 = np.where(
                self._u_vel_of_particle >= 0,
                self._grid.x_of_node[nodes_from_particle] + self._dx / 2,
                self._grid.x_of_node[nodes_from_particle] - self._dx / 2,
            )
            x_at_x1 = np.where(
                self._u_vel_of_particle >= 0,
                self._grid.x_of_node[nodes_from_particle] - self._dx / 2,
                self._grid.x_of_node[nodes_from_particle] + self._dx / 2,
            )
            y_at_y2 = np.where(
                self._v_vel_of_particle >= 0,
                self._grid.y_of_node[nodes_from_particle] + self._dy / 2,
                self._grid.y_of_node[nodes_from_particle] - self._dy / 2,
            )
            y_at_y1 = np.where(
                self._v_vel_of_particle >= 0,
                self._grid.y_of_node[nodes_from_particle] - self._dy / 2,
                self._grid.y_of_node[nodes_from_particle] + self._dy / 2,
            )

            # Getting velocity around the particle
            u_vel_at_x2 = np.where(
                link_at_x2 >= 0, self._vel_at_N[link_at_x2], self._vel_at_N[link_at_x1]
            )
            u_vel_at_x1 = np.where(
                link_at_x1 >= 0, self._vel_at_N[link_at_x1], self._vel_at_N[link_at_x2]
            )
            v_vel_at_y2 = np.where(
                link_at_y2 >= 0, self._vel_at_N[link_at_y2], self._vel_at_N[link_at_y1]
            )
            v_vel_at_y1 = np.where(
                link_at_y1 >= 0, self._vel_at_N[link_at_y1], self._vel_at_N[link_at_y2]
            )

            # Calculating gradients for path line tracing
            gradient_x_direction = (u_vel_at_x2 - u_vel_at_x1) / self._dx
            gradient_y_direction = (v_vel_at_y2 - v_vel_at_y1) / self._dy

            # Calculating entry velocity for each particle
            self._u_vel_of_particle = u_vel_at_x2 - gradient_x_direction * (
                x_at_x2 - self._x_of_particle
            )
            self._v_vel_of_particle = v_vel_at_y2 - gradient_y_direction * (
                y_at_y2 - self._y_of_particle
            )
            self._u_vel_of_particle = np.where(
                self._u_vel_of_particle < 1e-10, 0, self._u_vel_of_particle
            )
            self._v_vel_of_particle = np.where(
                self._v_vel_of_particle < 1e-10, 0, self._v_vel_of_particle
            )

            ### Calculation accoss x-direction
            # Avoiding divisions by zero
            tempCalc1 = np.where(
                self._u_vel_of_particle == 0, 9999, self._u_vel_of_particle
            )
            tempCalc2 = np.where(u_vel_at_x1 == 0, 9999, u_vel_at_x1)
            tempCalc3 = np.where(gradient_x_direction == 0, 9999, gradient_x_direction)
            TAUx = (1 / tempCalc3) * np.log(abs(tempCalc1 / tempCalc2))

            # Calculation when gradient is equal to zero
            tempCalc4 = abs((self._x_of_particle - x_at_x1) / tempCalc2)
            TAUx = np.where(gradient_x_direction == 0, tempCalc4, TAUx)

            # Calculation when:
            # a) Uxp/Ux1 = 1,
            # b) Uxp,Vyp = 0,
            # c) Ux1,Vy1 = 0, and
            # d) Uxp/Ux1, Vxp/Vy1 = -1
            tempCalc5 = self._u_vel_of_particle / tempCalc2
            TAUx = np.where(tempCalc5 == 1, tempCalc4, TAUx)
            TAUx = np.where(self._u_vel_of_particle == 0, remaining_time, TAUx)
            TAUx = np.where(u_vel_at_x1 == 0, remaining_time, TAUx)
            TAUx = np.where(tempCalc5 < 0, remaining_time, TAUx)
            TAUx = np.where(TAUx > self._dt, self._dt, TAUx)
            TAUx = np.where(TAUx < 0, 0, TAUx)

            ### Calculation across y-direction
            # Avoiding divisions by zero
            tempCalc1 = np.where(
                self._v_vel_of_particle == 0, 9999, self._v_vel_of_particle
            )
            tempCalc2 = np.where(v_vel_at_y1 == 0, 9999, v_vel_at_y1)
            tempCalc3 = np.where(gradient_y_direction == 0, 9999, gradient_y_direction)
            TAUy = (1 / tempCalc3) * np.log(abs(tempCalc1 / tempCalc2))

            # Calculation when gradient is equal to zero
            tempCalc4 = abs((self._y_of_particle - y_at_y1) / tempCalc2)
            TAUy = np.where(gradient_y_direction == 0, tempCalc4, TAUy)

            # Calculation when
            # a) Vyp/Vy1 = 1,
            # b) Uxp,Vyp = 0,
            # c) Ux1,Vy1 = 0, and
            # d) Uxp/Ux1, Vxp/Vy1 = -1
            tempCalc5 = self._v_vel_of_particle / tempCalc2
            TAUy = np.where(tempCalc5 == 1, tempCalc4, TAUy)
            TAUy = np.where(self._v_vel_of_particle == 0, remaining_time, TAUy)
            TAUy = np.where(v_vel_at_y1 == 0, remaining_time, TAUy)
            TAUy = np.where(tempCalc5 < 0, remaining_time, TAUy)
            TAUy = np.where(TAUy > self._dt, self._dt, TAUy)
            TAUy = np.where(TAUy < 0, 0, TAUy)

            # Obtaining TAU = min(TAUx, TAUy, (dt - sum_partial_times))
            TAUx = abs(TAUx)
            TAUy = abs(TAUy)
            TAU = np.array((TAUx, TAUy, remaining_time)).min(axis=0)
            # TAU  = np.where(TAU < 1e-10, 0, TAU)

            # Calculating exit point Xe, Ye
            tempCalc1 = np.where(gradient_x_direction == 0, 9999, gradient_x_direction)
            tempCalc2 = np.where(gradient_y_direction == 0, 9999, gradient_y_direction)

            # Exit point Xe (tempCalc3) and Ye (tempCalc4)
            tempCalc3 = x_at_x2 - (1 / tempCalc1) * (
                u_vel_at_x2
                - self._u_vel_of_particle / np.exp(gradient_x_direction * TAU)
            )
            tempCalc4 = y_at_y2 - (1 / tempCalc2) * (
                v_vel_at_y2
                - self._v_vel_of_particle / np.exp(gradient_y_direction * TAU)
            )

            tempCalc3 = np.where(
                gradient_x_direction == 0,
                self._x_of_particle - u_vel_at_x2 * TAU,
                tempCalc3,
            )
            tempCalc4 = np.where(
                gradient_y_direction == 0,
                self._y_of_particle - v_vel_at_y2 * TAU,
                tempCalc4,
            )

            tempCalc3 = np.where(
                self._u_vel_of_particle == 0, self._x_of_particle, tempCalc3
            )
            tempCalc4 = np.where(
                self._v_vel_of_particle == 0, self._y_of_particle, tempCalc4
            )

            self._x_at_exit_point = np.where(
                keep_tracing, tempCalc3, self._x_at_exit_point
            )
            self._y_at_exit_point = np.where(
                keep_tracing, tempCalc4, self._y_at_exit_point
            )

            # Updating sum of partial time-steps, TAU
            sum_partial_times = np.where(
                keep_tracing, sum_partial_times + TAU, self._dt
            )

            # Checking remaining_time == 0 (dt = sum_partial_times)
            remaining_time = np.where(
                remaining_time == 0, 0, self._dt - sum_partial_times
            )

            # Correcting entry velocity
            tempCalc1 = np.where(self._u_vel_of_particle == 0, 1, 0)
            tempCalc2 = np.where(self._v_vel_of_particle == 0, 1, 0)
            remaining_time = np.where((tempCalc1 * tempCalc2) == 1, 0, remaining_time)

            # Correction for static particles
            remaining_time = np.where(
                abs(self._x_of_particle - self._x_at_exit_point) < 1e-7,
                0,
                remaining_time,
            )
            remaining_time = np.where(
                abs(self._y_of_particle - self._y_at_exit_point) < 1e-7,
                0,
                remaining_time,
            )

            # Stop tracing if a particle hits the boundary
            # Keep tracing for all particles
            tempCalc1 = np.repeat(True, len(keep_tracing))
            # Particle hits the left edge?
            tempCalc2 = np.isin(
                self._x_at_exit_point, self._grid.x_of_node[self._nodes_at_left_edge]
            )
            # If above True, stop tracing for that particle
            tempCalc1 = np.where(tempCalc2, False, tempCalc1)
            # Particle hits the right edge?
            tempCalc2 = np.isin(
                self._x_at_exit_point, self._grid.x_of_node[self._nodes_at_right_edge]
            )
            # If above True, stop tracing for that particle
            tempCalc1 = np.where(tempCalc2, False, tempCalc1)
            # Particle hits the top edge?
            tempCalc2 = np.isin(
                self._y_at_exit_point, self._grid.y_of_node[self._nodes_at_top_edge]
            )
            # If above True, stop tracing for that particle
            tempCalc1 = np.where(tempCalc2, False, tempCalc1)
            # Particle hits the bottom edge?
            tempCalc2 = np.isin(
                self._y_at_exit_point, self._grid.y_of_node[self._nodes_at_bottom_edge]
            )
            # If above True, stop tracing for that particle
            tempCalc1 = np.where(tempCalc2, False, tempCalc1)
            # Where particles reach the boundary, remaining time is equal to zero
            # remaining_time = np.where(not tempCalc1, 0, remaining_time)
            remaining_time = np.where(~tempCalc1, 0, remaining_time)

            # Updating on particles that need to traced backwards
            keep_tracing = np.where(remaining_time > 0, True, False)

            # WHILE "np.any(remaining_time > 0)" END
        # DEF "path_line_tracing(self)" END

    def run_one_step(self):
        """Calculate water depth and water velocity for a time period dt."""

        # Getting velocity as U,V components
        self._u_vel = self._vel_at_N[self._horizontal_links]
        self._v_vel = self._vel_at_N[self._vertical_links]

        # Calculating Chezy coefficient
        self._chezy_at_nodes = self._h_at_N ** (1 / 6) / self._mannings_n
        self._chezy_at_links = self._h_at_N_at_links ** (1 / 6) / self._mannings_n

        # Computing V-velocity (vertical links) at U-velocity positions (horizontal links)
        tempCalc1 = self._grid.map_mean_of_horizontal_links_to_node(self._vel_at_N)
        self._u_vel_at_v_links = np.mean(
            tempCalc1[self._grid.nodes_at_link[self._vertical_links]], axis=1
        )

        # Computing U-velocity (horizontal links) at V-velocity positions (vertical links)
        tempCalc1 = self._grid.map_mean_of_vertical_links_to_node(self._vel_at_N)
        self._v_vel_at_u_links = np.mean(
            tempCalc1[self._grid.nodes_at_link[self._horizontal_links]], axis=1
        )

        """ Setting A-faces """
        # Setting A-faces
        self._a_links = np.zeros_like(self._vel_at_N)

        # Setting dry links equal to 1 to avoid divisions by zero
        tempCalc1 = np.where(self._wet_links, self._chezy_at_links, 1)

        # Computing A-faces
        self._a_links[self._horizontal_links] = (
            self._h_at_N_at_links[self._horizontal_links]
            + self._g
            * self._dt
            * (
                self._vel_at_N[self._horizontal_links] ** 2
                + self._v_vel_at_u_links**2
            )
            ** (1 / 2)
            / tempCalc1[self._horizontal_links]
        )
        self._a_links[self._vertical_links] = (
            self._h_at_N_at_links[self._vertical_links]
            + self._g
            * self._dt
            * (self._vel_at_N[self._vertical_links] ** 2 + self._u_vel_at_v_links**2)
            ** (1 / 2)
            / tempCalc1[self._vertical_links]
        )

        # Using only wet-link values, and setting dry links equal to 1 to avoid
        # divisions by zero
        self._a_links = np.where(self._wet_links, self._a_links, 1)

        """ Path line tracing
        U-velocity, x-direction, horizontal links
        """
        # Getting the initial particle location at each volume faces
        tempB1 = [i in self._horizontal_links for i in self._active_links]
        self._x_of_particle = self._grid.xy_of_link[:, 0][self._active_links][tempB1]
        self._y_of_particle = self._grid.xy_of_link[:, 1][self._active_links][tempB1]

        # Getting the initial particle velocity
        tempB2 = [i in self._active_links for i in self._horizontal_links]
        self._u_vel_of_particle = self._u_vel[tempB2]
        self._v_vel_of_particle = self._v_vel_at_u_links[tempB2]

        # Getting a first 'exit' point to begin the loop
        self._x_at_exit_point = self._x_of_particle
        self._y_at_exit_point = self._y_of_particle

        # Calculating path line backwards on time
        self.path_line_tracing()

        """ Bicuatradic interpolation
        U-velocity, x-direction, around (p,q) location
        """
        self._UsL = np.zeros_like(self._u_vel)

        # Getting V-velocity at U-links
        temp_Vvel = np.zeros_like(self._vel_at_N)
        temp_Vvel[self._horizontal_links] = self._v_vel_at_u_links

        # Getting links around the particle and defining downstream direction based on velocity
        nearest_link_to_particle = self.find_nearest_link(
            self._x_at_exit_point, self._y_at_exit_point, objective_links="horizontal"
        )
        adjacent_links_to_particle = self.find_adjacent_links_at_link(
            nearest_link_to_particle, objective_links="horizontal"
        )

        link_at_B2 = nearest_link_to_particle
        link_at_A2 = adjacent_links_to_particle[:, 1]  # 1: N, top
        link_at_C2 = adjacent_links_to_particle[:, 3]  # 3: S, bottom
        link_at_A2 = np.where(temp_Vvel[link_at_B2] >= 0, link_at_A2, link_at_C2)
        link_at_C2 = np.where(temp_Vvel[link_at_B2] >= 0, link_at_C2, link_at_A2)
        # We avoid "-1" links close to the boundary
        link_at_A2 = np.where(link_at_A2 >= 0, link_at_A2, link_at_B2)
        link_at_C2 = np.where(link_at_C2 >= 0, link_at_C2, link_at_B2)

        # Getting the surrounding links to every particle from closest link
        # 0: E, Right
        # 2: W, Left
        link_at_A1 = self.find_adjacent_links_at_link(
            link_at_A2, objective_links="horizontal"
        )[:, 0]
        link_at_A3 = self.find_adjacent_links_at_link(
            link_at_A2, objective_links="horizontal"
        )[:, 2]

        # Selecting downstream link based on velocity direction
        link_at_A1 = np.where(self._vel_at_N[link_at_B2] >= 0, link_at_A1, link_at_A3)
        link_at_A3 = np.where(self._vel_at_N[link_at_B2] >= 0, link_at_A3, link_at_A1)

        link_at_B1 = self.find_adjacent_links_at_link(
            link_at_B2, objective_links="horizontal"
        )[
            :, 0
        ]  # 0: E, Right
        link_at_B3 = self.find_adjacent_links_at_link(
            link_at_B2, objective_links="horizontal"
        )[
            :, 2
        ]  # 2: W, Left

        # Selecting downstream link based on velocity direction
        link_at_B1 = np.where(self._vel_at_N[link_at_B2] >= 0, link_at_B1, link_at_B3)
        link_at_B3 = np.where(self._vel_at_N[link_at_B2] >= 0, link_at_B3, link_at_B1)

        link_at_C1 = self.find_adjacent_links_at_link(
            link_at_C2, objective_links="horizontal"
        )[
            :, 0
        ]  # 0: E, Right
        link_at_C3 = self.find_adjacent_links_at_link(
            link_at_C2, objective_links="horizontal"
        )[
            :, 2
        ]  # 2: W, Left

        # Selecting downstream link based on velocity direction
        link_at_C1 = np.where(self._vel_at_N[link_at_B2] >= 0, link_at_C1, link_at_C3)
        link_at_C3 = np.where(self._vel_at_N[link_at_B2] >= 0, link_at_C3, link_at_C1)

        # Getting velocity around the particle
        vel_at_A1 = np.where(
            link_at_A1 >= 0, self._vel_at_N[link_at_A1], self._vel_at_N[link_at_A2]
        )
        vel_at_A2 = self._vel_at_N[link_at_A2]
        vel_at_A3 = np.where(
            link_at_A3 >= 0, self._vel_at_N[link_at_A3], self._vel_at_N[link_at_A2]
        )

        vel_at_B1 = np.where(
            link_at_B1 >= 0, self._vel_at_N[link_at_B1], self._vel_at_N[link_at_B2]
        )
        vel_at_B2 = self._vel_at_N[link_at_B2]
        vel_at_B3 = np.where(
            link_at_B3 >= 0, self._vel_at_N[link_at_B3], self._vel_at_N[link_at_B2]
        )

        vel_at_C1 = np.where(
            link_at_C1 >= 0, self._vel_at_N[link_at_C1], self._vel_at_N[link_at_C2]
        )
        vel_at_C2 = self._vel_at_N[link_at_C2]
        vel_at_C3 = np.where(
            link_at_C3 >= 0, self._vel_at_N[link_at_C3], self._vel_at_N[link_at_C2]
        )

        # Getting coordinates around the particle
        x_at_2 = self._grid.xy_of_link[link_at_B2][:, 0]
        x_at_1 = np.where(
            self._vel_at_N[link_at_B2] >= 0, x_at_2 + self._dx, x_at_2 - self._dx
        )
        x_at_3 = np.where(
            self._vel_at_N[link_at_B2] >= 0, x_at_2 - self._dx, x_at_2 + self._dx
        )

        y_at_B = self._grid.xy_of_link[link_at_B2][:, 1]
        y_at_A = np.where(
            temp_Vvel[link_at_B2] >= 0, y_at_B + self._dy, y_at_B - self._dy
        )
        y_at_C = np.where(
            temp_Vvel[link_at_B2] >= 0, y_at_B - self._dy, y_at_B + self._dy
        )

        # Calculating the weights W(i,j) for k around x-direction
        W1 = (
            (self._x_at_exit_point - x_at_2)
            * (self._x_at_exit_point - x_at_3)
            / ((x_at_1 - x_at_2) * (x_at_1 - x_at_3))
        )
        W2 = (
            (self._x_at_exit_point - x_at_1)
            * (self._x_at_exit_point - x_at_3)
            / ((x_at_2 - x_at_1) * (x_at_2 - x_at_3))
        )
        W3 = (
            (self._x_at_exit_point - x_at_1)
            * (self._x_at_exit_point - x_at_2)
            / ((x_at_3 - x_at_1) * (x_at_3 - x_at_2))
        )

        # Interpolation by row around 'x_at_exit_point'
        A = W1 * vel_at_A1 + W2 * vel_at_A2 + W3 * vel_at_A3
        B = W1 * vel_at_B1 + W2 * vel_at_B2 + W3 * vel_at_B3
        C = W1 * vel_at_C1 + W2 * vel_at_C2 + W3 * vel_at_C3

        # Calculating the weghts W(i,j) for l around y-direction
        W1 = (
            (self._y_at_exit_point - y_at_B)
            * (self._y_at_exit_point - y_at_C)
            / ((y_at_A - y_at_B) * (y_at_A - y_at_C))
        )
        W2 = (
            (self._y_at_exit_point - y_at_A)
            * (self._y_at_exit_point - y_at_C)
            / ((y_at_B - y_at_A) * (y_at_B - y_at_C))
        )
        W3 = (
            (self._y_at_exit_point - y_at_A)
            * (self._y_at_exit_point - y_at_B)
            / ((y_at_C - y_at_A) * (y_at_C - y_at_B))
        )

        # Calculating UsL by bicuadratic interpolation
        self._UsL[tempB2] = W1 * A + W2 * B + W3 * C

        """ Computing viscous terms
        U-located particles
        """
        # Central difference scheme around x- and y- direction for U-located particles
        self._Uvis = np.zeros_like(self._u_vel)

        tempCalc1 = (
            self._eddy_viscosity
            * self._dt
            * (vel_at_B3 - 2 * vel_at_B2 + vel_at_B1)
            / (self._dx**2)
        )
        tempCalc2 = (
            self._eddy_viscosity
            * self._dt
            * (vel_at_C2 - 2 * vel_at_B2 + vel_at_A2)
            / (self._dy**2)
        )

        self._Uvis[tempB2] = tempCalc1 + tempCalc2

        """ Path line tracing
        V-velocity, y-direction, vertical links
        """
        # Getting the initial particle location at each volume faces
        tempB1 = [j in self._vertical_links for j in self._active_links]
        self._x_of_particle = self._grid.xy_of_link[:, 0][self._active_links][tempB1]
        self._y_of_particle = self._grid.xy_of_link[:, 1][self._active_links][tempB1]

        # Getting the initial particle velocity
        tempB2 = [j in self._active_links for j in self._vertical_links]
        self._v_vel_of_particle = self._v_vel[tempB2]
        self._u_vel_of_particle = self._u_vel_at_v_links[tempB2]

        # Getting a first 'exit' point to begin the loop
        self._x_at_exit_point = self._x_of_particle
        self._y_at_exit_point = self._y_of_particle

        # Calculating path line backwards on time
        self.path_line_tracing()

        """ Bicuatradic interpolation
        V-velocity, y-direction, around (p,q) location
        """
        self._VsL = np.zeros_like(self._v_vel)

        # Getting V-velocity at U-links
        temp_Uvel = np.zeros_like(self._vel_at_N)
        temp_Uvel[self._vertical_links] = self._u_vel_at_v_links

        # Getting links around the particle and defining downstream direction based on velocity
        nearest_link_to_particle = self.find_nearest_link(
            self._x_at_exit_point, self._y_at_exit_point, objective_links="vertical"
        )
        adjacent_links_to_particle = self.find_adjacent_links_at_link(
            nearest_link_to_particle, objective_links="vertical"
        )

        link_at_B2 = nearest_link_to_particle
        link_at_A2 = adjacent_links_to_particle[:, 1]  # 1: N, top
        link_at_C2 = adjacent_links_to_particle[:, 3]  # 3: S, bottom
        link_at_A2 = np.where(self._vel_at_N[link_at_B2] >= 0, link_at_A2, link_at_C2)
        link_at_C2 = np.where(self._vel_at_N[link_at_B2] >= 0, link_at_C2, link_at_A2)
        link_at_A2 = np.where(link_at_A2 >= 0, link_at_A2, link_at_B2)
        link_at_C2 = np.where(link_at_C2 >= 0, link_at_C2, link_at_B2)
        # We avoid "-1" links close to the boundary

        # Getting the surrounding links to every particle from closest link
        link_at_A1 = self.find_adjacent_links_at_link(
            link_at_A2, objective_links="vertical"
        )[
            :, 0
        ]  # 0: E, Left
        link_at_A3 = self.find_adjacent_links_at_link(
            link_at_A2, objective_links="vertical"
        )[
            :, 2
        ]  # 2: W, Right

        # Selecting downstream link based on velocity direction
        link_at_A1 = np.where(temp_Uvel[link_at_B2] >= 0, link_at_A1, link_at_A3)
        link_at_A3 = np.where(temp_Uvel[link_at_B2] >= 0, link_at_A3, link_at_A1)

        link_at_B1 = self.find_adjacent_links_at_link(
            link_at_B2, objective_links="vertical"
        )[
            :, 0
        ]  # 0: E, Left
        link_at_B3 = self.find_adjacent_links_at_link(
            link_at_B2, objective_links="vertical"
        )[
            :, 2
        ]  # 2: W, Right

        # Selecting downstream link based on velocity direction
        link_at_B1 = np.where(temp_Uvel[link_at_B2] >= 0, link_at_B1, link_at_B3)
        link_at_B3 = np.where(temp_Uvel[link_at_B2] >= 0, link_at_B3, link_at_B1)

        link_at_C1 = self.find_adjacent_links_at_link(
            link_at_C2, objective_links="vertical"
        )[
            :, 0
        ]  # 0: E, Left
        link_at_C3 = self.find_adjacent_links_at_link(
            link_at_C2, objective_links="vertical"
        )[
            :, 2
        ]  # 2: W, Right

        # Selecting downstream link based on velocity direction
        link_at_C1 = np.where(temp_Uvel[link_at_B2] >= 0, link_at_C1, link_at_C3)
        link_at_C3 = np.where(temp_Uvel[link_at_B2] >= 0, link_at_C3, link_at_C1)

        # Getting velocity around the particle
        vel_at_A1 = np.where(
            link_at_A1 >= 0, self._vel_at_N[link_at_A1], self._vel_at_N[link_at_A2]
        )
        vel_at_A2 = self._vel_at_N[link_at_A2]
        vel_at_A3 = np.where(
            link_at_A3 >= 0, self._vel_at_N[link_at_A3], self._vel_at_N[link_at_A2]
        )

        vel_at_B1 = np.where(
            link_at_B1 >= 0, self._vel_at_N[link_at_B1], self._vel_at_N[link_at_B2]
        )
        vel_at_B2 = self._vel_at_N[link_at_B2]
        vel_at_B3 = np.where(
            link_at_B3 >= 0, self._vel_at_N[link_at_B3], self._vel_at_N[link_at_B2]
        )

        vel_at_C1 = np.where(
            link_at_C1 >= 0, self._vel_at_N[link_at_C1], self._vel_at_N[link_at_C2]
        )
        vel_at_C2 = self._vel_at_N[link_at_C2]
        vel_at_C3 = np.where(
            link_at_C3 >= 0, self._vel_at_N[link_at_C3], self._vel_at_N[link_at_C2]
        )

        # Getting coordinates around the particle
        x_at_2 = self._grid.xy_of_link[link_at_B2][:, 0]
        x_at_1 = np.where(
            temp_Uvel[link_at_B2] >= 0, x_at_2 + self._dx, x_at_2 - self._dx
        )
        x_at_3 = np.where(
            temp_Uvel[link_at_B2] >= 0, x_at_2 - self._dx, x_at_2 + self._dx
        )

        y_at_B = self._grid.xy_of_link[link_at_B2][:, 1]
        y_at_A = np.where(
            self._vel_at_N[link_at_B2] >= 0, y_at_B + self._dy, y_at_B - self._dy
        )
        y_at_C = np.where(
            self._vel_at_N[link_at_B2] >= 0, y_at_B - self._dy, y_at_B + self._dy
        )

        # Calculating the weights W(i,j) for k around x-direction
        W1 = (
            (self._x_at_exit_point - x_at_2)
            * (self._x_at_exit_point - x_at_3)
            / ((x_at_1 - x_at_2) * (x_at_1 - x_at_3))
        )
        W2 = (
            (self._x_at_exit_point - x_at_1)
            * (self._x_at_exit_point - x_at_3)
            / ((x_at_2 - x_at_1) * (x_at_2 - x_at_3))
        )
        W3 = (
            (self._x_at_exit_point - x_at_1)
            * (self._x_at_exit_point - x_at_2)
            / ((x_at_3 - x_at_1) * (x_at_3 - x_at_2))
        )

        # Interpolation by row around 'y_at_exit_point'
        A = W1 * vel_at_A1 + W2 * vel_at_A2 + W3 * vel_at_A3
        B = W1 * vel_at_B1 + W2 * vel_at_B2 + W3 * vel_at_B3
        C = W1 * vel_at_C1 + W2 * vel_at_C2 + W3 * vel_at_C3

        # Calculating the weghts W(i,j) for l around y-direction
        W1 = (
            (self._y_at_exit_point - y_at_B)
            * (self._y_at_exit_point - y_at_C)
            / ((y_at_A - y_at_B) * (y_at_A - y_at_C))
        )
        W2 = (
            (self._y_at_exit_point - y_at_A)
            * (self._y_at_exit_point - y_at_C)
            / ((y_at_B - y_at_A) * (y_at_B - y_at_C))
        )
        W3 = (
            (self._y_at_exit_point - y_at_A)
            * (self._y_at_exit_point - y_at_B)
            / ((y_at_C - y_at_A) * (y_at_C - y_at_B))
        )

        # Calculating VsL by bicuadratic interpolation
        self._VsL[tempB2] = W1 * A + W2 * B + W3 * C

        """ Computing viscous terms
        V-located particles
        """
        # Central difference scheme around x- and y- direction for V-located particles
        self._Vvis = np.zeros_like(self._v_vel)

        tempCalc1 = (
            self._eddy_viscosity
            * self._dt
            * (vel_at_B3 - 2 * vel_at_B2 + vel_at_B1)
            / (self._dx**2)
        )
        tempCalc2 = (
            self._eddy_viscosity
            * self._dt
            * (vel_at_C2 - 2 * vel_at_B2 + vel_at_A2)
            / (self._dy**2)
        )

        self._Vvis[tempB2] = tempCalc1 + tempCalc2

        """ Computing advective terms F(U,V)
        """
        # Computing advective terms (FU, FV)
        self._f_vel = np.zeros_like(self._vel_at_N)

        # Adding semi-lagrangian and viscous terms
        tempCalc1 = self._UsL + self._Uvis
        tempCalc2 = self._VsL + self._Vvis

        # Including the results according with links directions
        self._f_vel[self._horizontal_links] = tempCalc1
        self._f_vel[self._vertical_links] = tempCalc2

        """ Setting G-faces
        """
        # Setting G-faces
        self._g_links = np.zeros_like(self._vel_at_N)

        # Computing G-faces
        self._g_links = self._h_at_N_at_links * self._f_vel - self._h_at_N_at_links * (
            1 - self._theta
        ) * self._g * self._dt / self._dx * self._grid.calc_diff_at_link(self._eta_at_N)

        # Using only wet-link values, and setting dry links equal to 0 to avoid
        # using wrong values
        self._g_links = np.where(self._wet_links, self._g_links, 0)
        self._g_links = np.where(
            self._grid.status_at_link == 4, 0, self._g_links
        )  # link is active (0), fixed (2), inactive (4)

        """ Solving semi-implicit scheme with PCG method
        """
        # Building the system of equations 'A*x=b'
        A = np.zeros(
            (self._number_of_nodes, self._number_of_nodes)
        )  # Full 'A' matrix with all nodes on it
        b = np.zeros(self._number_of_nodes)  # Full 'b' vector with all nodes on it

        # Getting surrounding locations for core nodes
        adjacent_nodes = self._grid.adjacent_nodes_at_node[
            self._core_nodes
        ]  # East, North, West, South
        adjacent_links = self._grid.links_at_node[
            self._core_nodes
        ]  # East, North, West, South
        nodes_location = np.append(
            adjacent_nodes, np.array([self._core_nodes]).T, axis=1
        )  # East, North, West, South, Center

        # Boolean to differentiate between core and boundary nodes
        tempB1 = np.isin(nodes_location, self._core_nodes)
        # Core node if tempB1 == True, boundary node if tempB1 == False
        tempB2 = ~tempB1
        # Boundary node if tempB2 == True, core node if tempB2 == False

        ## Building 'b' for the right-hand side of the system
        # Calculating 'delta' to build 'b'
        tempCalc1 = (
            self._h_at_N_at_links[adjacent_links] * self._vel_at_N[adjacent_links]
        )
        tempCalc2 = (
            self._eta_at_N[self._core_nodes]
            - (1 - self._theta)
            * self._dt
            / self._dx
            * (tempCalc1[:, 0] - tempCalc1[:, 2])
            - (1 - self._theta)
            * self._dt
            / self._dy
            * (tempCalc1[:, 1] - tempCalc1[:, 3])
        )

        tempCalc1 = (
            self._h_at_N_at_links[adjacent_links]
            * self._g_links[adjacent_links]
            / self._a_links[adjacent_links]
        )
        b[self._core_nodes] = (
            tempCalc2
            - self._theta * self._dt / self._dx * (tempCalc1[:, 0] - tempCalc1[:, 2])
            - self._theta * self._dt / self._dy * (tempCalc1[:, 1] - tempCalc1[:, 3])
        )

        ## Building 'A' for the left-hand side of the system
        # Calculating coefficients for the system of equations
        tempCalc1 = (
            self._h_at_N_at_links[adjacent_links] ** 2 / self._a_links[adjacent_links]
        )
        coefficients = [
            -tempCalc1[:, 0] * (self._g * self._theta * self._dt / self._dx) ** 2,
            -tempCalc1[:, 1] * (self._g * self._theta * self._dt / self._dy) ** 2,
            -tempCalc1[:, 2] * (self._g * self._theta * self._dt / self._dx) ** 2,
            -tempCalc1[:, 3] * (self._g * self._theta * self._dt / self._dy) ** 2,
            1
            + (tempCalc1[:, 0] + tempCalc1[:, 2])
            * (self._g * self._theta * self._dt / self._dx) ** 2
            + (tempCalc1[:, 1] + tempCalc1[:, 3])
            * (self._g * self._theta * self._dt / self._dy) ** 2,
        ]
        coefficients = np.array(coefficients).T

        ## For loop iterates through every row of 'nodes_location' to:
        # a) find the node number (row in A) and choose an equation, and
        # b) find the columns associated with adjacent nodes
        for row in range(nodes_location.shape[0]):
            # Getting the current node
            current_node = nodes_location[row, 4]
            # Getting the row associated with the current node
            current_rows_in_A = np.array([current_node])
            # Getting the columns associated with surrounding nodes
            current_cols_in_A = nodes_location[row, :]

            # Filling A matrix with coefficients associated with surrounding nodes
            A[np.ix_(current_rows_in_A, current_cols_in_A)] = (
                coefficients[row, :] * tempB1[row, :]
            )

            # Adding known terms (boundary nodes) to the right-hand side of the equation
            b[current_rows_in_A] = b[current_rows_in_A] - sum(
                self._eta_at_N[nodes_location[row, :]]
                * coefficients[row, :]
                * tempB2[row, :]
            )
        # for END

        # Extracting only core nodes to be solved
        left_hand_side = A[np.ix_(self._core_nodes, self._core_nodes)]
        right_hand_side = b[self._core_nodes]

        # Applying PCG method to 'LHS*eta=RHS' using np.diag() as a preconditioner for 'LHS'
        # Preconditioned conjugated gradient output flag:
        # 0  : successful exit
        # >0 : convergence to tolerance not achieved, number of iterations
        # <0 : illegal input or breakdown
        # Alternative preconditiner: Li = np.linalg.cholesky(left_hand_side)
        Mi = np.diag(np.diag(left_hand_side))
        pcg_results = sp.sparse.linalg.cg(
            left_hand_side,
            right_hand_side,
            M=Mi,
            tol=self._pcg_tolerance,
            maxiter=self._pcg_max_iterations,
            atol=0,
        )

        # Getting the new water surface elevation
        self._eta = np.zeros_like(self._eta_at_N)
        self._eta[self._core_nodes] = pcg_results[0]

        """ Boundary conditions
        Radiation Boundary Conditions of Roed & Smedstad (1984) applied on open boundaries
        Water surface elevation
        """
        ## Updating the new WSE ('eta') with the fixed nodes values
        if self._fixed_nodes_exist is True:
            self._eta[self._fixed_entry_nodes] = (
                self._entry_nodes_h_values - self._z[self._fixed_entry_nodes]
            )

        ## Getting the 1-line-upstream nodes from boundary nodes
        tempCalc1 = self._grid.active_adjacent_nodes_at_node[self._open_boundary_nodes]
        open_boundary_nodes_1_backwards = np.extract(tempCalc1 >= 0, tempCalc1)

        ## Getting the 2-line-upstream nodes from boundary nodes
        # Looking for these nodes
        tempCalc1 = np.tile(self._open_boundary_nodes, (4, 1)).T
        tempCalc2 = self._grid.active_adjacent_nodes_at_node[
            open_boundary_nodes_1_backwards
        ]  # Surrounding nodes to 1-line-upstream

        # Getting the node positions to extract them from all the surrounding nodes
        tempB1 = np.tile([0, 1, 2, 3], (len(self._open_boundary_nodes), 1))
        tempB2 = tempB1[tempCalc1 == tempCalc2]

        # It gives me the node indices to extract
        # folowing the face direction
        tempB1 = np.where(tempB2 == 0, 2, tempB2)
        tempB1 = np.where(tempB2 == 1, 3, tempB1)
        tempB1 = np.where(tempB2 == 2, 0, tempB1)
        tempB1 = np.where(tempB2 == 3, 1, tempB1)

        open_boundary_nodes_2_backwards = tempCalc2[
            [range(tempCalc2.shape[0])], tempB1
        ][0]

        ## Getting WSE at different locations
        eta_at_N_at_B = self._eta_at_N[self._open_boundary_nodes]
        eta_at_N_at_B_1 = self._eta_at_N[open_boundary_nodes_1_backwards]
        eta_at_N_1_at_B_1 = self._eta_at_N_1[open_boundary_nodes_1_backwards]
        eta_at_N_1_at_B_2 = self._eta_at_N_1[open_boundary_nodes_2_backwards]

        ## Computing boundary condition
        tempCalc1 = eta_at_N_at_B_1 - eta_at_N_1_at_B_1
        tempCalc2 = eta_at_N_1_at_B_1 - eta_at_N_1_at_B_2
        tempCalc3 = np.where(tempCalc2 == 0, 1, tempCalc2)

        Ce = np.where(tempCalc2 == 0, 0, tempCalc1 / tempCalc3 * (-self._dx / self._dt))
        Ce = np.where(tempCalc2 == 0, 0, Ce)

        # eta[open_boundary_nodes] = tempCalc1/tempCalc2
        self._eta[self._open_boundary_nodes] = np.where(
            Ce >= 0, eta_at_N_at_B_1, eta_at_N_at_B
        )
        self._eta = np.where(
            abs(self._eta) > abs(self._z), -self._z, self._eta
        )  # Correcting WSE below topographic elevation

        self._eta_at_links = self._grid.map_mean_of_link_nodes_to_link(self._eta)

        # Corner nodes treatment
        self._eta[self._corner_nodes] = np.mean(
            self._eta[self._adjacent_nodes_at_corner_nodes], axis=1
        )

        """ Updating water velocity
        """
        # tempB1 :  Considering only wet cells
        # tempB2 :  Cells with elevation below the water surface elevation
        tempB1 = np.where(
            abs(self._eta[self._grid.node_at_link_head])
            <= abs(self._z[self._grid.node_at_link_head] - self._threshold_depth),
            1,
            0,
        )
        tempB2 = np.where(
            abs(self._z[self._grid.node_at_link_head])
            > abs(self._eta[self._grid.node_at_link_tail]),
            1,
            0,
        )
        tempB1 = tempB1 + tempB2
        tempB1 = np.where(tempB1 > 1, 1, 0)

        # Updating water velocity
        tempCalc1 = (
            self._theta
            * self._g
            * self._dt
            / self._dx
            * (self._grid.calc_diff_at_link(self._eta))
            * self._h_at_N_at_links
            / self._a_links
        )
        self._vel = self._g_links / self._a_links - tempCalc1 * tempB1

        # Only updating velocity on wet cells
        self._vel = np.where(self._wet_links, self._vel, 0)

        """ Boundary conditions
        Radiation Boundary Conditions of Roed & Smedstad (1984) applied on open boundaries
        Water velocity
        """
        ## Updating the new Velocity with the fixed links values
        if self._fixed_links_exist is True:
            self._vel[self._fixed_entry_links] = self._entry_links_vel_values

        ## Getting the boundary links
        tempB1 = [i in self._open_boundary_links for i in self._active_links]
        open_boundary_active_links = self._grid.active_links[tempB1]

        ## Getting the 1-line-upstream links from boundary links
        tempCalc1 = np.tile(open_boundary_active_links, (4, 1)).T
        tempCalc2 = self._grid.links_at_node[open_boundary_nodes_1_backwards]

        # Getting the link positions to extract them from all the surrounding nodes
        tempB1 = np.tile([0, 1, 2, 3], (len(self._open_boundary_nodes), 1))
        tempB2 = tempB1[tempCalc1 == tempCalc2]

        # It gives me the link indices to extract
        # folowing the face direction
        tempB1 = np.where(tempB2 == 0, 2, tempB2)
        tempB1 = np.where(tempB2 == 1, 3, tempB1)
        # tempB1 is where the target link is located
        tempB1 = np.where(tempB2 == 2, 0, tempB1)
        # tempB2 is where the upstream link is located
        tempB1 = np.where(tempB2 == 3, 1, tempB1)

        open_boundary_active_links_1_backwards = tempCalc2[
            [range(tempCalc2.shape[0])], tempB1
        ][0]

        ### Getting the 2-line-upstream links from boundary nodes
        tempCalc1 = np.tile(open_boundary_active_links_1_backwards, (4, 1)).T
        tempCalc2 = self._grid.links_at_node[open_boundary_nodes_2_backwards]

        # Getting the link positions to extract them from all the surrounding nodes
        tempB1 = np.tile([0, 1, 2, 3], (len(self._open_boundary_nodes), 1))
        tempB2 = tempB1[(tempCalc1 == tempCalc2)]

        # It gives me the link indices to extract
        # folowing the face direction
        tempB1 = np.where(tempB2 == 0, 2, tempB2)
        tempB1 = np.where(tempB2 == 1, 3, tempB1)
        # tempB1 is where the target link is located
        tempB1 = np.where(tempB2 == 2, 0, tempB1)
        # tempB2 is where the upstream link is located
        tempB1 = np.where(tempB2 == 3, 1, tempB1)

        open_boundary_active_links_2_backwards = tempCalc2[
            [range(tempCalc2.shape[0])], tempB1
        ][0]

        ## Getting water velocity at different locations
        vel_at_N_at_B = self._vel_at_N[open_boundary_active_links]
        vel_at_N_at_B_1 = self._vel_at_N[open_boundary_active_links_1_backwards]
        vel_at_N_1_at_B_1 = self._vel_at_N_1[open_boundary_active_links_1_backwards]
        vel_at_N_1_at_B_2 = self._vel_at_N_1[open_boundary_active_links_2_backwards]

        ## Computing boundary condition
        tempCalc1 = vel_at_N_at_B_1 - vel_at_N_1_at_B_1
        tempCalc2 = vel_at_N_1_at_B_1 - vel_at_N_1_at_B_2
        tempCalc3 = np.where(tempCalc2 == 0, 1, tempCalc2)

        Ce = np.where(tempCalc2 == 0, 0, tempCalc1 / tempCalc3 * (-self._dx / self._dt))
        Ce = np.where(tempCalc2 == 0, 0, Ce)

        self._vel[open_boundary_active_links] = np.where(
            Ce >= 0, vel_at_N_at_B_1, vel_at_N_at_B
        )

        """ Updating water depth at links
        """
        # Using only values where the WSE is above the topographic elevation
        tempB1 = np.where(abs(self._eta) <= abs(self._z - self._threshold_depth), 1, 0)

        # Updating water depth at links
        tempCalc1 = (
            self._z_at_links + self._eta[self._grid.node_at_link_head]
        ) * tempB1[self._grid.node_at_link_head]
        tempCalc2 = (
            self._z_at_links + self._eta[self._grid.node_at_link_tail]
        ) * tempB1[self._grid.node_at_link_tail]
        tempCalc3 = np.zeros_like(self._h_at_N_at_links)

        self._h_at_links = np.array((tempCalc1, tempCalc2, tempCalc3)).max(axis=0)

        # Applying boundary condition at links
        self._h_at_links[self._open_boundary_links] = (
            self._z_at_links[self._open_boundary_links]
            + self._eta_at_links[self._open_boundary_links]
        )

        # Wet cells threshold
        self._h_at_links = np.where(
            self._h_at_links < self._threshold_depth, 0, self._h_at_links
        )

        # Updating wet links
        self._wet_links = np.where(
            self._h_at_links >= self._threshold_depth, True, False
        )
        self._vel = self._vel * self._wet_links

        """ Calculating average water depth at nodes
        """
        # If a node is dry, using only surrounding links such that 'WSE' is above 'z'
        # If a node is wet, using all surrounding links even if 'WSE' is below 'z' (jumps)

        # Checking surrounding wet links
        surrounding_links = self._grid.links_at_node[self._core_nodes]

        # Checking whether the core node is wet (T) or dry (F)
        tempB1 = abs(self._eta[self._core_nodes]) < abs(
            self._z[self._core_nodes] - self._threshold_depth
        )

        # Checking whether surrounding links are wet (T) or dry (F)
        tempB2 = self._wet_links[surrounding_links]

        # Checking whether surrounding 'WSE' links are above (T) or below (F) 'z' at nodes
        tempB3 = (
            abs(self._eta_at_links[surrounding_links])
            < abs(self._z[self._core_nodes] - self._threshold_depth)[:, None]
        )

        # Getting the number of wet links around each core node, satisfying tempB2,
        # and avoiding divisions by zero
        tempCalc2 = np.sum(tempB2 * 1, axis=1)
        tempCalc2 = np.where(tempCalc2 == 0, -9999, tempCalc2)

        # Getting the number of wet links around each core node, satisfying tempB3,
        # and avoiding divisions by zero
        tempCalc3 = np.sum(tempB2 * tempB3 * 1, axis=1)
        tempCalc3 = np.where(tempCalc3 == 0, -9999, tempCalc3)

        # Updating water depth
        # h = h_at_N - rmg.calc_net_flux_at_node(h_at_links*vel) # (influx if negative)
        self._h[self._core_nodes] = np.where(
            tempCalc3 > 0,
            np.sum(self._h_at_links[surrounding_links] * tempB2 * tempB3, axis=1)
            / tempCalc3,
            0,
        )  # Dry nodes, tempB1 == False

        ### Updating boundary nodes
        if self._fixed_nodes_exist is True:
            self._h[self._fixed_entry_nodes] = self._entry_nodes_h_values
        self._h[self._open_boundary_nodes] = (
            self._eta[self._open_boundary_nodes] + self._z[self._open_boundary_nodes]
        )
        self._h = np.where(self._h < self._threshold_depth, 0, self._h)

        # Corner nodes treatment
        self._h[self._corner_nodes] = np.mean(
            self._h[self._adjacent_nodes_at_corner_nodes], axis=1
        )

        # Updating wet nodes
        self._wet_nodes = np.where(self._h >= self._threshold_depth, True, False)

        """ Storing values in the grid
        """
        self._grid.at_node["surface_water__depth"] = self._h
        self._grid.at_link["surface_water__velocity"] = self._vel
        self._grid.at_node["surface_water__elevation"] = self._eta + (
            self._max_elevation + self._additional_z
        )

        self._grid.at_link["surface_water__velocity_at_N-1"] = self._vel_at_N
        self._grid.at_node["surface_water__elevation_at_N-1"] = self._eta_at_N + (
            self._max_elevation + self._additional_z
        )
        self._grid.at_node["surface_water__elevation_at_N-2"] = self._eta_at_N_1 + (
            self._max_elevation + self._additional_z
        )

        """ Storing values at previous time steps
        """
        self._eta_at_N = self._eta.copy()
        self._eta_at_N_1 = self._eta_at_N.copy()
        self._eta_at_N_2 = self._eta_at_N_1.copy()

        self._h_at_N = self._h.copy()
        self._h_at_N_at_links = self._h_at_links.copy()

        self._vel_at_N = self._vel.copy()
        self._vel_at_N_1 = self._vel_at_N.copy()
