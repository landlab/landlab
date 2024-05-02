#!/usr/bin/env python
"""Rock uplift along a normal fault.

Landlab component that implements rock uplift by a normal fault. Note
that this component does not make any attempt to advect topography
laterally.
"""

import numpy as np

from landlab import Component
from landlab import FieldError

TWO_PI = 2.0 * np.pi


class NormalFault(Component):
    """NormalFault implements relative rock motion due to a normal fault.

    The fault can have an arbitrary trace given by two points (x1, y1) and
    (x2, y2) in the `fault_trace` input parameter. These value of these points
    is in model-space coordinates and is not based on node id values or number
    of rows and columns.

    This NormalFault component permits two primary methods for enacting fault
    motion.

    1. **run_one_step**: The throw rate is provided through the
       ``fault_throw_rate_through_time`` parameter. This rate can be constant or
       arbitrary. See the NormalFault tutorial in the landlab tutorials repository
       for an extensive example. In this case, the NormalFault component will keep
       track of the cumulative amount of model-run-time and set the rate based on
       interpolating the provided rate-time history. *NOTE: this means that the
       model run timesteps must align with the time-rate relationship provided to
       NormalFault*. Improving this is on the developers todo list but is of low
       priority.

    2. **run_one_earthquake**: A single uplift event of size dz can be
       specified by this method. If NormalFault is used in this way, any
       specifications provided in the ``fault_throw_rate_through_time`` keyword
       argument will be ignored.

    Note that the NormalFault component does not prevent a user from combining
    the **run_one_step** and **run_one_earthquake** methods. It is encumbent
    upon the user, however, to ensure that these two methods are used in
    combination correctly for their specific use case.

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    None Listed
    """

    _name = "NormalFault"

    _unit_agnostic = True

    _info = {
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        }
    }

    def __init__(
        self,
        grid,
        faulted_surface="topographic__elevation",
        fault_throw_rate_through_time=(("time", [0]), ("rate", [0.001])),
        fault_dip_angle=90.0,
        fault_trace=(("x1", 0), ("y1", 0), ("x2", 1), ("y2", 1)),
        include_boundaries=False,
    ):
        """Instantiation of a NormalFault.

        Parameters
        ----------
        grid : ModelGrid
        faulted_surface : str or list of str
            Surface that is modified by the NormalFault component. Must be a
            field name or a list of field names if the fault should uplift more
            than one field. Default value is `topographic__elevation`.
            If the faulted surface does not yet exist, it will be ingored. The
            ``run_one_step`` method will check to see an ignored field has been
            added and if it has been, it will modify it.
        fault_throw_rate_through_time : dict, optional
            Dictionary that specifies the time varying throw rate on the fault.
            Expected format is:
            `fault_throw_rate_through_time = {'time': array, 'rate': array}`
            Default value is a constant rate of 0.001 (units not specified).
            This is acomplished by providing the dictionary
            `{'time': [0], 'rate': [0.001]}`. NormalFault uses numpy.interp
            to interpolate the time and rate pattern to the current model time.
            This function uses the first value for all values less than the
            first, and the last value for all values greater than the last, and
            thus providing only one number results in all times getting a rate
            of that value.
        fault_dip_angle : float, optional
            Dip angle of the fault in degrees.  Default value is 90 degrees.
        fault_trace : dictionary, optional
            Dictionary that specifies the coordinates of two locations on the
            fault trace. Expected format is

            .. code-block:: python

                fault_trace = {"x1": float, "y1": float, "x2": float, "y2": float}

            where the vector from `(x1, y1)` to `(x2, y2)` defines the
            strike of the fault trace. The orientation of the fault dip relative
            to the strike follows the right hand rule.
            Default is for the fault to strike NE.
        include_boundaries : boolean, optional
            Flag to indicate if model grid boundaries should be uplifted. If
            set to `True` uplifted model grid boundaries will be set to the
            average value of their upstream nodes. Default value is False.

        Examples
        --------
        Create a grid on which we will run the NormalFault component.

        >>> from landlab import RasterModelGrid
        >>> from landlab.components import NormalFault
        >>> grid = RasterModelGrid((6, 6), xy_spacing=10)

        Add an elevation field.

        >>> z = grid.add_zeros("topographic__elevation", at="node")

        Set the parameter values for the NormalFault component.

        >>> param_dict = {
        ...     "faulted_surface": "topographic__elevation",
        ...     "fault_dip_angle": 90.0,
        ...     "fault_throw_rate_through_time": {
        ...         "time": [0, 9, 10],
        ...         "rate": [0, 0, 0.05],
        ...     },
        ...     "fault_trace": {"y1": 0, "x1": 0, "y2": 30, "x2": 60},
        ...     "include_boundaries": False,
        ... }

        Instantiate a NormalFault component.

        >>> nf = NormalFault(grid, **param_dict)
        >>> nf.faulted_nodes.reshape(grid.shape)
        array([[False, False, False, False, False, False],
               [False,  True, False, False, False, False],
               [False,  True,  True,  True, False, False],
               [False,  True,  True,  True,  True, False],
               [False,  True,  True,  True,  True, False],
               [False, False, False, False, False, False]])

        As we can see, only a subset of the nodes have been identified as
        *faulted nodes*. Because we have set include_boundaries' to False none
        of the boundary nodes are faulted nodes.

        Next we will run the NormalFault for 30 1-year timesteps.

        >>> dt = 1.0
        >>> for i in range(30):
        ...     nf.run_one_step(dt)
        ...
        >>> z.reshape(grid.shape)
        array([[0., 0., 0., 0., 0., 0.],
               [0., 1., 0., 0., 0., 0.],
               [0., 1., 1., 1., 0., 0.],
               [0., 1., 1., 1., 1., 0.],
               [0., 1., 1., 1., 1., 0.],
               [0., 0., 0., 0., 0., 0.]])

        This results in uplift of the faulted nodes, as we would expect.

        If the user knows how much uplift (dz) they want to occur in an event,
        they can use the **run_one_earthquake** function with a specified dz.
        In this case fault_throw_rate_through_time will be ignored.

        >>> nf.run_one_earthquake(dz=100)
        >>> z.reshape(grid.shape)
        array([[  0.,   0.,   0.,   0.,   0.,   0.],
               [  0., 101.,   0.,   0.,   0.,   0.],
               [  0., 101., 101., 101.,   0.,   0.],
               [  0., 101., 101., 101., 101.,   0.],
               [  0., 101., 101., 101., 101.,   0.],
               [  0.,   0.,   0.,   0.,   0.,   0.]])

        Next, we make a very simple landscape model. We need a few components
        and we will set include_boundaries to True.

        >>> from landlab.components import FastscapeEroder, FlowAccumulator
        >>> grid = RasterModelGrid((6, 6), xy_spacing=10)
        >>> z = grid.add_zeros("topographic__elevation", at="node")
        >>> param_dict = {
        ...     "faulted_surface": "topographic__elevation",
        ...     "fault_dip_angle": 90.0,
        ...     "fault_throw_rate_through_time": {
        ...         "time": [0, 900, 1000],
        ...         "rate": [0, 0, 0.05],
        ...     },
        ...     "fault_trace": {"y1": 0, "x1": 0, "y2": 30, "x2": 60},
        ...     "include_boundaries": True,
        ... }

        >>> nf = NormalFault(grid, **param_dict)
        >>> fr = FlowAccumulator(grid)
        >>> fs = FastscapeEroder(grid, K_sp=0.01)

        Run this model for 300 100-year timesteps.

        >>> dt = 100.0
        >>> for i in range(300):
        ...     nf.run_one_step(dt)
        ...     fr.run_one_step()
        ...     fs.run_one_step(dt)
        ...
        >>> z.reshape(grid.shape).round(decimals=2)
        array([[ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
               [ 5.  ,  5.  ,  0.  ,  0.  ,  0.  ,  0.  ],
               [ 7.39,  7.38,  2.38,  2.89,  0.  ,  0.  ],
               [ 9.36, 11.43,  5.51,  6.42,  3.54,  3.54],
               [15.06, 15.75, 10.6 , 11.42,  8.54,  8.54],
               [15.06, 15.06, 10.7 , 11.42,  8.54,  8.54]])

        The faulted nodes have been uplifted and eroded! Note that here the
        boundary nodes are also uplifted.

        NormalFault keeps track of internal time.

        For example, if a user wanted to only run NormalFault every tenth
        timestep (or some more seismogenically reasonable set of times).

        >>> grid = RasterModelGrid((6, 6), xy_spacing=10)
        >>> z = grid.add_zeros("topographic__elevation", at="node")
        >>> nf = NormalFault(grid, **param_dict)
        >>> fr = FlowAccumulator(grid)
        >>> fs = FastscapeEroder(grid, K_sp=0.01)
        >>> model_time = 0.0
        >>> dt = 100.0
        >>> for i in range(300):
        ...     if i % 10 == 0:
        ...         nf.run_one_step(dt * 10)
        ...     fr.run_one_step()
        ...     fs.run_one_step(dt)
        ...     model_time += dt
        ...
        >>> model_time
        30000.0
        >>> nf.current_time
        30000.0
        """
        fault_throw_rate_through_time = dict(fault_throw_rate_through_time)
        fault_trace = dict(fault_trace)

        super().__init__(grid)

        # save a reference to the grid

        # get the surface to be faulted
        self._surfaces = {}
        self._not_yet_instantiated = []
        if isinstance(faulted_surface, list):
            # if faulted surface is a list, then itterate through multiple
            # surfaces and save
            for surf in faulted_surface:
                try:
                    self._surfaces[surf] = grid.at_node[surf]
                except FieldError:
                    self._not_yet_instantiated.append(surf)
        else:
            self._surfaces[faulted_surface] = grid.at_node[faulted_surface]

        if fault_dip_angle > 90.0:
            raise ValueError(
                "NormaFault fault_dip_angle must be less than 90 " "degrees."
            )

        # get the fault throw parameter values from the parameter dictionary
        self._throw_time = np.array(fault_throw_rate_through_time["time"])
        self._throw_rate = np.array(fault_throw_rate_through_time["rate"])
        self._fault_dip = np.deg2rad(fault_dip_angle)
        self._uplift = self._throw_rate * np.sin(self._fault_dip)

        # Identify in current boundaries will be included
        self._include_boundaries = include_boundaries

        # Instantiate record of current time.
        self._current_time = 0.0

        # get the fault trace dictionary and use to to calculate where the
        # faulted nodes are located.
        self._fault_trace = fault_trace
        dx = self._fault_trace["x2"] - self._fault_trace["x1"]
        dy = self._fault_trace["y2"] - self._fault_trace["y1"]
        self._fault_azimuth = np.mod(np.arctan2(dy, dx), TWO_PI)
        self._fault_anti_azimuth = self._fault_azimuth + np.pi
        # deal with the edge case in which dx == 0
        if dx == 0:
            self._dy_over_dx = 0.0
            self._fault_trace_y_intercept = 0.0
            self._fault_trace_x_intercept = self._fault_trace["x2"]
        else:
            self._dy_over_dx = dy / dx
            self._fault_trace_y_intercept = self._fault_trace["y1"] - (
                self._dy_over_dx * self._fault_trace["x1"]
            )
            self._fault_trace_x_intercept = 0.0

        # set the considered nodes based on whether the boundaries will be
        # included in the faulted terrain.
        if self._include_boundaries:
            potential_nodes = np.arange(self._grid.size("node"))
        else:
            potential_nodes = self._grid.core_nodes

        # select those nodes that are on the correct side of the fault
        dx_pn = self._grid.x_of_node[potential_nodes] - self._fault_trace_x_intercept
        dy_pn = self._grid.y_of_node[potential_nodes] - self._fault_trace_y_intercept
        potential_angles = np.mod(np.arctan2(dy_pn, dx_pn), TWO_PI)
        if self._fault_anti_azimuth <= TWO_PI:
            faulted_node_ids = potential_nodes[
                (
                    (potential_angles > self._fault_azimuth)
                    & (potential_angles <= (self._fault_anti_azimuth))
                )
            ]
        else:
            faulted_node_ids = potential_nodes[
                (
                    (potential_angles > self._fault_azimuth)
                    | (potential_angles <= np.mod(self._fault_anti_azimuth, TWO_PI))
                )
            ]

        # save a n-node array of boolean identifing faulted nodes.
        self._faulted_nodes = np.zeros(self._grid.size("node"), dtype=bool)
        self._faulted_nodes[faulted_node_ids] = True

    @property
    def faulted_nodes(self):
        """At node array indicating which nodes are on the upthrown block."""
        return self._faulted_nodes

    def _check_surfaces(self):
        if len(self._not_yet_instantiated) > 0:
            still_not_instantiated = []
            for surf in self._not_yet_instantiated:
                if surf in self._grid.at_node:
                    self._surfaces[surf] = self._grid.at_node[surf]
                else:
                    still_not_instantiated.append(surf)
            self._not_yet_instantiated = still_not_instantiated

    def run_one_earthquake(self, dz):
        """Run one earthquake with uplift of magnitude ``dz``."""
        self._check_surfaces()

        # save z before uplift only if using include boundaries.
        if self._include_boundaries:
            surfs_before_uplift = {}
            for surf_name in self._surfaces:
                surfs_before_uplift[surf_name] = self._surfaces[surf_name].copy()

        # uplift the faulted_nodes
        for surf_name in self._surfaces:
            self._surfaces[surf_name][self._faulted_nodes] += dz

        # if faulted nodes includes boundaries we must do some extra work because
        # landlab components will typically not erode these boundaries. This means
        # they will be uplifted but not eroded.

        if self._include_boundaries:
            #  here our goal is to set faulted boundaries to average of open
            # node faulted neighbors

            # create boolean of the faulted boundary nodes
            faulted_boundaries = self._faulted_nodes.copy()
            faulted_boundaries[self._grid.core_nodes] = False

            core_nodes = np.zeros(self._grid.size("node"), dtype=bool)
            core_nodes[self._grid.core_nodes] = True

            neighbor_is_core = core_nodes[self._grid.adjacent_nodes_at_node]
            neighbor_is_faulted = self._faulted_nodes[self._grid.adjacent_nodes_at_node]

            neighbor_for_averaging = neighbor_is_faulted & neighbor_is_core

            # Identify nodes that have at least one adjacent node that is both
            # faulted and not a boundary node.
            # average the pre-uplift topography on those adjacent nodes and assign
            # to the boundary node.
            # here we use the pre-uplift elevation because other steps in the model
            # may diminish this topography.

            averaged = neighbor_for_averaging[faulted_boundaries].sum(axis=1) == 1
            if any(averaged):
                averaged_nodes = np.where(faulted_boundaries)[0][np.where(averaged)[0]]
                for surf_name in self._surfaces:
                    elevations_to_average = surfs_before_uplift[surf_name][
                        self._grid.adjacent_nodes_at_node
                    ]
                    elevations_to_average[self._grid.adjacent_nodes_at_node == -1] = (
                        np.nan
                    )
                    elevations_to_average[~neighbor_for_averaging] = np.nan
                    self._surfaces[surf_name][averaged_nodes] = np.nanmean(
                        elevations_to_average[averaged_nodes], axis=1
                    )

            # identify any boundary nodes that are not being averaged. This will
            # happen at the corners on RasterModelGrids. Average over adjacent
            # nodes that are faulted. These nodes will be boundary nodes.
            # here we use the current topography as we will have just updated the
            # adjacent nodes in the prior block.
            if any(~averaged):
                un_averaged_nodes = np.where(faulted_boundaries)[0][
                    np.where(~averaged)[0]
                ]
                for surf_name in self._surfaces:
                    elevations_to_average = self._surfaces[surf_name][
                        self._grid.adjacent_nodes_at_node
                    ]
                    elevations_to_average[self._grid.adjacent_nodes_at_node == -1] = (
                        np.nan
                    )
                    elevations_to_average[~neighbor_is_faulted] = np.nan
                    self._surfaces[surf_name][un_averaged_nodes] = np.nanmean(
                        elevations_to_average[un_averaged_nodes], axis=1
                    )

    def run_one_step(self, dt):
        """Run_one_step method for NormalFault.

        Parameters
        ----------
        dt : float
            Time increment used to advance the NormalFault component.
        current_time : float, optional
            If NormalFault is not being advanced by dt every timestep with all
            components, its internal time may be incorrect, this can be remedied
            by providing a value for current time. Default value is None which
            results in the internal timekeeping not being changed.
        """
        # calculate the current uplift rate
        current_uplift_rate = np.interp(
            self._current_time, self._throw_time, self._throw_rate
        )

        # run one earthquake of size current_uplift_rate * dt
        self.run_one_earthquake(current_uplift_rate * dt)

        # increment time
        self._current_time += dt
