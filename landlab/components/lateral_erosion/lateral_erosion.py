#! /usr/env/python
# -*- coding: utf-8 -*-
"""
April 4, 2019 Starting to rehab lateral erosion
ALangston


"""

import numpy as np

from landlab import Component, FieldError, RasterModelGrid
from landlab.components.flow_accum import FlowAccumulator
from landlab.components.flow_director import FlowDirectorD8
from landlab.components.flow_routing import DepressionFinderAndRouter
from landlab.components.lateral_erosion.node_finder import Node_Finder


class LateralEroder(Component):
    """
    Laterally erode neighbor node through fluvial erosion.

    Landlab component that finds a neighbor node to laterally erode and calculates lateral erosion.

    Construction:
        LateralEroder(grid, latero_mech="UC", alph=0.8, Kv=None, Kl_ratio=1.0, solver="basic", inlet_on=False, inlet_node=None, inlet_area=None, qsinlet=None)

    Parameteters
    ------------
    grid : ModelGrid
        A Landlab grid object
    latero_mech : string, optional (defaults to UC)
        Lateral erosion algorithm, choices are "UC" for undercutting-slump model and "TB" for total block erosion
    alph : float, optional (defaults to 0.8)
        Parameter describing potential for deposition, dimensionless
    Kv : float, node array, or field name
        Bedrock erodibility in vertical direction, 1/years
    Kl_ratio : float, optional (defaults to 1.0)
        Ratio of lateral to vertical bedrock erodibility, dimensionless
    solver : string
        Solver options:
            (1) 'basic' (default): explicit forward-time extrapolation.
                Simple but will become unstable if time step is too large or
                if bedrock erodibility is vry high.
            (2) 'adaptive': subdivides global time step as needed to
                prevent slopes from reversing.
    inlet_node : integer, optional
        Node location of inlet (source of water and sediment)
    inlet_area : float, optional
        Drainage area at inlet node, must be specified if inlet node is "on", m^2
    qsinlet : float, optional
        Sediment flux supplied at inlet, optional. m3/year

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY
    >>> from landlab.components import FlowAccumulator, LateralEroder
    >>> np.random.seed(2010)

    Define grid and initial topography

    * 5x4 grid with baselevel in the lower left corner
    * All other boundary nodes closed
    * Initial topography is plane tilted up to the upper right with noise

    >>> nr = 5
    >>> nc = 4
    >>> dx=10
    >>> mg = RasterModelGrid((nr, nc), xy_spacing = dx)
    >>> mg.set_status_at_node_on_edges(right=CLOSED_BOUNDARY, top=CLOSED_BOUNDARY, \
                                  left=CLOSED_BOUNDARY, bottom=CLOSED_BOUNDARY)
    >>> mg.status_at_node[1] = FIXED_VALUE_BOUNDARY
    >>> mg.add_zeros('node', 'topographic__elevation')
    array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.])
    >>> rand_noise=np.array([ 0.00436992,  0.03225985,  0.03107455,  0.00461312,  0.03771756,
    ...    0.02491226,  0.09613959,  0.07792969,  0.08707156,  0.03080568,
    ...    0.01242658,  0.08827382,  0.04475065,  0.07391732,  0.08221057,
    ...    0.02909259,  0.03499337,  0.09423741,  0.01883171,  0.09967794])
    >>> mg.at_node['topographic__elevation'] += (mg.node_y / 10. + mg.node_x / 10. + rand_noise)
    >>> U=0.001
    >>> dt=100

    Instantiate flow accumulation and lateral eroder and run each for one step

    >>> fa = FlowAccumulator(mg,
    ...                     surface='topographic__elevation',
    ...                     flow_director='FlowDirectorD8',
    ...                     runoff_rate=None,
    ...                     depression_finder=None)
    >>> latero = LateralEroder(mg,latero_mech="UC", Kv=0.001, Kl_ratio=1.5)

    Run one step of flow accumulation and lateral erosion to get the dzlat array
    needed for the next part of the test.

    >>> fa.run_one_step()
    >>> (mg, dzlat,)=latero.run_one_step(mg,dt,)

    Evolve the landscape until the first occurence of lateral erosion. Save arrays
    volume of lateral erosion and topographic elevation before and after the first
    occurence of lateral erosion

    >>> while min(dzlat)==0.:
    ...         oldlatvol=mg['node'][ 'volume__lateral_erosion'].copy()
    ...         oldelev=mg['node']['topographic__elevation'].copy()
    ...         fa.run_one_step()
    ...         (mg, dzlat,)=latero.run_one_step(mg,dt,)
    ...         newlatvol=mg['node'][ 'volume__lateral_erosion']
    ...         newelev=mg['node']['topographic__elevation']
    ...         mg.at_node['topographic__elevation'][mg.core_nodes] += U*dt

    Before lateral erosion occurs, volume__lateral_erosion has values at nodes 6 and 10.

    >>> np.around(oldlatvol, decimals=0) # doctest: +NORMALIZE_WHITESPACE
    array([  0.,   0.,   0.,   0.,
             0.,   0.,  79.,   0.,
             0.,   0.,  24.,   0.,
             0.,   0.,   0.,   0.,
             0.,   0.,   0.,   0.])


    After lateral erosion occurs at node 6, volume__lateral_erosion is reset to 0

    >>> np.around(newlatvol, decimals=0) # doctest: +NORMALIZE_WHITESPACE
    array([  0.,   0.,   0.,   0.,
             0.,   0.,   0.,   0.,
             0.,   0.,  24.,   0.,
             0.,   0.,   0.,   0.,
             0.,   0.,   0.,   0.])


    After lateral erosion at node 6, elevation at node 6 is reduced by -1.41
    (the height stored in dzlat[6]).

    >>> np.around(oldelev, decimals=2) # doctest: +NORMALIZE_WHITESPACE
    array([ 0.  ,  1.03,  2.03,  3.  ,
           1.04,  1.77,  2.45,  4.08,
           2.09,  2.65,  3.18,  5.09,
           3.04,  3.65,  4.07,  6.03,
           4.03,  5.09,  6.02,  7.1 ])

    >>> np.around(newelev, decimals=2) # doctest: +NORMALIZE_WHITESPACE
    array([ 0.  ,  1.03,  2.03,  3.  ,
           1.04,  1.77,  1.03,  4.08,
           2.09,  2.65,  3.18,  5.09,
           3.04,  3.65,  4.07,  6.03,
           4.03,  5.09,  6.02,  7.1 ])

    >>> np.around(dzlat, decimals=2) # doctest: +NORMALIZE_WHITESPACE
    array([ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  , -1.41,  0.  ,  0.  ,
            0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,
            0.  ,  0.  ])
    """

    _name = "LateralEroder"

    _input_var_names = (
        "topographic__elevation",
        "drainage_area",
        "flow__receiver_node",
        "flow__upstream_node_order",
        "topographic__steepest_slope",
    )

    _output_var_names = ("topographic__elevation", "dzdt", "dzlat", "vollat", "qs_in")
    _var_units = {
        "topographic__elevation": "m",
        "drainage_area": "m2",
        "flow__receiver_node": "-",
        "flow__upstream_node_order": "-",
        "topographic__steepest_slope": "-",
        "dzdt": "m/y",
        "dzlat": "m/y",
        "vollat": "m3",
        "qs_in": "m3/y",
    }

    def __init__(
        self,
        grid,
        latero_mech="UC",
        alph=0.8,
        Kv=None,
        Kl_ratio=1.0,
        solver="basic",
        inlet_on=False,
        inlet_node=None,
        inlet_area=None,
        qsinlet=0.0,
    ):  # input_stream,
        self._grid = grid
        # Create fields needed for this component if not already existing
        if "volume__lateral_erosion" in grid.at_node:
            self.vol_lat = grid.at_node["volume__lateral_erosion"]
        else:
            self.vol_lat = grid.add_zeros("node", "volume__lateral_erosion")
        if "qs_in" in grid.at_node:
            self.qs_in = grid.at_node["qs_in"]
        else:
            self.qs_in = grid.add_zeros("node", "qs_in")

        # you can specify the type of lateral erosion model you want to use. But if you don't
        # the default is the undercutting-slump model
        assert latero_mech in ("UC", "TB")
        if latero_mech == "TB":
            self._TB = True
            self._UC = False
        else:
            self._UC = True
            self._TB = False
        # option use adaptive time stepping. Default is fixed dt supplied by user
        if solver == "basic":
            self.run_one_step = self.run_one_step_basic
        elif solver == "adaptive":
            self.run_one_step = self.run_one_step_adaptive
            self.frac = 0.3  # for time step calculations
        self.alph = alph
        self.Kv = Kv  # can be overwritten with spatially variable
        self.inlet_on = False  # will be overwritten below if inlet area is provided
        self.Klr = float(Kl_ratio)  # default ratio of Kv/Kl is 1. Can be overwritten

        self.dzdt = grid.add_zeros("node", "dzdt")  # elevation change rate (M/Y)
        # optional inputs
        self.inlet_on = inlet_on
        if inlet_on is True:
            if inlet_node is None:
                raise ValueError("inlet_on is true, but no inlet_node is provided.")
            else:
                self.inlet_node = inlet_node
                if inlet_area is None:
                    raise ValueError(
                        "Inlet area must be provided if inlet node is active. "
                        + "No inlet area was found."
                    )
                self.inlet_area = inlet_area
                # runoff is an array with values of the area of each node (dx**2)
                runoffinlet = np.ones(grid.number_of_nodes) * grid.dx ** 2
                # Change the runoff at the inlet node to node area + inlet node
                runoffinlet[inlet_node] += inlet_area
                _ = grid.add_field(
                    "node", "water__unit_flux_in", runoffinlet, noclobber=False
                )
                # set qsinlet at inlet node. This doesn't have to be provided, defaults
                # to 0.
                self.qsinlet = qsinlet
                self.qs_in[self.inlet_node] = self.qsinlet
        # below, adding flag calling for Kv to be specified.
        if self.Kv is None:
            raise ValueError(
                "Kv must be set as a float, node array, or "
                + "field name. It was None."
            )
        # handling Kv for floats (inwhich case it populates an array N_nodes long) or
        # for arrays of Kv. Checks that length of Kv array is good.
        if isinstance(Kv, (float, int)):
            self.Kv = np.ones(self.grid.number_of_nodes) * float(Kv)
        else:
            self.Kv = np.asarray(Kv, dtype=float)
            if len(self.Kv) != self.grid.number_of_nodes:
                raise TypeError("Supplied value of Kv is not n_nodes long")

    def run_one_step_basic(
        self, grid, dt=None, Klr=None, inlet_area_ts=None, qsinlet_ts=None, **kwds
    ):
        if Klr is None:  # Added10/9 to allow changing rainrate (indirectly this way.)
            Klr = self.Klr

        UC = self._UC
        TB = self._TB
        inlet_on = self.inlet_on  # this is a true/false flag
        Kv = self.Kv
        qs_in = self.qs_in
        dzdt = self.dzdt
        alph = self.alph
        self.dt = dt
        vol_lat = self.grid.at_node["volume__lateral_erosion"]
        kw = 10.0
        F = 0.02
        # May 2, runoff calculated below (in m/s) is important for calculating
        # discharge and water depth correctly. renamed runoffms to prevent
        # confusion with other uses of runoff
        runoffms = (Klr * F / kw) ** 2
        # Kl is calculated from ratio of lateral to vertical K parameters
        Kl = Kv * Klr
        z = grid.at_node["topographic__elevation"]
        # clear qsin for next loop
        qs_in = grid.add_zeros("node", "qs_in", noclobber=False)
        qs = grid.add_zeros("node", "qs", noclobber=False)
        lat_nodes = np.zeros(grid.number_of_nodes, dtype=int)
        dzlat = np.zeros(grid.number_of_nodes)
        dzver = np.zeros(grid.number_of_nodes)
        vol_lat_dt = np.zeros(grid.number_of_nodes)
        if inlet_on is True:
            inlet_node = self.inlet_node
            # if a value is passed with qsinlet_ts, qsinlet has changed with this timestep,
            # so reset qsinlet to qsinlet_ts
            if qsinlet_ts is not None:
                qsinlet = qsinlet_ts
                qs_in[inlet_node] = qsinlet
            # if nothing is passed with qsinlet_ts, qsinlet remains the same from
            # initialized parameters
            else:  # qsinlet_ts==None:
                qsinlet = self.qsinlet
                qs_in[inlet_node] = qsinlet
            if inlet_area_ts is not None:
                inlet_area = inlet_area_ts
                runoffinlet = np.ones(grid.number_of_nodes) * grid.dx ** 2
                # Change the runoff at the inlet node to node area + inlet node
                runoffinlet[inlet_node] += inlet_area
                _ = grid.add_field(
                    "node", "water__unit_flux_in", runoffinlet, noclobber=False
                )
                # if inlet area has changed with time (so we have a new inlet area here)
                fa = FlowAccumulator(
                    grid,
                    surface="topographic__elevation",
                    flow_director="FlowDirectorD8",
                    runoff_rate=None,
                    depression_finder="DepressionFinderAndRouter",
                    router="D8",
                )
                (da, q) = fa.accumulate_flow()
                q = grid.at_node["surface_water__discharge"]
                # this is the drainage area that I need for code below with an inlet set
                # by spatially varible runoff.
                da = q / grid.dx ** 2
            else:  # inletarea not changing with time
                q = grid.at_node["surface_water__discharge"]
                # this is the drainage area that I need for code below with an inlet set
                # by spatially varible runoff.
                da = q / grid.dx ** 2
        #                print ("da", da[inlet_node])
        # if inlet flag is not on, proceed as normal.
        else:
            # renamed this drainage area set by flow router
            da = grid.at_node["drainage_area"]
        # flow__upstream_node_order is node array contianing downstream to
        # upstream order list of node ids
        s = grid.at_node["flow__upstream_node_order"]
        max_slopes = grid.at_node["topographic__steepest_slope"]
        flowdirs = grid.at_node["flow__receiver_node"]

        # make a list l, where node status is interior (signified by label 0) in s
        interior_s = s[np.where((grid.status_at_node[s] == 0))[0]]
        dwnst_nodes = interior_s.copy()
        # reverse list so we go from upstream to down stream
        dwnst_nodes = dwnst_nodes[::-1]
        max_slopes[:] = max_slopes.clip(0)
        for i in dwnst_nodes:
            # calc deposition and erosion
            dep = alph * qs_in[i] / da[i]
            ero = -Kv[i] * da[i] ** (0.5) * max_slopes[i]
            dzver[i] = dep + ero
            # potential lateral erosion initially set to 0
            petlat = 0.0
            # water depth in meters, needed for lateral erosion calc
            wd = 0.4 * (da[i] * runoffms) ** 0.35

            # Choose lateral node for node i. If node i flows downstream, continue.
            # if node i is the first cell at the top of the drainage network, don't go
            # into this loop because in this case, node i won't have a "donor" node
            if i in flowdirs:
                # Node_finder picks the lateral node to erode based on angle
                # between segments between three nodes
                [lat_node, inv_rad_curv] = Node_Finder(grid, i, flowdirs, da)
                # node_finder returns the lateral node ID and the radius of curvature
                lat_nodes[i] = lat_node
                # if the lateral node is not 0 or -1 continue. lateral node may be
                # 0 or -1 if a boundary node was chosen as a lateral node. then
                # radius of curavature is also 0 so there is no lateral erosion
                if lat_node > 0:
                    # if the elevation of the lateral node is higher than primary node,
                    # calculate a new potential lateral erosion (L/T), which is negative
                    if z[lat_node] > z[i]:
                        petlat = -Kl[i] * da[i] * max_slopes[i] * inv_rad_curv
                        # the calculated potential lateral erosion is mutiplied by the length of the node
                        # and the bank height, then added to an array, vol_lat_dt, for volume eroded
                        # laterally  *per timestep* at each node. This vol_lat_dt is reset to zero for
                        # each timestep loop. vol_lat_dt is added to itself in case more than one primary
                        # nodes are laterally eroding this lat_node
                        # volume of lateral erosion per timestep
                        vol_lat_dt[lat_node] += abs(petlat) * grid.dx * wd

            # send sediment downstream. sediment eroded from vertical incision
            # and lateral erosion is sent downstream
            qs_in[flowdirs[i]] += (
                qs_in[i] - (dzver[i] * grid.dx ** 2) - (petlat * grid.dx * wd)
            )  # qsin to next node
        qs[:] = qs_in - (dzver * grid.dx ** 2)
        dzdt[:] = dzver * dt
        vol_lat[:] += vol_lat_dt * dt
        # this loop determines if enough lateral erosion has happened to change
        # the height of the neighbor node.
        for i in dwnst_nodes:
            lat_node = lat_nodes[i]
            wd = 0.4 * (da[i] * runoffms) ** 0.35
            if lat_node > 0:  # greater than zero now bc inactive neighbors are value -1
                if z[lat_node] > z[i]:
                    # vol_diff is the volume that must be eroded from lat_node so that its
                    # elevation is the same as node downstream of primary node
                    # UC model: this would represent undercutting (the water height at
                    # node i), slumping, and instant removal.
                    if UC == 1:
                        voldiff = (z[i] + wd - z[flowdirs[i]]) * grid.dx ** 2
                    # TB model: entire lat node must be eroded before lateral erosion
                    # occurs
                    if TB == 1:
                        voldiff = (z[lat_node] - z[flowdirs[i]]) * grid.dx ** 2
                    # if the total volume eroded from lat_node is greater than the volume
                    # needed to be removed to make node equal elevation,
                    # then instantaneously remove this height from lat node. already has
                    # timestep in it
                    if vol_lat[lat_node] >= voldiff:
                        dzlat[lat_node] = z[flowdirs[i]] - z[lat_node]  # -0.001
                        # after the lateral node is eroded, reset its volume eroded to
                        # zero
                        vol_lat[lat_node] = 0.0
        # combine vertical and lateral erosion
        dz = dzdt + dzlat
        # change height of landscape
        z[:] = dz + z
        grid["node"]["topographic__elevation"][grid.core_nodes] = z[grid.core_nodes]
        return grid, dzlat

    # %%
    def run_one_step_adaptive(
        self, grid, dt=None, Klr=None, inlet_area_ts=None, qsinlet_ts=None, **kwds
    ):

        if Klr is None:  # Added10/9 to allow changing rainrate (indirectly this way.)
            Klr = self.Klr
        UC = self._UC
        TB = self._TB
        inlet_on = self.inlet_on  # this is a true/false flag
        Kv = self.Kv
        frac = self.frac
        qs_in = self.qs_in
        dzdt = self.dzdt
        alph = self.alph
        self.dt = dt
        vol_lat = self.grid.at_node["volume__lateral_erosion"]
        kw = 10.0
        F = 0.02
        runoffms = (Klr * F / kw) ** 2
        Kl = Kv * Klr
        z = grid.at_node["topographic__elevation"]
        # clear qsin for next loop
        qs_in = grid.add_zeros("node", "qs_in", noclobber=False)
        qs = grid.add_zeros("node", "qs", noclobber=False)
        lat_nodes = np.zeros(grid.number_of_nodes, dtype=int)
        dzlat = np.zeros(grid.number_of_nodes)
        dzver = np.zeros(grid.number_of_nodes)
        vol_lat_dt = np.zeros(grid.number_of_nodes)

        if inlet_on is True:
            # define inlet_node
            inlet_node = self.inlet_node
            # if a value is passed with qsinlet_ts, qsinlet has changed with this timestep,
            # so reset qsinlet to qsinlet_ts
            if qsinlet_ts is not None:
                qsinlet = qsinlet_ts
                qs_in[inlet_node] = qsinlet
            # if nothing is passed with qsinlet_ts, qsinlet remains the same from
            # initialized parameters
            else:  # qsinlet_ts==None:
                qsinlet = self.qsinlet
                qs_in[inlet_node] = qsinlet

            if inlet_area_ts is not None:
                inlet_area = inlet_area_ts
                runoffinlet = np.ones(grid.number_of_nodes) * grid.dx ** 2
                # Change the runoff at the inlet node to node area + inlet node
                runoffinlet[inlet_node] += inlet_area
                _ = grid.add_field(
                    "node", "water__unit_flux_in", runoffinlet, noclobber=False
                )
                # if inlet area has changed with time (so we have a new inlet area here)
                fa = FlowAccumulator(
                    grid,
                    surface="topographic__elevation",
                    flow_director="FlowDirectorD8",
                    runoff_rate=None,
                    depression_finder="DepressionFinderAndRouter",
                    router="D8",
                )
                (da, q) = fa.accumulate_flow()
                q = grid.at_node["surface_water__discharge"]
                # this is the drainage area that I need for code below with an inlet set
                # by spatially varible runoff.
                da = q / grid.dx ** 2
            else:
                q = grid.at_node["surface_water__discharge"]
                # this is the drainage area that I need for code below with an inlet set
                # by spatially varible runoff.
                da = q / grid.dx ** 2
        # if inlet flag is not on, proceed as normal.
        else:
            # renamed this drainage area set by flow router
            da = grid.at_node["drainage_area"]
        s = grid.at_node["flow__upstream_node_order"]
        max_slopes = grid.at_node["topographic__steepest_slope"]
        flowdirs = grid.at_node["flow__receiver_node"]
        interior_s = s[np.where((grid.status_at_node[s] == 0))[0]]
        dwnst_nodes = interior_s.copy()
        # reverse list so we go from upstream to down stream
        dwnst_nodes = dwnst_nodes[::-1]
        # local time
        time = 0
        globdt = dt

        while time < globdt:
            max_slopes[:] = max_slopes.clip(0)
            # here calculate dzdt for each node, with initial time step
            for i in dwnst_nodes:
                dep = alph * qs_in[i] / da[i]
                ero = -Kv[i] * da[i] ** (0.5) * max_slopes[i]
                dzver[i] = dep + ero
                petlat = 0.0
                # water depth in meters, needed for lateral erosion calc
                wd = 0.4 * (da[i] * runoffms) ** 0.35

                if i in flowdirs:
                    # Node_finder picks the lateral node to erode based on angle
                    # between segments between three nodes
                    [lat_node, inv_rad_curv] = Node_Finder(grid, i, flowdirs, da)
                    # node_finder returns the lateral node ID and the radius of
                    # curvature
                    lat_nodes[i] = lat_node
                    # if the lateral node is not 0 or -1 continue.
                    if lat_node > 0:
                        # if the elevation of the lateral node is higher than primary node,
                        # calculate a new potential lateral erosion (L/T), which is
                        # negative
                        if z[lat_node] > z[i]:
                            petlat = -Kl[i] * da[i] * max_slopes[i] * inv_rad_curv
                            # the calculated potential lateral erosion is mutiplied by the length of the node
                            # and the bank height, then added to an array, vol_lat_dt, for volume eroded
                            # laterally  *per timestep* at each node. This vol_lat_dt is reset to zero for
                            # each timestep loop. vol_lat_dt is added to itself more than one primary nodes are
                            # laterally eroding this lat_node
                            # volume of lateral erosion per timestep
                            vol_lat_dt[lat_node] += abs(petlat) * grid.dx * wd
                # send sediment downstream. sediment eroded from vertical incision
                # and lateral erosion is sent downstream
                qs_in[flowdirs[i]] += (
                    qs_in[i] - (dzver[i] * grid.dx ** 2) - (petlat * grid.dx * wd)
                )  # qsin to next node
            # summing qs for this entire timestep
            qs[:] += qs_in - (dzver * grid.dx ** 2)
            dzdt[:] = dzver
            # Do a time-step check
            # If the downstream node is eroding at a slower rate than the
            # upstream node, there is a possibility of flow direction reversal,
            # or at least a flattening of the landscape.
            # Limit dt so that this flattening or reversal doesn't happen.
            # How close you allow these two points to get to eachother is
            # determined by the variable frac.
            # dtn is an arbitrarily large number to begin with, but will be adapted as we step through
            # the nodes
            dtn = dt * 50  # starting minimum timestep for this round
            for i in dwnst_nodes:
                # are points converging? ie, downstream eroding slower than upstream
                dzdtdif = dzdt[flowdirs[i]] - dzdt[i]
                # if points converging, find time to zero slope
                if dzdtdif > 1.0e-5 and max_slopes[i] > 1e-5:
                    # time to flat between points
                    dtflat = (z[i] - z[flowdirs[i]]) / dzdtdif
                    # if time to flat is smaller than dt, take the lower value
                    if dtflat < dtn:
                        dtn = dtflat
                    #                        assert dtn>0, "dtn <0 at dtflat"
                    # if dzdtdif*dtflat will make upstream lower than downstream, find
                    # time to flat
                    if dzdtdif * dtflat > (z[i] - z[flowdirs[i]]):
                        dtn = (z[i] - z[flowdirs[i]]) / dzdtdif
            dtn *= frac
            # new minimum timestep for this round of nodes
            dt = min(abs(dtn), dt)
            assert dt > 0.0, "timesteps less than 0."

            # vol_lat is the total volume eroded from the lateral nodes through
            # the entire model run. So vol_lat is itself plus vol_lat_dt (for current loop)
            # times stable timestep size
            vol_lat[:] += vol_lat_dt * dt
            # this loop determines if enough lateral erosion has happened to change
            # the height of the neighbor node.
            for i in dwnst_nodes:
                lat_node = lat_nodes[i]
                wd = 0.4 * (da[i] * runoffms) ** 0.35
                if (
                    lat_node > 0
                ):  # greater than zero now bc inactive neighbors are value -1
                    if z[lat_node] > z[i]:
                        # vol_diff is the volume that must be eroded from lat_node so that its
                        # elevation is the same as node downstream of primary node
                        # UC model: this would represent undercutting (the water height
                        # at node i), slumping, and instant removal.
                        if UC == 1:
                            voldiff = (z[i] + wd - z[flowdirs[i]]) * grid.dx ** 2
                        # TB model: entire lat node must be eroded before lateral
                        # erosion occurs
                        if TB == 1:
                            voldiff = (z[lat_node] - z[flowdirs[i]]) * grid.dx ** 2
                        # if the total volume eroded from lat_node is greater than the volume
                        # needed to be removed to make node equal elevation,
                        # then instantaneously remove this height from lat node. already
                        # has timestep in it
                        if vol_lat[lat_node] >= voldiff:
                            dzlat[lat_node] = z[flowdirs[i]] - z[lat_node]  # -0.001
                            # after the lateral node is eroded, reset its volume eroded
                            # to zero
                            vol_lat[lat_node] = 0.0

            # multiply dzver(changed to dzdt above) by timestep size and combine with lateral erosion
            # dzlat, which is already a length for the chosen time step
            dz = dzdt * dt + dzlat
            # change height of landscape
            z[:] = dz + z
            grid["node"]["topographic__elevation"][grid.core_nodes] = z[grid.core_nodes]
            # update elapsed time
            time = dt + time
            # check to see that you are within 0.01% of the storm duration, if so
            # done, if not continue

            if time > 0.9999 * globdt:
                time = globdt

            else:
                dt = globdt - time
                qs_in = grid.zeros(centering="node")
                # recalculate flow directions
                fa = FlowAccumulator(
                    grid,
                    surface="topographic__elevation",
                    flow_director="FlowDirectorD8",
                    runoff_rate=None,
                    depression_finder="DepressionFinderAndRouter",
                    router="D8",
                )
                (da, q) = fa.accumulate_flow()
                if inlet_on:
                    #                   #if inlet on, reset drainage area and qsin to reflect inlet conditions
                    # this is the drainage area that I need for code below with an inlet
                    # set by spatially varible runoff.
                    da = q / grid.dx ** 2
                    qs_in[inlet_node] = qsinlet
                else:
                    # otherwise, drainage area is just drainage area. *** could remove the
                    # below line to speed up model. It's not really necessary.
                    # renamed this drainage area set by flow router
                    da = grid.at_node["drainage_area"]
                s = grid.at_node["flow__upstream_node_order"]
                max_slopes = grid.at_node["topographic__steepest_slope"]
                q = grid.at_node["surface_water__discharge"]
                flowdirs = grid.at_node["flow__receiver_node"]
                interior_s = s[np.where((grid.status_at_node[s] == 0))[0]]
                dwnst_nodes = interior_s.copy()
                dwnst_nodes = dwnst_nodes[::-1]

                lat_nodes = np.zeros(grid.number_of_nodes, dtype=int)
                dzlat = np.zeros(grid.number_of_nodes)
                vol_lat_dt = np.zeros(grid.number_of_nodes)
                dzver = np.zeros(grid.number_of_nodes)

        return grid, dzlat
