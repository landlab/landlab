#!/usr/bin/env python

import numpy as np

from landlab import Component
from landlab import RasterModelGrid


class KinematicWaveRengers(Component):
    """
    This code is based on an overland flow model by Francis Rengers and
    colleagues, after Julien et al., 1995. It uses an explicit face-centered
    solution to a depth-varying Manning's equation, broadly following, e.g.,
    Mugler et al., 2011.
    It was implemented in Landlab by DEJH, March '16. Please cite
    Rengers et al., 2016, Model Predictions of Water Runoff in Steep
    Catchments after Wildfire, WRR.

    Note: You will not have a good day if you have pits present in your topo
    before routing flow with this component. Fill pits before instantiating
    the component (or call :func:`update_topographic_params` once you have
    filled after instantiation).

    Note this module assumes that the topography DOES NOT change during the
    run. If it does, call :func:`update_topographic_params` to update the
    component to the new topo.

    Boundary condition control can be... interesting with this component.
    Be sure to close boundaries you do not wish water to leave - or enter! -
    through. To allow free water discharge from the grid edge it is
    recommended to use fixed gradient boundary conditions at the open edges.
    The component will then set the fixed gradient as equal to the underlying
    topographic gradient throughout the run.

    It is also possible to fix the water depth at the open edge, but this
    is not really recommended.

    Construction::

        KinematicWaveRengers(grid, mannings_n=0.03, critical_flow_depth=0.003,
                             mannings_epsilon=0.33333333, dt_max=0.3,
                             max_courant=0.2, min_surface_water_depth=1.e-8)

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.
    mannings_n : float
        A value to use for Manning's n in the Manning discharge equation.
    critical_flow_depth : float (m)
        An index flow depth for the depth-varying Manning's equation,
        controlling the depth at which the effective Manning's n begins to
        increase.
    mannings_epsilon : float
        An exponent for the depth-varying Manning's equation, controlling the
        rate of increase of effective Manning's n at small flow depths.
    dt_max : float or None (s)
        The largest permitted internal timestep for the component. If the
        Courant criterion produces a more restrictive condition, that will be
        used instead.
    max_courant : float
        The maximum permitted Courant number for the courant stability
        criterion.
    min_surface_water_depth : float (m)
        A water depth below which surface water thickness may never fall, to
        ensure model stabilty.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import KinematicWaveRengers
    >>> mg = RasterModelGrid((5, 10), 10.0)
    >>> mg.status_at_node[mg.nodes_at_left_edge] = mg.BC_NODE_IS_FIXED_GRADIENT
    >>> mg.status_at_node[mg.nodes_at_top_edge] = mg.BC_NODE_IS_CLOSED
    >>> mg.status_at_node[mg.nodes_at_bottom_edge] = mg.BC_NODE_IS_CLOSED
    >>> mg.status_at_node[mg.nodes_at_right_edge] = mg.BC_NODE_IS_CLOSED
    >>> _ = mg.add_field("node", "topographic__elevation", 0.05 * mg.node_x)
    >>> _ = mg.add_empty("node", "surface_water__depth")
    >>> mg.at_node["surface_water__depth"].fill(1.0e-8)
    >>> dt = 60.0  # 1 min intervals
    >>> rain_intensities = (1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5)
    >>> kw = KinematicWaveRengers(mg)
    >>> for i in rain_intensities:
    ...     kw.run_one_step(dt, rainfall_intensity=i)
    ...
    >>> mg.at_node["surface_water__depth"]
    array([1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
           1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
           1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
           1.00000000e-08, 2.95578314e-03, 2.95578314e-03,
           2.90945761e-03, 2.82912876e-03, 2.70127141e-03,
           2.51202011e-03, 2.24617193e-03, 1.88032853e-03,
           1.35451064e-03, 1.00000000e-08, 2.95578314e-03,
           2.95578314e-03, 2.90945761e-03, 2.82912876e-03,
           2.70127141e-03, 2.51202011e-03, 2.24617193e-03,
           1.88032853e-03, 1.35451064e-03, 1.00000000e-08,
           2.95578314e-03, 2.95578314e-03, 2.90945761e-03,
           2.82912876e-03, 2.70127141e-03, 2.51202011e-03,
           2.24617193e-03, 1.88032853e-03, 1.35451064e-03,
           1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
           1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
           1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
           1.00000000e-08, 1.00000000e-08])
    """

    _name = "KinematicWaveRengers"

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
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "surface_water__discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": "Volumetric discharge of surface water",
        },
        "surface_water__velocity": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/s",
            "mapping": "node",
            "doc": "Speed of water flow above the surface",
        },
    }

    def __init__(
        self,
        grid,
        mannings_n=0.03,
        critical_flow_depth=0.003,
        mannings_epsilon=0.33333333,
        dt_max=0.3,
        max_courant=0.2,
        min_surface_water_depth=1.0e-8,
    ):
        """Initialize the kinematic wave approximation overland flow component."""
        super().__init__(grid)

        if not isinstance(self.grid, RasterModelGrid):
            ValueError("KinematicWaveRengers: grid must be regular")

        if np.isclose(dt_max, 0.0):
            raise ValueError("KinematicWaveRengers: dt must be > 0.0")

        active = np.nonzero(self.grid.status_at_node != self.grid.BC_NODE_IS_CLOSED)
        self._h = self.grid.at_node["surface_water__depth"]
        self._active = active
        self._hc = critical_flow_depth
        self._n = mannings_n
        self._negepsilon = -mannings_epsilon

        self._dt_max = dt_max
        self._min_surface_water_depth = min_surface_water_depth
        self._active_depths = self.grid.at_node["surface_water__depth"][active]
        all_grads = self.grid.calc_grad_at_link("topographic__elevation")
        hoz_grads = self.grid.map_mean_of_horizontal_active_links_to_node(all_grads)
        vert_grads = self.grid.map_mean_of_vertical_active_links_to_node(all_grads)
        self._hozslopept5 = np.fabs(hoz_grads[active]) ** 0.5
        self._vertslopept5 = np.fabs(vert_grads[active]) ** 0.5
        self._velx = self.grid.zeros("node", dtype=float)
        self._vely = self.grid.zeros("node", dtype=float)
        self._qy = np.zeros(grid.number_of_nodes + 1, dtype=float)
        self._qx = np.zeros(grid.number_of_nodes + 1, dtype=float)
        self._poshozgrads = hoz_grads > 0.0
        self._posvertgrads = vert_grads > 0.0
        if np.isclose(self.grid.dx, self.grid.dy):
            self._equaldims = True
            self._courant_prefactor = max_courant * self.grid.dx
        else:
            self._equaldims = False
            self._courant_prefactor = max_courant * self.grid.dx * self.grid.dy
        self._neighbors = self.grid.adjacent_nodes_at_node.copy()
        self._neighbors[self._neighbors == self.grid.BAD_INDEX] = -1
        self._water_balance = []
        self._actives_BCs = (
            self.grid.status_at_node[active] == self.grid.BC_NODE_IS_FIXED_VALUE
        )
        self._actives_BCs_water_depth = self._h[active][self._actives_BCs]
        fixed_grad_nodes = self.grid.fixed_gradient_boundary_nodes.copy()
        fixed_grad_anchors = self.grid.fixed_gradient_boundary_node_anchor_node

        # ^add this value to the anchor nodes to update the BCs
        # these also need to be mapped to active_IDs:
        blank_nodes = self.grid.zeros("node", dtype=bool)
        blank_nodes[fixed_grad_nodes] = True
        self._fixed_grad_nodes_active = np.where(blank_nodes[active])[0]
        blank_nodes.fill(False)
        blank_nodes[fixed_grad_anchors] = True
        self._fixed_grad_anchors_active = np.where(blank_nodes[active])[0]

        # create outputs
        self.initialize_output_fields()

    def run_one_step(
        self,
        dt,
        rainfall_intensity=0.00001,
        update_topography=False,
        track_min_depth=False,
    ):
        """Update fields with current hydrologic conditions.

        Parameters
        ----------
        rain_intensity : float or array (m/s)
            The rainfall intensity across the grid (water input rate at each
            node).
        update_topography : bool
            Set to true if the topography of the grid evolves during the run.
        track_min_depth : bool
            At *very* low rainfall inputs, there is a possibility this
            component could allow creation of small amounts of water mass.
            Set to true to track this mass, and use the :func:`water_balance`
            property to investigate its evolution through time.
        """
        elapsed_time_in_dt = 0.0  # this is only since the start of the timestep
        active = self._active
        self._hnew = self._h[active]
        _hnew = self._hnew
        if update_topography:
            self.update_topographic_params()
        while elapsed_time_in_dt < dt:
            internal_dt = self.calc_grads_and_timesteps(
                update_topography, track_min_depth
            )
            remaining_dt = dt - elapsed_time_in_dt
            # now reduce timestep is needed if limited by total tstep length
            internal_dt = min(internal_dt, remaining_dt).clip(0.0)
            # this section uses our final-array-val-is-zero trick
            _qx_left = self._qx[self._neighbors[:, 2]].clip(min=0.0)
            _qx_right = self._qx[self._neighbors[:, 0]].clip(max=0.0)
            _qy_top = self._qy[self._neighbors[:, 1]].clip(min=0.0)
            _qy_bottom = self._qy[self._neighbors[:, 3]].clip(max=0.0)
            # FR's rainfall handling was here. We're going to assume that the
            # component is being driven by a "LL style" rainfall record, where
            # the provided rainfall_intensity is constant across the provide
            # dt. If it's not, it needs to be handled outside the component.

            # now add the rainfall input
            if type(rainfall_intensity) is not np.ndarray:
                _hnew += internal_dt * rainfall_intensity
            else:
                _hnew += internal_dt * rainfall_intensity[active]
            # set the BCs
            _hnew[self._actives_BCs] = self._actives_BCs_water_depth
            # flux it round
            _hnew -= internal_dt / self.grid.dx * np.fabs(self._qy[active])
            _hnew -= internal_dt / self.grid.dy * np.fabs(self._qx[active])
            _hnew += internal_dt / self.grid.dx * (_qy_top - _qy_bottom)[active]
            _hnew += internal_dt / self.grid.dy * (_qx_left - _qx_right)[active]
            _hnew[self._fixed_grad_nodes_active] = _hnew[
                self._fixed_grad_anchors_active
            ]
            # update the internal clock
            elapsed_time_in_dt += internal_dt

        # update the actual fields
        self._h[active] = _hnew
        self.grid.at_node["surface_water__velocity"][:] = np.sqrt(
            np.square(self._velx) + np.square(self._vely)
        )
        self.grid.at_node["surface_water__discharge"][:] = np.sqrt(
            np.square(self._qx[:-1]) + np.square(self._qy[:-1])
        )

    def calc_grads_and_timesteps(self, update_topography, track_min_depth):
        """
        Perform the first part of the calculation for the main run, mainly
        velocities and fluxes. The main objective of this part of the
        calculation is to derive the stable timestep for the run.

        Parameters
        ----------
        update_topography : bool
            If False, the underlying surface topography is assumed unchanged
            since the last run.
        track_min_depth : bool
            If True, the internal list _water_balance will be appended with
            the volumetric fractional change in mass balance during the run.

        Returns
        -------
        internal_dt : float
            The internally stable timestep that will be used on this loop.
        """
        active = self._active
        _hnew = self._hnew
        if update_topography:
            self.update_topographic_params()
        # assert the minimum water depth - this could introduce an element of
        # mass gain, but should remain minor
        _hnew.clip(self._min_surface_water_depth, out=_hnew)
        if track_min_depth:
            self._water_balance.append(
                (_hnew - self._h[active]).sum() / self._h[active].sum()
            )
        n = self._n * (_hnew / self._hc) ** self._negepsilon
        twothirds_hnewbyn = _hnew**0.66666666 / n
        self._vely[active] = twothirds_hnewbyn * self._vertslopept5
        self._velx[active] = twothirds_hnewbyn * self._hozslopept5
        self._vely[self._posvertgrads] *= -1.0
        self._velx[self._poshozgrads] *= -1.0
        self._qy[active] = self._vely[active] * _hnew  # m**2/s
        self._qx[active] = self._velx[active] * _hnew  # m**2/s
        max_vely = np.fabs(self._vely).max()
        max_velx = np.fabs(self._velx).max()
        if self._equaldims:
            courant_dt = self._courant_prefactor / (max_velx + max_vely)
        else:
            # note prefactor is NOT THE SAME as above in this case
            courant_dt = self._courant_prefactor / (
                self.grid.dy * max_velx + self.grid.dx * max_vely
            )
        if self._dt_max is not None:
            internal_dt = np.min((self._dt_max, courant_dt))
        else:
            internal_dt = courant_dt
        self._internal_dt = internal_dt

        return internal_dt

    def update_topographic_params(self):
        """
        If the topo changes during the run, change the held params used by
        :func:`run_one_step`.
        """
        active = np.where(self.grid.status_at_node != self.grid.BC_NODE_IS_CLOSED)[0]
        all_grads = self.grid.calculate_gradients_at_links("topographic__elevation")
        hoz_grads = self.grid.map_mean_of_horizontal_active_links_to_node(all_grads)
        vert_grads = self.grid.map_mean_of_vertical_active_links_to_node(all_grads)
        self._hozslopept5 = np.fabs(hoz_grads[active]) ** 0.5
        self._vertslopept5 = np.fabs(vert_grads[active]) ** 0.5
        self._poshozgrads = hoz_grads > 0.0
        self._posvertgrads = vert_grads > 0.0
        fixed_grad_nodes = self.grid.fixed_gradient_boundary_nodes
        fixed_grad_anchors = self.grid.fixed_gradient_boundary_node_anchor_node
        # ^add this value to the anchor nodes to update the BCs
        # these also need to be mapped to active_IDs:
        blank_nodes = self.grid.zeros("node", dtype=bool)
        blank_nodes[fixed_grad_nodes] = True
        self._fixed_grad_nodes_active = np.where(blank_nodes[active])[0]
        blank_nodes.fill(False)
        blank_nodes[fixed_grad_anchors] = True
        self._fixed_grad_anchors_active = np.where(blank_nodes[active])[0]
        # check is the grid topology has changed...
        if not np.all(np.equal(self._active, active)):
            self._active = active
            self._velx.fill(0.0)
            self._vely.fill(0.0)
            self._qy.fill(0.0)
            self._qx.fill(0.0)
            self._neighbors = self.grid.adjacent_nodes_at_node.copy()
            self._neighbors[self._neighbors == self.grid.BAD_INDEX] = -1
            self._actives_BCs = (
                self.grid.status_at_node[active] == self.grid.BC_NODE_IS_FIXED_VALUE
            )
            self._actives_BCs_water_depth = self._h[self._actives_BCs]

    @property
    def water_balance(self):
        """
        Return a list of the fractional gain/loss of water mass during the
        run, if it was tracked using the track_min_depth flag.
        """
        if self._water_balance == []:
            raise ValueError("No record of water balance was found!")
        else:
            return self._water_balance

    @property
    def internal_timestep(self):
        """
        Return the internal timestep last used by the kinematic wave component.
        """
        try:
            return self._internal_dt
        except AttributeError:
            # the component hasn't started running yet
            _ = self.calc_grads_and_timesteps(False, False)
            return self._internal_dt
