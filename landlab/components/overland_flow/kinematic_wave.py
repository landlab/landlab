#!/usr/bin/env python

import numpy as np

from landlab import Component, RasterModelGrid, CLOSED_BOUNDARY
from landlab import BAD_INDEX_VALUE, FIXED_VALUE_BOUNDARY
from landlab import FIXED_GRADIENT_BOUNDARY
from ...utils.decorators import use_file_name_or_kwds


class KinematicWave(Component):

    """
########### THINK MORE ABOUT BC SETTING - how does water leave?
    This code is based on an overland flow model by Francis Rengers and
    colleagues, after Julien et al., 1995. It uses an explicit face-centered
    solution to a depth-varying Manning's equation, broadly following, e.g.,
    Mugler et al., 2011.
    It was implemented in Landlab by DEJH, March '16. Please cite
    Rengers et al., in review, Model Predictions of Water Runoff in Steep
    Catchments after Wildfire, WRR.

    Note this module assumes that the topography DOES NOT change during the
    run. If it does, call :func:`update_topographic_params` to update the
    component to the new topo.

    Construction::

        KinematicWave(grid, mannings_n=0.03, critical_flow_depth=0.003,
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
    >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
    >>> mg = RasterModelGrid((5, 10), spacing=10.)
    >>> mg.status_at_node[mg.nodes_at_top_edge] = CLOSED_BOUNDARY
    >>> mg.status_at_node[mg.nodes_at_bottom_edge] = CLOSED_BOUNDARY
    >>> mg.status_at_node[mg.nodes_at_right_edge] = CLOSED_BOUNDARY
    >>> _ = mg.add_field('node', 'topographic__elevation', 0.05*mg.node_x)
    >>> _ = mg.add_empty('node', 'surface_water__depth')
    >>> mg.at_node['surface_water__depth'].fill(1.e-8)
    >>> dt = 60.  # 1 min intervals
    >>> rain_intensities = (1.e-5, 1.e-5, 1.e-5, 1.e-5, 1.e-5)
    >>> kw = KinematicWave(mg)
    >>> for i in rain_intensities:
    ...     kw.update_one_timestep(dt, rainfall_intensity=i)
    >>> mg.at_node['surface_water__depth']
    array([  1.00000000e-08,   1.00000000e-08,   1.00000000e-08,
             1.00000000e-08,   1.00000000e-08,   1.00000000e-08,
             1.00000000e-08,   1.00000000e-08,   1.00000000e-08,
             1.00000000e-08,   1.00000000e-08,   2.95578314e-03,
             2.90945761e-03,   2.82912876e-03,   2.70127141e-03,
             2.51202011e-03,   2.24617193e-03,   1.88032853e-03,
             1.35451064e-03,   1.00000000e-08,   1.00000000e-08,
             2.95578314e-03,   2.90945761e-03,   2.82912876e-03,
             2.70127141e-03,   2.51202011e-03,   2.24617193e-03,
             1.88032853e-03,   1.35451064e-03,   1.00000000e-08,
             1.00000000e-08,   2.95578314e-03,   2.90945761e-03,
             2.82912876e-03,   2.70127141e-03,   2.51202011e-03,
             2.24617193e-03,   1.88032853e-03,   1.35451064e-03,
             1.00000000e-08,   1.00000000e-08,   1.00000000e-08,
             1.00000000e-08,   1.00000000e-08,   1.00000000e-08,
             1.00000000e-08,   1.00000000e-08,   1.00000000e-08,
             1.00000000e-08,   1.00000000e-08])
    
    """

    _name = 'KinematicWave'

    _input_var_names = (
        'topographic__elevation',
        'surface_water__depth',
    )

    _output_var_names = (
        'surface_water__depth'
    )

    _var_units = {
        'topographic__elevation': 'm',
        'surface_water__depth': 'm'
    }

    _var_mapping = {
        'topographic__elevation': 'node',
        'surface_water__depth': 'node'
    }

    _var_doc = {
        'topographic__elevation': 'Land surface topographic elevation',
        'surface_water__depth': 'Depth of water above the surface'
    }

    @use_file_name_or_kwds
    def __init__(self, grid, mannings_n=0.03, critical_flow_depth=0.003,
                 mannings_epsilon=0.33333333, dt_max=0.3, max_courant=0.2,
                 min_surface_water_depth=1.e-8, **kwds):
        """Initialize the kinematic wave approximation overland flow component.
        """

        assert isinstance(grid, RasterModelGrid), 'grid must be regular'
        self._grid = grid
        active = np.where(self.grid.status_at_node != CLOSED_BOUNDARY)[0]
        assert not np.all(self.grid.status_at_node ==
                          FIXED_GRADIENT_BOUNDARY)
        self._h = self.grid.at_node['surface_water__depth']
        self._active = active
        self._hc = critical_flow_depth
        self._n = mannings_n
        self._negepsilon = -mannings_epsilon
        assert not np.isclose(dt_max, 0.)
        self.dt_max = dt_max
        self.min_surface_water_depth = min_surface_water_depth
        self._active_depths = self.grid.at_node[
            'surface_water__depth'][active]
        all_grads = self.grid.calculate_gradients_at_links(
            'topographic__elevation')
        hoz_grads = self.grid.map_mean_of_horizontal_active_links_to_node(
            all_grads)
        vert_grads = self.grid.map_mean_of_vertical_active_links_to_node(
            all_grads)
        self.hozslopept5 = np.fabs(hoz_grads[active])**0.5
        self.vertslopept5 = np.fabs(vert_grads[active])**0.5
        self.velx = self.grid.zeros('node', dtype=float)
        self.vely = self.grid.zeros('node', dtype=float)
        self.qy = np.zeros(grid.number_of_nodes+1, dtype=float)
        self.qx = np.zeros(grid.number_of_nodes+1, dtype=float)
        self.poshozgrads = hoz_grads > 0.
        self.posvertgrads = vert_grads > 0.
        if np.isclose(self.grid.dx, self.grid.dy):
            self.equaldims = True
            self.courant_prefactor = max_courant*self.grid.dx
        else:
            self.equaldims = False
            self.courant_prefactor = max_courant*self.grid.dx*self.grid.dy
        self._neighbors = self.grid.neighbors_at_node.copy()
        self._neighbors[self._neighbors == BAD_INDEX_VALUE] = -1
        self._water_balance = []
        self.actives_BCs = (self.grid.status_at_node[active] ==
                            FIXED_VALUE_BOUNDARY)
        self.actives_BCs_water_depth = self._h[active][self.actives_BCs]

    def update_one_timestep(self, dt, rainfall_intensity=0.00001,
                            update_topography=False, track_min_depth=False):
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
        elapsed_time_in_dt = 0.  # this is only since the start of the timestep
        active = self._active
        hnew = self._h[active]
        if update_topography:
            self.update_topographic_params()
        while elapsed_time_in_dt < dt:
            # assert the minimum water depth - this could introduce an
            # element of mass gain, but should remain minor
            hnew = hnew.clip(self.min_surface_water_depth)
            if track_min_depth:
                self._water_balance.append(
                    (hnew-self._h[active]).sum()/self._h[active].sum())
            n = self._n * (hnew/self._hc)**self._negepsilon
            twothirdshnewbyn = hnew**0.66666666 / n
            self.vely[active] = twothirdshnewbyn * self.vertslopept5
            self.velx[active] = twothirdshnewbyn * self.hozslopept5
            self.vely[self.posvertgrads] *= -1.
            self.velx[self.poshozgrads] *= -1.
            self.qy[active] = self.vely[active] * hnew  # m**2/s
            self.qx[active] = self.velx[active] * hnew  # m**2/s
            maxvely = np.fabs(self.vely).max()
            maxvelx = np.fabs(self.velx).max()
            if self.equaldims:
                courant_dt = self.courant_prefactor/(maxvelx + maxvely)
            else:
                # note prefactor is NOT THE SAME as above in this case
                courant_dt = self.courant_prefactor/(self.grid.dy*maxvelx +
                                                     self.grid.dx*maxvely)
            if self.dt_max is not None:
                internal_dt = np.min((self.dt_max, courant_dt))
            else:
                internal_dt = courant_dt
            remaining_dt = dt - elapsed_time_in_dt
            # now reduce timestep is needed if limited by total tstep length
            internal_dt = min(internal_dt, remaining_dt).clip(0.)
            # this section uses our final-array-val-is-zero trick
            qx_left = self.qx[self._neighbors[:, 2]].clip(
                min=0.)
            qx_right = self.qx[self._neighbors[:, 0]].clip(
                max=0.)
            qy_top = self.qy[self._neighbors[:, 1]].clip(
                min=0.)
            qy_bottom = self.qy[self._neighbors[:, 3]].clip(
                max=0.)
            # FR's rainfall handling was here. We're going to assume that the
            # component is being driven by a "LL style" rainfall record, where
            # the provided rainfall_intensity is constant across the provide
            # dt. If it's not, it needs to be handled outside the component.

            # now add the rainfall input
            if type(rainfall_intensity) is not np.ndarray:
                hnew += internal_dt * rainfall_intensity
            else:
                hnew += internal_dt * rainfall_intensity[active]
            # set the BCs
            hnew[self.actives_BCs] = self.actives_BCs_water_depth
            # flux it round
            hnew -= internal_dt/self.grid.dx*np.fabs(self.qy[active])
            hnew -= internal_dt/self.grid.dy*np.fabs(self.qx[active])
            hnew += internal_dt/self.grid.dx*(qy_top - qy_bottom)[active]
            hnew += internal_dt/self.grid.dy*(qx_left - qx_right)[active]
            # update the internal clock
            elapsed_time_in_dt += internal_dt

        # update the actual field
        self._h[active] = hnew

    def update_topographic_params(self):
        """
        If the topo changes during the run, change the held params used by
        :func:`update_one_timestep`.
        """
        active = np.where(self.grid.status_at_node != CLOSED_BOUNDARY)[0]
        all_grads = self.grid.calculate_gradients_at_links(
            'topographic__elevation')
        hoz_grads = self.grid.map_mean_of_horizontal_active_links_to_node(
            all_grads)
        vert_grads = self.grid.map_mean_of_vertical_active_links_to_node(
            all_grads)
        self.hozslopept5 = np.fabs(hoz_grads[active])**0.5
        self.vertslopept5 = np.fabs(vert_grads[active])**0.5
        self.poshozgrads = hoz_grads > 0.
        self.posvertgrads = vert_grads > 0.
        # check is the grid topology has changed...
        if not np.all(np.equal(self._active, active)):
            self._active = active
            self.velx.fill(0.)
            self.vely.fill(0.)
            self.qy.fill(0.)
            self.qx.fill(0.)
            self._neighbors = self.grid.neighbors_at_node.copy()
            self._neighbors[self._neighbors == BAD_INDEX_VALUE] = -1
            self.actives_BCs = (self.grid.status_at_node[active] ==
                                FIXED_VALUE_BOUNDARY)
            self.actives_BCs_water_depth = self._h[self.actives_BCs]

    @property
    def water_balance(self):
        """
        Return a list of the fractional gain/loss of water mass during the
        run, if it was tracked using the track_min_depth flag.
        """
        if self._water_balance == []:
            raise ValueError('No record of water balance was found!')
        else:
            return self._water_balance
