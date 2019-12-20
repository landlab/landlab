import numpy as np

from landlab.components.erosion_deposition.generalized_erosion_deposition import (
    DEFAULT_MINIMUM_TIME_STEP,
    _GeneralizedErosionDeposition,
)
from landlab.utils.return_array import return_array_at_node

from ..depression_finder.lake_mapper import _FLOODED
from .cfuncs import calculate_qs_in

ROOT2 = np.sqrt(2.0)  # syntactic sugar for precalculated square root of 2
TIME_STEP_FACTOR = 0.5  # factor used in simple subdivision solver


class Space(_GeneralizedErosionDeposition):
    """Stream Power with Alluvium Conservation and Entrainment (SPACE)

    See the publication:
    Shobe, C. M., Tucker, G. E., and Barnhart, K. R.: The SPACE 1.0 model: a
    Landlab component for 2-D calculation of sediment transport, bedrock
    erosion, and landscape evolution, Geosci. Model Dev., 10, 4577-4604,
    `https://doi.org/10.5194/gmd-10-4577-2017 <https://www.geosci-model-dev.net/10/4577/2017/>`_, 2017.

    Note: If timesteps are large enough that Es*dt (sediment erosion)
    exceeds sediment thickness H, the 'adaptive' solver is necessary to
    subdivide timesteps. Compare Es and H arrays to determine whether
    timesteps are appropriate or too large for the 'basic' solver.

    Parameters
    ----------
    grid : ModelGrid
        Landlab ModelGrid object
    K_sed : float, field name, or array
        Erodibility for sediment (units vary).
    K_br : float, field name, or array
        Erodibility for bedrock (units vary).
    F_f : float
        Fraction of permanently suspendable fines in bedrock [-].
    phi : float
        Sediment porosity [-].
    H_star : float
        Sediment thickness required for full entrainment [L].
    v_s : float
        Effective settling velocity for chosen grain size metric [L/T].
    m_sp : float
        Drainage area exponent (units vary)
    n_sp : float
        Slope exponent (units vary)
    sp_crit_sed : float, field name, or array
        Critical stream power to erode sediment [E/(TL^2)]
    sp_crit_br : float, field name, or array
        Critical stream power to erode rock [E/(TL^2)]
    discharge_field : float, field name, or array
        Discharge [L^2/T]. The default is to use the grid field
        'surface_water__discharge', which is simply drainage area
        multiplied by the default rainfall rate (1 m/yr). To use custom
        spatially/temporally varying rainfall, use 'water__unit_flux_in'
        to specify water input to the FlowAccumulator.
    erode_flooded_nodes : bool (optional)
        Whether erosion occurs in flooded nodes identified by a
        depression/lake mapper (e.g., DepressionFinderAndRouter). When set
        to false, the field *flood_status_code* must be present on the grid
        (this is created by the DepressionFinderAndRouter). Default True.
    solver : string
        Solver to use. Options at present include:
            (1) 'basic' (default): explicit forward-time extrapolation.
                Simple but will become unstable if time step is too large.
            (2) 'adaptive': subdivides global time step as needed to
                prevent slopes from reversing and alluvium from going
                negative.

    Examples
    ---------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import (FlowAccumulator,
    ...                                 DepressionFinderAndRouter,
    ...                                 Space,
    ...                                 FastscapeEroder)
    >>> np.random.seed(seed = 5000)

    Define grid and initial topography:

    *  5x5 grid with baselevel in the lower left corner
    *  All other boundary nodes closed
    *  Initial topography is plane tilted up to the upper right with
       noise

    >>> mg = RasterModelGrid((5, 5), xy_spacing=10.0)
    >>> _ = mg.add_zeros('topographic__elevation', at='node')
    >>> mg.at_node['topographic__elevation'] += (mg.node_y / 10. +
    ...     mg.node_x / 10. + np.random.rand(len(mg.node_y)) / 10.)
    >>> mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
    ...                                        left_is_closed=True,
    ...                                        right_is_closed=True,
    ...                                        top_is_closed=True)
    >>> mg.set_watershed_boundary_condition_outlet_id(
    ...     0, mg.at_node['topographic__elevation'], -9999.)
    >>> fsc_dt = 100.
    >>> space_dt = 100.

    Instantiate Fastscape eroder, flow router, and depression finder

    >>> fr = FlowAccumulator(mg, flow_director='D8')
    >>> df = DepressionFinderAndRouter(mg)
    >>> fsc = FastscapeEroder(
    ...     mg,
    ...     K_sp=.001,
    ...     m_sp=.5,
    ...     n_sp=1,
    ...     erode_flooded_nodes=False)

    Burn in an initial drainage network using the Fastscape eroder:

    >>> for x in range(100):
    ...     fr.run_one_step()
    ...     df.map_depressions()
    ...     fsc.run_one_step(dt=fsc_dt)
    ...     mg.at_node['topographic__elevation'][0] -= 0.001 # Uplift

    Add some soil to the drainage network:

    >>> _ = mg.add_zeros('soil__depth', at='node', dtype=float)
    >>> mg.at_node['soil__depth'] += 0.5
    >>> mg.at_node['topographic__elevation'] += mg.at_node['soil__depth']

    Instantiate the Space component:

    >>> ha = Space(
    ...     mg,
    ...     K_sed=0.00001,
    ...     K_br=0.00000000001,
    ...     F_f=0.5,
    ...     phi=0.1,
    ...     H_star=1.,
    ...     v_s=0.001,
    ...     m_sp=0.5,
    ...     n_sp = 1.0,
    ...     sp_crit_sed=0,
    ...     sp_crit_br=0,
    ...     erode_flooded_nodes=False)

    Now run the Space component for 2000 short timesteps:

    >>> for x in range(2000): #Space component loop
    ...     fr.run_one_step()
    ...     df.map_depressions()
    ...     ha.run_one_step(dt=space_dt)
    ...     mg.at_node['bedrock__elevation'][0] -= 2e-6 * space_dt

    Now we test to see if soil depth and topography are right:

    >>> np.around(mg.at_node['soil__depth'], decimals=3) # doctest: +NORMALIZE_WHITESPACE
    array([ 0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.495,  0.492,
            0.491,  0.5  ,  0.5  ,  0.492,  0.492,  0.49 ,  0.5  ,  0.5  ,
            0.491,  0.49 ,  0.484,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,
            0.5  ])

    >>> np.around(mg.at_node['topographic__elevation'], decimals=3) # doctest: +NORMALIZE_WHITESPACE
    array([ 0.423,  1.536,  2.573,  3.511,  4.561,  1.582,  0.424,  0.428,
            0.438,  5.51 ,  2.54 ,  0.428,  0.428,  0.438,  6.526,  3.559,
            0.438,  0.438,  0.45 ,  7.553,  4.559,  5.541,  6.57 ,  7.504,
            8.51 ])
    """

    _name = "Space"

    _info = {
        "flow__link_to_receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "ID of link downstream of each node, which carries the discharge",
        },
        "flow__receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        "flow__upstream_node_order": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array containing downstream-to-upstream ordered list of node IDs",
        },
        "sediment__flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m3/s",
            "mapping": "node",
            "doc": "Sediment flux (volume per unit time of sediment entering each node)",
        },
        "soil__depth": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of soil or weathered bedrock",
        },
        "surface_water__discharge": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": "Volumetric discharge of surface water",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "topographic__steepest_slope": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "The steepest *downhill* slope",
        },
    }

    _cite_as = """@Article{gmd-10-4577-2017,
                  AUTHOR = {Shobe, C. M. and Tucker, G. E. and Barnhart, K. R.},
                  TITLE = {The SPACE~1.0 model: a~Landlab component for 2-D calculation of sediment transport, bedrock erosion, and landscape evolution},
                  JOURNAL = {Geoscientific Model Development},
                  VOLUME = {10},
                  YEAR = {2017},
                  NUMBER = {12},
                  PAGES = {4577--4604},
                  URL = {https://www.geosci-model-dev.net/10/4577/2017/},
                  DOI = {10.5194/gmd-10-4577-2017}
                  }"""

    def __init__(
        self,
        grid,
        K_sed=0.02,
        K_br=0.02,
        F_f=0.0,
        phi=0.3,
        H_star=0.1,
        v_s=1.0,
        m_sp=0.5,
        n_sp=1.0,
        sp_crit_sed=0.0,
        sp_crit_br=0.0,
        discharge_field="surface_water__discharge",
        solver="basic",
        erode_flooded_nodes=True,
        dt_min=DEFAULT_MINIMUM_TIME_STEP,
    ):
        """Initialize the Space model."""
        if grid.at_node["flow__receiver_node"].size != grid.size("node"):
            msg = (
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that SPACE is compatible with "
                "route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )
            raise NotImplementedError(msg)

        super(Space, self).__init__(
            grid,
            m_sp=m_sp,
            n_sp=n_sp,
            phi=phi,
            F_f=F_f,
            v_s=v_s,
            dt_min=dt_min,
            discharge_field=discharge_field,
            erode_flooded_nodes=erode_flooded_nodes,
        )

        # space specific inits
        self._H_star = H_star
        self._soil__depth = grid.at_node["soil__depth"]

        if "bedrock__elevation" in grid.at_node:
            self._bedrock__elevation = grid.at_node["bedrock__elevation"]
        else:
            self._bedrock__elevation = grid.add_zeros(
                "bedrock__elevation", at="node", dtype=float
            )

            self._bedrock__elevation[:] = (
                self._topographic__elevation - self._soil__depth
            )

        self._Es = np.zeros(grid.number_of_nodes)
        self._Er = np.zeros(grid.number_of_nodes)

        # K's and critical values can be floats, grid fields, or arrays
        self._K_sed = return_array_at_node(grid, K_sed)
        self._K_br = return_array_at_node(grid, K_br)

        self._sp_crit_sed = return_array_at_node(grid, sp_crit_sed)
        self._sp_crit_br = return_array_at_node(grid, sp_crit_br)

        # Handle option for solver
        if solver == "basic":
            self.run_one_step = self.run_one_step_basic
        elif solver == "adaptive":
            self.run_one_step = self.run_with_adaptive_time_step_solver
            self._time_to_flat = np.zeros(grid.number_of_nodes)
            self._porosity_factor = 1.0 / (1.0 - self._phi)
        else:
            raise ValueError(
                "Parameter 'solver' must be one of: " + "'basic', 'adaptive'"
            )

    def _calc_erosion_rates(self):
        """Calculate erosion rates."""
        # if sp_crits are zero, then this colapses to correct all the time.
        omega_sed = self._K_sed * self._Q_to_the_m * np.power(self._slope, self._n_sp)
        omega_br = self._K_br * self._Q_to_the_m * np.power(self._slope, self._n_sp)

        omega_sed_over_sp_crit = np.divide(
            omega_sed,
            self._sp_crit_sed,
            out=np.zeros_like(omega_sed),
            where=self._sp_crit_sed != 0,
        )

        omega_br_over_sp_crit = np.divide(
            omega_br,
            self._sp_crit_br,
            out=np.zeros_like(omega_br),
            where=self._sp_crit_br != 0,
        )

        self._sed_erosion_term = omega_sed - self._sp_crit_sed * (
            1.0 - np.exp(-omega_sed_over_sp_crit)
        )
        self._br_erosion_term = omega_br - self._sp_crit_br * (
            1.0 - np.exp(-omega_br_over_sp_crit)
        )

        self._Es = self._sed_erosion_term * (
            1.0 - np.exp(-self._soil__depth / self._H_star)
        )
        self._Er = self._br_erosion_term * np.exp(-self._soil__depth / self._H_star)

    @property
    def Es(self):
        """Sediment erosion term."""
        return self._Es

    @property
    def Er(self):
        """Bedrock erosion term."""
        return self._Er

    @property
    def H(self):
        """Sediment thickness."""
        return self._H

    def run_one_step_basic(self, dt=1.0):
        """Calculate change in rock and alluvium thickness for a time period
        'dt'.

        Parameters
        ----------
        dt : float
            Model timestep [T]
        """
        # Choose a method for calculating erosion:
        self._calc_hydrology()
        self._calc_erosion_rates()

        if not self._erode_flooded_nodes:
            flood_status = self._grid.at_node["flood_status_code"]
            flooded_nodes = np.nonzero(flood_status == _FLOODED)[0]
        else:
            flooded_nodes = []

        flooded = np.full(self._grid.number_of_nodes, False, dtype=bool)
        flooded[flooded_nodes] = True

        self._qs_in[:] = 0

        # iterate top to bottom through the stack, calculate qs
        # cythonized version of calculating qs_in
        calculate_qs_in(
            np.flipud(self._stack),
            self._flow_receivers,
            self._cell_area_at_node,
            self._q,
            self._qs,
            self._qs_in,
            self._Es,
            self._Er,
            self._v_s,
            self._F_f,
        )

        self._depo_rate[self._q > 0] = self._qs[self._q > 0] * (
            self._v_s / self._q[self._q > 0]
        )

        # now, the analytical solution to soil thickness in time:
        # need to distinguish D=kqS from all other cases to save from blowup!

        # distinguish cases:
        blowup = self._depo_rate == self._K_sed * self._Q_to_the_m * self._slope

        # first, potential blowup case:
        # positive slopes, not flooded
        pos_not_flood = (self._q > 0) & (blowup) & (self._slope > 0) & (~flooded)
        self._soil__depth[pos_not_flood] = self._H_star * np.log(
            ((self._sed_erosion_term[pos_not_flood] / (1 - self._phi)) / self._H_star)
            * dt
            + np.exp(self._soil__depth[pos_not_flood] / self._H_star)
        )
        # positive slopes, flooded
        pos_flood = (self._q > 0) & (blowup) & (self._slope > 0) & (flooded)
        self._soil__depth[pos_flood] = (
            self._depo_rate[pos_flood] / (1 - self._phi)
        ) * dt

        # non-positive slopes, not flooded
        non_pos_not_flood = (self._q > 0) & (blowup) & (self._slope <= 0) & (~flooded)
        self._soil__depth[non_pos_not_flood] += (
            self._depo_rate[non_pos_not_flood] / (1 - self._phi) * dt
        )

        # more general case:
        pos_not_flood = (self._q > 0) & (~blowup) & (self._slope > 0) & (~flooded)

        self._soil__depth[pos_not_flood] = self._H_star * np.log(
            (
                1
                / (
                    (self._depo_rate[pos_not_flood] / (1 - self._phi))
                    / (self._sed_erosion_term[pos_not_flood] / (1 - self._phi))
                    - 1
                )
            )
            * (
                np.exp(
                    (
                        self._depo_rate[pos_not_flood] / (1 - self._phi)
                        - (self._sed_erosion_term[pos_not_flood] / (1 - self._phi))
                    )
                    * (dt / self._H_star)
                )
                * (
                    (
                        (
                            self._depo_rate[pos_not_flood]
                            / (1 - self._phi)
                            / (self._sed_erosion_term[pos_not_flood] / (1 - self._phi))
                        )
                        - 1
                    )
                    * np.exp(self._soil__depth[pos_not_flood] / self._H_star)
                    + 1
                )
                - 1
            )
        )

        # places where slope <= 0 but not flooded:
        neg_slope_not_flooded = (
            (self._q > 0) & (~blowup) & (self._slope <= 0) & (~flooded)
        )
        self._soil__depth[neg_slope_not_flooded] += (
            self._depo_rate[neg_slope_not_flooded] / (1 - self._phi) * dt
        )

        # flooded nodes:
        flooded_nodes = (self._q > 0) & (~blowup) & (flooded)
        self._soil__depth[flooded_nodes] += (
            self._depo_rate[flooded_nodes] / (1 - self._phi) * dt
        )

        # where discharge exists
        discharge_exists = self._q > 0
        self._bedrock__elevation[discharge_exists] += dt * (
            -self._br_erosion_term[discharge_exists]
            * (np.exp(-self._soil__depth[discharge_exists] / self._H_star))
        )

        # finally, determine topography by summing bedrock and soil
        cores = self._grid.core_nodes
        self._topographic__elevation[cores] = (
            self._bedrock__elevation[cores] + self._soil__depth[cores]
        )

    def run_with_adaptive_time_step_solver(self, dt=1.0):
        """Run step with CHILD-like solver that adjusts time steps to prevent
        slope flattening.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> import numpy as np

        >>> rg = RasterModelGrid((3, 4))
        >>> z = rg.add_zeros('topographic__elevation', at='node')
        >>> z[:] = 0.1 * rg.x_of_node
        >>> H = rg.add_zeros('soil__depth', at='node')
        >>> H += 0.1
        >>> br = rg.add_zeros('bedrock__elevation', at='node')
        >>> br[:] = z - H

        >>> fa = FlowAccumulator(rg, flow_director='FlowDirectorSteepest')
        >>> fa.run_one_step()
        >>> sp = Space(rg, K_sed=1.0, K_br=0.1,
        ...            F_f=0.5, phi=0.0, H_star=1., v_s=1.0,
        ...            m_sp=0.5, n_sp = 1.0, sp_crit_sed=0,
        ...            sp_crit_br=0, solver='adaptive')
        >>> sp.run_one_step(dt=10.0)

        >>> np.round(sp.Es[5:7], 4)
        array([ 0.0029,  0.0074])
        >>> np.round(sp.Er[5:7], 4)
        array([ 0.0032,  0.0085])
        >>> np.round(H[5:7], 3)
        array([ 0.088,  0.078])
        """

        # Initialize remaining_time, which records how much of the global time
        # step we have yet to use up.
        remaining_time = dt

        z = self._grid.at_node["topographic__elevation"]
        br = self._grid.at_node["bedrock__elevation"]
        H = self._grid.at_node["soil__depth"]
        r = self._flow_receivers
        time_to_flat = np.zeros(len(z))
        time_to_zero_alluv = np.zeros(len(z))
        dzdt = np.zeros(len(z))
        cores = self._grid.core_nodes

        first_iteration = True

        if not self._erode_flooded_nodes:
            flood_status = self._grid.at_node["flood_status_code"]
            flooded_nodes = np.nonzero(flood_status == _FLOODED)[0]
        else:
            flooded_nodes = []

        flooded = np.full(self._grid.number_of_nodes, False, dtype=bool)
        flooded[flooded_nodes] = True

        # Outer WHILE loop: keep going until time is used up
        while remaining_time > 0.0:

            # Update all the flow-link slopes.
            #
            # For the first iteration, we assume this has already been done
            # outside the component (e.g., by flow router), but we need to do
            # it ourselves on subsequent iterations.
            if not first_iteration:
                # update the link slopes
                self._update_flow_link_slopes()
                # update where nodes are flooded. This shouuldn't happen because
                # of the dynamic timestepper, but just in case, we update here.
                new_flooded_nodes = np.where(self._slope < 0)[0]
                flooded_nodes = np.asarray(
                    np.unique(np.concatenate((flooded_nodes, new_flooded_nodes))),
                    dtype=np.int64,
                )
            else:
                first_iteration = False

            # Calculate rates of entrainment
            self._calc_hydrology()
            self._calc_erosion_rates()

            # CORRECTION HERE?
            self._Es[flooded_nodes] = 0.0
            self._Er[flooded_nodes] = 0.0

            # Zero out sediment influx for new iteration
            self._qs_in[:] = 0

            calculate_qs_in(
                np.flipud(self._stack),
                self._flow_receivers,
                self._cell_area_at_node,
                self._q,
                self._qs,
                self._qs_in,
                self._Es,
                self._Er,
                self._v_s,
                self._F_f,
            )

            self._depo_rate[self._q > 0] = self._qs[self._q > 0] * (
                self._v_s / self._q[self._q > 0]
            )
            # TODO handle flooded nodes in the above fn

            # Now look at upstream-downstream node pairs, and recording the
            # time it would take for each pair to flatten. Take the minimum.
            dzdt[cores] = self._depo_rate[cores] - (self._Es[cores] + self._Er[cores])
            rocdif = dzdt - dzdt[r]
            zdif = z - z[r]
            time_to_flat[:] = remaining_time

            converging = np.where(rocdif < 0.0)[0]
            time_to_flat[converging] = -(
                TIME_STEP_FACTOR * zdif[converging] / rocdif[converging]
            )
            time_to_flat[np.where(zdif <= 0.0)[0]] = remaining_time

            # From this, find the maximum stable time step with regard to slope
            # evolution.
            dt_max1 = np.amin(time_to_flat)

            # Next we consider time to exhaust regolith
            time_to_zero_alluv[:] = remaining_time
            dHdt = self._porosity_factor * (self._depo_rate - self._Es)
            decreasing_H = np.where(dHdt < 0.0)[0]
            time_to_zero_alluv[decreasing_H] = -(
                TIME_STEP_FACTOR * H[decreasing_H] / dHdt[decreasing_H]
            )

            # Now find the smallest time that would lead to near-empty alluv
            dt_max2 = np.amin(time_to_zero_alluv)

            # Take the smaller of the limits
            dt_max = max(self._dt_min, min(dt_max1, dt_max2))

            # Now a vector operation: apply dzdt and dhdt to all nodes
            br[cores] -= self._Er[cores] * dt_max
            H[cores] += dHdt[cores] * dt_max
            z[cores] = br[cores] + H[cores]

            # Update remaining time and continue
            remaining_time -= dt_max
