"""Grid-based simulation of lateral erosion by channels in a drainage network.

Benjamin Campforts
"""

import numpy as np

from landlab import Component
from landlab import RasterModelGrid
from landlab.components.depression_finder.floodstatus import FloodStatus
from landlab.grid.nodestatus import NodeStatus
from landlab.utils.return_array import return_array_at_node

from .ext.calc_sequential_ero_depo import _sequential_ero_depo

ROOT2 = np.sqrt(2.0)  # syntactic sugar for precalculated square root of 2
TIME_STEP_FACTOR = 0.5  # factor used in simple subdivision solver


class SpaceLargeScaleEroder(Component):
    """Stream Power with Alluvium Conservation and Entrainment (SPACE) large scale eroder

    The SPACE_large_Scale_eroder is based on the SPACE component and is designed
    to be more robust against large time steps and coded in such a way that mass
    conservation is explicitly conserved during calculation.

    See the publication:
    Shobe, C. M., Tucker, G. E., and Barnhart, K. R.: The SPACE 1.0 model: a
    Landlab component for 2-D calculation of sediment transport, bedrock
    erosion, and landscape evolution, Geosci. Model Dev., 10, 4577-4604,
    `https://doi.org/10.5194/gmd-10-4577-2017 <https://www.geosci-model-dev.net/10/4577/2017/>`_, 2017.

    Unlike other some other fluvial erosion componets in Landlab, in this
    component (and :py:class:`~landlab.components.ErosionDeposition`) no
    erosion occurs in depressions or in areas with adverse slopes. There is no
    ability to pass a keyword argument ``erode_flooded_nodes``.

    If a depressions are handled (as indicated by the presence of the field
    "flood_status_code" at nodes), then deposition occurs throughout the
    depression and sediment is passed out of the depression. Where pits are
    encountered, then all sediment is deposited at that node only.

    Note: In the current version, we do not provide an adaptive time stepper.
    This will be addded in future versions of this component.

    For more explanation and examples,
    check out the correponding notebook of this component

    Examples
    ---------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import PriorityFloodFlowRouter, SpaceLargeScaleEroder
    >>> import matplotlib.pyplot as plt  # For plotting results; optional
    >>> from landlab import imshow_grid  # For plotting results; optional

    >>> num_rows = 20
    >>> num_columns = 20
    >>> node_spacing = 100.0
    >>> mg = RasterModelGrid((num_rows, num_columns), xy_spacing=node_spacing)
    >>> node_next_to_outlet = num_columns + 1
    >>> np.random.seed(seed=5000)
    >>> _ = mg.add_zeros("topographic__elevation", at="node")
    >>> _ = mg.add_zeros("soil__depth", at="node")
    >>> mg.at_node["soil__depth"][mg.core_nodes] = 2.0
    >>> _ = mg.add_zeros("bedrock__elevation", at="node")
    >>> mg.at_node["bedrock__elevation"] += (
    ...     mg.node_y / 10.0 + mg.node_x / 10.0 + np.random.rand(len(mg.node_y)) / 10.0
    ... )
    >>> mg.at_node["bedrock__elevation"][:] = mg.at_node["topographic__elevation"]
    >>> mg.at_node["topographic__elevation"][:] += mg.at_node["soil__depth"]
    >>> mg.set_closed_boundaries_at_grid_edges(
    ...     bottom_is_closed=True,
    ...     left_is_closed=True,
    ...     right_is_closed=True,
    ...     top_is_closed=True,
    ... )
    >>> mg.set_watershed_boundary_condition_outlet_id(
    ...     0, mg.at_node["topographic__elevation"], -9999.0
    ... )

    >>> fr = PriorityFloodFlowRouter(mg, flow_metric="D8", suppress_out=True)
    >>> sp = SpaceLargeScaleEroder(
    ...     mg,
    ...     K_sed=0.01,
    ...     K_br=0.001,
    ...     F_f=0.0,
    ...     phi=0.0,
    ...     H_star=1.0,
    ...     v_s=5.0,
    ...     m_sp=0.5,
    ...     n_sp=1.0,
    ...     sp_crit_sed=0,
    ...     sp_crit_br=0,
    ... )
    >>> timestep = 10.0
    >>> elapsed_time = 0.0
    >>> count = 0
    >>> run_time = 1e4
    >>> sed_flux = np.zeros(int(run_time // timestep))
    >>> while elapsed_time < run_time:
    ...     fr.run_one_step()
    ...     _ = sp.run_one_step(dt=timestep)
    ...     sed_flux[count] = mg.at_node["sediment__flux"][node_next_to_outlet]
    ...     elapsed_time += timestep
    ...     count += 1
    ...

    Plot the results.

    >>> fig = plt.figure()
    >>> plot = plt.subplot()
    >>> _ = imshow_grid(
    ...     mg,
    ...     "topographic__elevation",
    ...     plot_name="Sediment flux",
    ...     var_name="Sediment flux",
    ...     var_units=r"m$^3$/yr",
    ...     grid_units=("m", "m"),
    ...     cmap="terrain",
    ... )
    >>> _ = plt.figure()
    >>> _ = imshow_grid(
    ...     mg,
    ...     "sediment__flux",
    ...     plot_name="Sediment flux",
    ...     var_name="Sediment flux",
    ...     var_units=r"m$^3$/yr",
    ...     grid_units=("m", "m"),
    ...     cmap="terrain",
    ... )
    >>> fig = plt.figure()
    >>> sedfluxplot = plt.subplot()
    >>> _ = sedfluxplot.plot(
    ...     np.arange(len(sed_flux)) * timestep, sed_flux, color="k", linewidth=1.0
    ... )
    >>> _ = sedfluxplot.set_xlabel("Time [yr]")
    >>> _ = sedfluxplot.set_ylabel(r"Sediment flux [m$^3$/yr]")

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Shobe, C., Tucker, G., Barnhart, K. (2017). The SPACE 1.0 model: a Landlab
    component for 2-D calculation of sediment transport, bedrock erosion, and
    landscape evolution. Geoscientific Model Development  10(12), 4577 - 4604.
    https://dx.doi.org/10.5194/gmd-10-4577-2017

    **Additional References**

    None Listed

    """  # noqa: B950

    _name = "SpaceLargeScaleEroder"

    _unit_agnostic = True

    _info = {
        "flow__link_to_receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": True,
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
        "sediment__influx": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m3/s",
            "mapping": "node",
            "doc": "Sediment flux (volume per unit time of sediment entering each node)",
        },
        "sediment__outflux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m3/s",
            "mapping": "node",
            "doc": "Sediment flux (volume per unit time of sediment leaving each node)",
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
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "The steepest *downhill* slope",
        },
    }

    _cite_as = """
    @Article{gmd-10-4577-2017,
        AUTHOR = {Shobe, C. M. and Tucker, G. E. and Barnhart, K. R.},
        TITLE = {The SPACE~1.0 model: a~Landlab component for 2-D calculation
                 of sediment transport, bedrock erosion, and landscape evolution},
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
        v_s_lake=None,
        m_sp=0.5,
        n_sp=1.0,
        sp_crit_sed=0.0,
        sp_crit_br=0.0,
        discharge_field="surface_water__discharge",
        erode_flooded_nodes=False,
        thickness_lim=100,
    ):
        """Initialize the SpaceLargeScaleEroder model.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        K_sed : float, array of float, or str, optional
            Erodibility for sediment (units vary) as either a number or a field name.
        K_br : float, array of float, or str, optional
            Erodibility for bedrock (units vary) as either a number or a field name.
        F_f : float, optional
            Fraction of permanently suspendable fines in bedrock [-].
        phi : float, optional
            Sediment porosity [-].
        H_star : float, optional
            Sediment thickness required for full entrainment [L].
        v_s : float, optional
            Effective settling velocity for chosen grain size metric [L/T].
        v_s_lake : float, optional
            Effective settling velocity in lakes for chosen grain size metric [L/T].
        m_sp : float, optional
            Drainage area exponent (units vary).
        n_sp : float, optional
            Slope exponent (units vary).
        sp_crit_sed : float, array of float, or str, optional
            Critical stream power to erode sediment [E/(TL^2)].
        sp_crit_br : float, array of float, or str, optional
            Critical stream power to erode rock [E/(TL^2)]
        discharge_field : float, array of float, or str, optional
            Discharge [L^2/T]. The default is to use the grid field
            'surface_water__discharge', which is simply drainage area
            multiplied by the default rainfall rate (1 m/yr). To use custom
            spatially/temporally varying rainfall, use 'water__unit_flux_in'
            to specify water input to the FlowAccumulator.
        erode_flooded_nodes : bool, optional
            Whether erosion occurs in flooded nodes identified by a
            depression/lake mapper (e.g., DepressionFinderAndRouter). When set
            to false, the field *flood_status_code* must be present on the grid
            (this is created by the DepressionFinderAndRouter). Default True.
        """
        if grid.at_node["flow__receiver_node"].size != grid.size("node"):
            raise NotImplementedError(
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that SpaceLargeScaleEroder is compatible with "
                "route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )

        super().__init__(grid)

        self._soil__depth = grid.at_node["soil__depth"]
        self._topographic__elevation = grid.at_node["topographic__elevation"]

        if "bedrock__elevation" in grid.at_node:
            self._bedrock__elevation = grid.at_node["bedrock__elevation"]
        else:
            self._bedrock__elevation = grid.add_zeros(
                "bedrock__elevation", at="node", dtype=float
            )

            self._bedrock__elevation[:] = (
                self._topographic__elevation - self._soil__depth
            )

        # Check consistency of bedrock, soil and topogarphic elevation fields
        np.testing.assert_almost_equal(
            grid.at_node["bedrock__elevation"] + grid.at_node["soil__depth"],
            grid.at_node["topographic__elevation"],
            decimal=5,
            err_msg=(
                "The sum of bedrock elevation and topographic elevation should "
                "be equal"
            ),
        )

        # specific inits
        self._thickness_lim = thickness_lim
        self._H_star = H_star

        self._sed_erosion_term = np.zeros(grid.number_of_nodes)
        self._br_erosion_term = np.zeros(grid.number_of_nodes)
        self._Es = np.zeros(grid.number_of_nodes)
        self._Er = np.zeros(grid.number_of_nodes)

        # K's and critical values can be floats, grid fields, or arrays
        # use setters defined below
        self._K_sed = K_sed
        self._K_br = K_br

        self._sp_crit_sed = return_array_at_node(grid, sp_crit_sed)
        self._sp_crit_br = return_array_at_node(grid, sp_crit_br)

        self._erode_flooded_nodes = erode_flooded_nodes

        self._flow_receivers = grid.at_node["flow__receiver_node"]
        self._stack = grid.at_node["flow__upstream_node_order"]
        self._slope = grid.at_node["topographic__steepest_slope"]

        self.initialize_output_fields()

        self._qs = grid.at_node["sediment__outflux"]
        self._q = return_array_at_node(grid, discharge_field)

        # for backward compatibility (remove in 3.0.0+)
        grid.at_node["sediment__flux"] = grid.at_node["sediment__outflux"]

        self._Q_to_the_m = np.zeros(grid.number_of_nodes)
        self._S_to_the_n = np.zeros(grid.number_of_nodes)

        # store other constants
        self._m_sp = np.float64(m_sp)
        self._n_sp = np.float64(n_sp)
        self._phi = np.float64(phi)
        self._v_s = np.float64(v_s)

        if isinstance(grid, RasterModelGrid):
            self._link_lengths = grid.length_of_d8
        else:
            self._link_lengths = grid.length_of_link

        if v_s_lake is None:
            self._v_s_lake = np.float64(v_s)
        else:
            self._v_s_lake = np.float64(v_s_lake)
        self._F_f = np.float64(F_f)

        if phi >= 1.0:
            raise ValueError("Porosity must be < 1.0")

        if F_f > 1.0:
            raise ValueError("Fraction of fines must be <= 1.0")

        if phi < 0.0:
            raise ValueError("Porosity must be > 0.0")

        if F_f < 0.0:
            raise ValueError("Fraction of fines must be > 0.0")

    @property
    def K_br(self):
        """Erodibility of bedrock(units depend on m_sp)."""
        return self._K_br

    @K_br.setter
    def K_br(self, new_val):
        self._K_br = return_array_at_node(self._grid, new_val)

    @property
    def K_sed(self):
        """Erodibility of sediment(units depend on m_sp)."""
        return self._K_sed

    @K_sed.setter
    def K_sed(self, new_val):
        self._K_sed = return_array_at_node(self._grid, new_val)

    @property
    def fraction_fines(self):
        """Fraction of permanently suspendable fines in bedrock [-]."""
        return self._F_f

    @property
    def sediment_porosity(self):
        """Sediment porosity [-]."""
        return self._phi

    @property
    def settling_velocity(self):
        """Effective settling velocity for chosen grain size metric [L/T]."""
        return self._v_s

    @property
    def drainage_area_exp(self):
        """Drainage area exponent (units vary)."""
        return self._m_sp

    @property
    def slope_exp(self):
        """Slope exponent (units vary)."""
        return self._n_sp

    @property
    def Es(self):
        """Sediment erosion term."""
        return self._Es

    @property
    def Er(self):
        """Bedrock erosion term."""
        return self._Er

    @property
    def sediment_influx(self):
        """Volumetric sediment influx to each node."""
        return self.grid.at_node["sediment__influx"]

    def _calc_erosion_rates(self):
        """Calculate erosion rates."""

        H = self.grid.at_node["soil__depth"]

        # if sp_crits are zero, then this colapses to correct all the time.
        if np.isclose(self._n_sp, 1.0):
            S_to_the_n = self._slope
        else:
            S_to_the_n = np.power(self._slope, self._n_sp)
        omega_sed = self._K_sed * self._Q_to_the_m * S_to_the_n
        omega_br = self._K_br * self._Q_to_the_m * S_to_the_n

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
        ) / (
            1 - self._phi
        )  # convert from a volume to a mass flux.
        self._br_erosion_term = omega_br - self._sp_crit_br * (
            1.0 - np.exp(-omega_br_over_sp_crit)
        )

        self._Es = self._sed_erosion_term * (1.0 - np.exp(-H / self._H_star))
        self._Er = self._br_erosion_term * np.exp(-H / self._H_star)

        # if the soil layer becomes exceptionally thick (e.g. because of
        # landslide derived sediment deposition(,) the algorithm will become
        # unstable because np.exp(x) with x > 709 yeilds inf values.
        # Therefore soil depth is temporqlly topped of at 200m and the remaining
        # values are added back after the space component has run

        self._Es[H > self._thickness_lim] = self._sed_erosion_term[
            H > self._thickness_lim
        ]
        self._Er[H > self._thickness_lim] = 0

    def run_one_step_basic(self, dt=10):
        node_status = self.grid.status_at_node

        z = self.grid.at_node["topographic__elevation"]
        br = self.grid.at_node["bedrock__elevation"]
        H = self.grid.at_node["soil__depth"]
        link_to_rcvr = self.grid.at_node["flow__link_to_receiver_node"]
        area = self.grid.cell_area_at_node

        r = self.grid.at_node["flow__receiver_node"]
        stack_flip_ud = np.flipud(self.grid.at_node["flow__upstream_node_order"])
        # Select core nodes where qs >0
        stack_flip_ud_sel = stack_flip_ud[
            (node_status[stack_flip_ud] == NodeStatus.CORE)
            & (self._q[stack_flip_ud] > 0.0)
        ]
        slope = (z - z[r]) / self._link_lengths[link_to_rcvr]

        # Choose a method for calculating erosion:
        self._Q_to_the_m[:] = np.power(self._q, self._m_sp)
        self._calc_erosion_rates()

        if "flood_status_code" in self.grid.at_node:
            flood_status = self.grid.at_node["flood_status_code"]
            flooded_nodes = np.nonzero(flood_status == FloodStatus.FLOODED)[0]
        else:
            flooded_nodes = np.nonzero([slope < 0])[1]

        self._Es[flooded_nodes] = 0.0
        self._Er[flooded_nodes] = 0.0
        self._sed_erosion_term[flooded_nodes] = 0.0
        self._br_erosion_term[flooded_nodes] = 0.0

        self.sediment_influx[:] = 0

        K_sed_vector = np.broadcast_to(self._K_sed, self._q.shape)

        vol_SSY_riv = _sequential_ero_depo(
            stack_flip_ud_sel,
            r,
            area,
            self._q,
            self._qs,
            self.sediment_influx,
            self._Es,
            self._Er,
            self._Q_to_the_m,
            slope,
            H,
            br,
            self._sed_erosion_term,
            self._br_erosion_term,
            K_sed_vector,
            self._v_s,
            self._phi,
            self._F_f,
            self._H_star,
            dt,
            self._thickness_lim,
        )

        V_leaving_riv = np.sum(self.sediment_influx[self.grid.boundary_nodes]) * dt
        # Update topography
        cores = self._grid.core_nodes
        z[cores] = br[cores] + H[cores]

        return vol_SSY_riv, V_leaving_riv

    def run_one_step(self, dt):
        vol_SSY_riv, V_leaving_riv = self.run_one_step_basic(dt)
        return vol_SSY_riv, V_leaving_riv
