import warnings

import numpy as np
import scipy.constants

from landlab import Component
from landlab import MissingKeyError
from landlab.utils.decorators import make_return_array_immutable


class SedDepEroder(Component):
    """
    This module implements sediment flux dependent channel incision
    following::

        E = f(Qs, Qc) * ([a stream power-like term] - [an optional threshold]),

    where E is the bed erosion rate, Qs is the volumetric sediment flux
    into a node, and Qc is the volumetric sediment transport capacity at
    that node.

    This component is under active research and development; proceed with its
    use at your own risk.

    The details of the implementation are a function of the two key
    arguments, *sed_dependency_type* and *Qc*. The former controls the
    shape of the sediment dependent response function f(Qs, Qc), the
    latter controls the way in which sediment transport capacities are
    calculated (primarily, whether a full Meyer-Peter Muller approach is
    used, or whether simpler stream-power-like equations can be assumed).
    For Qc, 'power_law' broadly follows the assumptions in Gasparini et
    al. 2006, 2007; 'MPM' broadly follows those in Hobley et al., 2011.
    Note that a convex-up channel can result in many cases assuming MPM,
    unless parameters b and c are carefully tuned.

    If ``Qc == 'power_law'``::

        E  = K_sp * f(Qs, Qc) * A ** m_sp * S ** n_sp;
        Qc = K_t * A ** m_t * S ** n_t

    If ``Qc == 'MPM'``::

        shear_stress = fluid_density * g * depth * S
                     = fluid_density * g * (mannings_n/k_w) ** 0.6 * (
                       k_Q* A ** c_sp) ** (0.6 * (1. - b_sp)) * S ** 0.7,
                       for consistency with MPM

        E = K_sp * f(Qs, Qc) * (shear_stress ** a_sp - [threshold_sp])

        Qc = 8 * C_MPM * sqrt((sed_density-fluid_density)/fluid_density *
             g * D_char**3) * (shields_stress - threshold_shields)**1.5

        shields_stress = shear_stress / (g * (sed_density-fluid_density) *
                         D_char)

    If you choose Qc='MPM', you may provide thresholds for both channel
    incision and shields number, or alternatively set either or both of
    these threshold dynamically. The minimum shear stress can be made
    equivalent to the Shields number using *set_threshold_from_Dchar*,
    for full consistency with the MPM approach (i.e., the threshold
    becomes a function of the characteristic grain size on the bed). The
    Shields threshold itself can also be a weak function of slope if
    *slope_sensitive_threshold*, following Lamb et al. 2008,
    taustar_c = 0.15 * S ** 0.25.

    The component is able to handle flooded nodes, if created by a lake
    filler. It assumes the flow paths found in the fields already reflect
    any lake routing operations, and then requires the optional argument
    *flooded_depths* be passed to the run method. A flooded depression
    acts as a perfect sediment trap, and will be filled sequentially
    from the inflow points towards the outflow points.

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    Hobley, D. E. J., Sinclair, H. D., Mudd, S. M., and Cowie, P. A.: Field
    calibration of sediment ﬂux dependent river incision, J. Geophys. Res.,
    116, F04017, doi:10.1029/2010JF001935, 2011.

    """

    _name = "SedDepEroder"

    _unit_agnostic = False

    _info = {
        "channel__bed_shear_stress": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "Pa",
            "mapping": "node",
            "doc": (
                "Shear exerted on the bed of the channel, assuming all "
                "discharge travels along a single, self-formed channel"
            ),
        },
        "channel__depth": {
            "dtype": float,
            "intent": "out",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of the a single channel carrying all runoff through the node",
        },
        "channel__discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": (
                "Volumetric water flux of the a single channel carrying all "
                "runoff through the node"
            ),
        },
        "channel__width": {
            "dtype": float,
            "intent": "out",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Width of the a single channel carrying all runoff through the node",
        },
        "channel_sediment__relative_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": (
                "The fluvial_sediment_flux_into_node divided by the "
                "fluvial_sediment_transport_capacity"
            ),
        },
        "channel_sediment__volumetric_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": "Total volumetric fluvial sediment flux brought into the node from upstream",
        },
        "channel_sediment__volumetric_transport_capacity": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": (
                "Volumetric transport capacity of a channel carrying all runoff "
                "through the node, assuming the Meyer-Peter Muller transport equation"
            ),
        },
        "drainage_area": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m**2",
            "mapping": "node",
            "doc": "Upstream accumulated surface area contributing to the node's discharge",
        },
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

    def __init__(
        self,
        grid,
        K_sp=1.0e-6,
        g=scipy.constants.g,
        rock_density=2700,
        sediment_density=2700,
        fluid_density=1000,
        runoff_rate=1.0,
        sed_dependency_type="generalized_humped",
        kappa_hump=13.683,
        nu_hump=1.13,
        phi_hump=4.24,
        c_hump=0.00181,
        Qc="power_law",
        m_sp=0.5,
        n_sp=1.0,
        K_t=1.0e-4,
        m_t=1.5,
        n_t=1.0,
        # these params for Qc='MPM':
        C_MPM=1.0,
        a_sp=1.0,
        b_sp=0.5,
        c_sp=1.0,
        k_w=2.5,
        k_Q=2.5e-7,
        mannings_n=0.05,
        threshold_shear_stress=None,
        Dchar=0.05,
        set_threshold_from_Dchar=True,
        set_Dchar_from_threshold=False,
        threshold_Shields=0.05,
        slope_sensitive_threshold=False,
        # params for model numeric behavior:
        pseudoimplicit_repeats=5,
        return_stream_properties=False,
        # flooded node info
        flooded_depths=None,
    ):
        """Constructor for the class.

        Parameters
        ----------
        grid : a ModelGrid
            A grid.
        K_sp : float (time unit must be *years*)
            K in the stream power equation; the prefactor on the erosion
            equation (units vary with other parameters).
        g : float (m/s**2)
            Acceleration due to gravity.
        rock_density : float (Kg m**-3)
            Bulk intact rock density.
        sediment_density : float (Kg m**-3)
            Typical density of loose sediment on the bed.
        fluid_density : float (Kg m**-3)
            Density of the fluid.
        runoff_rate : float, array or field name (m/s)
            The rate of excess overland flow production at each node (i.e.,
            rainfall rate less infiltration).
        pseudoimplicit_repeats : int
            Number of loops to perform with the pseudoimplicit iterator,
            seeking a stable solution. Convergence is typically rapid.
        return_stream_properties : bool
            Whether to perform a few additional calculations in order to set
            the additional optional output fields, 'channel__width',
            'channel__depth', and 'channel__discharge' (default False).
        sed_dependency_type : {'generalized_humped', 'None', 'linear_decline',
                               'almost_parabolic'}
            The shape of the sediment flux function. For definitions, see
            Hobley et al., 2011. 'None' gives a constant value of 1.
            NB: 'parabolic' is currently not supported, due to numerical
            stability issues at channel heads.
        Qc : {'power_law', 'MPM'}
            Whether to use simple stream-power-like equations for both
            sediment transport capacity and erosion rate, or more complex
            forms based directly on the Meyer-Peter Muller equation and a
            shear stress based erosion model consistent with MPM (per
            Hobley et al., 2011).

        If ``sed_dependency_type == 'generalized_humped'``...

        kappa_hump : float
            Shape parameter for sediment flux function. Primarily controls
            function amplitude (i.e., scales the function to a maximum of 1).
            Default follows Leh valley values from Hobley et al., 2011.
        nu_hump : float
            Shape parameter for sediment flux function. Primarily controls
            rate of rise of the "tools" limb. Default follows Leh valley
            values from Hobley et al., 2011.
        phi_hump : float
            Shape parameter for sediment flux function. Primarily controls
            rate of fall of the "cover" limb. Default follows Leh valley
            values from Hobley et al., 2011.
        c_hump : float
            Shape parameter for sediment flux function. Primarily controls
            degree of function asymmetry. Default follows Leh valley values
            from Hobley et al., 2011.

        If ``Qc == 'power_law'``...

        m_sp : float
            Power on drainage area in the erosion equation.
        n_sp : float
            Power on slope in the erosion equation.
        K_t : float (time unit must be in *years*)
            Prefactor in the transport capacity equation.
        m_t : float
            Power on drainage area in the transport capacity equation.
        n_t : float
            Power on slope in the transport capacity equation.

        if ``Qc == 'MPM'``...

        C_MPM : float
            A prefactor on the MPM relation, allowing tuning to known sediment
            saturation conditions (leave as 1. in most cases).
        a_sp : float
            Power on shear stress to give erosion rate.
        b_sp : float
            Power on drainage area to give channel width.
        c_sp : float
            Power on drainage area to give discharge.
        k_w : float (unit variable with b_sp)
            Prefactor on A**b_sp to give channel width.
        k_Q : float (unit variable with c_sp, but time unit in *seconds*)
            Prefactor on A**c_sp to give discharge.
        mannings_n : float
            Manning's n for the channel.
        threshold_shear_stress : None or float (Pa)
            The threshold shear stress in the equation for erosion rate. If
            None, implies that *set_threshold_from_Dchar* is True, and this
            parameter will get set from the Dchar value and critical Shields
            number.
        Dchar :None, float, array, or field name (m)
            The characteristic grain size on the bed, that controls the
            relationship between critical Shields number and critical shear
            stress. If None, implies that *set_Dchar_from_threshold* is True,
            and this parameter will get set from the threshold_shear_stress
            value and critical Shields number.
        set_threshold_from_Dchar : bool
            If True (default), threshold_shear_stress will be set based on
            Dchar and threshold_Shields.
        set_Dchar_from_threshold : bool
            If True, Dchar will be set based on threshold_shear_stress and
            threshold_Shields. Default is False.
        threshold_Shields : None or float
            The threshold Shields number. If None, implies that
            *slope_sensitive_threshold* is True.
        slope_sensitive_threshold : bool
            If True, the threshold_Shields will be set according to
            0.15 * S ** 0.25, per Lamb et al., 2008 & Hobley et al., 2011.
        flooded_depths : array or field name (m)
            Depths of flooding at each node, zero where no lake. Note that the
            component will dynamically update this array as it fills nodes
            with sediment (...but does NOT update any other related lake
            fields).
        """
        super().__init__(grid)

        if "flow__receiver_node" in grid.at_node and grid.at_node[
            "flow__receiver_node"
        ].size != grid.size("node"):
            raise NotImplementedError(
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that SedDepEroder is compatible with "
                "route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )
        self._flooded_depths = flooded_depths
        self._pseudoimplicit_repeats = pseudoimplicit_repeats

        self._link_S_with_trailing_blank = np.zeros(grid.number_of_links + 1)
        # ^needs to be filled with values in execution
        self._countactive_links = np.zeros_like(
            self._link_S_with_trailing_blank, dtype=int
        )
        self._countactive_links[:-1] = 1

        self._K_unit_time = K_sp / 31557600.0
        # ^...because we work with dt in seconds
        # set gravity
        self._g = g
        self._rock_density = rock_density
        self._sed_density = sediment_density
        self._fluid_density = fluid_density
        self._relative_weight = (
            (self._sed_density - self._fluid_density) / self._fluid_density * self._g
        )
        # ^to accelerate MPM calcs
        self._rho_g = self._fluid_density * self._g
        self._type = sed_dependency_type
        assert self._type in (
            "generalized_humped",
            "None",
            "linear_decline",
            "almost_parabolic",
        )
        self._Qc = Qc
        assert self._Qc in ("MPM", "power_law")
        self._return_ch_props = return_stream_properties
        if return_stream_properties:
            assert self._Qc == "MPM", (
                "Qc must be 'MPM' to return stream " + "properties"
            )
        if isinstance(runoff_rate, (float, int)):
            self._runoff_rate = float(runoff_rate)
        elif isinstance(runoff_rate, str):
            self._runoff_rate = self._grid.at_node[runoff_rate]
        else:
            self._runoff_rate = np.array(runoff_rate)
            assert runoff_rate.size == self._grid.number_of_nodes

        if self._Qc == "MPM":
            if threshold_shear_stress is not None:
                self._thresh = threshold_shear_stress
                self._set_threshold = True
                # ^flag for sed_flux_dep_incision to see if the threshold was
                # manually set.
                # print("Found a shear stress threshold to use: ", self._thresh)
            else:
                warnings.warn("Found no incision threshold to use.", stacklevel=2)
                self._thresh = 0.0
                self._set_threshold = False
            self._a = a_sp
            self._b = b_sp
            self._c = c_sp

            self._k_Q = k_Q
            self._k_w = k_w
            self._mannings_n = mannings_n
            if mannings_n < 0.0 or mannings_n > 0.2:
                warnings.warn("Manning's n outside it's typical range", stacklevel=2)

            self._diffusivity_power_on_A = 0.9 * self._c * (1.0 - self._b)
            # ^i.e., q/D**(1/6)

            self._override_threshold = set_threshold_from_Dchar
            self._override_Dchar = set_Dchar_from_threshold
            if self._override_threshold:
                assert self._set_threshold is False, (
                    "If set_threshold_from_Dchar, threshold_Shields must be "
                    + "set to None"
                )
                assert self._override_Dchar is False
            if self._override_Dchar:
                assert self._override_threshold is False

            self._shields_crit = threshold_Shields
            self._lamb_flag = slope_sensitive_threshold
            if self._lamb_flag:
                assert self._shields_crit is None, (
                    "If slope_sensitive_threshold, threshold_Shields must "
                    + "be set to None"
                )

        elif self._Qc == "power_law":
            self._m = m_sp
            self._n = n_sp
            self._Kt = K_t / 31557600.0  # in sec
            self._mt = m_t
            self._nt = n_t

        # now conditional inputs
        if self._type == "generalized_humped":
            self._kappa = kappa_hump
            self._nu = nu_hump
            self._phi = phi_hump
            self._c = c_hump

        if self._Qc == "MPM":
            if Dchar is not None:
                if isinstance(Dchar, (int, float)):
                    self._Dchar_in = float(Dchar)
                elif isinstance(Dchar, str):
                    self._Dchar_in = self._grid.at_node[Dchar]
                else:
                    self._Dchar_in = np.array(Dchar)
                    assert self._Dchar_in.size == self._grid.number_of_nodes
                assert (
                    not self._override_Dchar
                ), "If set_Dchar_from_threshold, Dchar must be set to None"
            else:
                assert self._override_Dchar
                # remember the threshold getting set is already tau**a
                if not self._lamb_flag:
                    self._Dchar_in = (
                        self._thresh
                        / self._g
                        / (self._sed_density - self._fluid_density)
                        / self._shields_crit
                    )
                else:
                    self._Dchar_in = None
            self._C_MPM = C_MPM

            if self._override_threshold:
                # print("Overriding any supplied threshold...")
                try:
                    self._thresh = (
                        self._shields_crit
                        * self._g
                        * (self._sed_density - self._fluid_density)
                        * self._Dchar_in
                    )
                except AttributeError:
                    self._thresh = (
                        self._shields_crit
                        * self._g
                        * (self._sed_density - self._fluid_density)
                        * Dchar
                    )

            # new 11/12/14
            self._point6onelessb = 0.6 * (1.0 - self._b)
            self._shear_stress_prefactor = (
                self._fluid_density * self._g * (self._mannings_n / self._k_w) ** 0.6
            )

            if self._set_threshold is False or self._override_threshold:
                try:
                    self._shields_prefactor_to_shear = (
                        (self._sed_density - self._fluid_density)
                        * self._g
                        * self._Dchar_in
                    )
                except AttributeError:  # no Dchar
                    self._shields_prefactor_to_shear_noDchar = (
                        self._sed_density - self._fluid_density
                    ) * self._g

            twothirds = 2.0 / 3.0
            self._Qs_prefactor = (
                4.0
                * self._C_MPM**twothirds
                * self._fluid_density**twothirds
                / (self._sed_density - self._fluid_density) ** twothirds
                * self._g ** (twothirds / 2.0)
                * mannings_n**0.6
                * self._k_w ** (1.0 / 15.0)
                * self._k_Q ** (0.6 + self._b / 15.0)
                / self._sed_density**twothirds
            )
            self._Qs_thresh_prefactor = (
                4.0
                * (
                    self._C_MPM
                    * self._k_w
                    * self._k_Q**self._b
                    / self._fluid_density**0.5
                    / (self._sed_density - self._fluid_density)
                    / self._g
                    / self._sed_density
                )
                ** twothirds
            )
            # both these are divided by sed density to give a vol flux
            self._Qs_power_onA = self._c * (0.6 + self._b / 15.0)
            self._Qs_power_onAthresh = twothirds * self._b * self._c

        self._cell_areas = np.empty(grid.number_of_nodes)
        self._cell_areas.fill(np.mean(grid.area_of_cell))
        self._cell_areas[grid.node_at_cell] = grid.area_of_cell

        # set up the necessary fields:
        self.initialize_output_fields()
        if self._return_ch_props:
            self.initialize_optional_output_fields()

    def get_sed_flux_function(self, rel_sed_flux):
        """Get the sediment flux function.

        Parameters
        ----------
        rel_sed_flux
        """
        if self._type == "generalized_humped":
            """Returns K*f(qs,qc)"""
            sed_flux_fn = (
                self._kappa
                * (rel_sed_flux**self._nu + self._c)
                * np.exp(-self._phi * rel_sed_flux)
            )
        elif self._type == "linear_decline":
            sed_flux_fn = 1.0 - rel_sed_flux
        elif self._type == "parabolic":
            raise MissingKeyError(
                "Pure parabolic (where intersect at zero flux is exactly "
                + "zero) is currently not supported, sorry. Try "
                + "almost_parabolic instead?"
            )
            sed_flux_fn = 1.0 - 4.0 * (rel_sed_flux - 0.5) ** 2.0
        elif self._type == "almost_parabolic":
            sed_flux_fn = np.where(
                rel_sed_flux > 0.1,
                1.0 - 4.0 * (rel_sed_flux - 0.5) ** 2.0,
                2.6 * rel_sed_flux + 0.1,
            )
        elif self._type == "None":
            sed_flux_fn = 1.0
        else:
            raise MissingKeyError(
                "Provided sed flux sensitivity type in input file was not "
                + "recognised!"
            )
        return sed_flux_fn

    def get_sed_flux_function_pseudoimplicit(
        self, sed_in, trans_cap_vol_out, prefactor_for_volume, prefactor_for_dz
    ):
        """Get the pseudoimplicit sediment flux function.

        Parameters
        ----------
        sed_in
        trans_cap_vol_out
        prefactor_for_volume
        prefactor_for_dz
        """
        rel_sed_flux_in = sed_in / trans_cap_vol_out
        rel_sed_flux = rel_sed_flux_in

        if self._type == "generalized_humped":
            """Returns K*f(qs,qc)"""

            def sed_flux_fn_gen(rel_sed_flux_in):
                return (
                    self._kappa
                    * (rel_sed_flux_in**self._nu + self._c)
                    * np.exp(-self._phi * rel_sed_flux_in)
                )

        elif self._type == "linear_decline":

            def sed_flux_fn_gen(rel_sed_flux_in):
                return 1.0 - rel_sed_flux_in

        elif self._type == "parabolic":
            raise MissingKeyError(
                "Pure parabolic (where intersect at zero flux is exactly "
                + "zero) is currently not supported, sorry. Try "
                + "almost_parabolic instead?"
            )

            def sed_flux_fn_gen(rel_sed_flux_in):
                return 1.0 - 4.0 * (rel_sed_flux_in - 0.5) ** 2.0

        elif self._type == "almost_parabolic":

            def sed_flux_fn_gen(rel_sed_flux_in):
                return np.where(
                    rel_sed_flux_in > 0.1,
                    1.0 - 4.0 * (rel_sed_flux_in - 0.5) ** 2.0,
                    2.6 * rel_sed_flux_in + 0.1,
                )

        elif self._type == "None":

            def sed_flux_fn_gen(rel_sed_flux_in):
                return 1.0

        else:
            raise MissingKeyError(
                "Provided sed flux sensitivity type in input file was not "
                + "recognised!"
            )

        for _ in range(self._pseudoimplicit_repeats):
            sed_flux_fn = sed_flux_fn_gen(rel_sed_flux)
            sed_vol_added = prefactor_for_volume * sed_flux_fn
            rel_sed_flux = rel_sed_flux_in + sed_vol_added / trans_cap_vol_out
            # print rel_sed_flux
            if rel_sed_flux >= 1.0:
                rel_sed_flux = 1.0
                break
            if rel_sed_flux < 0.0:
                rel_sed_flux = 0.0
                break
        last_sed_flux_fn = sed_flux_fn
        sed_flux_fn = sed_flux_fn_gen(rel_sed_flux)
        # this error could alternatively be used to break the loop
        error_in_sed_flux_fn = sed_flux_fn - last_sed_flux_fn
        dz = prefactor_for_dz * sed_flux_fn
        sed_flux_out = rel_sed_flux * trans_cap_vol_out
        return dz, sed_flux_out, rel_sed_flux, error_in_sed_flux_fn

    def run_one_step(self, dt):
        """Run the component across one timestep increment, dt.

        Erosion occurs according to the sediment dependent rules specified
        during initialization. Method is fully equivalent to the :func:`erode`
        method.

        Parameters
        ----------
        dt : float (years, only!)
            Timestep for which to run the component.
        """

        grid = self._grid
        node_z = grid.at_node["topographic__elevation"]
        node_A = grid.at_node["drainage_area"]
        flow_receiver = grid.at_node["flow__receiver_node"]
        s_in = grid.at_node["flow__upstream_node_order"]
        node_S = grid.at_node["topographic__steepest_slope"]

        if isinstance(self._flooded_depths, str):
            flooded_depths = grid.at_node[self._flooded_depths]
            # also need a map of initial flooded conds:
            flooded_nodes = flooded_depths > 0.0
        elif isinstance(self._flooded_depths, np.ndarray):
            assert self._flooded_depths.size == self._grid.number_of_nodes
            flooded_nodes = self._flooded_depths > 0.0
            # need an *updateable* record of the pit depths
        else:
            # if None, handle in loop
            flooded_nodes = None
        steepest_link = "flow__link_to_receiver_node"
        link_length = np.empty(grid.number_of_nodes, dtype=float)
        link_length.fill(np.nan)
        draining_nodes = np.not_equal(grid.at_node[steepest_link], self._grid.BAD_INDEX)
        core_draining_nodes = np.intersect1d(
            np.where(draining_nodes)[0], grid.core_nodes, assume_unique=True
        )
        link_length[core_draining_nodes] = grid.length_of_d8[
            grid.at_node[steepest_link][core_draining_nodes]
        ]

        if self._Qc == "MPM":
            if self._Dchar_in is not None:
                self._Dchar = self._Dchar_in
            else:
                assert not self._set_threshold, (
                    "Something is seriously wrong with your model " + "initialization."
                )
                assert self._override_threshold, (
                    "You need to confirm to the module you intend it to "
                    + "internally calculate a shear stress threshold, "
                    + "with set_threshold_from_Dchar in the input file."
                )
                # we need to adjust the thresholds for the Shields number
                # & gs dynamically:
                variable_thresh = (
                    self._shields_crit
                    * self._g
                    * (self._sed_density - self._fluid_density)
                    * self._Dchar
                )
            if self._lamb_flag:
                variable_shields_crit = 0.15 * node_S**0.25
                try:
                    variable_thresh = (
                        variable_shields_crit * self._shields_prefactor_to_shear
                    )
                except AttributeError:
                    variable_thresh = (
                        variable_shields_crit
                        * self._shields_prefactor_to_shear_noDchar
                        * self._Dchar
                    )

            node_Q = self._k_Q * self._runoff_rate * node_A**self._c
            shear_stress_prefactor_timesAparts = (
                self._shear_stress_prefactor * node_Q**self._point6onelessb
            )
            try:
                transport_capacities_thresh = (
                    self._thresh
                    * self._Qs_thresh_prefactor
                    * self._runoff_rate ** (0.66667 * self._b)
                    * node_A**self._Qs_power_onAthresh
                )
            except AttributeError:
                transport_capacities_thresh = (
                    variable_thresh
                    * self._Qs_thresh_prefactor
                    * self._runoff_rate ** (0.66667 * self._b)
                    * node_A**self._Qs_power_onAthresh
                )

            transport_capacity_prefactor_withA = (
                self._Qs_prefactor
                * self._runoff_rate ** (0.6 + self._b / 15.0)
                * node_A**self._Qs_power_onA
            )

            internal_t = 0.0
            break_flag = False
            dt_secs = dt * 31557600.0
            counter = 0
            rel_sed_flux = np.empty_like(node_Q)
            # excess_vol_overhead = 0.

            while 1:
                # ^use the break flag, to improve computational efficiency for
                # runs which are very stable
                # we assume the drainage structure is forbidden to change
                # during the whole dt
                # note slopes will be *negative* at pits
                # track how many loops we perform:
                counter += 1
                downward_slopes = node_S.clip(0.0)
                # this removes the tendency to transfer material against
                # gradient, including in any lake depressions
                # we DON'T immediately zero trp capacity in the lake.
                # positive_slopes = np.greater(downward_slopes, 0.)
                slopes_tothe07 = downward_slopes**0.7
                transport_capacities_S = (
                    transport_capacity_prefactor_withA * slopes_tothe07
                )
                trp_diff = (transport_capacities_S - transport_capacities_thresh).clip(
                    0.0
                )
                transport_capacities = np.sqrt(trp_diff * trp_diff * trp_diff)
                shear_stress = shear_stress_prefactor_timesAparts * slopes_tothe07
                shear_tothe_a = shear_stress**self._a

                dt_this_step = dt_secs - internal_t
                # ^timestep adjustment is made AFTER the dz calc
                node_vol_capacities = transport_capacities * dt_this_step

                sed_into_node = np.zeros(grid.number_of_nodes, dtype=float)
                dz = np.zeros(grid.number_of_nodes, dtype=float)
                cell_areas = self._cell_areas
                try:
                    raise NameError
                    # ^tripped out deliberately for now; doesn't appear to
                    # accelerate much
                    weave.inline(
                        self._routing_code,
                        [
                            "len_s_in",
                            "sed_into_node",
                            "transport_capacities",
                            "dz",
                            "cell_areas",
                            "dt_this_step",
                            "flow__receiver_node",
                        ],
                    )
                except NameError:
                    for i in s_in[::-1]:  # work downstream
                        cell_area = cell_areas[i]
                        if flooded_nodes is not None:
                            flood_depth = flooded_depths[i]
                        else:
                            flood_depth = 0.0
                        sed_flux_into_this_node = sed_into_node[i]
                        node_capacity = transport_capacities[i]
                        # ^we work in volume flux, not volume per se here
                        node_vol_capacity = node_vol_capacities[i]
                        if flood_depth > 0.0:
                            node_vol_capacity = 0.0
                            # requires special case handling - as much sed as
                            # possible is dumped here, then the remainder
                            # passed on
                        if sed_flux_into_this_node < node_vol_capacity:
                            # ^note incision is forbidden at capacity
                            # flooded nodes never enter this branch
                            # #implementing the pseudoimplicit method:
                            try:
                                thresh = variable_thresh
                            except NameError:  # it doesn't exist
                                thresh = self._thresh
                            dz_prefactor = (
                                self._K_unit_time
                                * dt_this_step
                                * (shear_tothe_a[i] - thresh).clip(0.0)
                            )
                            vol_prefactor = dz_prefactor * cell_area
                            (
                                dz_here,
                                sed_flux_out,
                                rel_sed_flux_here,
                                error_in_sed_flux,
                            ) = self.get_sed_flux_function_pseudoimplicit(
                                sed_flux_into_this_node,
                                node_vol_capacity,
                                vol_prefactor,
                                dz_prefactor,
                            )
                            # note now dz_here may never create more sed than
                            # the out can transport...
                            assert sed_flux_out <= node_vol_capacity, (
                                "failed at node "
                                + str(s_in.size - i)
                                + " with rel sed flux "
                                + str(sed_flux_out / node_capacity)
                            )
                            rel_sed_flux[i] = rel_sed_flux_here
                            vol_pass = sed_flux_out
                        else:
                            rel_sed_flux[i] = 1.0
                            vol_dropped = sed_flux_into_this_node - node_vol_capacity
                            dz_here = -vol_dropped / cell_area
                            # with the pits, we aim to inhibit incision, but
                            # depo is OK. We have already zero'd any adverse
                            # grads, so sed can make it to the bottom of the
                            # pit but no further in a single step, which seems
                            # raeasonable. Pit should fill.
                            if flood_depth <= 0.0:
                                vol_pass = node_vol_capacity
                            else:
                                height_excess = -dz_here - flood_depth
                                # ...above water level
                                if height_excess <= 0.0:
                                    vol_pass = 0.0
                                    # dz_here is already correct
                                    flooded_depths[i] += dz_here
                                else:
                                    dz_here = -flood_depth
                                    vol_pass = height_excess * cell_area
                                    # ^bit cheeky?
                                    flooded_depths[i] = 0.0
                                    # note we must update flooded depths
                                    # transiently to conserve mass
                            # do we need to retain a small downhill slope?
                            # ...don't think so. Will resolve itself on next
                            # timestep.

                        dz[i] -= dz_here
                        sed_into_node[flow_receiver[i]] += vol_pass

                break_flag = True

                node_z[grid.core_nodes] += dz[grid.core_nodes]

                if break_flag:
                    break
                # do we need to reroute the flow/recalc the slopes here?
                # -> NO, slope is such a minor component of Diff we'll be OK
                # BUT could be important not for the stability, but for the
                # actual calc. So YES.
                node_S = np.zeros_like(node_S)
                node_S[core_draining_nodes] = (node_z - node_z[flow_receiver])[
                    core_draining_nodes
                ] / link_length[core_draining_nodes]
                internal_t += dt_this_step  # still in seconds, remember

        elif self._Qc == "power_law":
            transport_capacity_prefactor_withA = self._Kt * node_A**self._mt
            erosion_prefactor_withA = self._K_unit_time * node_A**self._m
            # ^doesn't include S**n*f(Qc/Qc)
            internal_t = 0.0
            break_flag = False
            dt_secs = dt * 31557600.0
            counter = 0
            rel_sed_flux = np.empty_like(node_A)
            while 1:
                counter += 1
                # print counter
                downward_slopes = node_S.clip(0.0)
                # positive_slopes = np.greater(downward_slopes, 0.)
                slopes_tothen = downward_slopes**self._n
                slopes_tothent = downward_slopes**self._nt
                transport_capacities = (
                    transport_capacity_prefactor_withA * slopes_tothent
                )
                erosion_prefactor_withS = (
                    erosion_prefactor_withA * slopes_tothen
                )  # no time, no fqs
                # shear_tothe_a = shear_stress**self._a

                dt_this_step = dt_secs - internal_t
                # ^timestep adjustment is made AFTER the dz calc
                node_vol_capacities = transport_capacities * dt_this_step

                sed_into_node = np.zeros(grid.number_of_nodes, dtype=float)
                dz = np.zeros(grid.number_of_nodes, dtype=float)
                cell_areas = self._cell_areas
                for i in s_in[::-1]:  # work downstream
                    cell_area = cell_areas[i]
                    if flooded_nodes is not None:
                        flood_depth = flooded_depths[i]
                    else:
                        flood_depth = 0.0
                    sed_flux_into_this_node = sed_into_node[i]
                    node_capacity = transport_capacities[i]
                    # ^we work in volume flux, not volume per se here
                    node_vol_capacity = node_vol_capacities[i]
                    if flood_depth > 0.0:
                        node_vol_capacity = 0.0
                    if sed_flux_into_this_node < node_vol_capacity:
                        # ^note incision is forbidden at capacity
                        dz_prefactor = dt_this_step * erosion_prefactor_withS[i]
                        vol_prefactor = dz_prefactor * cell_area
                        (
                            dz_here,
                            sed_flux_out,
                            rel_sed_flux_here,
                            error_in_sed_flux,
                        ) = self.get_sed_flux_function_pseudoimplicit(
                            sed_flux_into_this_node,
                            node_vol_capacity,
                            vol_prefactor,
                            dz_prefactor,
                        )
                        # note now dz_here may never create more sed than the
                        # out can transport...
                        assert sed_flux_out <= node_vol_capacity, (
                            "failed at node "
                            + str(s_in.size - i)
                            + " with rel sed flux "
                            + str(sed_flux_out / node_capacity)
                        )
                        rel_sed_flux[i] = rel_sed_flux_here
                        vol_pass = sed_flux_out
                    else:
                        rel_sed_flux[i] = 1.0
                        vol_dropped = sed_flux_into_this_node - node_vol_capacity
                        dz_here = -vol_dropped / cell_area
                        try:
                            isflooded = flooded_nodes[i]
                        except TypeError:  # was None
                            isflooded = False
                        if flood_depth <= 0.0 and not isflooded:
                            vol_pass = node_vol_capacity
                            # we want flooded nodes which have already been
                            # filled to enter the else statement
                        else:
                            height_excess = -dz_here - flood_depth
                            # ...above water level
                            if height_excess <= 0.0:
                                vol_pass = 0.0
                                # dz_here is already correct
                                flooded_depths[i] += dz_here
                            else:
                                dz_here = -flood_depth
                                vol_pass = height_excess * cell_area
                                # ^bit cheeky?
                                flooded_depths[i] = 0.0

                    dz[i] -= dz_here
                    sed_into_node[flow_receiver[i]] += vol_pass
                break_flag = True

                node_z[grid.core_nodes] += dz[grid.core_nodes]

                if break_flag:
                    break
                # do we need to reroute the flow/recalc the slopes here?
                # -> NO, slope is such a minor component of Diff we'll be OK
                # BUT could be important not for the stability, but for the
                # actual calc. So YES.
                node_S = np.zeros_like(node_S)
                # print link_length[core_draining_nodes]
                node_S[core_draining_nodes] = (node_z - node_z[flow_receiver])[
                    core_draining_nodes
                ] / link_length[core_draining_nodes]
                internal_t += dt_this_step  # still in seconds, remember

        if self._return_ch_props:
            # add the channel property field entries,
            # 'channel__width', 'channel__depth', and 'channel__discharge'
            W = self._k_w * node_Q**self._b
            H = shear_stress / self._rho_g / node_S  # ...sneaky!
            grid.at_node["channel__width"][:] = W
            grid.at_node["channel__depth"][:] = H
            grid.at_node["channel__discharge"][:] = node_Q
            grid.at_node["channel__bed_shear_stress"][:] = shear_stress

        grid.at_node["channel_sediment__volumetric_transport_capacity"][
            :
        ] = transport_capacities
        grid.at_node["channel_sediment__volumetric_flux"][:] = sed_into_node
        grid.at_node["channel_sediment__relative_flux"][:] = rel_sed_flux
        # elevs set automatically to the name used in the function call.
        self._iterations_in_dt = counter

        return grid, grid.at_node["topographic__elevation"]

    @property
    @make_return_array_immutable
    def characteristic_grainsize(self):
        """Return the characteristic grain size used by the component.

        Particularly useful if the set_Dchar_from_threshold flag was True
        at initialization.

        Returns
        -------
        Dchar : float or array

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator, SedDepEroder
        >>> mg1 = RasterModelGrid((3, 4))
        >>> z1 = mg1.add_zeros("node", "topographic__elevation")
        >>> fa1 = FlowAccumulator(mg1)
        >>> thresh_shields = np.arange(1, mg1.number_of_nodes + 1, dtype=float)
        >>> thresh_shields /= 100.0
        >>> sde1 = SedDepEroder(
        ...     mg1,
        ...     threshold_shear_stress=100.0,
        ...     Qc="MPM",
        ...     Dchar=None,
        ...     set_threshold_from_Dchar=False,
        ...     set_Dchar_from_threshold=True,
        ...     threshold_Shields=thresh_shields,
        ...     g=9.81,
        ... )
        >>> sde1.characteristic_grainsize.reshape(mg1.shape)
        array([[0.59962823, 0.29981412, 0.19987608, 0.14990706],
               [0.11992565, 0.09993804, 0.08566118, 0.07495353],
               [0.06662536, 0.05996282, 0.05451166, 0.04996902]])

        >>> mg2 = RasterModelGrid((3, 4))
        >>> z2 = mg2.add_zeros("node", "topographic__elevation")
        >>> fa2 = FlowAccumulator(mg2)
        >>> sde2 = SedDepEroder(
        ...     mg2,
        ...     threshold_shear_stress=100.0,
        ...     Qc="MPM",
        ...     Dchar=None,
        ...     set_threshold_from_Dchar=False,
        ...     set_Dchar_from_threshold=True,
        ...     threshold_Shields=None,
        ...     slope_sensitive_threshold=True,
        ...     g=9.81,
        ... )
        >>> S = mg2.at_node["topographic__steepest_slope"]
        >>> S[:] = 0.05  # thresh = 100 Pa @ 5pc slope
        >>> sde2.characteristic_grainsize.reshape(mg2.shape)
        array([[0.08453729, 0.08453729, 0.08453729, 0.08453729],
               [0.08453729, 0.08453729, 0.08453729, 0.08453729],
               [0.08453729, 0.08453729, 0.08453729, 0.08453729]])
        """
        # Dchar is None means self._lamb_flag, Dchar is spatially variable,
        # and not calculated until the main loop
        assert (
            self._Qc == "MPM"
        ), "Characteristic grainsize is only calculated if Qc == 'MPM'"
        if self._Dchar_in is not None:
            return self._Dchar_in
        else:
            taustarcrit = (
                0.15 * self._grid.at_node["topographic__steepest_slope"] ** 0.25
            )
            Dchar = self._thresh / (
                self._g * (self._sed_density - self._fluid_density) * taustarcrit
            )
            return Dchar
