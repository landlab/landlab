#!/usr/bin/env python3
"""
Model bedrock incision and gravel transport and abrasion in a network of rivers.

@author: gtucker
"""
import numpy as np
from landlab import Component
from landlab import HexModelGrid
from landlab.grid.diagonals import DiagonalsMixIn
from cfuncs import _calc_sediment_influx
from cfuncs import _calc_sediment_rate_of_change
# from cfuncs import _estimate_max_time_step_size_ext
from cfuncs import _calc_sed_abrs_rate
from cfuncs import _calc_bedrock_abrs_rate
from cfuncs import _get_classes_fractions
from cfuncs import _calc_pluck_rate

_DT_MAX = 1.0e-2
_ONE_SIXTH = 1.0 / 6.0
_SEVEN_SIXTHS = 7.0 / 6.0
_D8_CHAN_LENGTH_FACTOR = (
    1.0  # 0.5 * (1.0 + 2.0**0.5)  # D8 raster: average of straight and diagonal
)
_SEC_PER_YEAR = 365.25 * 24.0 * 3600.0
_SQUARED_SECS_IN_A_YEAR = _SEC_PER_YEAR**2
_EARTH_GRAV = 9.81
_FIVE_BY_THREE = 5.0 / 3.0
_WICKERT_ROUGHNESS_FACTOR = 0.17
_POINT_FIVE = 0.5
_THREE_BY_TWO = 3.0 / 2.0


def _calc_chan_width_fixed_width(coeff, expt, discharge, out=None):
    """Calculate channel width using empirical formula.

    b = coeff * Q**expt

    Parameters
    ----------
    coeff : float
        Empirical coefficient
    expt : float
        Empirical exponent
    discharge : array
        Discharge
    out : array, optional
        Array to store output

    Returns
    -------
    array
        Channel width

    Examples
    --------
    >>> _calc_chan_width_fixed_width(1.0, 0.5, np.array([1.0]))
    array([1.])
    >>> _calc_chan_width_fixed_width(6e-8, 0.5, np.array([100.0]))
    array([6.e-07])
    """

    if out is None:
        out = np.empty_like(discharge)
    out[:] = coeff * discharge**expt
    return out


def _calc_near_threshold_width(
    median_size, tau_star_c_median, discharge, slope, g_star, SG, epsilon, out=None
):
    """Utility function that calculates and returns channel width implied by
    discharge, slope, and grain diameter, using Wickert & Schildgen (2019)
    equation 16.

    Parameters
    ----------
    median_size : array (float)
        Median grain size at node
    tau_star_c_median : array (float)
        Normalized critical shear stress for the median size
    discharge : array (float)
        Discharge at node
    slope : array (float)
        Slope at node
    g_star : float
        Gravity coeffcient that takes into conversion of sec to years
    SG : float
        Specific gravity
    epsilon : float
        Near-threshold coefficient
    out : array, optional
        Array to store output

    Returns
    -------
    array
        Channel width

    Examples
    --------
    >>> g_star = 9769603575225600
    >>> SG = 1.65
    >>> epsilon = 0.2
    >>> tau_star_c_median = 0.045
    >>> discharge = np.array([1000000.0])
    >>> slope = 0.1
    >>> median_size = np.array([0.01])
    >>> epsilon = 0.2
    >>> _calc_near_threshold_width(median_size, tau_star_c_median, discharge, slope, g_star, SG, epsilon, out=None)
    array([6.59246163])
    """
    if out is None:
        out = np.empty_like(discharge)

    out[:] = (
        _WICKERT_ROUGHNESS_FACTOR
        * g_star ** (-0.5)
        * (SG) ** (-_FIVE_BY_THREE)
        * (1 + epsilon) ** (-_FIVE_BY_THREE)
        * tau_star_c_median ** (-_FIVE_BY_THREE)
        * np.divide(
            discharge * slope ** (_SEVEN_SIXTHS),
            median_size ** (1.5),
            where=median_size != 0,
            out=np.zeros_like(discharge),
        )
    )
    return out


def _calc_shear_stress_coef(rho_w, mannings_n, g=_EARTH_GRAV):
    """
    Calculate prefactor for:
    tau = rho*g*n^(3/5)(Q...)

    Assuming discHarge is given in m^3/year

    so rho*g*n^(3/5) is the prefactor.

    Examples
    --------
    >>> np.round(int(_calc_shear_stress_coef(1000, 0.05)*1000))
    np.int64(51)
    """

    return rho_w * g * mannings_n ** (0.6) * (1 / _SEC_PER_YEAR) ** 0.6


def _calc_shear_stress(shear_stress_coef, discharge, width, slope, out=None):
    """Calculate shear stress using Manning's equation.

    Parameters
    ----------
    shear_stress_coef : float
        Shear stress coefficient
    discharge : array (float)
        Discharge at node
    width : array (float)
        Channel width at node
    slope : array (float)
        Slope at node
    out : array, optional
        Array to store output

    Returns
    -------
    array
        Shear stress at node

     Examples
    --------
    >>> rho_w = 1000
    >>> mannings_n = 0.05
    >>> shear_stress_coeff  = _calc_shear_stress_coef(rho_w, mannings_n, g=_EARTH_GRAV)
    >>> discharge = np.array([1000000.0])
    >>> width = np.array([10])
    >>> slope = 0.01
    >>> round(_calc_shear_stress(shear_stress_coeff, discharge, width, slope)[0])
    2
    """

    if out is None:
        out = np.empty_like(discharge)
    out[:] = (
        shear_stress_coef
        * np.divide(discharge, width, where=width > 0, out=np.zeros_like(out)) ** (0.6)
        * slope ** (0.7)
    )
    return out


def _tau_crit_fixed_width(
    rho_sed, rho_w, grain_size, tau_star_crit, g=_EARTH_GRAV, out=None
):
    """Calculate transport rate using Meyer-Peter Mueller.
    Parameters
    ----------
    rho_sed : float
        Density of sediment
    rho_w : float
        Density of water
    grain_size : float
        Grain size
    tau_star_crit : array
        Dimensionless critical shear stress
    g : float, optional
        Acceleration due to gravity
    out : array, optional
        Array to store output

    >>> tau_star_crit = np.array([0.045])
    >>> rho_sed = 2650
    >>> rho_w = 1000
    >>> grain_size = np.array([0.1])
    >>> g = _EARTH_GRAV
    >>> _tau_crit_fixed_width(rho_sed, rho_w, grain_size, tau_star_crit, g)
    array([72.83925])
    """

    if out is None:
        out = np.empty_like(tau_star_crit)

    out[:] = tau_star_crit * (rho_sed - rho_w) * g * grain_size
    return out


# def _calc_plucking_rate_fixed_width(coeff_pluck, width_chan, width_val, tau, tau_crit, out=None):
#     """ see Overleaf document for description; SMM """
#
#     if out is None:
#         out = np.empty_like(tau_crit)
#
#     ratio = np.divide(width_chan,
#                       width_val,
#                       where=width_val > 0,
#                       out=np.zeros_like(width_chan))
#
#     out[:] = coeff_pluck * (tau - tau_crit) ** 1.5 * ratio
#     return out


class ExtendedGravelBedrockEroder(Component):
    """Drainage network evolution of rivers with gravel alluvium overlying bedrock.

    Model drainage network evolution for a network of rivers that have
    a layer of gravel alluvium overlying bedrock.

    :class:`~.GravelBedrockEroder` is designed to operate together with a flow-routing
    component such as :class:`~.FlowAccumulator`, so that each grid node has
    a defined flow direction toward one of its neighbor nodes. Each core node
    is assumed to contain one outgoing fluvial channel, and (depending on
    the drainage structure) zero, one, or more incoming channels. These channels are
    treated as effectively sub-grid-scale features that are embedded in valleys
    that have a width of one grid cell.

    As with the :class:`~.GravelRiverTransporter` component, the rate of gravel
    transport out of a given node is calculated as the product of bankfull discharge,
    channel gradient (to the 7/6 power), a dimensionless transport coefficient, and
    an intermittency factor that represents the fraction of time that bankfull
    flow occurs. The derivation of the transport law is given by Wickert &
    Schildgen (2019), and it derives from the assumption that channels are
    gravel-bedded and that they "instantaneously" adjust their width such that
    bankfull bed shear stress is just slightly higher than the threshold for
    grain motion. The substrate is assumed to consist entirely of gravel-size
    material with a given bulk porosity. The component calculates the loss of
    gravel-sized material to abrasion (i.e., conversion to finer sediment, which
    is not explicitly tracked) as a function of the volumetric transport rate,
    an abrasion coefficient with units of inverse length, and the local transport
    distance (for example, if a grid node is carrying a gravel load ``Qs`` to a
    neighboring node ``dx`` meters downstream, the rate of gravel loss in volume per
    time per area at the node will be ``beta * Qs * dx``, where ``beta`` is the abrasion
    coefficient).

    Sediment mass conservation is calculated across each entire
    grid cell. For example, if a cell has surface area ``A``, a total volume influx
    ``Qin``, and downstream transport rate ``Qs``, the resulting rate of change of
    alluvium thickness will be ``(Qin - Qs / (A * (1 - phi))``, plus gravel produced by
    plucking erosion of bedrock (``phi`` is porosity).

    Bedrock is eroded by a combination of abrasion and plucking. Abrasion per unit
    channel length is calculated as the product of volumetric sediment discharge
    and an abrasion coefficient. Sediment produced by abrasion is assumed to
    go into wash load that is removed from the model domain. Plucking is calculated
    using a discharge-slope expression, and a user-defined fraction of plucked
    material is added to the coarse alluvium.

    Parameters
    ----------
    grid : ModelGrid
        A Landlab model grid object
    intermittency_factor : float (default 0.01)
        Fraction of time that bankfull flow occurs
    transport_coefficient : float (default 0.041)
        if near-threshold:
            Dimensionless transport efficiency factor; see Wickert & Schildgen 2019
        if empirical width:
            See Meyer-Peter Mueller in Wong & Parker 2006
    abrasion_coefficient : float (default 0.0 1/m) *DEPRECATED*
        Abrasion coefficient with units of inverse length
    sediment_porosity : float (default 0.35)
        Bulk porosity of bed sediment
    depth_decay_scale : float (default 1.0)
        Scale for depth decay in bedrock exposure function
    plucking_coefficient : float or (n_core_nodes,) array of float (default 1.0e-4 1/m)
        Rate coefficient for bedrock erosion by plucking
        if near-threshold:
            See Gabel et al. 2024
        if empirical width:
            See [insert description here, see Overleaf doc for details; SMM]
    self._n_classes : int (default 1)
        Number of sediment abradability classes
    init_thickness_per_class : float or (n_core_nodes,) array of float (default 1 / n-classes)
        Starting thickness for each sediment fraction
    abrasion_coefficients : iterable containing floats (default 0.0 1/m)
        Abrasion coefficients; should be same length as number of sed classes
    bedrock_abrasion_coefficients : float
        Abrasion coefficient for bedrock
    coarse_fractions_from_plucking : float or (n_core_nodes,) array of float (default 1.0)
        Fraction(s) of plucked material that becomes part of gravel sediment load
    rock_abrasion_index : int (default 0)
        If multiple classes, specifies which contains the abrasion
        coefficient for bedrock
    --The following are necessary only if using the empirical width calculation--
    tau_star_crit : float (default 0.045)
        Dimensionless critical shear stress;
    grav_accel : float (default to Earth, 9.81)
        Acceleration due to gravity;
    width_coeff : float (default 2?)
        Empirical coefficient for channel width;
        NOTE: units, if exponent is 1/2: 1/(m^0.5*s^0.5) -- need to make this conversion to years, not seconds
    width_expt : float (default )
        Empirical exponent for channel width;
    rho_w : float (default 1000)
        Density of water;
    rho_sed : float (default 2650)
        Density of sediment;
    D_50 : float (default 0.01)
        Median grain size;


    Notes
    -----
    The doctest below demonstrates approximate equilibrium between uplift, transport,
    and sediment abrasion in a case with effectively unlimited sediment. The
    analytical solution is:

    sediment input by uplift = sediment outflux + sediment loss to abrasion

    In math,

    U A = kq I Q S^(7/6) + 0.5 b Qs dx

    S = (U A / (kq I Q (1 + 0.5 b dx))) ^ 6/7

    S = (1.0e-4 1e6 / (0.041 0.01 10.0 1e6 (1 + 0.5 0.0005 1000.0)))^(6 / 7)

    ~ 0.0342

    The sediment abrasion rate should be, in volume per time, as follows:

    Qsout + abrasion loss rate = generation rate

    Qsout + 0.5 Qsout b dx = U dx^2

    Qsout = U dx^2 / (1 + 0.5 b dx)

    Qsout = 0.0001 1e6 / (1 + 0.5 0.0005 1e3)

    Qsout = 1e2 / 1.25 = 80

    Elevation of the single core node = 0.0342 x 1,000 m ~ 34.2 m.
    However, because of VERY long time steps, the post-erosion elevation is 1 m
    lower, at 33.2 m (it will be uplifted by a meter at the start of each step).

    # Examples
    # --------
    # >>> from landlab import RasterModelGrid
    # >>> from landlab.components import FlowAccumulator
    # >>> grid = RasterModelGrid((3, 3), xy_spacing=1000.0)
    # >>> elev = grid.add_zeros("topographic__elevation", at="node")
    # >>> elev[4] = 1.0
    # >>> sed = grid.add_zeros("soil__depth", at="node")
    # >>> sed[4] = 1.0e6
    # >>> grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    # >>> grid.status_at_node[5] = grid.BC_NODE_IS_FIXED_VALUE
    # >>> fa = FlowAccumulator(grid, runoff_rate=10.0)
    # >>> fa.run_one_step()
    # >>> eroder = ExtendedGravelBedrockEroder(
    # ...     grid, sediment_porosity=0.0, abrasion_coefficients=[0.0005]
    # ... )
    # >>> rock_elev = grid.at_node["bedrock__elevation"]
    # >>> fa.run_one_step()
    # >>> dt = 10000.0
    # >>> for _ in range(200):
    # ...     rock_elev[grid.core_nodes] += 1.0e-4 * dt
    # ...     elev[grid.core_nodes] += 1.0e-4 * dt
    # ...     eroder.run_one_step(dt)
    # ...
    # >>> int(elev[4])
    33
    """

    _name = "ExtendedGravelBedrockEroder"

    _unit_agnostic = True

    _info = {
        "bedload_sediment__rate_of_loss_to_abrasion": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/y",
            "mapping": "node",
            "doc": "Rate of bedload sediment volume loss to abrasion per unit area",
        },
        "bedload_sediment__volume_influx": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/y",
            "mapping": "node",
            "doc": "Volumetric incoming streamwise bedload sediment transport rate",
        },
        "bedload_sediment__volume_outflux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/y",
            "mapping": "node",
            "doc": "Volumetric outgoing streamwise bedload sediment transport rate",
        },
        "bedrock__abrasion_rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/y",
            "mapping": "node",
            "doc": "rate of bedrock lowering by abrasion",
        },
        "bedrock__elevation": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "elevation of the bedrock surface",
        },
        "bedrock__exposure_fraction": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "fractional exposure of bedrock",
        },
        "bedrock__plucking_rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/y",
            "mapping": "node",
            "doc": "rate of bedrock lowering by plucking",
        },
        "bedrock__lowering_rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/y",
            "mapping": "node",
            "doc": "Rate of lowering of bedrock surface",
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
        "sediment__rate_of_change": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/y",
            "mapping": "node",
            "doc": "Time rate of change of sediment thickness",
        },
        "soil__depth": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of soil or weathered bedrock",
        },
        "surface_water__discharge": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m**3/y",
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
        "grains__weight": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "kg",
            "mapping": "node",
            "doc": "",
        },
    }

    def __init__(
        self,
        grid,
        intermittency_factor=0.01,
        transport_coefficient=0.041,
        sediment_porosity=0.35,
        depth_decay_scale=1.0,
        plucking_coefficient=1.0e-4,
        epsilon=0.2,  # from parker 1978
        abrasion_coefficients=0.0,
        bedrock_abrasion_coefficient=0.01,
        fractions_from_plucking=1.0,
        rock_abrasion_index=0,
        rho_sed=2650,
        rho_water=1000,
        fixed_width_flag=1,
        fixed_width_coeff=0.002,
        fixed_width_expt=0.5,
        mannings_n=0.05,
        tau_star_c_median=0.045,
        alpha=-0.68,
    ):

        super().__init__(grid)
        super().initialize_output_fields()

        # Get the number of classes from SoilGrading component
        # In this version the number of classes can represent EITHER
        # lithology OR grain sizes
        if np.ndim(self._grid.at_node["grains__weight"]) > 1:
            self._n_classes = np.shape(self._grid.at_node["grains__weight"])[1]
        else:
            self._n_classes = 1
        # Create 2D arrays for input variables
        # The component will also check the validity of
        # the user-defined variable
        self._abr_coefs = self._create_2D_array_for_input_var(
            abrasion_coefficients, "abrasion_coefficients"
        )
        self._fractions_from_plucking = self._create_2D_array_for_input_var(
            fractions_from_plucking, "fractions_from_plucking"
        )

        self._tau_star_c = self._create_2D_array_for_input_var(0, "tau_star_c")

        # Recognize whether the component deals with lithologies or grain sizes
        self._get_classes_identity()

        # Non-array parameters
        self._trans_coef = transport_coefficient
        self._intermittency_factor = intermittency_factor
        self._sediment_porosity = sediment_porosity
        self._porosity_factor = 1.0 / (1.0 - self._sediment_porosity)
        self._depth_decay_scale = depth_decay_scale
        self._rock_abrasion_index = rock_abrasion_index
        self._rho_sed = rho_sed
        self._rho_water = rho_water
        self._SG = np.divide(self._rho_sed - self._rho_water, self._rho_water)
        self._fixed_width_coeff = fixed_width_coeff
        self._fixed_width_expt = fixed_width_expt
        self._calc_weight_threshold_to_deliv()
        self._shear_stress_coef = _calc_shear_stress_coef(
            rho_w=rho_water, mannings_n=mannings_n
        )
        self._fixed_width_flag = fixed_width_flag
        self._calc_gravity_coefficient_star()
        self._tau_star_c_median = tau_star_c_median
        self._alpha = alpha
        self._epsilon = epsilon
        # Pointers to field
        self._elev = grid.at_node["topographic__elevation"]
        self._sed = grid.at_node["soil__depth"]
        if "bedrock__elevation" in grid.at_node:
            self._bedrock__elevation = grid.at_node["bedrock__elevation"]
        else:
            self._bedrock__elevation = grid.add_zeros(
                "bedrock__elevation", at="node", dtype=float
            )
            self._bedrock__elevation[:] = self._elev - self._sed
        self._discharge = grid.at_node["surface_water__discharge"]
        self._slope = grid.at_node["topographic__steepest_slope"]
        self._receiver_node = grid.at_node["flow__receiver_node"]
        self._receiver_link = grid.at_node["flow__link_to_receiver_node"]
        self._sediment_influx = grid.at_node["bedload_sediment__volume_influx"]
        self._sediment_outflux = grid.at_node["bedload_sediment__volume_outflux"]
        self._dHdt = grid.at_node["sediment__rate_of_change"]
        self._rock_lowering_rate = grid.at_node["bedrock__lowering_rate"]
        self._rock_exposure_fraction = grid.at_node["bedrock__exposure_fraction"]
        self._rock_abrasion_rate = grid.at_node["bedrock__abrasion_rate"]
        self._pluck_rate = grid.at_node["bedrock__plucking_rate"]
        self._setup_length_of_flow_link()

        # 1D arrays in dimensions of n_nodes
        self._implied_width = np.zeros_like(self.grid.nodes.flatten()).astype(float)
        self._channel_width = np.zeros_like(self.grid.nodes.flatten()).astype(float)
        self._plucking_coef = np.zeros_like(self.grid.nodes.flatten()).astype(float)
        self._tau = np.zeros_like(self.grid.nodes.flatten())

        if isinstance(plucking_coefficient, float) or isinstance(
            plucking_coefficient, int
        ):
            self._plucking_coef[:] = plucking_coefficient
        elif isinstance(plucking_coefficient, np.ndarray) or isinstance(
            plucking_coefficient, list
        ):
            try:
                self._plucking_coef[:] = plucking_coefficient
            except ValueError:
                print("Plucking coefficient array does not match the number of nodes")
        else:
            print("Plucking coefficient format is not valid")
        self._br_abr_coef = np.zeros_like(self.grid.nodes.flatten(), dtype=float)
        if isinstance(bedrock_abrasion_coefficient, float) or isinstance(
            bedrock_abrasion_coefficient, int
        ):
            self._br_abr_coef[:] = plucking_coefficient
        elif isinstance(bedrock_abrasion_coefficient, np.ndarray) or isinstance(
            bedrock_abrasion_coefficient, list
        ):
            try:
                self._br_abr_coef[:] = bedrock_abrasion_coefficient
            except ValueError:
                print(
                    "Bedrock abrasion coefficient array does not match the number of nodes"
                )
        else:
            print("Bedrock abrasion coefficient format is not valid")
        # 2D arrays in dimensions of n_nodes x n_classes
        self._thickness_by_class = np.zeros((grid.number_of_nodes, self._n_classes))
        self._sed_influxes = np.zeros((grid.number_of_nodes, self._n_classes))
        self._sed_outfluxes = np.zeros((grid.number_of_nodes, self._n_classes))
        self._sed_abr_rates = np.zeros((grid.number_of_nodes, self._n_classes))
        self._br_abrasion_coef = np.zeros((grid.number_of_nodes, self._n_classes))
        self._dHdt_by_class = np.zeros((grid.number_of_nodes, self._n_classes))
        self._tau_star_c = np.zeros((grid.number_of_nodes, self._n_classes))
        self._excess_stress = np.zeros((grid.number_of_nodes, self._n_classes))
        self._get_sediment_thickness_by_class()

    def _create_2D_array_for_input_var(self, input_var, var_name):
        """ ""
        This procedure receive an input variable that will be used
        by the component across nodes and classes (classes represent lithology/grain size).
        The procedure will verify that the input variable match the number of classes
        and will return a 2D array (with dimensions of n_nodes x n_classes).
        If the input variable is already a 2D array, its match to n_nodes and n_classes
        will be examined

        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components.soil_grading import SoilGrading
        >>> xy_spacing=100
        >>> grid = RasterModelGrid((3, 4), xy_spacing=xy_spacing)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> sed_depth = 1000
        >>> porosity=0.5
        >>> sed_weight = sed_depth * 2650 * (1-porosity) # weight per grid node area
        >>> grains_weight =[sed_weight, sed_weight, sed_weight]
        >>> grain_sizes = [0.001, 0.01, 0.05]
        >>> sg = SoilGrading(grid,
        ...         meansizes=grain_sizes,
        ...         grains_weight=grains_weight,
        ...         phi=porosity)
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = ExtendedGravelBedrockEroder(grid)
        >>> input_var =[1]
        >>> np.ndim(eroder._create_2D_array_for_input_var(input_var,'test'))
        2
        """

        if np.ndim(input_var) == 2:
            input_var_array = input_var
        elif isinstance(input_var, int) or isinstance(input_var, float):
            input_var_array = np.zeros(
                (np.size(self._grid.nodes.flatten()), self._n_classes)
            )
            input_var_array[:, :] = input_var
        elif np.ndim(input_var) <= 1:
            if isinstance(input_var, list):
                input_var = np.array(input_var)
            if isinstance(input_var, tuple):
                input_var = np.array(list(input_var))

            input_var_array = (
                np.ones((np.size(self._grid.nodes.flatten()), np.size(input_var)))
                * input_var[np.newaxis, :]
            )
        else:
            raise ValueError(f"{var_name} array format is invalid")

        # Make sure the array match the grid number of nodes and
        # the number of classes (lithology/grain sizes).
        if np.shape(input_var_array)[1] > 1:
            if np.shape(input_var_array)[1] != self._n_classes:
                raise ValueError(
                    f"{var_name} array dont match the number of lithology classes"
                )
        if np.shape(input_var_array)[0] != np.size(self._grid.nodes.flatten()):
            raise ValueError(
                f"The size of {var_name} array dont match the number of grid nodes"
            )

        return input_var_array

    def _get_classes_identity(self):
        """""
        The following procedure will identify if the classes represent lithology or grain sizes.
        A flag (self._classes_identity) that indicates whether the component
        deals with lithologies or grain sizes defined as:
        0 = single lithology and grain size
        1 = multiple lithologies
        2 = multiple grain sizes
        """ ""

        self._classes_identity = 0
        is_multi_abrasion = np.size(np.unique(self._abr_coefs)) > 1
        if is_multi_abrasion:
            self._classes_identity = 1
        elif self._n_classes > 1:
            self._classes_identity = 2

    def _calc_gravity_coefficient_star(self):
        """ ""
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components.soil_grading import SoilGrading
        >>> xy_spacing=100
        >>> grid = RasterModelGrid((3, 4), xy_spacing=xy_spacing)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> sed_depth = 1000
        >>> porosity=0.5
        >>> sed_weight = sed_depth * 2650 * (1-porosity)
        >>> grains_weight =[sed_weight, sed_weight, sed_weight]
        >>> grain_sizes = [0.001, 0.01, 0.05]
        >>> sg = SoilGrading(grid,
        ...         meansizes=grain_sizes,
        ...         grains_weight=grains_weight,
        ...         phi=porosity)
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = ExtendedGravelBedrockEroder(grid)
        >>> eroder._g_star
        9769603575225600.0
        """
        self._g_star = _EARTH_GRAV * _SQUARED_SECS_IN_A_YEAR

    def _setup_length_of_flow_link(self):
        """Set up a float or array containing length of the flow link from
        each node, which is needed for the abrasion rate calculations.

        Note: if raster, assumes grid.dx == grid.dy

        """
        if isinstance(self.grid, HexModelGrid):
            self._flow_link_length_over_cell_area = (
                self.grid.spacing / self.grid.area_of_cell[0]
            )
            self._flow_length_is_variable = False
            self._grid_has_diagonals = False
        elif isinstance(self.grid, DiagonalsMixIn):
            self._flow_length_is_variable = False
            self._grid_has_diagonals = True
            self._flow_link_length_over_cell_area = (
                _D8_CHAN_LENGTH_FACTOR * self.grid.dx / self.grid.area_of_cell[0]
            )
        else:
            self._flow_length_is_variable = True
            self._update_flow_link_length_over_cell_area()
            self._grid_has_diagonals = False

    def _get_sediment_thickness_by_class(self):
        """ "
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components.soil_grading import SoilGrading
        >>> xy_spacing=100
        >>> grid = RasterModelGrid((3, 4), xy_spacing=xy_spacing)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> sed_depth = 3
        >>> porosity=0
        >>> rho_sed = 2650
        >>> sed_weight = (sed_depth * rho_sed  * (1-porosity))/3 # weight per node area
        >>> grains_weight =[sed_weight, sed_weight, sed_weight]
        >>> grain_sizes = [0.001, 0.01, 0.05]
        >>> sg = SoilGrading(grid,
        ...         meansizes=grain_sizes,
        ...         grains_weight=grains_weight,
        ...         phi=porosity,
        ...         soil_density=rho_sed )
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = ExtendedGravelBedrockEroder(grid, sediment_porosity=porosity, rho_sed=rho_sed)
        >>> eroder._thickness_by_class[grid.core_nodes[0]]
        array([1., 1., 1.])
        """

        # Grains weight is the weight per grid cell area
        self._thickness_by_class = np.divide(
            self._grid.at_node["grains__weight"],
            (self._rho_sed * (1 - self._sediment_porosity)),
        )

    def _update_flow_link_length_over_cell_area(self):
        """Update the ratio of the length of link along which water flows out of
        each node to the area of the node's cell."""
        self._flow_link_length_over_cell_area = (
            self.grid.length_of_link[self._receiver_link[self.grid.core_nodes]]
            / self.grid.area_of_cell[self.grid.cell_at_node[self.grid.core_nodes]]
        )

    def calc_rock_exposure_fraction(self):
        """Update the bedrock exposure fraction.

        The result is stored in the ``bedrock__exposure_fraction`` field.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components.soil_grading import SoilGrading
        >>> xy_spacing=100
        >>> grid = RasterModelGrid((3, 4), xy_spacing=xy_spacing)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> sed_depth = 1000
        >>> porosity=0.5
        >>> sed_weight = sed_depth * xy_spacing * xy_spacing * 2650 * (1-porosity)
        >>> grains_weight =[sed_weight, sed_weight, sed_weight]
        >>> grain_sizes = [0.001, 0.01, 0.05]
        >>> sg = SoilGrading(grid,
        ...         meansizes=grain_sizes,
        ...         grains_weight=grains_weight,
        ...         phi=porosity)
        >>> grid.at_node['soil__depth'][6] = 0.0
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = ExtendedGravelBedrockEroder(grid)
        >>> eroder.calc_rock_exposure_fraction()
        >>> eroder._rock_exposure_fraction[5:7]
        array([0., 1.])
        >>> grid.at_node['soil__depth'][5] = 1.0  # exposure frac should be 1/e ~ 0.3679
        >>> grid.at_node['soil__depth'][6] = 2.0  # exposure frac should be 1/e^2 ~ 0.1353
        >>> eroder.calc_rock_exposure_fraction()
        >>> np.round(eroder._rock_exposure_fraction[5:7], 4)
        array([0.3679, 0.1353])
        """
        self._rock_exposure_fraction[:] = np.exp(-self._sed / self._depth_decay_scale)

    def _calc_tau_star_c(self):
        """Calculate tau_star_c based on Komar 1987

        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components.soil_grading import SoilGrading
        >>> xy_spacing=100
        >>> grid = RasterModelGrid((3, 4), xy_spacing=xy_spacing)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> sed_depth = 1000
        >>> porosity=0.5
        >>> sed_weight = sed_depth * xy_spacing * xy_spacing * 2650 * (1-porosity)
        >>> grains_weight =[sed_weight, sed_weight, sed_weight]
        >>> grain_sizes = [0.001, 0.01, 0.05]
        >>> sg = SoilGrading(grid,
        ...         meansizes=grain_sizes,
        ...         grains_weight=grains_weight,
        ...         phi=porosity)
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = ExtendedGravelBedrockEroder(grid)
        >>> eroder._calc_tau_star_c()
        >>> eroder._tau_star_c[grid.core_nodes[0]]
        array([0.21538354, 0.045     , 0.01506305])
        """

        fractions_sizes = self._grid.at_node["grains_classes__size"]
        median_size_at_node = self._grid.at_node["median_size__weight"][:, np.newaxis]
        tau_star_c = self._tau_star_c
        tau_star_c[:] = np.inf
        if self._n_classes > 1:
            tau_star_c[self.grid.core_nodes, :] = (
                self._tau_star_c_median
                * np.divide(
                    fractions_sizes[self.grid.core_nodes, :],
                    median_size_at_node[self.grid.core_nodes],
                )
                ** self._alpha
            )
        else:
            tau_star_c[self.grid.core_nodes, :] = (
                self._tau_star_c_median
                * np.divide(
                    fractions_sizes[self.grid.core_nodes, np.newaxis],
                    median_size_at_node[self.grid.core_nodes],
                )
                ** self._alpha
            )

    def _update_slopes(self):
        """Update self._slope.

        Result is stored in field ``topographic__steepest_slope``.
        """
        dz = np.maximum(self._elev - self._elev[self._receiver_node], 0.0)
        if self._flow_length_is_variable or self._grid_has_diagonals:
            if self._grid_has_diagonals:
                link_len = self.grid.length_of_d8
            else:
                link_len = self.grid.length_of_link
            self._slope[self.grid.core_nodes] = (
                dz[self.grid.core_nodes] / link_len[self.grid.core_nodes]
            )
        else:
            self._slope[self.grid.core_nodes] = (
                dz[self.grid.core_nodes] / self.grid.spacing
            )

    def calc_abrasion_rate(self):
        """Update the volume rate of bedload loss to abrasion, per unit area.

        Here we use the average of incoming and outgoing sediment flux to
        calculate the loss rate to abrasion in each sediment class.
        The result is stored in self._sed_abr_rates.

        The factor dx (node spacing) appears in the denominator to represent
        flow segment length (i.e., length of the link along which water is
        flowing in the cell) divided by cell area. This would need to be updated
        to handle non-raster and/or non-uniform grids.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components.soil_grading import SoilGrading
        >>> xy_spacing=100
        >>> grid = RasterModelGrid((3, 4), xy_spacing=xy_spacing)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> sed_depth = 1000
        >>> porosity=0.5
        >>> sed_weight = sed_depth * xy_spacing * xy_spacing * 2650 * (1-porosity)
        >>> grains_weight =[sed_weight, sed_weight, sed_weight]
        >>> grain_sizes = [0.001, 0.01, 0.05]
        >>> sg = SoilGrading(grid,
        ...         meansizes=grain_sizes,
        ...         grains_weight=grains_weight,
        ...         phi=porosity)
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = ExtendedGravelBedrockEroder(grid, abrasion_coefficients=[1])
        >>> eroder._sed_influxes[grid.core_nodes[0]]=1
        >>> eroder.calc_abrasion_rate()
        >>> eroder._sed_abr_rates[grid.core_nodes[0]]
        array([0.005, 0.005, 0.005])
        """

        num_sed_classes = self._n_classes
        num_core_nodes = np.size(self._grid.core_nodes)
        flow_link_length_over_cell_area = self._flow_link_length_over_cell_area
        core_nodes = self._grid.core_nodes
        sed_abr_rates = self._sed_abr_rates
        sed_abr_coeff = self._abr_coefs
        sed_influxes = self._sed_influxes
        sed_outfluxes = self._sed_outfluxes
        _calc_sed_abrs_rate(
            num_sed_classes,
            num_core_nodes,
            flow_link_length_over_cell_area,
            core_nodes,
            sed_abr_rates,
            sed_abr_coeff,
            sed_influxes,
            sed_outfluxes,
        )

    def calc_bedrock_abrasion_rate(self):
        """Update the rate of bedrock abrasion.

        Note: assumes _abrasion (of sediment) and _rock_exposure_fraction
        have already been updated. Like _abrasion, the rate is a length
        per time (equivalent to rate of lowering of the bedrock surface by
        abrasion). Result is stored in the field ``bedrock__abrasion_rate``.

        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components.soil_grading import SoilGrading
        >>> xy_spacing=100
        >>> grid = RasterModelGrid((3, 4), xy_spacing=xy_spacing)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> sed_depth = 1000
        >>> porosity=0.5
        >>> sed_weight = sed_depth * xy_spacing * xy_spacing * 2650 * (1-porosity)
        >>> grains_weight =[sed_weight, sed_weight, sed_weight]
        >>> grain_sizes = [0.001, 0.01, 0.05]
        >>> sg = SoilGrading(grid,
        ...         meansizes=grain_sizes,
        ...         grains_weight=grains_weight,
        ...         phi=porosity)
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = ExtendedGravelBedrockEroder(grid, abrasion_coefficients=[1])
        >>> eroder._sed_influxes[grid.core_nodes[0]]=1
        >>> eroder._br_abr_coef[grid.core_nodes[0]]=10**-4
        >>> eroder._rock_exposure_fraction[grid.core_nodes[0]]=0.1
        >>> eroder.calc_bedrock_abrasion_rate()
        >>> int(round(eroder._rock_abrasion_rate[grid.core_nodes[0]]*10**8))
        15
        """
        num_sed_classes = self._n_classes
        num_core_nodes = np.size(self._grid.core_nodes)
        flow_link_length_over_cell_area = self._flow_link_length_over_cell_area
        rock_exposure_fraction = self._rock_exposure_fraction
        core_nodes = self._grid.core_nodes
        bedrock_abrasion_rate = self._rock_abrasion_rate
        bedrock_abr_coeff = self._br_abr_coef
        sed_influxes = self._sed_influxes
        sed_outfluxes = self._sed_outfluxes
        _calc_bedrock_abrs_rate(
            num_sed_classes,
            num_core_nodes,
            flow_link_length_over_cell_area,
            rock_exposure_fraction,
            core_nodes,
            bedrock_abrasion_rate,
            bedrock_abr_coeff,
            sed_influxes,
            sed_outfluxes,
        )

    def _calc_tau_star(self):
        """Calculate tau star at node
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components.soil_grading import SoilGrading
        >>> xy_spacing=100
        >>> grid = RasterModelGrid((3, 4), xy_spacing=xy_spacing)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> sed_depth = 1000
        >>> porosity=0.5
        >>> sed_weight = sed_depth * xy_spacing * xy_spacing * 2650 * (1-porosity)
        >>> grains_weight =[sed_weight, sed_weight, sed_weight]
        >>> grain_sizes = [0.001, 0.01, 0.05]
        >>> sg = SoilGrading(grid,
        ...         meansizes=grain_sizes,
        ...         grains_weight=grains_weight,
        ...         phi=porosity)
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = ExtendedGravelBedrockEroder(grid, abrasion_coefficients=[1])
        >>> eroder._tau[grid.core_nodes[0]] = 10
        >>> eroder._calc_tau_star()
        >>> np.round(eroder._tau_star[grid.core_nodes[0]],3)
        array([0.618, 0.062, 0.012])


        """

        if self._n_classes > 1:
            self._tau_star = np.divide(
                self._tau[:, np.newaxis],
                (
                    (self._SG * self._rho_water)
                    * _EARTH_GRAV
                    * self._grid.at_node["grains_classes__size"]
                ),
            )
        else:
            self._tau_star = np.divide(
                self._tau,
                (
                    (self._SG * self._rho_water)
                    * _EARTH_GRAV
                    * self._grid.at_node["grains_classes__size"]
                ),
            )
            self._tau_star = self._tau_star[:, np.newaxis]

    def calc_transport_rate(self):
        """Calculate and return bed-load transport rate.

        Calculation uses Wickert-Schildgen approach, and provides
        volume per time rate. Transport rate is modulated by available
        sediment, using the exponential function ``(1 - exp(-H / Hs))``,
        so that transport rate approaches zero as sediment thickness
        approaches zero. Rate is a volume per time. The result is
        stored in the ``bedload_sediment__volume_outflux`` field.

        Examples
        """

        self._calc_width()
        self._tau[:] = _calc_shear_stress(
            shear_stress_coef=self._shear_stress_coef,
            discharge=self._grid.at_node["surface_water__discharge"],
            width=self._channel_width,
            slope=self._slope,
        )

        self._calc_tau_star()
        self._calc_tau_star_c()

        self._excess_stress[:] = self._tau_star - self._tau_star_c
        self._excess_stress[self._excess_stress < 0] = 0

        # Get sediment flux
        qs = self._calc_qs(excess_stress=self._excess_stress)
        Qs = self._channel_width[:, np.newaxis] * qs
        self._sed_outfluxes[:] = Qs * (
            1.0 - self._rock_exposure_fraction[:, np.newaxis]
        )
        self._sed_outfluxes[
            self._grid.at_node["grains__weight"] <= self._weight_threshold_to_deliv
        ] = 0
        if self._n_classes > 1:
            self._sediment_outflux[:] = np.sum(self._sed_outfluxes, axis=1)
        else:
            self._sediment_outflux[:] = self._sed_outfluxes[:, 0]

    def _calc_qs(self, excess_stress, out=None):

        if out is None:
            out = np.empty_like(self._grid.at_node["grains_classes__size"])

        if self._n_classes > 1:
            out[:] = (
                3.97
                * self._SG ** (0.5)
                * self._intermittency_factor
                * self._g_star ** (0.5)
                * (excess_stress) ** (1.5)
                * self._grid.at_node["grains_classes__size"] ** (1.5)
            )
        else:
            out[:] = (
                3.97
                * self._SG ** (0.5)
                * self._intermittency_factor
                * self._g_star ** (0.5)
                * (excess_stress[:, 0]) ** (1.5)
                * self._grid.at_node["grains_classes__size"] ** (1.5)
            )
            out = out[:, np.newaxis]
        return out

    def _calc_weight_threshold_to_deliv(self):
        """Calc minimal weight threshold to deliver"""
        # d_min = np.min(self._grid.at_node['grains_classes__size'])
        d_min = 0.1
        self._weight_threshold_to_deliv = (
            d_min * self._rho_sed * (1 - self._sediment_porosity)
        )

    def calc_bedrock_plucking_rate(self):
        """Update the rate of bedrock erosion by plucking.

        The rate is a volume per area per time [L/T], equivalent to the
        rate of lowering of the bedrock surface relative to the underlying
        material as a result of plucking. Result is stored in the field
        ``bedrock__plucking_rate``.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components.soil_grading import SoilGrading
        >>> xy_spacing=100
        >>> grid = RasterModelGrid((3, 4), xy_spacing=xy_spacing)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> sed_depth = 1000
        >>> porosity=0.5
        >>> sed_weight = sed_depth * xy_spacing * xy_spacing * 2650 * (1-porosity)
        >>> grains_weight =[sed_weight, sed_weight, sed_weight]
        >>> grain_sizes = [0.001, 0.01, 0.05]
        >>> sg = SoilGrading(grid,
        ...         meansizes=grain_sizes,
        ...         grains_weight=grains_weight,
        ...         phi=porosity)
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = ExtendedGravelBedrockEroder(grid, abrasion_coefficients=[1], plucking_coefficient=[1])
        >>> eroder._channel_width[:] = 1
        >>> eroder._intermittency_factor = 0.1
        >>> eroder._rock_exposure_fraction[:] = 1
        >>> eroder._excess_stress[:] = 1
        >>> eroder._flow_link_length_over_cell_area = 1
        >>> eroder.calc_bedrock_plucking_rate()
        >>> print(round(eroder._pluck_rate[grid.core_nodes[0]],1))
        0.1
        """

        if self._n_classes > 1:
            classes_fractions = self._get_classes_fractions()
        else:
            classes_fractions = np.ones_like(self._excess_stress)
        _calc_pluck_rate(
            self._n_classes,
            self.grid.number_of_core_nodes,
            self._intermittency_factor,
            self._flow_link_length_over_cell_area,
            self._grid.core_nodes,
            self._plucking_coef,
            self._channel_width,
            self._rock_exposure_fraction,
            classes_fractions,
            self._excess_stress,
            self._pluck_rate,
        )

    def _get_classes_fractions(self):

        out = np.zeros((self._grid.number_of_nodes, self._n_classes))
        classes_fractions = _get_classes_fractions(
            self._n_classes,
            self.grid.number_of_core_nodes,
            self._grid.core_nodes,
            self._grid.at_node["grains__weight"],
            out,
        )
        return classes_fractions

    def calc_sediment_influx(self):
        """
        Update the volume influx at each node for each sediment class.

        Results are stored:
        (1) Total flux in the field ``bedload_sediment__volume_influx``
        (2) Per size class in the _sed_influxes array
        """
        _calc_sediment_influx(
            self._n_classes,
            self.grid.number_of_core_nodes,
            self._sediment_influx,
            self._sed_influxes,
            self._sediment_outflux,
            self._sed_outfluxes,
            self.grid.core_nodes,
            self._receiver_node,
        )

    def calc_sediment_rate_of_change(self):
        """
        Calculate and store time rate of change of sediment thickness
        by sediment class. Result stored in self._dHdt_by_class.

        The doctest below illustrates some of the steps in calculating
        the rate of change of sediment thickness. calc_transport_rate()
        sets the _sediment_outflux at each core node. The value of
        outflux can be calculated as follows:

        outflux = transport coefficient x intermittency factor
                  x discharge x slope^(7/6) x (1 - rock exposure fraction)

        Using default parameters, assuming a runoff rate of unity, and
        given the negligible bedrock exposure fraction,
        for the upstream grid cell this works out to:

        outflux = 0.041 x 0.01 x (100 x 100) x 0.01^(7/6) x 1
                ~ 0.19 m3/y

        For the downstream node, the discharge is doubled, so the
        flux is doubled.

        The call to calc_sediment_influx() then passes the outflux to
        each downstream neighbor (the "receiver" node), in the array
        _sediment_influx. For purposes of the test, we are *not*
        calculating sediment loss to abrasion, which therefore remains
        at its initial value of zero for each node.

        The call to calc_sediment_rate_of_change() adds up the influx
        and subtracts the outflux for each core node, dividing by its
        cell area and multiplying by the porosity factor to get a rate
        of thickness change. For both nodes, we have a net outflux of
        ~0.19 m3/y, a porosity factor of 1 / (1 - 0.35) ~ 1.54, and a
        grid cell area of 10,000 m2, so the rate of change should be
        ~2.93 x 10^-5 m/y. This is stored in _dHdt_by_class; because
        we have only on sediment class, these are the values that should
        be assigned to the single class at each of the two core nodes.

        Examples
        --------
        # >>> import numpy as np
        # >>> from landlab import RasterModelGrid
        # >>> from landlab.components import FlowAccumulator
        # >>> grid = RasterModelGrid((3, 4), xy_spacing=100.0)
        # >>> elev = grid.add_zeros("topographic__elevation", at="node")
        # >>> elev[:] = 0.01 * grid.x_of_node
        # >>> sed = grid.add_zeros("soil__depth", at="node")
        # >>> sed[:] = 100.0
        # >>> grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
        # >>> grid.status_at_node[4] = grid.BC_NODE_IS_FIXED_VALUE
        # >>> fa = FlowAccumulator(grid)
        # >>> fa.run_one_step()
        # >>> eroder = ExtendedGravelBedrockEroder(grid)
        # >>> eroder.calc_transport_rate()
        # >>> eroder.calc_sediment_influx()
        # >>> eroder.calc_sediment_rate_of_change()
        # >>> np.round(eroder._sed_outfluxes[0, 4:7], 3)
        # array([ 0.   ,  0.038,  0.019])
        # >>> np.round(eroder._sed_influxes[0, 4:7], 3)
        # array([ 0.038,  0.019,  0.   ])
        # >>> np.round(eroder._dHdt_by_class[0, 5:7], 8)
        # array([ -2.93000000e-06,  -2.93000000e-06])
        """

        _calc_sediment_rate_of_change(
            self._n_classes,
            self.grid.number_of_core_nodes,
            self._porosity_factor,
            self.grid.area_of_cell[0],
            self.grid.core_nodes,
            self._fractions_from_plucking,
            self._dHdt,
            self._pluck_rate,
            self._dHdt_by_class,
            self._sed_influxes,
            self._sed_outfluxes,
            self._sed_abr_rates,
            self._br_abrasion_coef,
        )

    def update_rates(self):
        """Update rate of sediment thickness change, and rate of bedrock lowering by abrasion
        and plucking.

        Combined rate of rock lowering relative to underlying material is stored in the field
        ``bedrock__lowering_rate``.

        The doctest below evaluates the code against the following
        calculation. We have two core nodes, each with a gradient of
        0.01 m/m and a surface area of 10,000 m2.

        The sediment cover is 0.69315 m, which happens to give a
        cover fraction of ~0.5.

        We run just one time step, in which the transport and abrasion
        rates apply to an initial slope of 0.01 and a discharge of
        10,000 and 20,000 m3/y at the upstream and downstream nodes,
        respectively.

        Transport rate, upstream node:

        coefficient x intermittency x discharge x slope^(7/6) x cover
        0.041       x 0.01          x 10,000 m3/y x 0.004642 x 0.5
        ~0.0095 m3/y

        Downstream node: ~0.019 m3/y (because of 2x discharge)

        Sediment influxes: at the open boundary node, should equal the
        outflux from the downstream core node (0.019 m3/y); at the
        downstream core node, should equal outflux from the upstream
        node (0.095 m3/y), and at the upstream node should be zero.

        With 50% bedrock exposure, the bedrock plucking rate for the
        upstream node should be:

        pluck coef x intermittency x discharge x slope^(7/6) x exposure
        / width of grid cell

        1e-4 x 0.01 x 10,000 x 0.01^(7/6) x 0.5 / 100 ~ 2.32e-7

        For the downstream node, it should be twice this value.

        With an abrasion coefficient of 0.001 1/m (fast!), the
        sediment abrasion rate at the upstream node should derive
        from the average of influx and outflux:

        0.5 x (0.0095 + 0) m3/y x 0.001 1/m  x 100 m ~ 0.000475 m3/y
        = 4.75e-8 m/y lowering equivalent

        For the downstream node,

        0.5 x (0.019 + 0.0095) m3/y x 0.001 1/m  x 100 m ~ 0.001425 m3/y
        = 1.425e-7 m/y lowering equivalent

        The sediment rate of change should be the sum of 3 terms:
        (1) the difference between influx and outflux, multiplied by
        the porosity factor 1 / (1 - phi) and divided by cell area.
        That equates to:
        -1 / (1 - 0.35) x 0.0095 / 10000 ~ -1.46 x 10^-6
        (2) the abrasion lowering rate above.
        (3) the plucking rate times the coarse fraction derived from
        plucking:
        + 2.32e-7 m/y x 0.25 ~ 5.8e-8 upstream, and 1.16e-7 downstream.
        Their sum is (approximately, subject to rounding errors):
        ~ -1.46e-6 - 1.43e-7 + 1.16e-7 = -1.487e-6 m/y downstream,
        ~ -1.46e-6 - 4.8e-8 + 5.8e-8 = -1.45e-6 m/y upstream.

        This should be stored in both _dHdt_by_class (the rate of
        thickening per sediment class, of which there's only one in
        this case) and _dHdt (which totals up all the classes).

        This rate is then extrapolated for one time step of 1000 y,
        so that the thickness is reduced by about 1.5 mm. The resulting
        thickness should be its original value minus this amount:
        0.69315 - 1.45e-3 = 0.6917 m, 0.69315 - 1.5e-3 = 0.69165 m

        Rock elevation should reduced by the loss of sediment plus the loss
        of rock. Rock loss is (plucking rate + rock abrasion rate) x time step.
        Upstream: (2 - 0.69315) - (2.32e-7 + 0.5 * 4.76e-8) m x 1000 y = 1.3065942 m
        Downstream: (1 - 0.69315) - (2.32e-7 + 0.5 * 4.76e-8) m x 1000 y = 0.3063145 m

        Finally, total elevation is the sum of the new sediment thickness
        and rock elevation (again, approximate, subject to rounding):
        Upstream: 0.3063 + 0.6917 = 0.9980
        Downstream: 1.3066 + 0.6917 = 1.9983

        Examples
        --------
        # >>> import numpy as np
        # >>> from landlab import RasterModelGrid
        # >>> from landlab.components import FlowAccumulator
        # >>> grid = RasterModelGrid((3, 4), xy_spacing=100.0)
        # >>> elev = grid.add_zeros("topographic__elevation", at="node")
        # >>> elev[:] = 0.01 * grid.x_of_node
        # >>> sed = grid.add_zeros("soil__depth", at="node")
        # >>> sed[:] = 0.69315
        # >>> rock = grid.add_zeros("bedrock__elevation", at="node")
        # >>> rock[:] = elev - sed
        # >>> grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
        # >>> grid.status_at_node[4] = grid.BC_NODE_IS_FIXED_VALUE
        # >>> fa = FlowAccumulator(grid)
        # >>> fa.run_one_step()
        # >>> eroder = ExtendedGravelBedrockEroder(
        # ...     grid,
        # ...     abrasion_coefficients=[0.001],
        # ...     coarse_fractions_from_plucking=[0.25],
        # ... )
        # >>> eroder._bedrock__elevation[5:7]
        # array([ 0.30685,  1.30685])
        # >>> eroder.run_one_step(1000.0)
        # >>> np.round(eroder._sed_outfluxes[0, 5:7], 3)
        # array([ 0.019,  0.01 ])
        # >>> np.round(eroder._sed_influxes[0, 4:7], 3)
        # array([ 0.019,  0.01 ,  0.   ])
        # >>> np.round(eroder._pluck_rate[5:7], 10)
        # array([  4.64200000e-07,   2.32100000e-07])
        # >>> np.round(eroder._sed_abr_rates[0, 5:7], 9)
        # array([  1.43000000e-07,   4.80000000e-08])
        # >>> np.round(eroder._dHdt_by_class[0, 5:7], 8)
        # array([ -1.50000000e-06,  -1.45000000e-06])
        # >>> np.round(eroder._dHdt[5:7], 8)
        # array([ -1.50000000e-06,  -1.45000000e-06])
        # >>> np.round(sed[5:7], 4)
        # array([ 0.6916,  0.6917])
        # >>> np.round(eroder._bedrock__elevation[5:7], 5)
        # array([ 0.30631,  1.30659])
        # >>> np.round(elev[4:7], 4)
        array([ 0.    ,  0.998 ,  1.9983])
        """
        self._init_vars()
        self._update_slopes()
        self.calc_rock_exposure_fraction()
        self.calc_transport_rate()
        self.calc_sediment_influx()

        if self._flow_length_is_variable:
            self._update_flow_link_length_over_cell_area()
        self.calc_bedrock_plucking_rate()
        if np.amax(self._abr_coefs) > 0.0:
            self.calc_abrasion_rate()
        if np.amax(self._br_abr_coef) > 0.0:
            self.calc_bedrock_abrasion_rate()
        self.calc_sediment_rate_of_change()
        self._rock_lowering_rate[self.grid.core_nodes] = (
            self._pluck_rate[self.grid.core_nodes]
            + self._rock_abrasion_rate[self.grid.core_nodes]
        )

    def _update_rock_sed_and_elev(self, dt):
        """Update rock elevation, sediment thickness, and elevation
        using current rates of change extrapolated forward by time dt.
        """
        weights_at_node = self._grid.at_node["grains__weight"]

        # Update grains weight based on fluxes
        weight_dt_by_class = (
            self._dHdt_by_class * self._rho_sed * (1 - self._sediment_porosity) * dt
        )
        if self._n_classes > 1:
            weights_at_node[:] += weight_dt_by_class
            weights_at_node[weights_at_node <= 0] = 0  # Ensure non-negative weight

            # Update sediment thickness (sum of all grain size classes)
            self._sed[self.grid.core_nodes] = np.sum(
                self._grid.at_node["grains__weight"][self.grid.core_nodes], axis=1
            ) / (self._rho_sed * (1 - self._sediment_porosity))
        else:
            weights_at_node[:] += weight_dt_by_class[:, 0]
            weights_at_node[weights_at_node <= 0] = 0  # Ensure non-negative weight

            # Update sediment thickness (sum of all grain size classes)
            self._sed[self.grid.core_nodes] = self._grid.at_node["grains__weight"][
                self.grid.core_nodes
            ] / (self._rho_sed * (1 - self._sediment_porosity))
        # Update bedrock lowering
        self._bedrock__elevation[self.grid.core_nodes] -= (
            self._rock_lowering_rate[self.grid.core_nodes] * dt
        )

        # Update elevation
        self._elev[self.grid.core_nodes] = (
            self._bedrock__elevation[self.grid.core_nodes]
            + self._sed[self.grid.core_nodes]
        )

    def _calc_sed_transport_coeff(self):
        tr_coeff = (
            self._sediment_transport_rate_coeff
        )  # From Wong and Parker 2006 (they call it phi).
        epsilon = self._epsilon
        tau_star_c = self._tau_star_c
        kQs = self._kQs
        SG = self._SG
        kQs[:] = np.divide(
            0.17 * tr_coeff * epsilon ** (3 / 2),
            SG * (1 + epsilon) ** (5 / 3) * tau_star_c ** (1 / 6),
        )

    def _estimate_max_time_step_size(self, upper_limit_dt=1.0e6):
        """
        Estimate the maximum possible time-step size that avoids
        flattening any streamwise slope or exhausting sediment.

        The ``upper_limit_dt`` parameter handles the special case of
        a nonexistent upper limit, which only occurs when there are
        no nodes at which either sediment or slope gradient is
        declining. Value is arbitrary as long as it is >= the user-provided
        global time-step size (in :meth:`~.run_one_step`).

        Parameters
        ----------
        dt : float (default 1.0e6)
            Maximum time step size
        """

        dhdt_by_class = self._dHdt_by_class
        dh_by_class = self._grid.at_node["grains__weight"] / (
            self._rho_sed * (1 - self._sediment_porosity)
        )

        if self._n_classes > 1:
            sed_is_declining = np.logical_and(dhdt_by_class < 0.0, dh_by_class > 0.0)
        else:
            sed_is_declining = np.logical_and(
                dhdt_by_class[:, 0] < 0.0, dh_by_class > 0.0
            )

        if np.any(sed_is_declining):
            min_time_to_exhaust_sed = np.amin(
                dh_by_class[sed_is_declining] / np.abs(dhdt_by_class[sed_is_declining])
            )
        else:
            min_time_to_exhaust_sed = upper_limit_dt

        dzdt = self._dHdt - self._rock_lowering_rate
        rate_diff = dzdt[self._receiver_node] - dzdt
        height_above_rcvr = self._elev - self._elev[self._receiver_node]
        slope_is_declining = np.logical_and(rate_diff > 0.0, height_above_rcvr > 0.0)
        if np.any(slope_is_declining):
            min_time_to_flatten_slope = np.amin(
                height_above_rcvr[slope_is_declining] / rate_diff[slope_is_declining]
            )
        else:
            min_time_to_flatten_slope = upper_limit_dt
        min_dt = 0.5 * min(min_time_to_exhaust_sed, min_time_to_flatten_slope)

        return min_dt

    def _calc_width(self):
        """Calculate width
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components.soil_grading import SoilGrading
        >>> xy_spacing=100
        >>> grid = RasterModelGrid((3, 4), xy_spacing=xy_spacing)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> sed_depth = 1000
        >>> porosity=0.5
        >>> sed_weight = sed_depth * xy_spacing * xy_spacing * 2650 * (1-porosity)
        >>> grains_weight =[sed_weight, sed_weight, sed_weight]
        >>> grain_sizes = [0.001, 0.01, 0.05]
        >>> sg = SoilGrading(grid,
        ...         meansizes=grain_sizes,
        ...         grains_weight=grains_weight,
        ...         phi=porosity)
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> grid.at_node['surface_water__discharge'][grid.core_nodes[0]] = 10000000
        >>> eroder = ExtendedGravelBedrockEroder(grid, fixed_width_flag=True)
        >>> eroder._calc_width()
        >>> round(eroder._channel_width[grid.core_nodes[0]])
        6
        """

        discharge = self._grid.at_node["surface_water__discharge"]
        if self._fixed_width_flag:
            coeff = self._fixed_width_coeff
            expt = self._fixed_width_expt
            _calc_chan_width_fixed_width(
                coeff, expt, discharge, out=self._channel_width
            )
        else:

            # SOME EXPLANATION
            # MORE EFFCIENT WAY?
            # row, col = np.where(
            #     self.grid.at_node['median_size__weight'][:, np.newaxis] == self._grid.at_node['grains_classes__size'])
            # tau_star_c_median = self._tau_star_c[row, col]

            tau_star_c_median = self._tau_star_c_median
            median_size = self.grid.at_node["median_size__weight"]
            slope = self._slope
            g_star = self._g_star
            SG = self._SG
            epsilon = self._epsilon

            self._channel_width[:] = _calc_near_threshold_width(
                median_size, tau_star_c_median, discharge, slope, g_star, SG, epsilon
            )

    def _init_vars(self):
        """Initialize variables"""
        self._rock_abrasion_rate.fill(0.0)
        self._br_abr_coef.fill(0.0)
        self._sed_influxes.fill(0.0)
        self._sed_outfluxes.fill(0.0)
        self._rock_exposure_fraction.fill(0.0)
        self._excess_stress.fill(0.0)

    def run_one_step(self, global_dt):
        """Advance solution by time interval global_dt, subdividing
        into sub-steps as needed."""
        time_remaining = global_dt
        while time_remaining > 0.0:
            self.update_rates()
            max_dt = self._estimate_max_time_step_size()
            this_dt = min(max_dt, time_remaining)
            this_dt = max(this_dt, _DT_MAX)
            self._update_rock_sed_and_elev(this_dt)
            time_remaining -= this_dt


if __name__ == "__main__":
    import doctest

    doctest.testmod(name="calc_rock_exposure_fraction")
