from __future__ import print_function

import warnings

from six.moves import range

import numpy as np
from time import sleep
from landlab import ModelParameterDictionary, CLOSED_BOUNDARY, Component

from landlab.core.model_parameter_dictionary import MissingKeyError
from landlab.field.scalar_data_fields import FieldError
from landlab.grid.base import BAD_INDEX_VALUE
from landlab.utils.decorators import make_return_array_immutable


def sed_flux_fn_gen_genhump(rel_sed_flux_in, kappa, nu, c, phi, norm):
    """
    Returns f(qs,qc) assuming a generalized humped function per
    Hobley et al., 2011.
    """
    return (
        norm*kappa*(rel_sed_flux_in**nu + c)*np.exp(-phi*rel_sed_flux_in))


def sed_flux_fn_gen_lindecl(rel_sed_flux_in, kappa, nu, c, phi, norm):
    """
    Returns f(qs,qc) assuming a linear decline model (see e.g. Gasparini
    et al., 2006).
    kappa, nu, c, phi, & norm are all dummy variables here.
    """
    return 1.-rel_sed_flux_in


def sed_flux_fn_gen_almostparabolic(rel_sed_flux_in, kappa, nu, c, phi,
                                    norm):
    """
    Returns f(qs,qc) assuming the almost parabolic humped model (see
    Gasparini et al., 2006).
    kappa, nu, c, phi, & norm are all dummy variables here.
    """
    return np.where(rel_sed_flux_in > 0.1,
                    1. - 4.*(rel_sed_flux_in-0.5)**2.,
                    2.6*rel_sed_flux_in+0.1)


def sed_flux_fn_gen_lindecl(rel_sed_flux_in, kappa, nu, c, phi, norm):
    """
    Returns f(qs,qc) assuming a linear decline model (see e.g. Gasparini
    et al., 2006).
    kappa, nu, c, phi, & norm are all dummy variables here.
    """
    return 1.-rel_sed_flux_in


def sed_flux_fn_gen_const(rel_sed_flux_in, kappa, nu, c, phi, norm):
    """
    Returns 1, and thus no sed flux effects.
    kappa, nu, c, phi, & norm are all dummy variables here.
    """
    return 1.


def get_sed_flux_function_pseudoimplicit(sed_in_bydt,
                                         trans_cap_vol_out_bydt,
                                         prefactor_for_volume_bydt,
                                         prefactor_for_dz_bydt,
                                         sed_flux_fn_gen,
                                         kappa, nu, c, phi, norm,
                                         pseudoimplicit_repeats):
    """
    This function uses a pseudoimplicit method to calculate the sediment
    flux function for a node, and also returns dz/dt and the rate of
    sediment output from the node.

    Note that this method now operates in PER TIME units; this was not
    formerly the case.

    Parameters
    ----------
    sed_in_bydt : float
        Total rate of incoming sediment, sum(Q_s_in)/dt
    trans_cap_vol_out_bydt : float
        Volumetric transport capacity as a rate (i.e., m**3/s) on outgoing
        link
    prefactor_for_volume_bydt : float
        Equal to K*A**m*S**n * cell_area
    prefactor_for_dz_bydt : float
        Equal to K*A**m*S**n (both prefactors are passed for computational
        efficiency)
    sed_flux_fn_gen : function
        Function to calculate the sed flux function. Takes inputs
        rel_sed_flux_in, kappa, nu, c, phi, norm, where last 5 are dummy
        unless type is generalized_humped.
    kappa, nu, c, phi, norm : float
        Params for the sed flux function. Values if generalized_humped,
        zero otherwise.
    pseudoimplicit_repeats : int
        Maximum number of loops to perform with the pseudoimplicit
        iterator, seeking a stable solution. Convergence is typically
        rapid.

    Returns
    -------
    dzbydt : float
        Rate of change of substrate elevation
    sed_flux_out_bydt : float
        Q_s/dt on the outgoing link
    rel_sed_flux : float
        f(Q_s/Q_c)
    error_in_sed_flux_fn : float
        Measure of how well converged rel_sed_flux is

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import SedDepEroder
    >>> mg = RasterModelGrid((25, 25), 10.)
    >>> sde = SedDepEroder(mg, sed_dependency_type='almost_parabolic')
    >>> (dzbydt, sed_flux_out_bydt, rel_sed_flux, error_in_sed_flux_fn) = (
    ...   get_sed_flux_function_pseudoimplicit(1.e3, 1.e6, 1.e4, 10.,
    ...                                        sde.sed_flux_fn_gen,
    ...                                        0., 0., 0., 0., 0.,
    ...                                        50))
    >>> dzbydt
    1.0393380000000001
    >>> sed_flux_out_bydt
    1776.1689999999999
    >>> rel_sed_flux
    0.001776169
    >>> error_in_sed_flux_fn
    0.0013337999999999961
    """
    rel_sed_flux_in = sed_in_bydt/trans_cap_vol_out_bydt
    rel_sed_flux = rel_sed_flux_in

    last_sed_flux_fn = 2.  # arbitrary value out of calc range
    for i in range(pseudoimplicit_repeats):
        sed_flux_fn = sed_flux_fn_gen(
            rel_sed_flux, kappa, nu, c, phi, norm)
        sed_vol_added_bydt = prefactor_for_volume_bydt*sed_flux_fn
        prop_added = sed_vol_added_bydt/trans_cap_vol_out_bydt
        poss_new = prop_added + rel_sed_flux
        if poss_new > 1.:
            poss_new = 1.
        rel_sed_flux = 0.5 * (rel_sed_flux_in + poss_new)

        if rel_sed_flux < 0.:
            rel_sed_flux = 0.
        error_in_sed_flux_fn = abs(sed_flux_fn-last_sed_flux_fn)
        if error_in_sed_flux_fn < 0.01:
            break
        last_sed_flux_fn = sed_flux_fn
    # note that the method will silently terminate even if we still have
    # bad convergence. Note this is very rare.

    dzbydt = prefactor_for_dz_bydt*sed_flux_fn
    sed_flux_out_bydt = rel_sed_flux*trans_cap_vol_out_bydt
    return dzbydt, sed_flux_out_bydt, rel_sed_flux, error_in_sed_flux_fn


def iterate_downstream(s_in, cell_areas,
                       hillslope_sediment,
                       hillslope_sediment_flux,
                       river_volume_flux_into_node,
                       transport_capacities,
                       erosion_prefactor_withS,
                       rel_sed_flux,
                       is_it_TL,
                       vol_drop_rate,
                       flow_receiver, pseudoimplicit_repeats,
                       dzbydt, sed_flux_fn_gen,
                       kappa, nu, c, phi, norm,
                       flooded_nodes=None,
                       flooded_depths=None):
    """
    Iterates down a drainage network, redistributing sediment and solving
    the sediment flux dependent incision equations.

    Parameters
    ----------
    s_in : array
        The upstream node order
    cell_areas : array
        The areas of all cells in the grid
    hillslope_sediment : array
        Depth of sediment in the channel at each node.
    hillslope_sediment_flux : array
        The existing volume of sediment on the channel bed at a node,
        expressed as volume per unit time of the timestep. This turns
        the accumulated sediment on the bed into a virtual sediment supply
        from upstream, such that at the end of the step the same depth of
        sediment would be present at the node if no other transport
        occurred.
    river_volume_flux_into_node : array
        Total ""true" river flux coming into node from upstream.
    transport_capacities : array
        The bedload transport capacity at each node, expressed as a flux.
    erosion_prefactor_withS : array
        Equal to K * A**m * S**n at nodes
    rel_sed_flux : array
        The sediment flux as a function of the transport capacity.
    is_it_TL : boolean array
        Describes whether the sediment transported at the node is at
        capacity or not.
    vol_drop_rate : array
        Flux of sediment ending up on the bed during transport.
    flow_receiver : array
        The downstream node ID.
    pseudoimplicit_repeats : int
        Maximum number of loops to perform with the pseudoimplicit
        iterator, seeking a stable solution. Convergence is typically
        rapid.
    dzbydt : array
        The rate of change of *bedrock* surface elevation.
    sed_flux_fn_gen : function
        Function to calculate the sed flux function. Takes inputs
        rel_sed_flux_in, kappa, nu, c, phi, norm, where last 5 are dummy
        unless type is generalized_humped.
    kappa, nu, c, phi, norm : float
        Params for the sed flux function. Values if generalized_humped,
        zero otherwise.
    flooded_nodes : boolean nnodes array or None
        Logical array of flooded nodes, if any
    flooded_depths : float nnodes array or None
        Depth of standing water at node (0. where not flooded)
    """
    for i in s_in[::-1]:  # work downstream
        cell_area = cell_areas[i]
        if flooded_nodes is not None:
            flood_depth = flooded_depths[i]
        else:
            flood_depth = 0.
        sed_flux_into_this_node_bydt = (
            hillslope_sediment_flux[i] +
            river_volume_flux_into_node[i])
        node_capacity = transport_capacities[i]
        # ^we work in volume flux, not volume per se here
        if flood_depth - hillslope_sediment[i] > 0.:
            # ...node is truly flooded
            node_capacity = 0.
        if sed_flux_into_this_node_bydt < node_capacity:
            # ^note incision is forbidden at capacity
            dz_prefactor_bydt = erosion_prefactor_withS[i]
            vol_prefactor_bydt = dz_prefactor_bydt*cell_area
            (dzbydt_here, sed_flux_out_bydt, rel_sed_flux_here,
             error_in_sed_flux) = get_sed_flux_function_pseudoimplicit(
                    sed_flux_into_this_node_bydt,
                    node_capacity,
                    vol_prefactor_bydt, dz_prefactor_bydt,
                    sed_flux_fn_gen,
                    kappa, nu, c, phi, norm,
                    pseudoimplicit_repeats)
            # note now dz_here may never create more sed than the
            # out can transport...
            assert sed_flux_out_bydt <= node_capacity, (
                'failed at node '+str(s_in.size-i) +
                ' with rel sed flux '+str(
                    sed_flux_out_bydt/node_capacity))
            rel_sed_flux[i] = rel_sed_flux_here
            vol_pass_rate = sed_flux_out_bydt
        else:
            is_it_TL[i] = True
            rel_sed_flux[i] = 1.
            dzbydt_here = 0.
            try:
                isflooded = flooded_nodes[i]
            except TypeError:  # was None
                isflooded = False
            # flood handling adjusted as of 29/6/17 (untested) such
            # that sed can flux through as long as the sed depth
            # at the start of the loop was capable of filling the
            # hole in the underlying topo.
            # Under this scheme, there should be no such thing as a
            # lake hosted only in surface sediment
            if (flood_depth - hillslope_sediment[i] <= 0. and
                    not isflooded):
                vol_pass_rate = node_capacity
                # we want flooded nodes which have already been
                # filled to enter the else statement
            else:
                vol_pass_rate = 0.
            vol_drop_rate[i] = (
                sed_flux_into_this_node_bydt - vol_pass_rate)
            assert vol_drop_rate[i] >= 0.
        dzbydt[i] -= dzbydt_here
        # ^minus returns us to the correct sign convention
        river_volume_flux_into_node[
            flow_receiver[i]] += vol_pass_rate


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
    calculated.
    ***Note that this new implementation only permits a power_law formulation
    for Qc***
    For Qc, 'power_law' broadly follows the assumptions in Gasparini et
    al. 2006, 2007.

    If ``Qc == 'power_law'``::

        E  = K_sp * f(Qs, Qc) * A ** m_sp * S ** n_sp;
        Qc = K_t * A ** m_t * S ** n_t

    The component uses the field fluvial_sediment__depth as the record of the
    sediment on the bed at any given time. This may be set and/or freely
    modified by the user, but will be created by the component if not found.
    Tools-driven bedrock erosion is permitted only when this layer thickness
    is reduced to zero. The sediment recorded in fluvial_sediment__depth is
    considered loose, and freely transportable by clearwater flow.
    Note that the topography tracked by 'topographic__elevation' is the true
    surface topography, not the bedrock topography (...as is consistent with
    that field name). The bedrock topography can be found from the
    topographic__elevation less the channel_sediment__depth.

    The component is able to handle flooded nodes, if created by a lake
    filler. It assumes the flow paths found in the fields already reflect
    any lake routing operations, and then requires the optional argument
    *flooded_depths* be passed to the run method. A flooded depression
    acts as a perfect sediment trap, and will be filled sequentially
    from the inflow points towards the outflow points.

    Construction::

        SedDepEroder(grid, K_sp=1.e-6, g=9.81, rock_density=2700,
                     sediment_density=2700, fluid_density=1000,
                     runoff_rate=1.,
                     sed_dependency_type='generalized_humped',
                     kappa_hump=13.683, nu_hump=1.13, phi_hump=4.24,
                     c_hump=0.00181, Qc='power_law', m_sp=0.5, n_sp=1.,
                     K_t=1.e-4, m_t=1.5, n_t=1.,
                     pseudoimplicit_repeats=50,
                     return_stream_properties=False)

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
        Maximum number of loops to perform with the pseudoimplicit iterator,
        seeking a stable solution. Convergence is typically rapid. The
        component counts the total number of times it "maxed out" the
        loop to seek a stable solution (error in sed flux fn <1%) with the
        internal variable "_pseudoimplicit_aborts".
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

    Examples
    --------
    >>> from six.moves import range
    >>> import numpy as np
    >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
    >>> from landlab.components import FlowRouter, SedDepEroder

    >>> mg = RasterModelGrid((10, 3), 200.)
    >>> for edge in (mg.nodes_at_left_edge, mg.nodes_at_top_edge,
    ...              mg.nodes_at_right_edge):
    ...     mg.status_at_node[edge] = CLOSED_BOUNDARY

    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> th = mg.add_zeros('node', 'channel_sediment__depth')
    >>> th += 0.0007

    >>> fr = FlowRouter(mg)
    >>> sde = SedDepEroder(mg, K_sp=1.e-4,
    ...                    sed_dependency_type='almost_parabolic',
    ...                    Qc='power_law', K_t=1.e-4, external_sediment=True)

    >>> z[:] = mg.node_y/10000.

    >>> initz = z.copy()

    >>> dt = 100.
    >>> up = 0.05

    >>> for i in range(1):
    ...     fr.run_one_step()
    ...     sde.run_one_step(dt)

    Where TL conditions predominate, the incision is linked only to movement
    of the sediment layer:

    >>> TLs = mg.core_nodes[sde.is_it_TL[mg.core_nodes]]
    >>> np.allclose(th[TLs] + (initz - z)[TLs], 0.0007)
    True

    Otherwise, incision procedes in the naked nodes according to the sediment
    dependent bedrock incision rules:

    >>> incising_nodes = mg.core_nodes[
    ...     np.logical_not(sde.is_it_TL)[mg.core_nodes]]
    >>> np.all((initz - z)[incising_nodes] > 0.0007)
    True
    """

    _name = 'SedDepEroder'

    _input_var_names = (
        'topographic__elevation',
        'drainage_area',
        'flow__receiver_node',
        'flow__upstream_node_order',
        'topographic__steepest_slope',
        'flow__link_to_receiver_node',
        'flow__sink_flag',
        'channel_sediment__depth'
    )

    _output_var_names = (
        'topographic__elevation',
        'channel_sediment__depth',
        'channel__bed_shear_stress',
        'channel_sediment__volumetric_transport_capacity',
        'channel_sediment__volumetric_flux',
        'channel_sediment__relative_flux',
        'channel__discharge',
        'channel__width',  # optional
        'channel__depth',  # optional
    )

    _optional_var_names = (
        'channel__width',
        'channel__depth'
    )

    _var_units = {'topographic__elevation': 'm',
                  'drainage_area': 'm**2',
                  'flow__receiver_node': '-',
                  'topographic__steepest_slope': '-',
                  'flow__upstream_node_order': '-',
                  'flow__link_to_receiver_node': '-',
                  'flow__sink_flag': '-',
                  'channel_sediment__depth': 'm',
                  'channel__bed_shear_stress': 'Pa',
                  'channel_sediment__volumetric_transport_capacity': 'm**3/s',
                  'channel_sediment__volumetric_flux': 'm**3/s',
                  'channel_sediment__relative_flux': '-',
                  'channel__discharge': 'm**3/s',
                  'channel__width': 'm',
                  'channel__depth': 'm'
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'drainage_area': 'node',
                    'flow__receiver_node': 'node',
                    'topographic__steepest_slope': 'node',
                    'flow__upstream_node_order': 'node',
                    'flow__link_to_receiver_node': 'node',
                    'flow__sink_flag': 'node',
                    'channel_sediment__depth': 'node',
                    'channel__bed_shear_stress': 'node',
                    'channel_sediment__volumetric_transport_capacity': 'node',
                    'channel_sediment__volumetric_flux': 'node',
                    'channel_sediment__relative_flux': 'node',
                    'channel__discharge': 'node',
                    'channel__width': 'node',
                    'channel__depth': 'node'
                    }

    _var_type = {'topographic__elevation': float,
                 'drainage_area': float,
                 'flow__receiver_node': int,
                 'topographic__steepest_slope': float,
                 'flow__upstream_node_order': int,
                 'flow__link_to_receiver_node': int,
                 'flow__sink_flag': bool,                            # CHECK
                 'channel_sediment__depth': float,
                 'channel__bed_shear_stress': float,
                 'channel_sediment__volumetric_transport_capacity': float,
                 'channel_sediment__volumetric_flux': float,
                 'channel_sediment__relative_flux': float,
                 'channel__discharge': float,
                 'channel__width': float,
                 'channel__depth': float
                 }

    _var_doc = {
        'topographic__elevation': 'Land surface topographic elevation',
        'drainage_area':
            ("Upstream accumulated surface area contributing to the node's " +
             "discharge"),
        'flow__receiver_node':
            ('Node array of receivers (node that receives flow from current ' +
             'node)'),
        'topographic__steepest_slope':
            'Node array of steepest *downhill* slopes',
        'flow__upstream_node_order':
            ('Node array containing downstream-to-upstream ordered list of ' +
             'node IDs'),
        'flow__link_to_receiver_node':
            'ID of link downstream of each node, which carries the discharge',
        'flow__sink_flag': 'Boolean array, True at local lows',
        'channel_sediment__depth':
            ('Loose fluvial sediment at each node. Can be' +
             ' freely entrained by the flow, and must be to permit erosion. ' +
             'Note that the sediment is assumed to be distributed across the' +
             ' whole cell area.'),
        'channel__bed_shear_stress':
            ('Shear exerted on the bed of the channel, assuming all ' +
             'discharge travels along a single, self-formed channel'),
        'channel_sediment__volumetric_transport_capacity':
            ('Volumetric transport capacity of a channel carrying all runoff' +
             ' through the node, assuming the Meyer-Peter Muller transport ' +
             'equation'),
        'channel_sediment__volumetric_flux':
            ('Total volumetric fluvial sediment flux brought into the node ' +
             'from upstream'),
        'channel_sediment__relative_flux':
            ('The fluvial_sediment_flux_into_node divided by the fluvial_' +
             'sediment_transport_capacity'),
        'channel__discharge':
            ('Volumetric water flux of the a single channel carrying all ' +
             'runoff through the node'),
        'channel__width':
            ('Width of the a single channel carrying all runoff through the ' +
             'node'),
        'channel__depth':
            ('Depth of the a single channel carrying all runoff through the ' +
             'node')
    }

    def __init__(self, grid, K_sp=1.e-6, g=9.81,
                 rock_density=2700, sediment_density=2700, fluid_density=1000,
                 runoff_rate=1.,
                 sed_dependency_type='generalized_humped', kappa_hump=13.683,
                 nu_hump=1.13, phi_hump=4.24, c_hump=0.00181,
                 Qc='power_law', m_sp=0.5, n_sp=1., K_t=1.e-4, m_t=1.5, n_t=1.,
                 # params for model numeric behavior:
                 pseudoimplicit_repeats=50, return_stream_properties=False,
                 **kwds):
        """Constructor for the class."""
        self._grid = grid
        self.pseudoimplicit_repeats = pseudoimplicit_repeats

        self._K_unit_time = K_sp/31557600.
        # ^...because we work with dt in seconds
        # set gravity
        self.g = g
        self.rock_density = rock_density
        self.sed_density = sediment_density
        self.fluid_density = fluid_density
        self.relative_weight = (
            (self.sed_density-self.fluid_density)/self.fluid_density*self.g)
        # ^to accelerate MPM calcs
        self.rho_g = self.fluid_density*self.g
        self.type = sed_dependency_type
        assert self.type in ('generalized_humped', 'None', 'linear_decline',
                             'almost_parabolic')
        # now conditional inputs
        if self.type == 'generalized_humped':
            self.kappa = kappa_hump
            self.nu = nu_hump
            self.phi = phi_hump
            self.c = c_hump
            self.norm = None
        else:
            self.kappa = 0.
            self.nu = 0.
            self.phi = 0.
            self.c = 0.
            self.norm = 0.
        # set the sed flux fn for later on:
        self.set_sed_flux_fn_gen()

        self.Qc = Qc
        assert self.Qc in ('MPM', 'power_law')
        self.return_ch_props = return_stream_properties
        if return_stream_properties:
            assert self.Qc == 'MPM', ("Qc must be 'MPM' to return stream " +
                                      "properties")
        if type(runoff_rate) in (float, int):
            self.runoff_rate = float(runoff_rate)
        elif type(runoff_rate) is str:
            self.runoff_rate = self.grid.at_node[runoff_rate]
        else:
            self.runoff_rate = np.array(runoff_rate)
            assert runoff_rate.size == self.grid.number_of_nodes

        if self.Qc == 'MPM':
            raise TypeError('MPM is no longer a permitted value for Qc!')
        elif self.Qc == 'power_law':
            self._m = m_sp
            self._n = n_sp
            self._Kt = K_t/31557600.  # in sec
            self._mt = m_t
            self._nt = n_t

        self._hillslope_sediment_flux_wzeros = self.grid.zeros('node',
                                                               dtype=float)
        external_sediment = True  # now hardwired
        if external_sediment is True:
            self._ext_sed = True
            try:
                self._hillslope_sediment = self.grid.at_node[
                    'channel_sediment__depth']  # a field was present
            except FieldError:
                self._hillslope_sediment = self.grid.add_zeros(
                    'node', 'channel_sediment__depth')
        else:
            self._ext_sed = False

        self.cell_areas = np.empty(grid.number_of_nodes)
        self.cell_areas.fill(np.mean(grid.area_of_cell))
        self.cell_areas[grid.node_at_cell] = grid.area_of_cell
        self._pseudoimplicit_aborts = 0
        self._total_DL_calls = 0
        self._error_at_abort = []

        # set up the necessary fields:
        self.initialize_output_fields()
        if self.return_ch_props:
            self.initialize_optional_output_fields()

    def set_sed_flux_fn_gen(self):
        """
        Sets the property self.sed_flux_fn_gen that controls which sed flux
        function to use elsewhere in the component.
        """
        if self.type == 'generalized_humped':
            # work out the normalization param:
            max_val = 0.
            for i in range(0, 1, 0.001):
                sff = sed_flux_fn_gen_genhump(
                    i, self.kappa, self.nu, self.c, self.phi, 1.)
                max_val = max((sff, max_val))
            self.norm = 1./max_val
            self.sed_flux_fn_gen = sed_flux_fn_gen_genhump
        elif self.type == 'None':
            self.sed_flux_fn_gen = sed_flux_fn_gen_const
        elif self.type == 'linear_decline':
            self.sed_flux_fn_gen = sed_flux_fn_gen_lindecl
        elif self.type == 'almost_parabolic':
            self.sed_flux_fn_gen = sed_flux_fn_gen_almostparabolic

    def erode(self, dt, flooded_depths=None, **kwds):
        """Erode and deposit on the channel bed for a duration of *dt*.

        Erosion occurs according to the sediment dependent rules specified
        during initialization.

        Parameters
        ----------
        dt : float (years, only!)
            Timestep for which to run the component.
        flooded_depths : array or field name (m)
            Depths of flooding at each node, zero where no lake. Note that the
            component will dynamically update this array as it fills nodes
            with sediment (...but does NOT update any other related lake
            fields).
        """
        grid = self.grid
        node_z = grid.at_node['topographic__elevation']
        node_A = grid.at_node['drainage_area']
        flow_receiver = grid.at_node['flow__receiver_node']
        s_in = grid.at_node['flow__upstream_node_order']
        node_S = grid.at_node['topographic__steepest_slope']
        elev_less_sed = node_z - self._hillslope_sediment

        dt_secs = dt * 31557600.

        if type(flooded_depths) is str:
            flooded_depths = self.grid.at_node[flooded_depths]
            # also need a map of initial flooded conds:
            flooded_nodes = flooded_depths > 0.
        elif type(flooded_depths) is np.ndarray:
            assert flooded_depths.size == self.grid.number_of_nodes
            flooded_nodes = flooded_depths > 0.
            # need an *updateable* record of the pit depths
        else:
            # if None, handle in loop
            flooded_nodes = None
            flooded_depths = None
        steepest_link = 'flow__link_to_receiver_node'
        link_length = np.empty(grid.number_of_nodes, dtype=float)
        link_length.fill(np.nan)
        draining_nodes = np.not_equal(grid.at_node[steepest_link],
                                      BAD_INDEX_VALUE)
        core_draining_nodes = np.intersect1d(np.where(draining_nodes)[0],
                                             grid.core_nodes,
                                             assume_unique=True)
        link_length[core_draining_nodes] = grid._length_of_link_with_diagonals[
            grid.at_node[steepest_link][core_draining_nodes]]

        # calc fluxes from hillslopes into each node:
        # we're going to budget the sediment coming IN as part of the fluvial
        # budget, on the basis this flux will be higher than the out-flux in
        # a convergent flow node
        # ...worst that can happen is that the river can't move it, and dumps
        # it in the first node
        # Take care to add this sed to the cover, but not to actually move it.
        # turn depth into a supply flux:
        sed_in_cells = self._hillslope_sediment[self.grid.node_at_cell]
        flux_in_cells = sed_in_cells * self.grid.area_of_cell / dt_secs
        self._hillslope_sediment_flux_wzeros[
            self.grid.node_at_cell] = flux_in_cells
        self._voldroprate = self.grid.zeros('node', dtype=float)
        self._hillslope_sediment.fill(0.)
        # ^ this will get refilled below. For now, the sed has been "fully
        # mobilised" off the bed; at the end of the step it can resettle.

        if self.Qc == 'power_law':
            transport_capacity_prefactor_withA = self._Kt * node_A**self._mt
            erosion_prefactor_withA = self._K_unit_time * node_A**self._m
            # ^doesn't include S**n*f(Qc/Qc)
            downward_slopes = node_S.clip(0.)
            # for the stability condition:
            if np.isclose(self._n, 1.):
                wave_stab_cond_denominator = erosion_prefactor_withA
            else:
                wave_stab_cond_denominator = (
                    erosion_prefactor_withA * downward_slopes**(self._n - 1.))
                # this is a bit cheeky; we're basically assuming the slope
                # won't change much as the tstep evolves
                # justifiable as the -1 suppresses the sensitivity to S
                # for now, we're going to be conservative and just exclude the
                # role of the value of the sed flux fn in modifying the calc;
                # i.e., we assume f(Qs,Qc) = 1 in the stability calc.

            max_tstep_wave = 0.7 * np.nanmin(link_length /
                                             wave_stab_cond_denominator)
            # ^adding additional scaling per CHILD, tho CHILD uses 0.2
            self.wave_denom = wave_stab_cond_denominator
            self.link_length = link_length
            self.max_tstep_wave = max_tstep_wave

            t_elapsed_internal = 0.
            break_flag = False
            rel_sed_flux = np.empty_like(node_A)

            dzbydt = np.zeros(grid.number_of_nodes, dtype=float)
            self._loopcounter = 0
            while 1:
                flood_node = None  # int if timestep is limited by flooding
                conv_factor = 0.8
                flood_tstep = dt_secs
                slopes_tothen = downward_slopes**self._n
                slopes_tothent = downward_slopes**self._nt
                transport_capacities = (transport_capacity_prefactor_withA *
                                        slopes_tothent)
                erosion_prefactor_withS = (
                    erosion_prefactor_withA * slopes_tothen)  # no time, no fqs

                river_volume_flux_into_node = np.zeros(grid.number_of_nodes,
                                                       dtype=float)
                dzbydt.fill(0.)
                cell_areas = self.cell_areas

                self._is_it_TL = np.zeros(
                    self.grid.number_of_nodes, dtype=bool)

                iterate_downstream(s_in, cell_areas,
                                   self._hillslope_sediment,
                                   self._hillslope_sediment_flux_wzeros,
                                   river_volume_flux_into_node,
                                   transport_capacities,
                                   erosion_prefactor_withS,
                                   rel_sed_flux, self._is_it_TL,
                                   self._voldroprate, flow_receiver,
                                   self.pseudoimplicit_repeats,
                                   dzbydt, self.sed_flux_fn_gen,
                                   self.kappa, self.nu, self.c,
                                   self.phi, self.norm,
                                   flooded_nodes=flooded_nodes,
                                   flooded_depths=flooded_depths)

                # now perform a CHILD-like convergence-based stability test:
                ratediff = dzbydt[flow_receiver] - dzbydt
                # if this is +ve, the nodes are converging
                downstr_vert_diff = node_z - node_z[flow_receiver]
                # ^+ve when dstr is lower
                botharepositive = np.logical_and(ratediff > 0.,
                                                 downstr_vert_diff > 0.)
                try:
                    t_to_converge = np.amin(
                        downstr_vert_diff[botharepositive] /
                        ratediff[botharepositive])
                except ValueError:  # no node pair converges
                    t_to_converge = dt_secs
                t_to_converge *= conv_factor
                # ^arbitrary safety factor; CHILD uses 0.3
                # check this is a more restrictive condition than Courant:
                t_to_converge = min((t_to_converge, max_tstep_wave))
                if t_to_converge < 3600. and flood_node is not None:
                    t_to_converge = 3600.  # forbid tsteps < 1hr; a bit hacky
                # without this, it's possible for the component to get stuck in
                # a loop, presumably when the gradient is "supposed" to level
                # out. We make exception got nodes that need to be filled in
                # "just so"
                # the new handling of flooded nodes as of 25/10/16 should make
                # this redundant
                this_tstep = min((t_to_converge, dt_secs))
                t_elapsed_internal += this_tstep
                if t_elapsed_internal >= dt_secs:
                    break_flag = True
                    t_to_converge = dt_secs - t_elapsed_internal + this_tstep
                    self.t_to_converge = t_to_converge
                    this_tstep -= t_elapsed_internal - dt_secs

                # back-calc the sed budget in the nodes, as appropriate:
                if self._ext_sed:
                    self._hillslope_sediment[self.grid.node_at_cell] = (
                        self._voldroprate[self.grid.node_at_cell] /
                        self.grid.area_of_cell * this_tstep)
                # modify elevs; both sed & rock:
                # note dzbydt applies only to the ROCK surface
                elev_less_sed[grid.core_nodes] += (
                    dzbydt[grid.core_nodes] * this_tstep)
                node_z[grid.core_nodes] = (
                    elev_less_sed[grid.core_nodes] +
                    self._hillslope_sediment[grid.core_nodes])

                if break_flag:
                    break
                else:
                    self._loopcounter += 1
                # do we need to reroute the flow/recalc the slopes here?
                # -> NO, slope is such a minor component of Diff we'll be OK
                # BUT could be important not for the stability, but for the
                # actual calc. So YES.
                node_S = np.zeros_like(node_S)
                # print link_length[core_draining_nodes]
                node_S[core_draining_nodes] = (node_z-node_z[flow_receiver])[
                    core_draining_nodes]/link_length[core_draining_nodes]
                downward_slopes = node_S.clip(0.)
        else:
            raise TypeError  # should never trigger

        if self.return_ch_props:
            # add the channel property field entries,
            # 'channel__width', 'channel__depth', and 'channel__discharge'
            W = self.k_w*node_Q**self._b
            H = shear_stress/self.rho_g/node_S  # ...sneaky!
            grid.at_node['channel__width'][:] = W
            grid.at_node['channel__depth'][:] = H
            grid.at_node['channel__discharge'][:] = node_Q
            grid.at_node['channel__bed_shear_stress'][:] = shear_stress

        grid.at_node['channel_sediment__volumetric_transport_capacity'][
            :] = transport_capacities
        grid.at_node['channel_sediment__volumetric_flux'][
            :] = river_volume_flux_into_node
        # ^note this excludes the hillslope fluxes now, i.e., it's the
        # incoming sed flux to the node. Including the local flux in would
        # make little sense as it could easily exceed capacity.
        grid.at_node['channel_sediment__relative_flux'][:] = rel_sed_flux
        # elevs set automatically to the name used in the function call.

        return grid, grid.at_node['topographic__elevation']

    def run_one_step(self, dt, flooded_depths=None, **kwds):
        """Run the component across one timestep increment, dt.

        Erosion occurs according to the sediment dependent rules specified
        during initialization. Method is fully equivalent to the :func:`erode`
        method.

        Parameters
        ----------
        dt : float (years, only!)
            Timestep for which to run the component.
        flooded_depths : array or field name (m)
            Depths of flooding at each node, zero where no lake. Note that the
            component will dynamically update this array as it fills nodes
            with sediment (...but does NOT update any other related lake
            fields).
        """
        self.erode(dt=dt, flooded_depths=flooded_depths, **kwds)

    @property
    def is_it_TL(self):
        """Return a map of where erosion is purely transport-limited.
        """
        return self._is_it_TL

