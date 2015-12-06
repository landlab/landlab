from __future__ import print_function

import numpy as np
import inspect
from landlab import ModelParameterDictionary, CLOSED_BOUNDARY
from landlab import RasterModelGrid
from time import sleep
from landlab.utils import structured_grid as sgrid
import pylab

from landlab.core.model_parameter_dictionary import MissingKeyError
from landlab.field.scalar_data_fields import FieldError
from landlab.grid.base import BAD_INDEX_VALUE


class TransportLimitedEroder(object):
    """
    This component implements transport limited erosion for a landscape in which
    flow directions are fully convergent. i.e., all nodes in the landscape have
    a single, uniquely defined downstream node.

    The module can in principle take multiple transport laws, but at the moment
    only Meyer-Peter Muller (MPM) is implemented.

    This module is defined to interface neatly with the flow_routing module.

    Assumes grid does not deform or change during run.

    Note it is vital to ensure all your units match. t is assumed to be in
    years. Any length value is assumed to be in meters (including both dx
    and the uplift rate...!)

    DEJH Sept 14.
    (Only tested for raster grid so far)
    """

    def __init__(self, grid, params):
        self.initialize(grid, params)

    def initialize(self, grid, params_file):
        '''
        params_file is the name of the text file containing the parameters
        needed for this stream power component.

        ***Parameters for input file***
        OBLIGATORY:
            * Qc -> String. Controls how to set the carrying capacity.
                Either 'MPM', or a string giving the name of the model field
                where capacity values are stored on nodes.
                At the moment, only 'MPM' is permitted as a way to set the
                capacity automatically, but expansion would be trivial.
                If 'from_array', the module will attempt to set the capacity
                Note capacities must be specified as volume flux.
            *

            ...Then, assuming you set Qc=='MPM':
            * b_sp, c_sp -> Floats. These are the powers on discharge and
                drainage area in the equations used to control channel width and
                basin hydrology, respectively:
                        W = k_w * Q**b_sp
                        Q = k_Q * A**c_sp
                These parameters are used to constrain flow depth, and may be
                omitted if use_W or use_Q are set.
            *k_Q, k_w, mannings_n -> floats. These are the prefactors on the
                basin hydrology and channel width-discharge relations, and n
                from the Manning's equation, respectively. These are
                needed to allow calculation of shear stresses and hence carrying
                capacities from the local slope and drainage area alone.
                Don't know what to set these values to? k_w=2.5, k_Q=2.5e-7,
                mannings_n=0.05 give vaguely plausible numbers with b=0.5,
                c = 1.(e.g., for a drainage area ~350km2, like Boulder Creek
                at Boulder, => depth~1.3m, width~23m, shear stress ~O(200Pa)
                for an "annual-ish" flood). [If you want to continue playing
                with calibration, the ?50yr return time 2013 floods produced
                depths ~2.3m with Q~200m3/s]
            *Dchar -> float.  The characteristic grain diameter in meters
                (==D50 in most cases) used to calculate Shields numbers
                in the channel. If you want to define Dchar values at each node,
                don't set, and use the Dchar_if_used argument in erode()
                instead.

        OPTIONS:
            *rock_density -> in kg/m3 (defaults to 2700)
            *sediment_density -> in kg/m3 (defaults to 2700)
            *fluid_density -> in most cases water density, in kg/m3 (defaults to 1000)
            *g -> acceleration due to gravity, in m/s**2 (defaults to 9.81)

            *threshold_shields -> +ve float; the threshold taustar_crit.
                Defaults to 0.047, or if 'slope_sensitive_threshold' is set True,
                becomes a weak function of local slope following Lamb et al
                (2008):
                    threshold_shields=0.15*S**0.25
            *slope_sensitive_threshold -> bool, defaults to 'False'.
                If true, threshold_shields is set according to the Lamb
                equation,
                An exception will be raised if threshold_shields is
                also set.
            *dt -> +ve float. If set, this is the fixed timestep for this
                component. Can be overridden easily as a parameter in erode().
                If not set (default), this parameter MUST be set in erode().
            *use_W -> Bool; if True, component will look for node-centered data
                describing channel width in grid.at_node['channel_width'], and
                use it to implement incision ~ stream power per unit width.
                Defaults to False. NOT YET IMPLEMENTED
            *use_Q -> Bool. Overrides the basin hydrology relation, using an
                local water discharge value assumed already calculated and
                stored in grid.at_node['discharge']. NOT YET IMPLEMENTED
            *C_MPM -> float. Defaults to 1. Allows tuning of the MPM prefactor,
                which is calculated as
                    Qc = 8.*C_MPM*(taustar - taustarcrit)**1.5
                In almost all cases, tuning depth_equation_prefactor' is
                preferred to tuning this parameter.
            *return_stream_properties -> bool (default False).
                If True, this component will save the calculations for
                'channel_width', 'channel_depth', and 'channel_discharge' in
                those grid fields. (Requires some
                additional math, so is suppressed for speed by default).

        '''
        # this is the fraction we allow any given slope in the grid to evolve
        # by in one go (suppresses numerical instabilities)
        self.fraction_gradient_change = 0.25
        self.grid = grid
        # needs to be filled with values in execution
        self.link_S_with_trailing_blank = np.zeros(grid.number_of_links + 1)
        self.count_active_links = np.zeros_like(
            self.link_S_with_trailing_blank, dtype=int)
        self.count_active_links[:-1] = 1
        inputs = ModelParameterDictionary(params_file)
        try:
            self.g = inputs.read_float('g')
        except MissingKeyError:
            self.g = 9.81
        try:
            self.rock_density = inputs.read_float('rock_density')
        except MissingKeyError:
            self.rock_density = 2700.
        try:
            self.sed_density = inputs.read_float('sediment_density')
        except MissingKeyError:
            self.sed_density = 2700.
        try:
            self.fluid_density = inputs.read_float('fluid_density')
        except MissingKeyError:
            self.fluid_density = 1000.
        self.rho_g = self.fluid_density * self.g

        try:
            self.Qc = inputs.read_string('Qc')
        except MissingKeyError:
            raise MissingKeyError("Qc must be 'MPM' or a grid field name!")
        else:
            if self.Qc == 'MPM':
                self.calc_cap_flag = True
            else:
                self.calc_cap_flag = False
        try:
            self.return_ch_props = inputs.read_bool('return_stream_properties')
        except MissingKeyError:
            self.return_ch_props = False

        try:
            self.lamb_flag = inputs.read_bool('slope_sensitive_threshold')
        except:
            self.lamb_flag = False
        try:
            self.shields_crit = inputs.read_float('threshold_shields')
            # flag for sed_flux_dep_incision to see if the threshold was
            # manually set.
            self.set_threshold = True
            print("Found a threshold to use: ", self.shields_crit)
            assert self.lamb_flag == False
        except MissingKeyError:
            if not self.lamb_flag:
                self.shields_crit = 0.047
            self.set_threshold = False
        try:
            self.tstep = inputs.read_float('dt')
        except MissingKeyError:
            pass
        try:
            self.use_W = inputs.read_bool('use_W')
        except MissingKeyError:
            self.use_W = False
        try:
            self.use_Q = inputs.read_bool('use_Q')
        except MissingKeyError:
            self.use_Q = False
        try:
            self.return_capacity = inputs.read_bool('return_capacity')
        except MissingKeyError:
            self.return_capacity = False

        try:
            self._b = inputs.read_float('b_sp')
        except MissingKeyError:
            if self.use_W:
                self._b = 0.
            else:
                if self.calc_cap_flag:
                    raise NameError('b was not set')
        try:
            self._c = inputs.read_float('c_sp')
        except MissingKeyError:
            if self.use_Q:
                self._c = 1.
            else:
                if self.calc_cap_flag:
                    raise NameError('c was not set')
        try:
            self.Dchar_in = inputs.read_float('Dchar')
        except MissingKeyError:
            pass

        # assume Manning's equation to set the power on A for shear stress:
        self.shear_area_power = 0.6 * self._c * (1. - self._b)

        self.k_Q = inputs.read_float('k_Q')
        self.k_w = inputs.read_float('k_w')
        mannings_n = inputs.read_float('mannings_n')
        self.mannings_n = mannings_n
        if mannings_n < 0. or mannings_n > 0.2:
            print("***STOP. LOOK. THINK. You appear to have set Manning's n outside its typical range. Did you mean it? Proceeding...***")
            sleep(2)

        try:
            self.C_MPM = inputs.read_float('C_MPM')
        except MissingKeyError:
            self.C_MPM = 1.
        self.diffusivity_power_on_A = 0.9 * self._c * \
            (1. - self._b)  # i.e., q/D**(1/6)

        # new for v3:
        # set thresh in shear stress if poss at this stage:
        try:  # fails if no Dchar provided, or shields crit is being set dynamically from slope
            self.thresh = self.shields_crit * \
                (self.sed_density - self.fluid_density) * self.g * self.Dchar_in
        except AttributeError:
            try:
                self.shields_prefactor_to_shear = (
                    self.sed_density - self.fluid_density) * self.g * self.Dchar_in
            except AttributeError:  # no Dchar
                self.shields_prefactor_to_shear_noDchar = (
                    self.sed_density - self.fluid_density) * self.g
        twothirds = 2. / 3.
        self.Qs_prefactor = 4. * self.C_MPM**twothirds * self.fluid_density**twothirds / (self.sed_density - self.fluid_density)**twothirds * self.g**(
            twothirds / 2.) * mannings_n**0.6 * self.k_w**(1. / 15.) * self.k_Q**(0.6 + self._b / 15.) / self.sed_density**twothirds
        self.Qs_thresh_prefactor = 4. * (self.C_MPM * self.k_w * self.k_Q**self._b / self.fluid_density**0.5 / (
            self.sed_density - self.fluid_density) / self.g / self.sed_density)**twothirds
        # both these are divided by sed density to give a vol flux
        self.Qs_power_onA = self._c * (0.6 + self._b / 15.)
        self.Qs_power_onAthresh = twothirds * self._b * self._c

        if RasterModelGrid in inspect.getmro(grid.__class__):
            self.cell_areas = grid.dx * grid.dy
        else:
            self.cell_areas = np.empty(grid.number_of_nodes)
            self.cell_areas.fill(np.mean(grid.area_of_cell))
            self.cell_areas[grid.node_at_cell] = grid.area_of_cell
        self.bad_neighbor_mask = np.equal(
            grid.get_active_neighbors_at_node(bad_index=-1), -1)

        self.routing_code = """
            double sed_flux_into_this_node;
            double sed_flux_out_of_this_node;
            double flux_excess;
            for (int i=len_s_in; i>0; i--) {
                sed_flux_into_this_node = sed_into_node[i];
                sed_flux_out_of_this_node = transport_capacities[i];
                flux_excess = sed_flux_into_this_node - sed_flux_out_of_this_node;
                dz[i] = flux_excess/cell_areas*dt_this_step;
                sed_into_node[flow_receiver[i]] += sed_flux_out_of_this_node;
            }
            """

    def erode(self, grid, dt=None, node_elevs='topographic__elevation',
              node_drainage_areas='drainage_area',
              node_receiving_flow='flow_receiver',
              node_order_upstream='upstream_node_order',
              node_slope='topographic__steepest_slope',
              steepest_link='links_to_flow_receiver',
              runoff_rate_if_used=None,
              # W_if_used=None, Q_if_used=None,
              stability_condition='loose',
              Dchar_if_used=None, io=None):
        """
        This method calculates the fluvial sediment transport capacity at all
        nodes, then erodes or deposits sediment as required.

        *grid* & *dt* are the grid object and timestep (float) respectively.

        *node_elevs* tells the component where to look for the node elevations.
        Pass another string to override which grid field the component looks
        at, or pass a nnodes-long array of elevation values directly instead.

        *node_drainage_areas* tells the component where to look for the drainage
        area values. Change to another string to override which grid field the
        component looks at, or pass a nnodes-long array of drainage areas values
        directly instead.

        *node_receiving flow* tells the component where to look for the node
        ids which receive flow from each node. This is an output from the
        flow_routing module.
        Change to another string to override which grid field the
        component looks at, or pass a nnodes-long array of IDs
        directly instead.

        Alternatively, set *link_slopes* (and *link_node_mapping*) if this data
        is only available at links. 'planet_surface__derivative_of_elevation'
        is the default field name for link slopes. Override this name by
        setting the variable as the appropriate string, or override use of
        grid fields altogether by passing an array. *link_node_mapping* controls
        how the component maps these link values onto the arrays. We assume
        there is always a 1:1 mapping (pass the values already projected onto
        the nodes using slopes_at_nodes if not). Other components, e.g.,
        flow_routing.route_flow_dn, may provide the necessary outputs to make
        the mapping easier: e.g., just pass 'links_to_flow_reciever' from that
        module (the default name). If the component cannot find an existing
        mapping through this parameter, it will derive one on the fly, at
        considerable cost of speed (see on-screen reports).

        *slopes_from_elevs* allows the module to create gradients internally
        from elevations rather than have them provided. Set to True to force
        the component to look for the data in grid.at_node['topographic__elevation'];
        set to 'name_of_field' to override this name, or pass an nnode-array
        to use those values as elevations instead. Using this option is
        considerably slower than any of the alternatives, as it also has to
        calculate the link_node_mapping from stratch each time.

        In both these cases, at present the mapping is to use the maximum
        slope of --any-- link attached to the node as the representative node
        slope. This is primarily for speed, but may be a good idea to modify
        later.

        *runoff_rate_if_used* is the runoff rate in m/yr, if used (take care
        with the units...!). Either a float, or an nnodes-long array
        (NB: array functionality is untested).
        If specified, becomes a multiplicative modifier on k_Q, above.
        Ensure you adjust the precipitation time series so that the
        flood series you get makes sense! If not set, the precip rate
        is assumed to be rolled into the k_Q term already.

        *W_if_used* and *Q_if_used* must be provided if you set use_W and use_Q
        respectively in the component initialization. They can be either field
        names or nnodes arrays as in the other cases.

        *stability_condition* controls how we limit the internal timestep of
        the model to improve solution stability. 'loose' (default) employs a
        condition that prevents changes in the drainage structure of the
        existing channel network - slopes may not reverse - and this should
        prove adequate for most uses. 'tight' uses a considerably stricter Lax/
        Von Neumann criterion, but will be considerably slower.

        *Dchar_if_used* must be set as a grid field string or nnoodes-long array
        if 'Dchar' as a float was not provided in the input file. (If it was,
        this will be overridden).

        SETS: (as fields on the grid)
        ***Note the time units are SECONDS in these fields***
        - 'topographic__elevation' (m), the elevations (or your name)
        - 'fluvial_sediment_transport_capacity' (m**3/s), the volumetric
            transport capacities at each node
        - 'fluvial_sediment_flux_into_node' (m**3/s), the total volumetric sed
            flux entering each node
        If your stability_condition was 'tight':
        - 'effective_fluvial_diffusivity'
        If return_stream_properties was True:
        - 'channel_width' (m)
        - 'channel_depth' (m)
        - 'channel_discharge' (m**3/s)
        - 'channel_bed_shear_stress' (Pa)
        - The number of internal loops required per supplied time step is
            stored as a property of the class instance, self.iterations_in_dt.
            This can be useful for tracking computational load imposed by this
            component.

        RETURNS:
        grid, elevations

        """

        if runoff_rate_if_used != None:
            runoff_rate = runoff_rate_if_used
            assert type(runoff_rate) in (int, float, np.ndarray)
        else:
            runoff_rate = 1.

        if dt == None:
            dt = self.tstep
        try:
            self.Dchar = self.Dchar_in
        except AttributeError:
            try:
                self.Dchar = grid.at_node[Dchar_if_used]
            except FieldError:
                assert type(Dchar_if_used) == np.ndarray
                self.Dchar = Dchar_if_used

        if type(node_elevs) == str:
            node_z = grid.at_node[node_elevs]
        else:
            node_z = node_elevs

        if type(node_drainage_areas) == str:
            node_A = grid.at_node[node_drainage_areas]
        else:
            node_A = node_drainage_areas

        if type(node_receiving_flow) == str:
            flow_receiver = grid.at_node[node_receiving_flow]
        else:
            flow_receiver = node_receiving_flow

        # new V3:
        if type(node_order_upstream) == str:
            s_in = grid.at_node[node_order_upstream]
        else:
            s_in = node_order_upstream

        if type(node_slope) == str:
            node_S = grid.at_node[node_slope]
        else:
            node_S = node_slope

        if self.lamb_flag:
            variable_shields_crit = 0.15 * node_S**0.25
            try:
                variable_thresh = variable_shields_crit * self.shields_prefactor_to_shear
            except AttributeError:
                variable_thresh = variable_shields_crit * \
                    self.shields_prefactor_to_shear_noDchar * self.Dchar

        if type(steepest_link) == str:
            link_length = np.empty(grid.number_of_nodes, dtype=float)
            link_length.fill(np.nan)
            draining_nodes = np.not_equal(
                grid.at_node[steepest_link], BAD_INDEX_VALUE)
            core_draining_nodes = np.intersect1d(
                np.where(draining_nodes)[0], grid.core_nodes, assume_unique=True)
            link_length[core_draining_nodes] = grid.link_length[
                grid.at_node[steepest_link][core_draining_nodes]]
            # link_length=grid.dx
        else:
            link_length = grid.link_length[steepest_link]
        square_link_length = np.square(link_length)  # nans propagate forward

        try:
            transport_capacities_thresh = self.thresh * self.Qs_thresh_prefactor * \
                runoff_rate**(0.66667 * self._b) * \
                node_A**self.Qs_power_onAthresh
        except AttributeError:
            transport_capacities_thresh = variable_thresh * self.Qs_thresh_prefactor * \
                runoff_rate**(0.66667 * self._b) * \
                node_A**self.Qs_power_onAthresh

        transport_capacity_prefactor_withA = self.Qs_prefactor * \
            runoff_rate**(0.6 + self._b / 15.) * node_A**self.Qs_power_onA

        internal_t = 0.
        break_flag = False
        dt_secs = dt * 31557600.
        counter = 0

        while 1:  # use the break flag, to improve computational efficiency for runs which are very stable
            # we assume the drainage structure is forbidden to change during the whole dt
            # print "loop..."
            # note slopes will be *negative* at pits
            # track how many loops we perform:
            counter += 1
            downward_slopes = node_S.clip(0.)
            #positive_slopes = np.greater(downward_slopes, 0.)
            transport_capacities_S = transport_capacity_prefactor_withA * \
                (downward_slopes)**0.7
            trp_diff = (transport_capacities_S -
                        transport_capacities_thresh).clip(0.)
            transport_capacities = np.sqrt(trp_diff * trp_diff * trp_diff)

            if stability_condition == 'tight':
                mock_diffusivities = np.zeros_like(
                    transport_capacities, dtype=float)
                mock_diffusivities = transport_capacities / downward_slopes
                # we're relaxing the condition fivefold here, as the true
                # VonNeumann condition is VERY restrictive
                tstep_each_node = 10. * square_link_length / mock_diffusivities
                # if no node exceeds crit, tstep_each_node will just be nans
                # and infs
                # in seconds, nanmin avoids the pit nodes
                delta_t_internal = np.nanmin(tstep_each_node)
                if delta_t_internal == np.inf:  # no node exceeds crit
                    # nothing happened, so let the loop complete, awaiting more
                    # uplift
                    delta_t_internal = dt_secs
                if internal_t + delta_t_internal >= dt_secs:
                    dt_this_step = dt_secs - internal_t  # now in seconds
                    break_flag = True
                else:
                    # a min tstep was found (seconds). We terminate the loop
                    dt_this_step = delta_t_internal
            else:  # loose, gradient based method
                # and the adjustment is made AFTER the dz calc
                dt_this_step = dt_secs - internal_t

            sed_into_node = np.zeros(grid.number_of_nodes, dtype=float)
            dz = np.zeros(grid.number_of_nodes, dtype=float)
            len_s_in = s_in.size
            cell_areas = self.cell_areas

            for i in s_in[::-1]:  # work downstream
                sed_flux_into_this_node = sed_into_node[i]
                # we work in volume flux, not volume per se here
                sed_flux_out_of_this_node = transport_capacities[i]
                flux_excess = sed_flux_into_this_node - \
                    sed_flux_out_of_this_node  # gets deposited
                dz[i] = flux_excess / cell_areas * dt_this_step
                sed_into_node[flow_receiver[i]] += sed_flux_out_of_this_node

            if stability_condition == 'loose':
                elev_diff = node_z - node_z[flow_receiver]
                delta_dz = dz[flow_receiver] - dz
                # note the condition is that gradient may not change by >X%,
                # not must be >0
                node_flattening = self.fraction_gradient_change * elev_diff - delta_dz
                # note all these things are zero for a pit node
                most_flattened_nodes = np.argmin(
                    node_flattening[grid.core_nodes])
                # get it back to node number, not core_node number
                most_flattened_nodes = np.take(
                    grid.core_nodes, most_flattened_nodes)
                most_flattened_val = np.take(
                    node_flattening, most_flattened_nodes)
                if most_flattened_val >= 0.:
                    break_flag = True  # all nodes are stable
                else:  # a fraction < 1
                    dt_fraction = self.fraction_gradient_change * \
                        np.take(elev_diff, most_flattened_nodes) / \
                        np.take(delta_dz, most_flattened_nodes)
                    # print dt_fraction
                    # correct those elevs
                    dz *= dt_fraction
                    dt_this_step *= dt_fraction

            # print np.amax(dz), np.amin(dz)

            node_z[grid.core_nodes] += dz[grid.core_nodes]

            if break_flag:
                break
            # do we need to reroute the flow/recalc the slopes here? -> NO, slope is such a minor component of Diff we'll be OK
            # BUT could be important not for the stability, but for the actual
            # calc. So YES.
            node_S = np.zeros_like(node_S)
            # print link_length[core_draining_nodes]
            node_S[core_draining_nodes] = (
                node_z - node_z[flow_receiver])[core_draining_nodes] / link_length[core_draining_nodes]
            internal_t += dt_this_step  # still in seconds, remember

        self.grid = grid

        active_nodes = np.where(grid.status_at_node != CLOSED_BOUNDARY)[0]
        if io:
            try:
                io[active_nodes] += node_z[active_nodes]
            except TypeError:
                if type(io) == str:
                    elev_name = io
            else:
                return grid, io

        else:
            elev_name = node_elevs

        if self.return_ch_props:
            # add the channel property field entries,
            #'channel_width', 'channel_depth', and 'channel_discharge'
            Q = self.k_Q * runoff_rate * node_A**self._c
            W = self.k_w * Q**self._b
            H = Q**(0.6 * (1. - self._b)) * \
                (self.mannings_n / self.k_w)**0.6 * node_S**-0.3
            tau = self.fluid_density * self.g * H * node_S
            grid.at_node['channel_width'] = W
            grid.at_node['channel_depth'] = H
            grid.at_node['channel_discharge'] = Q
            grid.at_node['channel_bed_shear_stress'] = tau

        grid.at_node[
            'fluvial_sediment_transport_capacity'] = transport_capacities
        grid.at_node['fluvial_sediment_flux_into_node'] = sed_into_node
        # elevs set automatically to the name used in the function call.
        if stability_condition == 'tight':
            grid.at_node['effective_fluvial_diffusivity'] = mock_diffusivities
        self.iterations_in_dt = counter

        return grid, grid.at_node[elev_name]
